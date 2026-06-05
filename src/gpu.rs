//! GPU-backed pairwise distance computation via wgpu.
//!
//! Entry points are [`try_pairwise_rf`], [`try_pairwise_wrf`], [`try_pairwise_kf`].
//! Each returns `None` when no suitable GPU adapter is found or the problem is too
//! small to justify the transfer overhead, causing the caller to fall back to the
//! CPU path transparently.
//!
//! The GPU context ([`GpuContext`]) is created lazily on the first call and reused
//! for the lifetime of the process. Whether to use it at all is decided by the
//! caller via the `use_gpu` argument on the public API — there are no environment
//! overrides.

use std::sync::OnceLock;
use wgpu::util::DeviceExt;

use crate::distances::{dense_length_rows, split_bit_layout};
use crate::snapshot::Snapshots;

/// Minimum number of trees to justify the GPU transfer overhead.
pub(crate) const GPU_THRESHOLD: usize = 64;

/// Conservative upper bound on the work of a *single* pairwise dispatch, in
/// `n² × inner` units (`inner` = packed words for RF, split count for WRF/KF).
///
/// The whole n×n matrix is computed in one `queue.submit`. On a headless compute
/// GPU a single dispatch that runs too long trips the driver watchdog, which loses
/// the device — and wgpu *panics* on submit to a lost device, so it can't be caught
/// and recovered. Above this bound we fall back to the CPU (always correct, just
/// slower). Empirical: RF at 5000²×378 ≈ 9.5e9 completes in ~2.9 s on a Tesla P100
/// and is safe; the WRF dense path at 5000²×2400 ≈ 6e10 loses the device. The bound
/// sits above the proven-safe point with margin. A future chunked/banded dispatch
/// (bounded rows per submit) would remove this cap entirely.
const MAX_GPU_DISPATCH_WORK: u64 = 1.5e10 as u64;

// ── Global lazy device ────────────────────────────────────────────────────────

static GPU: OnceLock<Option<GpuContext>> = OnceLock::new();

fn get_gpu() -> Option<&'static GpuContext> {
    // The OnceLock guarantees this runs at most once, so the fallback warning is
    // emitted only on the first request and never repeats. Reached only when the
    // caller actually requested the GPU and the problem cleared GPU_THRESHOLD.
    GPU.get_or_init(|| {
        let ctx = GpuContext::try_new();
        if ctx.is_none() {
            eprintln!(
                "rapidtrees: GPU requested but no compatible GPU adapter was found; \
                 falling back to the CPU. On Linux this usually means the Vulkan \
                 loader/ICD is missing (wgpu has no CUDA backend)."
            );
        }
        ctx
    })
    .as_ref()
}

/// Human-readable label for the active GPU adapter (e.g. `Tesla P100-SXM2-16GB
/// (Vulkan, DiscreteGpu)`), or `None` if no compatible GPU is available.
///
/// Forces the lazy GPU context to initialise, so the first call performs adapter
/// selection (and emits the same one-time diagnostics as a real distance call).
/// Used by callers that want to confirm up front whether the GPU will actually be
/// used instead of inferring it from timings.
pub(crate) fn adapter_label() -> Option<String> {
    get_gpu().map(|ctx| ctx.label.clone())
}

// ── GpuContext ────────────────────────────────────────────────────────────────

pub(crate) struct GpuContext {
    device: wgpu::Device,
    queue: wgpu::Queue,
    rf_pipeline: wgpu::ComputePipeline,
    wrf_pipeline: wgpu::ComputePipeline,
    kf_pipeline: wgpu::ComputePipeline,
    /// Human-readable adapter description, e.g. `Tesla P100-SXM2-16GB (Vulkan, DiscreteGpu)`.
    label: String,
    /// Device limits granted at creation — used to reject buffers that would
    /// exceed `max_storage_buffer_binding_size` / `max_buffer_size` before dispatch.
    limits: wgpu::Limits,
}

impl GpuContext {
    fn try_new() -> Option<Self> {
        pollster::block_on(Self::init_async())
    }

    async fn init_async() -> Option<Self> {
        let instance = wgpu::Instance::new(wgpu::InstanceDescriptor {
            backends: wgpu::Backends::all(),
            ..Default::default()
        });

        let adapter = match instance
            .request_adapter(&wgpu::RequestAdapterOptions {
                power_preference: wgpu::PowerPreference::HighPerformance,
                force_fallback_adapter: false,
                compatible_surface: None,
            })
            .await
        {
            Some(adapter) => adapter,
            None => {
                // List what wgpu *did* enumerate so a failure here distinguishes
                // "no Vulkan device at all" from "device present but not selectable".
                for info in instance.enumerate_adapters(wgpu::Backends::all()) {
                    let info = info.get_info();
                    eprintln!(
                        "rapidtrees: wgpu saw adapter {} ({:?} backend, {:?}) but none was selected",
                        info.name, info.backend, info.device_type
                    );
                }
                return None;
            }
        };

        let info = adapter.get_info();
        // Reject software/CPU adapters — they're slower than the native CPU path.
        if info.device_type == wgpu::DeviceType::Cpu {
            eprintln!(
                "rapidtrees: only a software adapter ({}, {:?} backend) was found; \
                 using the CPU path instead",
                info.name, info.backend
            );
            return None;
        }

        // `downlevel_defaults()` caps `max_storage_buffer_binding_size` at 128 MiB.
        // Our packed-row and n×n output buffers blow past that on large inputs, which
        // wgpu rejects with a validation error mid-dispatch — and the readback then
        // hands back a garbage (all-zero) matrix that looks like a successful GPU run.
        // Request the adapter's real limits instead so large buffers are allowed where
        // the hardware supports them (NVIDIA Vulkan reports multi-GiB binding sizes);
        // we still guard every buffer against these limits before dispatch below.
        let (device, queue) = match adapter
            .request_device(
                &wgpu::DeviceDescriptor {
                    label: None,
                    required_features: wgpu::Features::empty(),
                    required_limits: adapter.limits(),
                    ..Default::default()
                },
                None,
            )
            .await
        {
            Ok(pair) => pair,
            Err(e) => {
                eprintln!(
                    "rapidtrees: GPU adapter {} found but request_device failed: {e}",
                    info.name
                );
                return None;
            }
        };

        // wgpu routes device errors (validation, device-lost, out-of-memory) to this
        // handler, which by default logs via the `log` crate. Callers that install no
        // logger would never see them, so the compute/readback would fail silently and
        // fall back to the CPU with no explanation. Print straight to stderr instead.
        device.on_uncaptured_error(Box::new(|err| {
            eprintln!("rapidtrees: wgpu device error: {err}");
        }));

        eprintln!(
            "rapidtrees: GPU active — {} ({:?} backend, {:?})",
            info.name, info.backend, info.device_type
        );

        let rf_pipeline = make_pipeline(&device, include_str!("../shaders/rf.wgsl"), "rf");
        let wrf_pipeline = make_pipeline(&device, include_str!("../shaders/wrf.wgsl"), "wrf");
        let kf_pipeline = make_pipeline(&device, include_str!("../shaders/kf.wgsl"), "kf");

        let limits = device.limits();
        Some(Self {
            device,
            queue,
            rf_pipeline,
            wrf_pipeline,
            kf_pipeline,
            label: format!("{} ({:?}, {:?})", info.name, info.backend, info.device_type),
            limits,
        })
    }
}

fn make_pipeline(device: &wgpu::Device, src: &str, label: &str) -> wgpu::ComputePipeline {
    let module = device.create_shader_module(wgpu::ShaderModuleDescriptor {
        label: Some(label),
        source: wgpu::ShaderSource::Wgsl(src.into()),
    });
    device.create_compute_pipeline(&wgpu::ComputePipelineDescriptor {
        label: Some(label),
        layout: None,
        module: &module,
        entry_point: "main",
        compilation_options: wgpu::PipelineCompilationOptions::default(),
        cache: None,
    })
}

// ── public(crate) entry points ────────────────────────────────────────────────

/// Try to compute all pairwise RF distances on the GPU.
///
/// Returns `None` if no GPU is available, the problem is below the size threshold,
/// or any GPU operation fails — letting the caller fall back to the CPU path.
pub(crate) fn try_pairwise_rf(snaps: &Snapshots) -> Option<Vec<usize>> {
    let n = snaps.snapshots.len();
    if n < GPU_THRESHOLD {
        return None;
    }

    let ctx = get_gpu()?;

    // Re-pack into u32 words for WGSL (the CPU path uses u64 words instead).
    let (bit_slot, kept_count, everywhere) = split_bit_layout(snaps);
    let words_u32 = kept_count.div_ceil(32);
    if words_u32 == 0 {
        // No non-universal splits remain — e.g. every tree is identical, so all
        // splits are universal. Universal splits cancel in RF, so every pairwise
        // distance is 0. Returning here also avoids handing wgpu a zero-length
        // storage buffer, which it rejects as a validation error.
        return Some(vec![0usize; n * n]);
    }

    let mut packed = vec![0u32; n * words_u32];
    for (row, snap) in packed.chunks_mut(words_u32).zip(&snaps.snapshots) {
        for &id in &snap.split_ids {
            let slot = bit_slot[id as usize];
            if slot != u32::MAX {
                let slot = slot as usize;
                row[slot >> 5] |= 1u32 << (slot & 31);
            }
        }
    }

    let kept_per_tree: Vec<u32> = snaps
        .snapshots
        .iter()
        .map(|snap| (snap.split_ids.len() - everywhere) as u32)
        .collect();

    // `packed` and the n×n output are each bound as a single storage buffer, so each
    // must satisfy the device's per-binding limit; otherwise the dispatch fails
    // validation and returns garbage. Fall back to the CPU when either won't fit.
    let packed_bytes = (n * words_u32 * 4) as u64;
    let output_bytes = (n * n * 4) as u64;
    if !ctx.storage_binding_fits(packed_bytes) || !ctx.storage_binding_fits(output_bytes) {
        return None;
    }

    // Guard against a single dispatch large enough to trip the GPU watchdog (which
    // would lose the device and panic on submit). Fall back to the CPU instead.
    if (n as u64) * (n as u64) * (words_u32 as u64) > MAX_GPU_DISPATCH_WORK {
        return None;
    }

    ctx.run_rf(&packed, &kept_per_tree, n, words_u32)
}

/// Try to compute all pairwise WRF distances on the GPU.
pub(crate) fn try_pairwise_wrf(snaps: &Snapshots) -> Option<Vec<f64>> {
    try_pairwise_length_metric(snaps, |ctx| &ctx.wrf_pipeline)
}

/// Try to compute all pairwise KF distances on the GPU.
pub(crate) fn try_pairwise_kf(snaps: &Snapshots) -> Option<Vec<f64>> {
    try_pairwise_length_metric(snaps, |ctx| &ctx.kf_pipeline)
}

/// Shared dense-length GPU path for WRF and KF; `select` picks which compute
/// pipeline to run (the two metrics differ only in the per-split accumulator).
fn try_pairwise_length_metric(
    snaps: &Snapshots,
    select: impl FnOnce(&GpuContext) -> &wgpu::ComputePipeline,
) -> Option<Vec<f64>> {
    let n = snaps.snapshots.len();
    if n < GPU_THRESHOLD {
        return None;
    }

    let ctx = get_gpu()?;
    let (rows_f64, n_splits) = dense_length_rows(snaps);
    if n_splits == 0 {
        return Some(vec![0.0f64; n * n]);
    }

    // The dense length rows and the n×n output are each a single storage binding;
    // each must fit the device limit or the dispatch fails. Fall back to CPU if not.
    let rows_bytes = (n * n_splits * 4) as u64;
    let output_bytes = (n * n * 4) as u64;
    if !ctx.storage_binding_fits(rows_bytes) || !ctx.storage_binding_fits(output_bytes) {
        return None;
    }

    // The dense path is O(n² × n_splits) in one dispatch — far heavier than RF — and
    // is what loses the device on large cells. Fall back to the CPU above the guard.
    if (n as u64) * (n as u64) * (n_splits as u64) > MAX_GPU_DISPATCH_WORK {
        return None;
    }

    let rows_f32: Vec<f32> = rows_f64.iter().map(|&v| v as f32).collect();
    ctx.run_length_metric(select(ctx), &rows_f32, n, n_splits)
        .map(|v| v.into_iter().map(|x| x as f64).collect())
}

// ── GpuContext compute helpers ────────────────────────────────────────────────

impl GpuContext {
    /// Whether a single storage buffer of `bytes` can be bound on this device.
    /// A buffer beyond `max_storage_buffer_binding_size` (or `max_buffer_size`)
    /// fails validation at bind-group creation, so the caller falls back to the CPU.
    fn storage_binding_fits(&self, bytes: u64) -> bool {
        bytes <= self.limits.max_storage_buffer_binding_size as u64
            && bytes <= self.limits.max_buffer_size
    }

    fn run_rf(&self, packed: &[u32], kept: &[u32], n: usize, words: usize) -> Option<Vec<usize>> {
        let params: [u32; 4] = [n as u32, words as u32, 0, 0];

        // Capture any validation error (bind-group/dispatch) raised below instead of
        // only logging it via on_uncaptured_error. A captured error means the readback
        // would be garbage, so we return None and let the caller fall back to the CPU.
        self.device.push_error_scope(wgpu::ErrorFilter::Validation);

        // `packed` is never empty here: try_pairwise_rf returns early when
        // words_u32 == 0, so this buffer always has at least one u32.
        let packed_buf = self
            .device
            .create_buffer_init(&wgpu::util::BufferInitDescriptor {
                label: Some("rf_packed"),
                contents: cast_u8(packed),
                usage: wgpu::BufferUsages::STORAGE,
            });
        let kept_buf = self
            .device
            .create_buffer_init(&wgpu::util::BufferInitDescriptor {
                label: Some("rf_kept"),
                contents: cast_u8(kept),
                usage: wgpu::BufferUsages::STORAGE,
            });
        let params_buf = self
            .device
            .create_buffer_init(&wgpu::util::BufferInitDescriptor {
                label: Some("rf_params"),
                contents: cast_u8(&params),
                usage: wgpu::BufferUsages::UNIFORM,
            });
        let output_size = (n * n * std::mem::size_of::<u32>()) as u64;
        let output_buf = self.device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("rf_output"),
            size: output_size,
            usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_SRC,
            mapped_at_creation: false,
        });
        let staging_buf = self.device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("rf_staging"),
            size: output_size,
            usage: wgpu::BufferUsages::MAP_READ | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });

        let layout = self.rf_pipeline.get_bind_group_layout(0);
        let bind_group = self.device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: None,
            layout: &layout,
            entries: &[
                wgpu::BindGroupEntry {
                    binding: 0,
                    resource: packed_buf.as_entire_binding(),
                },
                wgpu::BindGroupEntry {
                    binding: 1,
                    resource: kept_buf.as_entire_binding(),
                },
                wgpu::BindGroupEntry {
                    binding: 2,
                    resource: params_buf.as_entire_binding(),
                },
                wgpu::BindGroupEntry {
                    binding: 3,
                    resource: output_buf.as_entire_binding(),
                },
            ],
        });

        let workgroups = n.div_ceil(16) as u32;
        let mut encoder = self
            .device
            .create_command_encoder(&wgpu::CommandEncoderDescriptor { label: None });
        {
            let mut pass = encoder.begin_compute_pass(&wgpu::ComputePassDescriptor {
                label: None,
                ..Default::default()
            });
            pass.set_pipeline(&self.rf_pipeline);
            pass.set_bind_group(0, &bind_group, &[]);
            pass.dispatch_workgroups(workgroups, workgroups, 1);
        }
        encoder.copy_buffer_to_buffer(&output_buf, 0, &staging_buf, 0, output_size);
        self.queue.submit(std::iter::once(encoder.finish()));

        // If anything in the scope failed validation, the matrix is unreliable.
        if pollster::block_on(self.device.pop_error_scope()).is_some() {
            return None;
        }

        let bytes = readback_bytes(&self.device, &staging_buf, n * n * 4)?;
        Some(
            bytes
                .chunks_exact(4)
                .map(|b| u32::from_ne_bytes([b[0], b[1], b[2], b[3]]) as usize)
                .collect(),
        )
    }

    fn run_length_metric(
        &self,
        pipeline: &wgpu::ComputePipeline,
        rows: &[f32],
        n: usize,
        n_splits: usize,
    ) -> Option<Vec<f32>> {
        let params: [u32; 4] = [n as u32, n_splits as u32, 0, 0];

        // See run_rf: capture validation errors so a failed dispatch falls back to the
        // CPU instead of returning a garbage matrix.
        self.device.push_error_scope(wgpu::ErrorFilter::Validation);

        // `rows` is never empty here: try_pairwise_length_metric returns early when
        // n_splits == 0, so this buffer always holds n × n_splits f32 values.
        let rows_buf = self
            .device
            .create_buffer_init(&wgpu::util::BufferInitDescriptor {
                label: Some("len_rows"),
                contents: cast_u8(rows),
                usage: wgpu::BufferUsages::STORAGE,
            });
        let params_buf = self
            .device
            .create_buffer_init(&wgpu::util::BufferInitDescriptor {
                label: Some("len_params"),
                contents: cast_u8(&params),
                usage: wgpu::BufferUsages::UNIFORM,
            });
        let output_size = (n * n * std::mem::size_of::<f32>()) as u64;
        let output_buf = self.device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("len_output"),
            size: output_size,
            usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_SRC,
            mapped_at_creation: false,
        });
        let staging_buf = self.device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("len_staging"),
            size: output_size,
            usage: wgpu::BufferUsages::MAP_READ | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });

        let layout = pipeline.get_bind_group_layout(0);
        let bind_group = self.device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: None,
            layout: &layout,
            entries: &[
                wgpu::BindGroupEntry {
                    binding: 0,
                    resource: rows_buf.as_entire_binding(),
                },
                wgpu::BindGroupEntry {
                    binding: 1,
                    resource: params_buf.as_entire_binding(),
                },
                wgpu::BindGroupEntry {
                    binding: 2,
                    resource: output_buf.as_entire_binding(),
                },
            ],
        });

        let workgroups = n.div_ceil(16) as u32;
        let mut encoder = self
            .device
            .create_command_encoder(&wgpu::CommandEncoderDescriptor { label: None });
        {
            let mut pass = encoder.begin_compute_pass(&wgpu::ComputePassDescriptor {
                label: None,
                ..Default::default()
            });
            pass.set_pipeline(pipeline);
            pass.set_bind_group(0, &bind_group, &[]);
            pass.dispatch_workgroups(workgroups, workgroups, 1);
        }
        encoder.copy_buffer_to_buffer(&output_buf, 0, &staging_buf, 0, output_size);
        self.queue.submit(std::iter::once(encoder.finish()));

        // If anything in the scope failed validation, the matrix is unreliable.
        if pollster::block_on(self.device.pop_error_scope()).is_some() {
            return None;
        }

        let bytes = readback_bytes(&self.device, &staging_buf, n * n * 4)?;
        Some(
            bytes
                .chunks_exact(4)
                .map(|b| f32::from_ne_bytes([b[0], b[1], b[2], b[3]]))
                .collect(),
        )
    }
}

// ── Readback helpers ──────────────────────────────────────────────────────────

/// Map `buf`, wait for the GPU, and copy out exactly `byte_len` bytes. Returns
/// `None` if the map fails or the readback length doesn't match, letting the
/// caller fall back to the CPU path. Callers reinterpret the bytes themselves.
fn readback_bytes(device: &wgpu::Device, buf: &wgpu::Buffer, byte_len: usize) -> Option<Vec<u8>> {
    let slice = buf.slice(..);
    let (tx, rx) = std::sync::mpsc::channel();
    slice.map_async(wgpu::MapMode::Read, move |r| {
        let _ = tx.send(r);
    });
    device.poll(wgpu::Maintain::Wait);
    rx.recv().ok()?.ok()?;
    let data = slice.get_mapped_range();
    let bytes = data.to_vec();
    drop(data);
    buf.unmap();
    (bytes.len() == byte_len).then_some(bytes)
}

// ── Byte casting ──────────────────────────────────────────────────────────────

fn cast_u8<T>(data: &[T]) -> &[u8] {
    // Safety: T is a plain numeric type (u32, f32). Reinterpreting its bytes as u8
    // is always valid — no padding, no references, no UB.
    unsafe { std::slice::from_raw_parts(data.as_ptr() as *const u8, std::mem::size_of_val(data)) }
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::snapshot::Snapshots;

    const T0: &str = "(A:0.1,(B:0.1,(H:0.1,(D:0.1,(J:0.1,(((G:0.1,E:0.1):0.1,(F:0.1,I:0.1):0.1):0.1,C:0.1):0.1):0.1):0.1):0.1):0.1);";
    const T1: &str = "(A:0.1,(B:0.1,(D:0.1,((J:0.1,H:0.1):0.1,(((G:0.1,E:0.1):0.1,(F:0.1,I:0.1):0.1):0.1,C:0.1):0.1):0.1):0.1):0.1);";
    const T2: &str = "(A:0.1,(B:0.1,(D:0.1,(H:0.1,(J:0.1,(((G:0.1,E:0.1):0.1,(F:0.1,I:0.1):0.1):0.1,C:0.1):0.1):0.1):0.1):0.1):0.1);";

    fn three_snaps() -> Snapshots {
        Snapshots::from_newicks(&[T0, T1, T2], false).unwrap()
    }

    #[test]
    fn gpu_try_new_does_not_panic() {
        // GpuContext::try_new() must never panic — returning None is fine.
        let _ = GpuContext::try_new();
    }

    #[test]
    fn try_pairwise_rf_below_threshold_returns_none() {
        // 3 trees < GPU_THRESHOLD — must always return None regardless of GPU presence.
        let snaps = three_snaps();
        assert!(try_pairwise_rf(&snaps).is_none());
    }

    #[test]
    fn try_pairwise_wrf_below_threshold_returns_none() {
        let snaps = three_snaps();
        assert!(try_pairwise_wrf(&snaps).is_none());
    }

    #[test]
    fn try_pairwise_kf_below_threshold_returns_none() {
        let snaps = three_snaps();
        assert!(try_pairwise_kf(&snaps).is_none());
    }

    #[test]
    fn gpu_rf_matches_cpu_when_available() {
        // If a GPU adapter is available, RF results must be bit-identical to the CPU path.
        // If no GPU, test is a no-op.
        use crate::distances::pairwise_rf_packed;

        // Build a set above the threshold using repeated trees.
        let mut trees = vec![T0; GPU_THRESHOLD];
        trees.push(T1);
        trees.push(T2);
        let snaps = Snapshots::from_newicks(&trees, false).unwrap();

        let gpu_result = match try_pairwise_rf(&snaps) {
            Some(r) => r,
            None => return, // no GPU on this machine — skip
        };
        let cpu_result = pairwise_rf_packed(&snaps, None);
        assert_eq!(gpu_result, cpu_result, "GPU RF disagrees with CPU RF");
    }

    #[test]
    fn gpu_wrf_within_tolerance_when_available() {
        use crate::distances::pairwise_wrf_dense;

        let mut trees = vec![T0; GPU_THRESHOLD];
        trees.push(T1);
        trees.push(T2);
        let snaps = Snapshots::from_newicks(&trees, false).unwrap();

        let gpu_result = match try_pairwise_wrf(&snaps) {
            Some(r) => r,
            None => return,
        };
        let cpu_result = pairwise_wrf_dense(&snaps, None);
        let n = snaps.snapshots.len();
        for i in 0..n {
            for j in 0..n {
                let diff = (gpu_result[i * n + j] - cpu_result[i * n + j]).abs();
                assert!(
                    diff < 1e-4,
                    "GPU WRF[{i}][{j}] = {} vs CPU = {}; diff = {}",
                    gpu_result[i * n + j],
                    cpu_result[i * n + j],
                    diff,
                );
            }
        }
    }

    #[test]
    fn gpu_kf_within_tolerance_when_available() {
        use crate::distances::pairwise_kf_dense;

        let mut trees = vec![T0; GPU_THRESHOLD];
        trees.push(T1);
        trees.push(T2);
        let snaps = Snapshots::from_newicks(&trees, false).unwrap();

        let gpu_result = match try_pairwise_kf(&snaps) {
            Some(r) => r,
            None => return,
        };
        let cpu_result = pairwise_kf_dense(&snaps, None);
        let n = snaps.snapshots.len();
        for i in 0..n {
            for j in 0..n {
                let diff = (gpu_result[i * n + j] - cpu_result[i * n + j]).abs();
                assert!(
                    diff < 1e-4,
                    "GPU KF[{i}][{j}] = {} vs CPU = {}; diff = {}",
                    gpu_result[i * n + j],
                    cpu_result[i * n + j],
                    diff,
                );
            }
        }
    }
}
