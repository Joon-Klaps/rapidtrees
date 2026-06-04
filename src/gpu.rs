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

/// Upper bound on total GPU buffer allocation (1 GiB); above this we fall back to CPU.
const MAX_GPU_BYTES: u64 = 1 << 30;

// ── Global lazy device ────────────────────────────────────────────────────────

static GPU: OnceLock<Option<GpuContext>> = OnceLock::new();

fn get_gpu() -> Option<&'static GpuContext> {
    GPU.get_or_init(GpuContext::try_new).as_ref()
}

// ── GpuContext ────────────────────────────────────────────────────────────────

pub(crate) struct GpuContext {
    device: wgpu::Device,
    queue: wgpu::Queue,
    rf_pipeline: wgpu::ComputePipeline,
    wrf_pipeline: wgpu::ComputePipeline,
    kf_pipeline: wgpu::ComputePipeline,
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

        let adapter = instance
            .request_adapter(&wgpu::RequestAdapterOptions {
                power_preference: wgpu::PowerPreference::HighPerformance,
                force_fallback_adapter: false,
                compatible_surface: None,
            })
            .await?;

        // Reject software/CPU adapters — they're slower than the native CPU path.
        if adapter.get_info().device_type == wgpu::DeviceType::Cpu {
            return None;
        }

        let (device, queue) = adapter
            .request_device(
                &wgpu::DeviceDescriptor {
                    label: None,
                    required_features: wgpu::Features::empty(),
                    required_limits: wgpu::Limits::downlevel_defaults(),
                    ..Default::default()
                },
                None,
            )
            .await
            .ok()?;

        let rf_pipeline = make_pipeline(&device, include_str!("../shaders/rf.wgsl"), "rf");
        let wrf_pipeline = make_pipeline(&device, include_str!("../shaders/wrf.wgsl"), "wrf");
        let kf_pipeline = make_pipeline(&device, include_str!("../shaders/kf.wgsl"), "kf");

        Some(Self {
            device,
            queue,
            rf_pipeline,
            wrf_pipeline,
            kf_pipeline,
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
    if kept_count + everywhere == 0 {
        return Some(vec![0usize; n * n]);
    }
    let words_u32 = kept_count.div_ceil(32);

    let mut packed = vec![0u32; n * words_u32];
    if words_u32 > 0 {
        for (row, snap) in packed.chunks_mut(words_u32).zip(&snaps.snapshots) {
            for &id in &snap.split_ids {
                let slot = bit_slot[id as usize];
                if slot != u32::MAX {
                    let slot = slot as usize;
                    row[slot >> 5] |= 1u32 << (slot & 31);
                }
            }
        }
    }

    let kept_per_tree: Vec<u32> = snaps
        .snapshots
        .iter()
        .map(|snap| (snap.split_ids.len() - everywhere) as u32)
        .collect();

    // Estimate buffer sizes and bail if they exceed a safe limit.
    let packed_bytes = (n * words_u32 * 4) as u64;
    let output_bytes = (n * n * 4) as u64;
    if packed_bytes + output_bytes > MAX_GPU_BYTES {
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

    let total_bytes = (n * n_splits * 4 + n * n * 4) as u64;
    if total_bytes > MAX_GPU_BYTES {
        return None;
    }

    let rows_f32: Vec<f32> = rows_f64.iter().map(|&v| v as f32).collect();
    ctx.run_length_metric(select(ctx), &rows_f32, n, n_splits)
        .map(|v| v.into_iter().map(|x| x as f64).collect())
}

// ── GpuContext compute helpers ────────────────────────────────────────────────

impl GpuContext {
    fn run_rf(&self, packed: &[u32], kept: &[u32], n: usize, words: usize) -> Option<Vec<usize>> {
        let params: [u32; 4] = [n as u32, words as u32, 0, 0];

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
