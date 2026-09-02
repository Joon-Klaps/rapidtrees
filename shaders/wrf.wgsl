// Weighted Robinson-Foulds distance: f32 L1 over shared-candidate columns.
//
// One thread per (i, j) cell:
//   WRF(i, j) = Σ_k |rows[i][k] - rows[j][k]| + unique[i] + unique[j]
//
// `rows` is an n × stride matrix covering only the splits that at least two
// trees hold; a split held by one tree alone can never overlap, so it is left
// out of the matrix and its length is folded into that tree's `unique` term
// instead (|a - 0| = a for non-negative lengths). That filter is what keeps the
// buffer inside WebGPU's 128 MiB storage-binding limit at realistic tree
// counts — it cuts ~48 000 split columns to ~2 600 on diverse posteriors.
//
// Equivalent to the CPU's `Σlᵢ + Σlⱼ - 2·Σ min(lᵢ, lⱼ)` form, which assumes
// non-negative branch lengths. Missing lengths are 0.0. The direct |a - b|
// form is used rather than that sum-minus-overlap identity because the latter
// cancels catastrophically in f32 on near-identical trees.
//
// Buffers come from `rapidtrees::gpu_layout::weighted_rows`.

struct Params {
    n: u32,
    stride: u32,
    _pad0: u32,
    _pad1: u32,
}

@group(0) @binding(0) var<storage, read>       rows:   array<f32>;
@group(0) @binding(1) var<storage, read>       unique: array<f32>;
@group(0) @binding(2) var<uniform>             params: Params;
@group(0) @binding(3) var<storage, read_write> matrix: array<f32>;

@compute @workgroup_size(16, 16, 1)
fn main(@builtin(global_invocation_id) gid: vec3<u32>) {
    let i = gid.x;
    let j = gid.y;
    if i >= params.n || j >= params.n { return; }
    if i == j {
        matrix[i * params.n + j] = 0.0;
        return;
    }
    var dist: f32 = 0.0;
    let base_i = i * params.stride;
    let base_j = j * params.stride;
    for (var k: u32 = 0u; k < params.stride; k = k + 1u) {
        dist = dist + abs(rows[base_i + k] - rows[base_j + k]);
    }
    matrix[i * params.n + j] = dist + unique[i] + unique[j];
}
