// Kuhner-Felsenstein (branch score) distance: f32 Euclidean over shared columns.
//
// One thread per (i, j) cell:
//   KF(i, j) = sqrt(Σ_k (rows[i][k] - rows[j][k])² + unique[i] + unique[j])
//
// `rows` covers only splits at least two trees hold; a split held by one tree
// alone contributes (a - 0)² = a², which is folded into that tree's `unique`
// term. `unique` is therefore already squared — see
// `rapidtrees::gpu_layout::weighted_rows`, which applies the same `term` the
// CPU path uses.
//
// Uses the direct squared-difference form, NOT the Gram reformulation
// (‖a‖² + ‖b‖² - 2a·b): in f32 that bracket cancels to a large negative on
// near-identical trees, and sqrt of it is NaN. Direct form holds ~1e-5
// relative error.
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
    var sum_sq: f32 = 0.0;
    let base_i = i * params.stride;
    let base_j = j * params.stride;
    for (var k: u32 = 0u; k < params.stride; k = k + 1u) {
        let d = rows[base_i + k] - rows[base_j + k];
        sum_sq = sum_sq + d * d;
    }
    matrix[i * params.n + j] = sqrt(sum_sq + unique[i] + unique[j]);
}
