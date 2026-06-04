// Kuhner-Felsenstein (branch score) distance via direct f32 squared differences.
//
// Each thread (i, j) computes KF(i, j) = sqrt(Σ (rows[i][k] - rows[j][k])²).
// Uses the direct form — NOT the Gram reformulation (norm²+norm²-2·dot) — to avoid
// catastrophic cancellation on near-identical trees where the bracket rounds to
// a large negative number in f32. Direct form is safe at ~1e-5 relative error.

struct Params {
    n: u32,
    n_splits: u32,
    _pad0: u32,
    _pad1: u32,
}

@group(0) @binding(0) var<storage, read>       rows:   array<f32>;
@group(0) @binding(1) var<uniform>             params: Params;
@group(0) @binding(2) var<storage, read_write> matrix: array<f32>;

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
    let base_i = i * params.n_splits;
    let base_j = j * params.n_splits;
    for (var k: u32 = 0u; k < params.n_splits; k = k + 1u) {
        let d = rows[base_i + k] - rows[base_j + k];
        sum_sq = sum_sq + d * d;
    }
    matrix[i * params.n + j] = sqrt(sum_sq);
}
