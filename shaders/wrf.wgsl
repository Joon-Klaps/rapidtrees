// Weighted Robinson-Foulds distance via f32 L1 on dense length rows.
//
// Each thread (i, j) computes WRF(i, j) = Σ |rows[i][k] - rows[j][k]| over all splits.
// Branch lengths are stored as f32 (downcast from the f64 CPU representation).
// Missing splits have length 0.0. No cancellation — precision is safe.

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
    var dist: f32 = 0.0;
    let base_i = i * params.n_splits;
    let base_j = j * params.n_splits;
    for (var k: u32 = 0u; k < params.n_splits; k = k + 1u) {
        dist = dist + abs(rows[base_i + k] - rows[base_j + k]);
    }
    matrix[i * params.n + j] = dist;
}
