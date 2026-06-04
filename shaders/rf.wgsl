// Robinson-Foulds distance via bit-packed popcount.
//
// Each thread (i, j) computes RF(i, j) = kept[i] + kept[j] - 2 * popcount(row_i & row_j).
// The packed matrix stores one bit per kept (non-universal) split; each row is `words` u32s.
// Universal splits (present in every tree) are dropped before packing — they cancel in RF.

struct Params {
    n: u32,
    words: u32,  // u32 words per row
    _pad0: u32,
    _pad1: u32,
}

@group(0) @binding(0) var<storage, read>       packed: array<u32>;
@group(0) @binding(1) var<storage, read>       kept:   array<u32>;
@group(0) @binding(2) var<uniform>             params: Params;
@group(0) @binding(3) var<storage, read_write> matrix: array<u32>;

@compute @workgroup_size(16, 16, 1)
fn main(@builtin(global_invocation_id) gid: vec3<u32>) {
    let i = gid.x;
    let j = gid.y;
    if i >= params.n || j >= params.n { return; }
    if i == j {
        matrix[i * params.n + j] = 0u;
        return;
    }
    var n_common: u32 = 0u;
    let base_i = i * params.words;
    let base_j = j * params.words;
    for (var k: u32 = 0u; k < params.words; k = k + 1u) {
        n_common = n_common + countOneBits(packed[base_i + k] & packed[base_j + k]);
    }
    matrix[i * params.n + j] = kept[i] + kept[j] - 2u * n_common;
}
