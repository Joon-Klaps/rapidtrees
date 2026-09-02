// Robinson-Foulds distance via bit-packed popcount.
//
// One thread per (i, j) cell:
//   RF(i, j) = kept[i] + kept[j] - 2 * popcount(row_i & row_j)
//
// `packed` holds one bit per *kept* split, `words` u32s per row. Splits held by
// every tree are dropped before packing: they add equally to both `kept` terms
// and to the shared count, so they cancel exactly. WGSL has no u64, so rows are
// packed 32 bits to a word here where the CPU path in `distances.rs` uses 64.
//
// Buffers come from `rapidtrees::gpu_layout::split_bit_rows`.

struct Params {
    n: u32,
    words: u32,
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
    var shared_count: u32 = 0u;
    let base_i = i * params.words;
    let base_j = j * params.words;
    for (var k: u32 = 0u; k < params.words; k = k + 1u) {
        shared_count = shared_count + countOneBits(packed[base_i + k] & packed[base_j + k]);
    }
    matrix[i * params.n + j] = kept[i] + kept[j] - 2u * shared_count;
}
