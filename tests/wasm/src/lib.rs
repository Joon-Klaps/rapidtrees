//! A `rapidtrees` build that runs on `wasm32-unknown-unknown` without `wasm-bindgen`.
//!
//! The host (`run.ts`) copies a BEAST `.trees` file into linear memory with [`alloc`], calls [`run`], and reads the three distance matrices back through [`matrix`]. That path goes through `rapidtrees::io::read_trees`, the reader every binding shares, with the options the CLI defaults to, and then the pairwise kernels. A rayon pool, a `std::fs` call or a `std::time::Instant` that crept into that path compiles fine for wasm32 but traps here, which is the failure `cargo check` cannot see.
//!
//! The ABI is plain `extern "C"` on purpose: rapidtrees must not depend on `wasm-bindgen` (see CLAUDE.md), and a module with no imports instantiates in Node with nothing but `WebAssembly.instantiate`.

use rapidtrees::io::{ReadOptions, ReadTreesError, read_trees};
use std::sync::Mutex;

/// Matrices from the last [`run`]: RF, weighted RF, KF, each row-major `n * n`.
static MATRICES: Mutex<Vec<Vec<f64>>> = Mutex::new(Vec::new());

/// Reserve `len` bytes for the host to write the input into. The buffer is leaked: one file per instance.
#[unsafe(no_mangle)]
pub extern "C" fn alloc(len: usize) -> *mut u8 {
    Vec::<u8>::with_capacity(len).leak().as_mut_ptr()
}

/// Parse the BEAST file at `ptr..ptr + len` and compute all three pairwise matrices.
///
/// Returns the number of trees, or `-1` if the input is not UTF-8, `-2` if it holds fewer than two trees and `-3` if the trees fail to parse.
///
/// # Safety
///
/// `ptr` must come from [`alloc`] with at least `len` bytes, all of them written by the host.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn run(ptr: *const u8, len: usize) -> i32 {
    // SAFETY: guaranteed by the caller, per the contract above.
    let bytes = unsafe { std::slice::from_raw_parts(ptr, len) };
    let snaps = match read_trees(bytes, &ReadOptions::default()) {
        Ok((_, snaps)) => snaps,
        Err(ReadTreesError::Io(_)) => return -1,
        Err(ReadTreesError::NoTrees(_) | ReadTreesError::TooFewTrees { .. }) => return -2,
        Err(_) => return -3,
    };
    let rf = snaps.pairwise_rf(None).into_iter().map(f64::from).collect();
    let matrices = vec![rf, snaps.pairwise_wrf(None), snaps.pairwise_kf(None)];
    *MATRICES.lock().unwrap_or_else(|e| e.into_inner()) = matrices;
    i32::try_from(snaps.len()).unwrap_or(i32::MAX)
}

/// Pointer to matrix `metric` (0 = RF, 1 = weighted RF, 2 = KF) from the last [`run`], or null.
#[unsafe(no_mangle)]
pub extern "C" fn matrix(metric: usize) -> *const f64 {
    MATRICES
        .lock()
        .unwrap_or_else(|e| e.into_inner())
        .get(metric)
        .map_or(std::ptr::null(), |m| m.as_ptr())
}
