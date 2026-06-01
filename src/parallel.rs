//! One small place where "parallel vs sequential" is decided.
//!
//! Native builds (CLI / Python) parallelise the hot loops with [`rayon`]. The
//! WebAssembly build cannot: `wasm32-unknown-unknown` has no threads unless the
//! page is cross-origin isolated (COOP/COEP headers + `SharedArrayBuffer`),
//! which a static GitHub Pages site can't provide. So the `wasm` build compiles
//! without the `rayon` feature and these helpers fall back to plain iterators.
//!
//! Centralising the `#[cfg]` here keeps every call site in `distances.rs` and
//! `snapshot.rs` identical and backend-agnostic — they just call
//! [`for_each_row_mut`] / [`try_map_indexed`] with no conditional compilation.

/// Apply `f` to each `row_len`-sized mutable chunk ("row") of `data`, passing
/// the row's index. Runs across rayon's thread pool when the `rayon` feature is
/// enabled, otherwise sequentially.
#[cfg(feature = "rayon")]
pub(crate) fn for_each_row_mut<T, F>(data: &mut [T], row_len: usize, f: F)
where
    T: Send,
    F: Fn(usize, &mut [T]) + Sync,
{
    use rayon::prelude::*;
    data.par_chunks_mut(row_len)
        .enumerate()
        .for_each(|(i, row)| f(i, row));
}

/// Sequential fallback used when the `rayon` feature is disabled (e.g. wasm).
#[cfg(not(feature = "rayon"))]
pub(crate) fn for_each_row_mut<T, F>(data: &mut [T], row_len: usize, f: F)
where
    F: Fn(usize, &mut [T]),
{
    data.chunks_mut(row_len)
        .enumerate()
        .for_each(|(i, row)| f(i, row));
}

/// Map `f` over `items` (with each item's index), collecting into a `Vec` and
/// short-circuiting on the first `Err`. Runs in parallel via rayon when the
/// `rayon` feature is enabled, otherwise sequentially.
#[cfg(feature = "rayon")]
pub(crate) fn try_map_indexed<T, R, E, F>(items: &[T], f: F) -> Result<Vec<R>, E>
where
    T: Sync,
    R: Send,
    E: Send,
    F: Fn(usize, &T) -> Result<R, E> + Sync,
{
    use rayon::prelude::*;
    items
        .par_iter()
        .enumerate()
        .map(|(i, item)| f(i, item))
        .collect()
}

/// Sequential fallback used when the `rayon` feature is disabled (e.g. wasm).
#[cfg(not(feature = "rayon"))]
pub(crate) fn try_map_indexed<T, R, E, F>(items: &[T], f: F) -> Result<Vec<R>, E>
where
    F: Fn(usize, &T) -> Result<R, E>,
{
    items
        .iter()
        .enumerate()
        .map(|(i, item)| f(i, item))
        .collect()
}
