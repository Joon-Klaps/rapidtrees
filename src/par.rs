//! Parallelism shim.
//!
//! Every parallel loop in this crate goes through `use crate::par::*` rather
//! than `use rayon::prelude::*`. With the `parallel` feature on (the default)
//! that *is* rayon's prelude and nothing changes. With it off, the same method
//! names resolve to sequential standard-library iterators.
//!
//! This exists for wasm32. `rayon-core` spawns workers through
//! `std::thread::Builder` with no wasm special-casing, so on
//! `wasm32-unknown-unknown` the pool fails to build and the first `par_iter`
//! panics — the crate compiles but cannot run. Real threads in the browser need
//! `wasm-bindgen-rayon`, which in turn needs `SharedArrayBuffer` and a
//! cross-origin-isolated page. Being able to build without rayon at all means
//! the browser build can ship before that machinery is in place.
//!
//! Rayon's parallel iterators deliberately mirror the standard iterator API, so
//! the call sites are identical either way: `enumerate`, `zip`, `map`,
//! `for_each` and `collect` all exist on both. The sequential path is also
//! strictly more permissive — rayon requires closures to be `Send + Sync`,
//! plain iterators do not.

#[cfg(feature = "parallel")]
pub use rayon::prelude::*;

#[cfg(not(feature = "parallel"))]
pub use sequential::*;

#[cfg(not(feature = "parallel"))]
mod sequential {
    use std::slice::{ChunksMut, Iter};

    /// Sequential stand-in for `rayon::prelude::ParallelSliceMut`.
    pub trait ParallelSliceMut<T> {
        fn par_chunks_mut(&mut self, chunk_size: usize) -> ChunksMut<'_, T>;
    }

    impl<T> ParallelSliceMut<T> for [T] {
        #[inline]
        fn par_chunks_mut(&mut self, chunk_size: usize) -> ChunksMut<'_, T> {
            self.chunks_mut(chunk_size)
        }
    }

    /// Sequential stand-in for `rayon::prelude::IntoParallelRefIterator`.
    pub trait ParallelSlice<T> {
        fn par_iter(&self) -> Iter<'_, T>;
    }

    impl<T> ParallelSlice<T> for [T] {
        #[inline]
        fn par_iter(&self) -> Iter<'_, T> {
            self.iter()
        }
    }
}
