//! Parallelism shim.
//!
//! Every parallel loop in this crate goes through `use crate::par::*` rather
//! than `use rayon::prelude::*`. With the `parallel` feature on (the default)
//! that *is* rayon's prelude and nothing changes. With it off, the same method
//! names resolve to sequential standard-library iterators.
//!
//! This exists for wasm32, but the reason is narrower than "rayon does not work
//! in a browser". Rayon compiles for `wasm32-unknown-unknown` today, and every
//! parallel loop in this crate runs threaded there. What does not work is the
//! *default pool*: `rayon-core` spawns workers through `std::thread::Builder`,
//! which `wasm32-unknown-unknown` does not implement, so with `parallel` on and
//! nothing else done the first `par_iter` panics.
//!
//! Building the pool is therefore the caller's job, and deliberately not this
//! crate's. A browser worker is `new Worker(...)` plus handing it
//! `wasm_bindgen::memory()` so it shares the `SharedArrayBuffer` — irreducibly
//! JS interop, which is why it lives in `treetracer-core` (via
//! `wasm-bindgen-rayon`) rather than here. That crate builds rayon's *global*
//! pool, after which these call sites pick it up with no API between the two.
//!
//! Consumers therefore have three configurations, two of which work:
//!
//! | `parallel` | Pool built by caller | Result |
//! | --- | --- | --- |
//! | off | — | the sequential shim below; builds anywhere, no setup |
//! | on | yes | threaded; needs `+atomics`, `-Z build-std`, and a cross-origin-isolated page |
//! | on | no | **panics** on the first parallel loop |
//!
//! The threaded build needs nightly — `atomics` is still an unstable
//! `-C target-feature`. See the `check-wasm-threads` pixi task for the exact
//! invocation CI holds this to.
//!
//! There is deliberately no `thread_count()` helper here: `rayon`'s
//! `current_num_threads` *initialises* the global registry, so calling it to
//! ask whether a pool exists is precisely what triggers the panic above. The
//! consumer knows its own pool width from the `init_thread_pool` it called.
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
