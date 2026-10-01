//! Filling the symmetric output matrix, and small flat-array helpers.

use crate::par::*;
use std::sync::atomic::{AtomicUsize, Ordering};

/// Fill a symmetric `n × n` matrix a band of `band` rows at a time, one rayon
/// task per band.
///
/// `fill_band(i0, rows)` writes `rows[r * n + j]` for every row `i0 + r` of the
/// band and every `j > i0 + r`; with `band == 1` that is `rows[j]` of row `i0`.
/// `rows` arrives holding `T::default()` throughout, so a caller may accumulate
/// into it first. The last band may be shorter. The diagonal stays at
/// `T::default()` and the lower triangle is mirrored from the upper.
/// `progress` is bumped by each band's pair count as that band finishes.
pub(super) fn fill_symmetric<T, F>(
    n: usize,
    band: usize,
    progress: Option<&AtomicUsize>,
    fill_band: F,
) -> Vec<T>
where
    T: Copy + Default + Send,
    F: Fn(usize, &mut [T]) + Sync,
{
    debug_assert!(band > 0, "a band holds at least one row");
    let mut matrix = vec![T::default(); n * n];

    matrix
        .par_chunks_mut((n * band).max(1))
        .enumerate()
        .for_each(|(b, rows)| {
            let i0 = b * band;
            fill_band(i0, rows);
            if let Some(counter) = progress {
                let pairs = (i0..i0 + rows.len() / n).map(|i| n - i - 1).sum();
                counter.fetch_add(pairs, Ordering::Relaxed);
            }
        });

    // Mirror into the lower triangle in tiles: a plain row-read/column-write
    // sweep puts every write on its own cache line at large `n`.
    pub(super) const TILE: usize = 64;
    for i0 in (0..n).step_by(TILE) {
        for j0 in (i0..n).step_by(TILE) {
            for i in i0..(i0 + TILE).min(n) {
                for j in j0.max(i + 1)..(j0 + TILE).min(n) {
                    matrix[j * n + i] = matrix[i * n + j];
                }
            }
        }
    }

    matrix
}

/// Borrow row `i` of a flat, row-major matrix whose rows are `stride` wide.
#[inline]
pub(super) fn row_slice<T>(flat: &[T], i: usize, stride: usize) -> &[T] {
    &flat[i * stride..][..stride]
}

/// Replace each slot with the sum of the slots before it and return the total.
pub(super) fn exclusive_prefix_sum<'a, T>(slots: impl IntoIterator<Item = &'a mut T>) -> T
where
    T: 'a + Copy + Default + std::ops::Add<Output = T>,
{
    let mut total = T::default();
    for slot in slots {
        (*slot, total) = (total, total + *slot);
    }
    total
}
