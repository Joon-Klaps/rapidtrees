//! Runtime CPU dispatch for the pairwise row loops.
//!
//! A crate cannot choose the CPU it is compiled for. `cargo install` and the
//! PyPI wheels target the architecture's baseline, which on x86-64 predates
//! POPCNT and AVX2: RF's popcount compiles to an SSE2 bit trick and every
//! sweep runs at most two lanes wide. `-C target-cpu=native` fixes a build
//! made on the machine that runs it, and breaks a wheel on any CPU older than
//! the one that built it.
//!
//! So the one loop that matters, a row of the pairwise matrix, is compiled once
//! per x86-64 microarchitecture level, and the best level the running CPU
//! supports is picked once per process. Every level computes the same numbers:
//! the loops fix their own summation order, and Rust never fuses a multiply
//! and an add unless told to, so a wider vector changes the speed and not the
//! result.
//!
//! aarch64 needs none of this, since NEON, which carries the popcount, is part
//! of its baseline. wasm32 has no runtime detection and keeps its baseline.
//!
//! Set `RAPIDTREES_CPU` to a level's name (`baseline`, `x86-64-v2`,
//! `x86-64-v3`, `x86-64-v4`) to cap the level, for instance to measure what
//! dispatch is worth on one machine. It can only lower the level, never raise
//! it past what the CPU has; an unrecognised value is ignored.

use std::sync::OnceLock;

/// An instruction-set level the row loops are compiled for, in increasing
/// order of capability.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub(crate) enum Level {
    /// Whatever the target was compiled for; the only level off x86-64.
    Baseline,
    /// x86-64-v2: POPCNT and SSE4.2, every x86-64 CPU since about 2009.
    V2,
    /// x86-64-v3: AVX2, BMI1/2, FMA, LZCNT and MOVBE; Haswell (2013) and Zen on.
    V3,
    /// x86-64-v4 (AVX-512 F/BW/CD/DQ/VL) plus AVX512-VPOPCNTDQ, which RF's
    /// popcount wants; Ice Lake (2019) and Zen 4 on.
    V4,
}

impl Level {
    const ALL: [Level; 4] = [Level::Baseline, Level::V2, Level::V3, Level::V4];

    /// The name `RAPIDTREES_CPU` accepts and [`crate::cpu_level`] reports.
    pub(crate) fn name(self) -> &'static str {
        match self {
            Level::Baseline => "baseline",
            Level::V2 => "x86-64-v2",
            Level::V3 => "x86-64-v3",
            Level::V4 => "x86-64-v4",
        }
    }

    /// The level every row loop runs at: the best this CPU supports, capped by
    /// `RAPIDTREES_CPU` when that names a level. Decided once per process.
    pub(crate) fn current() -> Level {
        static CURRENT: OnceLock<Level> = OnceLock::new();
        *CURRENT.get_or_init(|| {
            let cap = std::env::var("RAPIDTREES_CPU").ok().and_then(|name| {
                Level::ALL
                    .into_iter()
                    .find(|level| level.name() == name.trim())
            });
            let best = Level::detected();
            cap.map_or(best, |cap| cap.min(best))
        })
    }

    /// The best level this CPU supports, detected once per process.
    pub(crate) fn detected() -> Level {
        static DETECTED: OnceLock<Level> = OnceLock::new();
        *DETECTED.get_or_init(detect)
    }

    /// Every level this CPU can run, lowest first.
    #[cfg(test)]
    pub(crate) fn supported() -> impl Iterator<Item = Level> {
        let best = Level::detected();
        Level::ALL.into_iter().filter(move |&level| level <= best)
    }
}

#[cfg(target_arch = "x86_64")]
fn detect() -> Level {
    let v2 = is_x86_feature_detected!("popcnt")
        && is_x86_feature_detected!("sse3")
        && is_x86_feature_detected!("ssse3")
        && is_x86_feature_detected!("sse4.1")
        && is_x86_feature_detected!("sse4.2");
    let v3 = v2
        && is_x86_feature_detected!("avx")
        && is_x86_feature_detected!("avx2")
        && is_x86_feature_detected!("bmi1")
        && is_x86_feature_detected!("bmi2")
        && is_x86_feature_detected!("fma")
        && is_x86_feature_detected!("lzcnt")
        && is_x86_feature_detected!("movbe");
    let v4 = v3
        && is_x86_feature_detected!("avx512f")
        && is_x86_feature_detected!("avx512bw")
        && is_x86_feature_detected!("avx512cd")
        && is_x86_feature_detected!("avx512dq")
        && is_x86_feature_detected!("avx512vl")
        && is_x86_feature_detected!("avx512vpopcntdq");
    match (v2, v3, v4) {
        (_, _, true) => Level::V4,
        (_, true, _) => Level::V3,
        (true, _, _) => Level::V2,
        _ => Level::Baseline,
    }
}

#[cfg(not(target_arch = "x86_64"))]
fn detect() -> Level {
    Level::Baseline
}

/// One row of a pairwise matrix: the unit of work handed to a core.
///
/// Implementations must mark [`RowKernel::fill_row`] `#[inline(always)]`.
/// That is what puts the whole loop inside each level's `#[target_feature]`
/// copy in [`fill_row`]: a call that is not inlined runs at the baseline
/// whatever its caller was compiled for, and a closure cannot be forced inline
/// on stable Rust, which is why this is a trait rather than a closure.
pub(crate) trait RowKernel: Sync {
    /// One matrix cell.
    type Cell: Copy + Default + Send;

    /// Write `row[j]` for every `j > i`.
    fn fill_row(&self, i: usize, row: &mut [Self::Cell]);
}

/// Run `kernel.fill_row(i, row)` compiled for `level`, or for the best level
/// the CPU has if `level` asks for more.
#[inline]
pub(crate) fn fill_row<K: RowKernel>(level: Level, kernel: &K, i: usize, row: &mut [K::Cell]) {
    #[cfg(target_arch = "x86_64")]
    {
        // SAFETY: each copy is only entered when `Level::detected` found every
        // feature it enables on this CPU, since `level` is capped by it here.
        match level.min(Level::detected()) {
            Level::V4 => return unsafe { x86::fill_row_v4(kernel, i, row) },
            Level::V3 => return unsafe { x86::fill_row_v3(kernel, i, row) },
            Level::V2 => return unsafe { x86::fill_row_v2(kernel, i, row) },
            Level::Baseline => {}
        }
    }
    #[cfg(not(target_arch = "x86_64"))]
    let _ = level;
    kernel.fill_row(i, row)
}

/// The same row loop compiled three more times, once per level above the
/// baseline. Each body is the inlined [`RowKernel::fill_row`].
#[cfg(target_arch = "x86_64")]
mod x86 {
    use super::RowKernel;

    /// # Safety
    /// The CPU must support every feature enabled here (x86-64-v2).
    #[target_feature(enable = "popcnt,sse3,ssse3,sse4.1,sse4.2")]
    pub(super) unsafe fn fill_row_v2<K: RowKernel>(kernel: &K, i: usize, row: &mut [K::Cell]) {
        kernel.fill_row(i, row)
    }

    /// # Safety
    /// The CPU must support every feature enabled here (x86-64-v3).
    #[target_feature(enable = "popcnt,sse3,ssse3,sse4.1,sse4.2,avx,avx2,bmi1,bmi2,fma,lzcnt,movbe")]
    pub(super) unsafe fn fill_row_v3<K: RowKernel>(kernel: &K, i: usize, row: &mut [K::Cell]) {
        kernel.fill_row(i, row)
    }

    /// # Safety
    /// The CPU must support every feature enabled here (x86-64-v4 plus
    /// AVX512-VPOPCNTDQ).
    #[target_feature(
        enable = "popcnt,sse3,ssse3,sse4.1,sse4.2,avx,avx2,bmi1,bmi2,fma,lzcnt,movbe,avx512f,avx512bw,avx512cd,avx512dq,avx512vl,avx512vpopcntdq"
    )]
    pub(super) unsafe fn fill_row_v4<K: RowKernel>(kernel: &K, i: usize, row: &mut [K::Cell]) {
        kernel.fill_row(i, row)
    }
}

#[cfg(test)]
mod tests {
    use super::Level;

    #[test]
    fn levels_are_ordered_and_named() {
        assert!(Level::Baseline < Level::V2 && Level::V2 < Level::V3 && Level::V3 < Level::V4);
        let names: Vec<&str> = Level::ALL.iter().map(|level| level.name()).collect();
        assert_eq!(names, ["baseline", "x86-64-v2", "x86-64-v3", "x86-64-v4"]);
    }

    #[test]
    fn current_never_exceeds_detected() {
        assert!(Level::current() <= Level::detected());
        assert_eq!(Level::supported().last(), Some(Level::detected()));
    }

    #[cfg(not(target_arch = "x86_64"))]
    #[test]
    fn only_the_baseline_off_x86_64() {
        assert_eq!(Level::detected(), Level::Baseline);
    }
}
