//! Bipartition deduplication: many trees' splits collapse into one table of
//! unique splits, and each tree becomes a list of `u32` IDs into it.
//!
//! This is what lets pairwise distance work on integers: a comparison is a
//! single integer compare, and a large run's working set fits in L2 rather
//! than DRAM.
//!
//! The fingerprint space is cut into [`SHARDS`] shards by the key's low bits,
//! each with its own table, so a chunk of trees is interned by every shard at
//! once instead of by one thread. Each shard walks the chunk in tree order, so
//! its local IDs are first-seen order within the shard whatever the thread
//! count or chunk size, and [`Interner::finish`] numbers them densely, shard by
//! shard. The IDs therefore depend on the input alone.
//!
//! Pendant edges never reach a table. Every tree has all `n` of them, so the
//! pendant of leaf `k` is split `k` outright, and a tree's IDs open with its
//! `n` pendants in leaf order, followed by its internal splits in part order.
//! A path that keeps neither branch lengths nor bipartitions (plain RF) does
//! not store them at all: they are held by every tree and cancel out of RF.

use super::Retain;
use super::Snapshots;
use super::build::Snapshot;
use super::clades::CladeTable;
use super::fingerprint::Fingerprint;
use crate::par::*;
use hashbrown::HashTable;

/// One tree's bipartitions in interned form (private implementation detail of [`Snapshots`]).
///
/// `split_ids` holds the tree's pendant edges first, in leaf order (when they
/// are kept), then its internal splits in part order: the kernels and the
/// export index by ID, so nothing downstream needs it sorted. `lengths[i]` is
/// the branch length of `split_ids[i]` (parallel arrays).
#[derive(Debug, Clone)]
pub(crate) struct InternSnap {
    pub(crate) split_ids: Vec<u32>,
    pub(crate) lengths: Vec<f64>,
}

const SHARD_BITS: u32 = 8;
const SHARDS: usize = 1 << SHARD_BITS;

/// Local IDs a shard can hand out: while a chunk is interned, an ID is
/// `local << SHARD_BITS | shard`, which has to fit a `u32`. All shards together
/// hold as many splits as a `u32` split ID can name.
const LOCAL_IDS: u32 = 1 << (u32::BITS - SHARD_BITS);

/// The shard of `key`. A canonical key is the smaller of `h` and `h ^ total`,
/// which biases its high bits but leaves the low ones uniform.
#[inline]
fn shard_of(key: Fingerprint) -> usize {
    key as usize & (SHARDS - 1)
}

/// The table hash of `key`: the fingerprint is already uniformly random, so
/// the bits above the shard bits serve as they are.
#[inline]
fn hash_of(key: Fingerprint) -> u64 {
    (key >> SHARD_BITS) as u64
}

/// What interning a part needs to know about the run.
#[derive(Clone, Copy)]
struct Run {
    num_leaves: u32,
    rooted: bool,
    /// Whether new splits record their leaf set, see [`Retain::bipartitions`].
    bipartitions: bool,
}

/// One shard of the split table.
///
/// The fingerprint and smaller-side size of each local ID live in separate
/// columns: a `(u128, u32)` tuple is 32 bytes, 12 of them padding.
#[derive(Default)]
struct Shard {
    table: HashTable<u32>,
    fps: Vec<Fingerprint>,
    sizes: Vec<u32>,
    /// How many trees hold each local ID.
    counts: Vec<u32>,
    clades: CladeTable,
}

impl Shard {
    /// Local ID of part `k` of `snap`, registering the split if it is new.
    ///
    /// A split is matched on its 128-bit fingerprint *and* the cardinality of
    /// its smaller side, which is equal for both sides of a bipartition and
    /// rules out a slice of the collision space for one integer compare. Only
    /// a split the shard has never seen records a leaf set, so that
    /// `O(subtree)` walk runs once per unique split rather than once per
    /// occurrence, and not at all when [`Retain::bipartitions`] is off.
    fn intern(&mut self, snap: &Snapshot, k: usize, run: Run) -> Result<u32, String> {
        let part = &snap.parts[k];
        let (key, size) = (part.key, part.size.min(run.num_leaves - part.size));
        let hash = hash_of(key);
        let Self {
            table,
            fps,
            sizes,
            counts,
            clades,
        } = self;
        if let Some(&id) = table.find(hash, |&id| {
            fps[id as usize] == key && sizes[id as usize] == size
        }) {
            counts[id as usize] += 1;
            return Ok(id);
        }

        let id = u32::try_from(fps.len())
            .ok()
            .filter(|&id| id < LOCAL_IDS)
            .ok_or_else(|| "Too many distinct splits for 32-bit split IDs.".to_string())?;
        // Grow by half rather than double: the shards fill evenly, so doubling
        // would leave every one of them oversized by the same large margin.
        // The floor stays small because it is paid [`SHARDS`] times over: a
        // floor of 1024 cost 6 MiB per collection however few splits it held.
        if fps.len() == fps.capacity() {
            let extra = (fps.capacity() / 2).max(16);
            fps.reserve_exact(extra);
            sizes.reserve_exact(extra);
            counts.reserve_exact(extra);
        }
        fps.push(key);
        sizes.push(size);
        counts.push(1);
        if run.bipartitions {
            snap.push_canonical(part, run.rooted, clades);
        }
        table.insert_unique(hash, id, |&id| hash_of(fps[id as usize]));
        Ok(id)
    }
}

/// One tree's internal (non-pendant) parts, grouped by shard.
///
/// `internal[j]` is the part index of the tree's `j`-th internal part, so `j`
/// is that part's slot among the tree's internal IDs. Shard `s` owns the slots
/// `order[starts[s]..starts[s + 1]]`, in part order.
struct Route {
    internal: Vec<u32>,
    order: Vec<u32>,
    starts: Vec<u32>,
}

impl Route {
    /// A counting sort of `snap`'s internal parts on their shard.
    fn of(snap: &Snapshot) -> Self {
        let internal: Vec<u32> = (0..snap.parts.len() as u32)
            .filter(|&k| snap.parts[k as usize].size != 1)
            .collect();
        let shard = |k: u32| shard_of(snap.parts[k as usize].key);
        let mut starts = vec![0u32; SHARDS + 1];
        for &k in &internal {
            starts[shard(k) + 1] += 1;
        }
        for s in 0..SHARDS {
            starts[s + 1] += starts[s];
        }
        let mut next = starts.clone();
        let mut order = vec![0u32; internal.len()];
        for (slot, &k) in internal.iter().enumerate() {
            let at = &mut next[shard(k)];
            order[*at as usize] = slot as u32;
            *at += 1;
        }
        Self {
            internal,
            order,
            starts,
        }
    }

    /// The internal slots whose splits belong to shard `s`.
    #[inline]
    fn shard(&self, s: usize) -> &[u32] {
        &self.order[self.starts[s] as usize..self.starts[s + 1] as usize]
    }
}

/// Incremental, sharded bipartition interner.
///
/// Snapshots are folded in a chunk at a time via [`Interner::push_chunk`],
/// which consumes them, so their parts are freed as soon as they are interned.
/// This keeps construction peak memory near the deduplicated footprint instead
/// of holding every tree's raw parts at once.
pub(super) struct Interner {
    shards: Vec<Shard>,
    snapshots: Vec<InternSnap>,
    words: usize,
    num_leaves: usize,
    rooted: bool,
    retain: Retain,
    /// Whether each tree stores its pendant edges, as IDs `0..n`. Only branch
    /// lengths and the bipartition table need them.
    keep_pendants: bool,
}

impl Interner {
    pub(super) fn new(
        n_trees: usize,
        words: usize,
        num_leaves: usize,
        rooted: bool,
        retain: Retain,
    ) -> Self {
        Self {
            shards: (0..SHARDS).map(|_| Shard::default()).collect(),
            snapshots: Vec::with_capacity(n_trees),
            words,
            num_leaves,
            rooted,
            retain,
            keep_pendants: retain.lengths || retain.bipartitions,
        }
    }

    /// Where a tree's internal IDs start: after its pendants, if it keeps them.
    fn base(&self) -> usize {
        if self.keep_pendants {
            self.num_leaves
        } else {
            0
        }
    }

    /// Intern a chunk of raw snapshots, given in tree order.
    ///
    /// Three passes, each parallel: every tree groups its internal parts by
    /// shard; every shard interns its parts of every tree, in tree order, into
    /// local IDs; every tree collects its IDs back from the shards, after its
    /// pendants. When [`Retain::lengths`] is off, branch lengths are dropped
    /// (RF-only paths never read them).
    ///
    /// # Errors
    /// Returns `Err` when the splits no longer fit 32-bit split IDs.
    pub(super) fn push_chunk(&mut self, raw: Vec<Snapshot>) -> Result<(), String> {
        let (retain, base, keep_pendants) = (self.retain, self.base(), self.keep_pendants);
        let run = Run {
            num_leaves: self.num_leaves as u32,
            rooted: self.rooted,
            bipartitions: retain.bipartitions,
        };
        let routes: Vec<Route> = raw.par_iter().map(Route::of).collect();

        // Per shard, its local IDs for the chunk in route order, and where
        // each tree's run of them starts.
        let locals: Vec<(Vec<u32>, Vec<usize>)> = self
            .shards
            .par_iter_mut()
            .enumerate()
            .map(|(s, shard)| {
                let mut ids = Vec::new();
                let mut starts = Vec::with_capacity(raw.len() + 1);
                for (snap, route) in raw.iter().zip(&routes) {
                    starts.push(ids.len());
                    for &slot in route.shard(s) {
                        let k = route.internal[slot as usize] as usize;
                        ids.push(shard.intern(snap, k, run)?);
                    }
                }
                starts.push(ids.len());
                Ok((ids, starts))
            })
            .collect::<Result<_, String>>()?;

        let interned: Vec<InternSnap> = raw
            .par_iter()
            .zip(&routes)
            .enumerate()
            .map(|(t, (snap, route))| {
                // The pendant of leaf `k` is split `k`, so a tree's pendant IDs
                // are `0..n` whatever its shape.
                let mut split_ids: Vec<u32> = (0..base as u32).collect();
                split_ids.resize(base + route.internal.len(), 0);
                for (s, (ids, starts)) in locals.iter().enumerate() {
                    let run = &ids[starts[t]..starts[t + 1]];
                    for (&slot, &local) in route.shard(s).iter().zip(run) {
                        split_ids[base + slot as usize] = (local << SHARD_BITS) | s as u32;
                    }
                }
                let mut lengths = Vec::new();
                if retain.lengths {
                    lengths.resize(split_ids.len(), 0.0);
                    for part in snap.parts.iter().filter(|part| part.size == 1) {
                        lengths[snap.leaf_order[part.first as usize] as usize] = part.length;
                    }
                    for (slot, &k) in route.internal.iter().enumerate() {
                        lengths[base + slot] = snap.parts[k as usize].length;
                    }
                }
                debug_assert!(keep_pendants || !retain.lengths);
                InternSnap { split_ids, lengths }
            })
            .collect();
        self.snapshots.extend(interned);
        Ok(())
    }

    /// Number every shard's local IDs densely, shard by shard after the
    /// pendants, and assemble the collection.
    pub(super) fn finish(self, leaf_names: Vec<String>) -> Snapshots {
        let base = self.base();
        let mut offsets = [0u32; SHARDS];
        let mut total = base as u32;
        for (offset, shard) in offsets.iter_mut().zip(&self.shards) {
            *offset = total;
            total += shard.fps.len() as u32;
        }

        let mut snapshots = self.snapshots;
        snapshots.par_iter_mut().for_each(|snap| {
            for id in &mut snap.split_ids[base..] {
                *id = offsets[*id as usize & (SHARDS - 1)] + (*id >> SHARD_BITS);
            }
        });

        // Every tree holds every pendant once, and pendant `k` names leaf `k`.
        let mut counts = Vec::with_capacity(total as usize);
        counts.resize(base, snapshots.len() as u32);
        // Sized up front: grown by doubling, the merged table would end up to
        // twice its size while the shards' tables are still alive.
        let pendants = if self.retain.bipartitions {
            self.num_leaves
        } else {
            0
        };
        let mut clades = CladeTable::with_capacity(
            total as usize,
            pendants
                + self
                    .shards
                    .iter()
                    .map(|shard| shard.clades.n_leaves())
                    .sum::<usize>(),
        );
        for leaf in 0..pendants as u32 {
            clades.push(std::iter::once(leaf));
        }
        for shard in self.shards {
            counts.extend_from_slice(&shard.counts);
            clades.append(shard.clades);
        }

        Snapshots {
            snapshots,
            clades,
            split_counts: counts,
            words_per_bitset: self.words,
            leaf_names,
        }
    }
}
