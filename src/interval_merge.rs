//! Generic sort + single-sweep interval merge with per-locus accumulation.
//!
//! Equivalent to `bedtools merge -c -o count`, but the caller folds arbitrary
//! per-interval payload (`meta`) into each merged locus's accumulator (`acc`).
//! Written generically so any command can accumulate its own stats.
//!
//! Half-open coordinates: `[start, end)`. Intervals that *touch*
//! (`next.start == cur.end`) as well as those that *overlap* (`next.start < cur.end`)
//! are merged; a gap (`next.start > cur.end`) starts a new locus.

use std::cmp::Ordering;

/// One interval to merge. `meta` carries whatever per-interval payload the caller
/// wants folded into the locus (strand, outcome, score, …).
pub(crate) struct Ivl<M> {
    pub chrom: String,
    pub start: u64, // 0-based, half-open [start, end)
    pub end: u64,
    pub meta: M,
}

/// A merged locus: the union interval plus the caller's accumulated stats `A`.
pub(crate) struct Locus<A> {
    pub chrom: String,
    pub start: u64,
    pub end: u64,
    pub acc: A,
}

/// Sort by `(chrom, start)`; sweep once, merging overlapping-or-touching intervals
/// and folding each absorbed interval's `meta` into the locus accumulator. Returns
/// loci in `(chrom, start)` order. O(k log k) time, O(k) memory.
///
/// The count is "free": the merge condition *is* the overlap condition, so any
/// interval absorbed into a locus overlaps it by construction — `fold` runs at the
/// moment of absorption, so per-locus stats are complete when the locus is emitted.
pub(crate) fn merge_and_count<M, A>(
    mut ivls: Vec<Ivl<M>>,
    new_acc: impl Fn() -> A,
    fold: impl Fn(&mut A, &M),
) -> Vec<Locus<A>> {
    if ivls.is_empty() {
        return Vec::new();
    }
    ivls.sort_by(|x, y| match x.chrom.cmp(&y.chrom) {
        Ordering::Equal => x.start.cmp(&y.start),
        other => other,
    });

    let mut out: Vec<Locus<A>> = Vec::new();
    let mut iter = ivls.into_iter();
    let first = iter.next().unwrap();
    let mut cur = Locus {
        chrom: first.chrom,
        start: first.start,
        end: first.end,
        acc: new_acc(),
    };
    fold(&mut cur.acc, &first.meta);

    for iv in iter {
        if iv.chrom == cur.chrom && iv.start <= cur.end {
            if iv.end > cur.end {
                cur.end = iv.end;
            }
            fold(&mut cur.acc, &iv.meta);
        } else {
            out.push(cur);
            cur = Locus {
                chrom: iv.chrom,
                start: iv.start,
                end: iv.end,
                acc: new_acc(),
            };
            fold(&mut cur.acc, &iv.meta);
        }
    }
    out.push(cur);
    out
}

// ───────────────────────────────────────────────────────────────────────────
// Tests
// ───────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    fn iv(chrom: &str, start: u64, end: u64) -> Ivl<()> {
        Ivl { chrom: chrom.to_string(), start, end, meta: () }
    }

    /// Merge with a plain count accumulator → (chrom, start, end, n).
    fn merged(ivls: Vec<Ivl<()>>) -> Vec<(String, u64, u64, u64)> {
        merge_and_count(ivls, || 0u64, |a, _| *a += 1)
            .into_iter()
            .map(|l| (l.chrom, l.start, l.end, l.acc))
            .collect()
    }

    #[test]
    fn empty() {
        assert!(merged(vec![]).is_empty());
    }

    #[test]
    fn single() {
        assert_eq!(merged(vec![iv("chr1", 10, 20)]), vec![("chr1".into(), 10, 20, 1)]);
    }

    #[test]
    fn overlap_merges() {
        assert_eq!(
            merged(vec![iv("chr1", 10, 20), iv("chr1", 15, 25)]),
            vec![("chr1".into(), 10, 25, 2)]
        );
    }

    #[test]
    fn exactly_touching_merges() {
        // half-open: [10,20) and [20,30) touch at 20 → merge
        assert_eq!(
            merged(vec![iv("chr1", 10, 20), iv("chr1", 20, 30)]),
            vec![("chr1".into(), 10, 30, 2)]
        );
    }

    #[test]
    fn one_bp_gap_does_not_merge() {
        // [10,20) then [21,30): 21 > 20 → gap
        assert_eq!(
            merged(vec![iv("chr1", 10, 20), iv("chr1", 21, 30)]),
            vec![("chr1".into(), 10, 20, 1), ("chr1".into(), 21, 30, 1)]
        );
    }

    #[test]
    fn contained_interval_merges_and_keeps_max_end() {
        assert_eq!(
            merged(vec![iv("chr1", 10, 100), iv("chr1", 20, 30)]),
            vec![("chr1".into(), 10, 100, 2)]
        );
    }

    #[test]
    fn multi_chrom_boundary() {
        // input out of chrom order → sorted, no cross-chrom merge
        assert_eq!(
            merged(vec![iv("chr2", 10, 20), iv("chr1", 10, 20)]),
            vec![("chr1".into(), 10, 20, 1), ("chr2".into(), 10, 20, 1)]
        );
    }

    #[test]
    fn unsorted_input_is_sorted_then_merged() {
        assert_eq!(
            merged(vec![iv("chr1", 400, 450), iv("chr1", 100, 150), iv("chr1", 120, 180)]),
            vec![("chr1".into(), 100, 180, 2), ("chr1".into(), 400, 450, 1)]
        );
    }

    #[test]
    fn fold_accumulates_strand_breakdown() {
        // Worked trace from the spec §6.4 (strand meta).
        let ivls = vec![
            Ivl { chrom: "chr7".into(), start: 100, end: 150, meta: '+' },
            Ivl { chrom: "chr7".into(), start: 120, end: 180, meta: '+' },
            Ivl { chrom: "chr7".into(), start: 175, end: 200, meta: '-' },
            Ivl { chrom: "chr7".into(), start: 400, end: 450, meta: '+' },
        ];
        let loci = merge_and_count(
            ivls,
            || (0u64, 0u64), // (n_plus, n_minus)
            |a, &s| match s {
                '+' => a.0 += 1,
                '-' => a.1 += 1,
                _ => {}
            },
        );
        assert_eq!(loci.len(), 2);
        assert_eq!((loci[0].chrom.as_str(), loci[0].start, loci[0].end, loci[0].acc), ("chr7", 100, 200, (2, 1)));
        assert_eq!((loci[1].chrom.as_str(), loci[1].start, loci[1].end, loci[1].acc), ("chr7", 400, 450, (1, 0)));
    }
}
