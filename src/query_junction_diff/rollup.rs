//! Per-side rollup of unmatched genomic junctions: how many distinct reads
//! report each exact `(chrom, start, end, strand)` locus unsupported by the
//! other side, plus (pass 2, see `mod.rs`) how many reads *total*, across the
//! whole comparison table, carry that same junction on this side.
//!
//! Pass 2 also collects every distinct junction seen on each side
//! (`UniqueJunctions`), for the summary's unique / shared / only-A / only-B
//! junction counts.

use std::collections::{HashMap, HashSet};
use std::io::Write;

use anyhow::Result;

use crate::io_utils::open_output;

use super::reconstruct::JunctionRecord;

/// An exact genomic locus + strand.
type UnmatchedKey = (String, u64, u64, char);
pub(super) type UnmatchedAcc = HashMap<UnmatchedKey, u64>;
/// Same key shape, used for pass 2's total-occurrence tally.
pub(super) type TotalsAcc = HashMap<UnmatchedKey, u64>;

/// Increment the rollup for one read's `matched_in_genomic == false` junctions
/// on this side. Dedupes within this read's own record list first, so a read
/// that happens to report the same unmatched locus more than once still only
/// increments the read-count by 1 for that key (the column is
/// "distinct reads", not a raw occurrence tally).
pub(super) fn accumulate_unmatched(acc: &mut UnmatchedAcc, records: &[JunctionRecord]) {
    let mut seen_this_read: HashSet<UnmatchedKey> = HashSet::new();
    for r in records {
        if r.matched_in_genomic {
            continue;
        }
        let key = (r.genomic.0.clone(), r.genomic.1, r.genomic.2, r.strand);
        if seen_this_read.insert(key.clone()) {
            *acc.entry(key).or_insert(0) += 1;
        }
    }
}

/// Pass 2: for one comparison-table row's one side, take its parsed
/// genomic junctions (from `parse_genomic_junction_str`), reattach `chrom`
/// (dropped at serialization — see that function's doc comment), pair each
/// junction with this row's own (whole-read) `strand`, and increment the
/// running total for any junction that's a key in `totals` (i.e. one of the
/// junctions already flagged by pass 1's `UnmatchedAcc`). Junctions not
/// already in `totals` are ignored — this only ever counts occurrences of
/// already-flagged junctions, never discovers new ones. Needs no `cs`-tag
/// reconstruction: pure genomic set membership, unaffected by the `-`-strand
/// query↔genomic pairing issue `reconstruct.rs` exists to solve, since this
/// never needs a query position.
pub(super) fn accumulate_side_total(chrom: &str, strand: char, junctions: &[(u64, u64)], totals: &mut TotalsAcc) {
    for &(start, end) in junctions {
        let key = (chrom.to_string(), start, end, strand);
        if let Some(count) = totals.get_mut(&key) {
            *count += 1;
        }
    }
}

/// Every distinct `(chrom, start, end, strand)` genomic junction seen on one
/// side, across every row of the comparison table — same key as the
/// unmatched tables (`strand` is the alignment's strand), so the counts line
/// up with those files.
///
/// Two-level (chrom → set of `(start, end, strand)`) rather than one set keyed
/// on the full tuple: the chrom string is looked up once per row instead of
/// allocated once per junction, and each entry stays a fixed-size 24 bytes.
/// Memory scales with the number of *distinct* junctions, not reads.
#[derive(Default)]
pub(super) struct UniqueJunctions(HashMap<String, HashSet<(u64, u64, char)>>);

impl UniqueJunctions {
    /// Add one row's one side: its chrom, alignment strand, and parsed
    /// genomic junctions.
    pub(super) fn observe(&mut self, chrom: &str, strand: char, junctions: &[(u64, u64)]) {
        if junctions.is_empty() {
            return;
        }
        if !self.0.contains_key(chrom) {
            self.0.insert(chrom.to_string(), HashSet::new());
        }
        let set = self.0.get_mut(chrom).expect("inserted above");
        for &(start, end) in junctions {
            set.insert((start, end, strand));
        }
    }

    pub(super) fn len(&self) -> u64 {
        self.0.values().map(|s| s.len() as u64).sum()
    }

    /// How many of this side's distinct junctions also appear in `other`.
    fn n_shared_with(&self, other: &UniqueJunctions) -> u64 {
        let mut n = 0;
        for (chrom, set) in &self.0 {
            if let Some(other_set) = other.0.get(chrom) {
                // Walk the smaller set, probe the larger.
                let (small, large) = if set.len() <= other_set.len() { (set, other_set) } else { (other_set, set) };
                n += small.iter().filter(|k| large.contains(k)).count() as u64;
            }
        }
        n
    }
}

/// Distinct-junction counts for the summary: per side, and split into shared
/// vs only-one-side across the whole table.
pub(super) struct UniqueJunctionCounts {
    pub n_a: u64,
    pub n_b: u64,
    pub n_shared: u64,
    pub n_only_a: u64,
    pub n_only_b: u64,
}

pub(super) fn unique_junction_counts(a: &UniqueJunctions, b: &UniqueJunctions) -> UniqueJunctionCounts {
    let n_a = a.len();
    let n_b = b.len();
    let n_shared = a.n_shared_with(b);
    UniqueJunctionCounts { n_a, n_b, n_shared, n_only_a: n_a - n_shared, n_only_b: n_b - n_shared }
}

/// Write the unmatched-junction table with both counts: `acc`'s
/// distinct-differing-reads count (unchanged from before pass 2 existed),
/// and `totals`'s whole-table occurrence count.
pub(super) fn write_unmatched_table(path: &str, acc: &UnmatchedAcc, totals: &TotalsAcc) -> Result<()> {
    let mut rows: Vec<&UnmatchedKey> = acc.keys().collect();
    rows.sort();
    let mut w = open_output(Some(path))?;
    writeln!(
        w,
        "chrom\tgenomic_start\tgenomic_end\tstrand\tn_reads_query_junctions_different\tn_reads_with_junction_total"
    )?;
    for key @ (chrom, start, end, strand) in rows {
        let n = acc[key];
        let total = totals.get(key).copied().unwrap_or(0);
        writeln!(w, "{chrom}\t{start}\t{end}\t{strand}\t{n}\t{total}")?;
    }
    w.flush()?;
    Ok(())
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn accumulate_side_total_counts_a_junction_of_interest() {
        let mut totals: TotalsAcc = HashMap::new();
        totals.insert(("chr1".to_string(), 100, 250, '+'), 0);

        accumulate_side_total("chr1", '+', &[(100, 250)], &mut totals);
        assert_eq!(totals[&("chr1".to_string(), 100, 250, '+')], 1);

        accumulate_side_total("chr1", '+', &[(100, 250)], &mut totals);
        assert_eq!(totals[&("chr1".to_string(), 100, 250, '+')], 2);
    }

    #[test]
    fn accumulate_side_total_ignores_junctions_not_of_interest() {
        let mut totals: TotalsAcc = HashMap::new();
        totals.insert(("chr1".to_string(), 100, 250, '+'), 0);

        accumulate_side_total("chr1", '+', &[(400, 800)], &mut totals);
        assert_eq!(totals[&("chr1".to_string(), 100, 250, '+')], 0);
    }

    #[test]
    fn accumulate_side_total_requires_matching_strand() {
        let mut totals: TotalsAcc = HashMap::new();
        totals.insert(("chr1".to_string(), 100, 250, '+'), 0);

        accumulate_side_total("chr1", '-', &[(100, 250)], &mut totals);
        assert_eq!(totals[&("chr1".to_string(), 100, 250, '+')], 0);
    }

    #[test]
    fn accumulate_side_total_requires_matching_chrom() {
        let mut totals: TotalsAcc = HashMap::new();
        totals.insert(("chr1".to_string(), 100, 250, '+'), 0);

        accumulate_side_total("chr2", '+', &[(100, 250)], &mut totals);
        assert_eq!(totals[&("chr1".to_string(), 100, 250, '+')], 0);
    }

    #[test]
    fn unique_junctions_dedupes_across_rows() {
        let mut u = UniqueJunctions::default();
        u.observe("chr1", '+', &[(100, 250), (400, 800)]);
        u.observe("chr1", '+', &[(100, 250)]);
        u.observe("chr1", '+', &[]);
        assert_eq!(u.len(), 2);
    }

    #[test]
    fn unique_junctions_keys_on_chrom_and_strand() {
        let mut u = UniqueJunctions::default();
        u.observe("chr1", '+', &[(100, 250)]);
        u.observe("chr1", '-', &[(100, 250)]);
        u.observe("chr2", '+', &[(100, 250)]);
        assert_eq!(u.len(), 3);
    }

    #[test]
    fn unique_junction_counts_splits_shared_and_only() {
        let mut a = UniqueJunctions::default();
        a.observe("chr1", '+', &[(100, 250), (400, 800)]);
        a.observe("chr2", '-', &[(10, 20)]);
        let mut b = UniqueJunctions::default();
        b.observe("chr1", '+', &[(100, 250), (900, 950)]);
        b.observe("chr1", '-', &[(400, 800)]); // same coords, other strand: not shared
        b.observe("chr3", '+', &[(5, 6)]);

        let c = unique_junction_counts(&a, &b);
        assert_eq!((c.n_a, c.n_b), (3, 4));
        assert_eq!(c.n_shared, 1);
        assert_eq!((c.n_only_a, c.n_only_b), (2, 3));
    }

    #[test]
    fn write_unmatched_table_includes_both_columns() {
        let mut acc: UnmatchedAcc = HashMap::new();
        acc.insert(("chr22".to_string(), 100, 250, '+'), 3);
        let mut totals: TotalsAcc = HashMap::new();
        totals.insert(("chr22".to_string(), 100, 250, '+'), 12);

        let path = std::env::temp_dir().join("maligno_test_write_unmatched_table.tsv");
        let path_str = path.to_str().unwrap();
        write_unmatched_table(path_str, &acc, &totals).unwrap();
        let content = std::fs::read_to_string(&path).unwrap();
        let _ = std::fs::remove_file(&path);

        let lines: Vec<&str> = content.lines().collect();
        assert_eq!(
            lines[0],
            "chrom\tgenomic_start\tgenomic_end\tstrand\tn_reads_query_junctions_different\tn_reads_with_junction_total"
        );
        assert_eq!(lines[1], "chr22\t100\t250\t+\t3\t12");
    }

    #[test]
    fn write_unmatched_table_defaults_total_to_zero_when_missing() {
        // Defensive: if `totals` somehow lacks a key `acc` has (inputs out of
        // sync), don't panic — report 0 rather than fabricating a number.
        let mut acc: UnmatchedAcc = HashMap::new();
        acc.insert(("chr22".to_string(), 100, 250, '+'), 1);
        let totals: TotalsAcc = HashMap::new();

        let path = std::env::temp_dir().join("maligno_test_write_unmatched_table_missing_total.tsv");
        let path_str = path.to_str().unwrap();
        write_unmatched_table(path_str, &acc, &totals).unwrap();
        let content = std::fs::read_to_string(&path).unwrap();
        let _ = std::fs::remove_file(&path);

        assert!(content.lines().nth(1).unwrap().ends_with("\t1\t0"));
    }
}
