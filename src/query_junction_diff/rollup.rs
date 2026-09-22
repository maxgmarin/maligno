//! Per-side rollup of unmatched genomic junctions: how many distinct reads
//! report each exact `(chrom, start, end, strand)` locus unsupported by the
//! other side, plus (pass 2, see `mod.rs`) how many reads *total*, across the
//! whole comparison table, carry that same junction on this side.

use std::collections::{HashMap, HashSet};
use std::io::Write;

use anyhow::Result;

use crate::io_utils::open_output;
use crate::junction::parse_genomic_junction_str;

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

/// Pass 2: for one comparison-table row's one side, parse its
/// genomic-junctions column, reattach `chrom` (dropped at serialization —
/// see `parse_genomic_junction_str`'s doc comment), pair each junction with
/// this row's own (whole-read) `strand`, and increment the running total for
/// any junction that's a key in `totals` (i.e. one of the junctions already
/// flagged by pass 1's `UnmatchedAcc`). Junctions not already in `totals` are
/// ignored — this only ever counts occurrences of already-flagged junctions,
/// never discovers new ones. Needs no `cs`-tag reconstruction: pure genomic
/// set membership, unaffected by the `-`-strand query↔genomic pairing issue
/// `reconstruct.rs` exists to solve, since this never needs a query position.
pub(super) fn accumulate_side_total(chrom: &str, strand: char, genomic_junctions_str: &str, totals: &mut TotalsAcc) {
    for (start, end) in parse_genomic_junction_str(genomic_junctions_str) {
        let key = (chrom.to_string(), start, end, strand);
        if let Some(count) = totals.get_mut(&key) {
            *count += 1;
        }
    }
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

        accumulate_side_total("chr1", '+', "((100, 250),)", &mut totals);
        assert_eq!(totals[&("chr1".to_string(), 100, 250, '+')], 1);

        accumulate_side_total("chr1", '+', "((100, 250),)", &mut totals);
        assert_eq!(totals[&("chr1".to_string(), 100, 250, '+')], 2);
    }

    #[test]
    fn accumulate_side_total_ignores_junctions_not_of_interest() {
        let mut totals: TotalsAcc = HashMap::new();
        totals.insert(("chr1".to_string(), 100, 250, '+'), 0);

        accumulate_side_total("chr1", '+', "((400, 800),)", &mut totals);
        assert_eq!(totals[&("chr1".to_string(), 100, 250, '+')], 0);
    }

    #[test]
    fn accumulate_side_total_requires_matching_strand() {
        let mut totals: TotalsAcc = HashMap::new();
        totals.insert(("chr1".to_string(), 100, 250, '+'), 0);

        accumulate_side_total("chr1", '-', "((100, 250),)", &mut totals);
        assert_eq!(totals[&("chr1".to_string(), 100, 250, '+')], 0);
    }

    #[test]
    fn accumulate_side_total_requires_matching_chrom() {
        let mut totals: TotalsAcc = HashMap::new();
        totals.insert(("chr1".to_string(), 100, 250, '+'), 0);

        accumulate_side_total("chr2", '+', "((100, 250),)", &mut totals);
        assert_eq!(totals[&("chr1".to_string(), 100, 250, '+')], 0);
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
