//! Per-side rollup of unmatched genomic junctions: how many distinct reads
//! report each exact `(chrom, start, end, strand)` locus unsupported by the
//! other side.

use std::collections::{HashMap, HashSet};
use std::io::Write;

use anyhow::Result;

use crate::io_utils::open_output;

use super::reconstruct::JunctionRecord;

/// An exact genomic locus + strand.
type UnmatchedKey = (String, u64, u64, char);
pub(super) type UnmatchedAcc = HashMap<UnmatchedKey, u64>;

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

pub(super) fn write_unmatched_table(path: &str, acc: &UnmatchedAcc) -> Result<()> {
    let mut rows: Vec<&UnmatchedKey> = acc.keys().collect();
    rows.sort();
    let mut w = open_output(Some(path))?;
    writeln!(w, "chrom\tgenomic_start\tgenomic_end\tstrand\tn_reads_query_junctions_different")?;
    for key @ (chrom, start, end, strand) in rows {
        let n = acc[key];
        writeln!(w, "{chrom}\t{start}\t{end}\t{strand}\t{n}")?;
    }
    w.flush()?;
    Ok(())
}
