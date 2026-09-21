//! Per-read, per-side splice-junction reconstruction and cross-matching.
//!
//! See the parent module's doc comment for why this re-derives junctions
//! from each side's `cs` tag instead of trusting the comparison table's own
//! `junctions`/`genomic_junctions` columns.

use std::collections::HashSet;
use std::io::Write;

use anyhow::Result;

use crate::cs_parser::parse_cs;

/// One reconstructed splice junction on one side of one read, paired and
/// cross-matched against the other side.
pub(super) struct JunctionRecord {
    pub junction_index: u32,
    pub query_pos: i64,
    pub genomic: (String, u64, u64),
    pub strand: char,
    pub matched_in_query: bool,
    pub matched_in_genomic: bool,
    pub other_side_aligned: bool,
}

/// Reconstruct one side's final, query-ascending, paired `(query_pos,
/// genomic)` junction list from its `cs` tag — applying the same
/// offset/strand conversion as `record.rs`'s `AlnInfo::from_paf`, but keeping
/// the query and genomic vectors correctly paired throughout via a single
/// stable-sort permutation applied to both together (see module doc comment
/// in `mod.rs`).
pub(super) fn build_side_junctions(
    cs: &str,
    strand: char,
    query_start: u64,
    query_end: u64,
    target_start: u64,
    chrom: &str,
    read_len: u64,
) -> Vec<(i64, (String, u64, u64))> {
    if cs.is_empty() {
        return Vec::new();
    }
    let stats = parse_cs(cs);
    let n = stats.raw_junctions.len();
    if n == 0 {
        return Vec::new();
    }

    // Same offset formula as record.rs:AlnInfo::from_paf.
    let offset: u64 = if strand == '+' { query_start } else { read_len.saturating_sub(query_end) };

    let mut query_vals: Vec<i64> = stats.raw_junctions.iter().map(|&q| (offset + q) as i64).collect();
    if strand == '-' {
        let ql = read_len as i64;
        query_vals = query_vals.into_iter().map(|j| ql - j).collect();
    }

    let genomic_vals: Vec<(String, u64, u64)> = stats
        .raw_genomic_junctions
        .iter()
        .map(|&(s, e)| (chrom.to_string(), target_start + s, target_start + e))
        .collect();

    // Sort BOTH vectors by one shared, stable permutation (derived from the
    // query values) rather than sorting `query_vals` alone — this is what
    // keeps query_pos[i]/genomic[i] paired to the same intron after sorting,
    // unlike record.rs's independent `junctions.sort_unstable()`.
    let mut perm: Vec<usize> = (0..n).collect();
    perm.sort_by_key(|&i| query_vals[i]);

    perm.into_iter().map(|i| (query_vals[i], genomic_vals[i].clone())).collect()
}

/// Build this side's membership sets for O(1) per-element lookup by the other
/// side. `junction.rs`'s `junction_set_stats`/`genomic_junction_set_stats`
/// only return aggregate counts, not per-element membership, so this is a
/// small dedicated helper rather than a reuse of those functions.
pub(super) fn side_position_sets(
    pairs: &[(i64, (String, u64, u64))],
) -> (HashSet<i64>, HashSet<(String, u64, u64)>) {
    let qset: HashSet<i64> = pairs.iter().map(|(q, _)| *q).collect();
    let gset: HashSet<(String, u64, u64)> = pairs.iter().map(|(_, g)| g.clone()).collect();
    (qset, gset)
}

pub(super) fn build_junction_records(
    pairs: Vec<(i64, (String, u64, u64))>,
    strand: char,
    other_query_set: &HashSet<i64>,
    other_genomic_set: &HashSet<(String, u64, u64)>,
    other_side_aligned: bool,
) -> Vec<JunctionRecord> {
    pairs
        .into_iter()
        .enumerate()
        .map(|(i, (q, g))| JunctionRecord {
            junction_index: (i + 1) as u32,
            query_pos: q,
            matched_in_query: other_side_aligned && other_query_set.contains(&q),
            matched_in_genomic: other_side_aligned && other_genomic_set.contains(&g),
            genomic: g,
            strand,
            other_side_aligned,
        })
        .collect()
}

pub(super) fn write_per_read_row(w: &mut dyn Write, read_name: &str, side: char, rec: &JunctionRecord) -> Result<()> {
    writeln!(
        w,
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
        read_name,
        side,
        rec.junction_index,
        rec.strand,
        rec.query_pos,
        rec.genomic.0,
        rec.genomic.1,
        rec.genomic.2,
        rec.matched_in_query as u8,
        rec.matched_in_genomic as u8,
        rec.other_side_aligned as u8,
    )?;
    Ok(())
}
