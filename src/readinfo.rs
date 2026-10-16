use std::collections::HashMap;
use std::io::{BufRead, Write};

use anyhow::{bail, Result};

use crate::columns::*;
use crate::io_utils::{escape_tsv_field, fmt_float, open_input, open_output};

// ── Column indices in alninfo rows ───────────────────────────────────────────
// Resolved by name against `ALNINFO_COLS` at compile time.
const COL_QUERY_NAME: usize = index_of(&ALNINFO_COLS, QUERY_NAME);
const COL_QUERY_LEN: usize = index_of(&ALNINFO_COLS, QUERY_LEN);
const COL_QUERY_START: usize = index_of(&ALNINFO_COLS, QUERY_START);
const COL_QUERY_END: usize = index_of(&ALNINFO_COLS, QUERY_END);
const COL_STRAND: usize = index_of(&ALNINFO_COLS, STRAND);
const COL_QUERY_TRANSCRIPT_STRAND: usize = index_of(&ALNINFO_COLS, QUERY_TRANSCRIPT_STRAND);
const COL_GENOMIC_INTRON_STRAND: usize = index_of(&ALNINFO_COLS, GENOMIC_INTRON_STRAND);
const COL_MQ: usize = index_of(&ALNINFO_COLS, MAPQ);
const COL_TARGET_NAME: usize = index_of(&ALNINFO_COLS, TARGET_NAME);
const COL_TARGET_START: usize = index_of(&ALNINFO_COLS, TARGET_START);
const COL_TARGET_END: usize = index_of(&ALNINFO_COLS, TARGET_END);
const COL_MS: usize = index_of(&ALNINFO_COLS, MS_SCORE);
const COL_AS: usize = index_of(&ALNINFO_COLS, AS_SCORE);
const COL_CS: usize = index_of(&ALNINFO_COLS, CS);
const COL_N_MATCH_EVENTS: usize = index_of(&ALNINFO_COLS, N_MATCH_EVENTS);
const COL_N_MATCH_BASES: usize = index_of(&ALNINFO_COLS, N_MATCH_BASES);
const COL_N_SUB_EVENTS: usize = index_of(&ALNINFO_COLS, N_SUBSTITUTION_EVENTS);
const COL_N_SUB_BASES: usize = index_of(&ALNINFO_COLS, N_SUBSTITUTION_BASES);
const COL_N_INS_EVENTS: usize = index_of(&ALNINFO_COLS, N_INSERTION_EVENTS);
const COL_N_INS_BASES: usize = index_of(&ALNINFO_COLS, N_INSERTION_BASES);
const COL_N_DEL_EVENTS: usize = index_of(&ALNINFO_COLS, N_DELETION_EVENTS);
const COL_N_DEL_BASES: usize = index_of(&ALNINFO_COLS, N_DELETION_BASES);
const COL_N_SPLICE_EVENTS: usize = index_of(&ALNINFO_COLS, N_INTRON_EVENTS);
const COL_N_SPLICE_BASES: usize = index_of(&ALNINFO_COLS, N_INTRON_BASES);
const COL_N_SC_START: usize = index_of(&ALNINFO_COLS, N_SOFTCLIP_BASES_START);
const COL_N_SC_END: usize = index_of(&ALNINFO_COLS, N_SOFTCLIP_BASES_END);
const COL_N_SC_EVENTS: usize = index_of(&ALNINFO_COLS, N_SOFTCLIP_EVENTS);
const COL_JUNCTIONS: usize = index_of(&ALNINFO_COLS, QUERY_JUNCTIONS);
const COL_ALN_IDENTITY: usize = index_of(&ALNINFO_COLS, ALN_IDENTITY);
const COL_QUERY_ALN_LEN: usize = index_of(&ALNINFO_COLS, QUERY_ALN_LEN);
const COL_QUERY_ALN_COV: usize = index_of(&ALNINFO_COLS, QUERY_ALN_COV);
const COL_GENOMIC_JUNCTIONS: usize = index_of(&ALNINFO_COLS, GENOMIC_JUNCTIONS);
const COL_TP_TAG: usize = index_of(&ALNINFO_COLS, TP_TAG);

/// The readinfo header line (no trailing newline), from `READINFO_COLS`.
pub fn readinfo_header() -> String {
    header_line(&READINFO_COLS)
}

// ── ReadInfo row ─────────────────────────────────────────────────────────────

/// Per-read summary row: one row of the readinfo table (`toolkit paf2tables
/// --readinfo`, `compare --emit-readinfo`).
/// Fields are typed for use in the `compare` subcommand; `raw_fields` carries
/// all non-key columns pre-formatted for TSV pass-through.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct ReadInfoRow {
    pub read_name: String,
    pub read_len: u64,
    pub target_chr: String,
    pub strand: char,
    /// Representative alignment's transcript strand, relative to the read.
    pub query_transcript_strand: char,
    /// Representative alignment's intron strand on the reference.
    pub genomic_intron_strand: char,
    pub mapq_best: u32,
    pub aln_score: i64,
    pub ms: i64,
    pub query_aln_cov: f64,
    pub query_aln_len: u64,
    pub seqid: f64,
    pub junctions: String,
    pub num_aln: u64,
    pub num_aln_maxscore: u64,
    pub n_match_events: u64,
    pub n_match_bases: u64,
    pub n_sub_events: u64,
    pub n_sub_bases: u64,
    pub n_ins_events: u64,
    pub n_ins_bases: u64,
    pub n_del_events: u64,
    pub n_del_bases: u64,
    pub n_splice_events: u64,
    pub n_splice_bases: u64,
    pub n_sc_start: u64,
    pub n_sc_end: u64,
    pub n_sc_events: u64,
    pub cs: String,
    /// Best alignment's genomic_junctions string (Python tuple-of-tuples form).
    pub genomic_junctions: String,
    /// Best alignment's query_start / query_end (0-based half-open, BED-style, on the read).
    pub query_start: u64,
    pub query_end: u64,
    /// Best alignment's target_start / target_end (0-based half-open, BED-style, on the reference).
    pub target_start: u64,
    pub target_end: u64,
    /// Best alignment's PAF `tp:A` tag (alignment type: 'P' primary, 'S'
    /// secondary, 'I' inversion of the primary). '*' when absent or unmapped.
    pub tp_tag: char,
    /// Mapped alignment rows whose `tp:A` is `P` or `I` (SAM primary + supplementary).
    pub num_aln_tpp: u64,
    /// Mapped alignment rows whose `tp:A` is `S` or `i` (secondary).
    pub num_aln_tps: u64,
    /// All fields except query_name and query_len, pre-formatted for pass-through,
    /// in `READINFO_COLS` order (pinned by `raw_fields_follow_the_readinfo_columns`).
    pub raw_fields: Vec<String>,
}

impl ReadInfoRow {
    pub fn write<W: Write + ?Sized>(&self, w: &mut W) -> std::io::Result<()> {
        write!(w, "{}\t{}", self.read_name, self.read_len)?;
        for f in &self.raw_fields {
            let escaped = escape_tsv_field(f);
            write!(w, "\t{escaped}")?;
        }
        writeln!(w)
    }
}

// ── Alninfo row (lightweight, strings only) ───────────────────────────────────

pub(crate) struct AlnRow {
    pub(crate) fields: Vec<String>,
    pub(crate) ms: i64,
    pub(crate) aln_score: i64,
    pub(crate) is_aligned: bool, // target_name != "*"
    pub(crate) mapq: u32,        // parsed once; used as tertiary sort key
    /// Final sort key, from `tp:A` (see [`tp_rank`]). Only breaks ties; it is
    /// not a key of `n_alns_tied_top_score`.
    pub(crate) tp_rank: u8,
}

/// Rank of an alignment row as the final selection key, highest first:
/// primary (`tp:A` `P`/`I`) 3, untagged 2, secondary (`S`/`i`) 1, unmapped 0.
/// So among alignments tied on ms, AS and MAPQ, the aligner's own primary
/// represents the read, and a mapped row always beats an unmapped placeholder.
fn tp_rank(is_aligned: bool, tp: Option<char>) -> u8 {
    match (is_aligned, tp) {
        (false, _) => 0,
        (true, Some('P' | 'I')) => 3,
        (true, Some('S' | 'i')) => 1,
        (true, _) => 2,
    }
}

impl AlnRow {
    /// query_name (read ID) of this alignment row — used to detect group
    /// boundaries when reading PAF straight into per-read groups.
    pub(crate) fn query_name(&self) -> &str {
        &self.fields[COL_QUERY_NAME]
    }
}

pub(crate) fn parse_aln_row(line: &str) -> Option<AlnRow> {
    let fields: Vec<String> = line.split('\t').map(str::to_owned).collect();
    if fields.len() <= COL_GENOMIC_JUNCTIONS {
        return None;
    }
    let ms: i64 = fields[COL_MS].parse().unwrap_or(0);
    let aln_score: i64 = fields[COL_AS].parse().unwrap_or(0);
    let is_aligned = fields[COL_TARGET_NAME] != "*";
    let mapq: u32 = fields[COL_MQ].parse().unwrap_or(0);
    let tp_rank = tp_rank(is_aligned, fields.get(COL_TP_TAG).and_then(|s| s.chars().next()));
    Some(AlnRow {
        fields,
        ms,
        aln_score,
        is_aligned,
        mapq,
        tp_rank,
    })
}

/// Collapse a group of alignment rows into a single per-read summary.
///
/// Sort order: ms desc → AS desc → MQ desc → primary before secondary
/// ([`tp_rank`]). The highest-scoring, highest-confidence alignment becomes the
/// representative, and among full ties the aligner's primary; every per-alignment
/// column (target, strand, scores, identity, coverage, junctions, cs,
/// coordinates, event counts) comes from that single row. Only the counts
/// (n_alns, n_alns_primsupp, n_alns_secondary, n_alns_tied_top_score) are computed over the
/// full group.
pub(crate) fn collapse_group(rows: &mut [AlnRow]) -> ReadInfoRow {
    rows.sort_by(|a, b| {
        b.ms
            .cmp(&a.ms)
            .then_with(|| b.aln_score.cmp(&a.aln_score))
            .then_with(|| b.mapq.cmp(&a.mapq))
            .then_with(|| b.tp_rank.cmp(&a.tp_rank))
    });

    let best = &rows[0];
    let bf = &best.fields;

    let read_name = bf[COL_QUERY_NAME].clone();
    let read_len: u64 = bf[COL_QUERY_LEN].parse().unwrap_or(0);
    let target_chr = bf[COL_TARGET_NAME].clone();
    let strand: char = bf[COL_STRAND].chars().next().unwrap_or('*');
    let query_transcript_strand: char = bf[COL_QUERY_TRANSCRIPT_STRAND].chars().next().unwrap_or('*');
    let genomic_intron_strand: char = bf[COL_GENOMIC_INTRON_STRAND].chars().next().unwrap_or('*');
    let mapq_best: u32 = best.mapq;
    // Best alignment's query span on the read and target span on the reference
    // (both 0-based half-open per PAF).
    let query_start: u64  = bf[COL_QUERY_START].parse().unwrap_or(0);
    let query_end: u64    = bf[COL_QUERY_END].parse().unwrap_or(0);
    let target_start: u64 = bf[COL_TARGET_START].parse().unwrap_or(0);
    let target_end: u64   = bf[COL_TARGET_END].parse().unwrap_or(0);

    // Per-alignment stats — all from the representative (best) alignment.
    let n_match_events: u64 = bf[COL_N_MATCH_EVENTS].parse().unwrap_or(0);
    let n_match_bases: u64 = bf[COL_N_MATCH_BASES].parse().unwrap_or(0);
    let n_sub_events: u64 = bf[COL_N_SUB_EVENTS].parse().unwrap_or(0);
    let n_sub_bases: u64 = bf[COL_N_SUB_BASES].parse().unwrap_or(0);
    let n_ins_events: u64 = bf[COL_N_INS_EVENTS].parse().unwrap_or(0);
    let n_ins_bases: u64 = bf[COL_N_INS_BASES].parse().unwrap_or(0);
    let n_del_events: u64 = bf[COL_N_DEL_EVENTS].parse().unwrap_or(0);
    let n_del_bases: u64 = bf[COL_N_DEL_BASES].parse().unwrap_or(0);
    let n_splice_events: u64 = bf[COL_N_SPLICE_EVENTS].parse().unwrap_or(0);
    let n_splice_bases: u64 = bf[COL_N_SPLICE_BASES].parse().unwrap_or(0);
    let n_sc_start: u64 = bf[COL_N_SC_START].parse().unwrap_or(0);
    let n_sc_end: u64 = bf[COL_N_SC_END].parse().unwrap_or(0);
    let n_sc_events: u64 = bf[COL_N_SC_EVENTS].parse().unwrap_or(0);
    let junctions = bf[COL_JUNCTIONS].clone();
    let cs = bf[COL_CS].clone();
    let genomic_junctions = bf[COL_GENOMIC_JUNCTIONS].clone();
    let tp_tag: char = bf.get(COL_TP_TAG).and_then(|s| s.chars().next()).unwrap_or('*');
    let aln_score: i64 = best.aln_score;
    let ms: i64 = best.ms;
    let query_aln_cov: f64 = bf[COL_QUERY_ALN_COV].parse().unwrap_or(f64::NAN);
    let query_aln_len: u64 = bf[COL_QUERY_ALN_LEN].parse().unwrap_or(0);
    let seqid: f64 = bf[COL_ALN_IDENTITY].parse().unwrap_or(f64::NAN);

    // Counts over all rows in the group.
    let mut num_aln: u64 = 0;
    let mut num_aln_tpp: u64 = 0;
    let mut num_aln_tps: u64 = 0;

    for row in rows.iter() {
        let rf = &row.fields;
        if row.is_aligned {
            num_aln += 1;
            // Split mapped rows by PAF `tp:A`; a missing tag counts toward neither.
            match rf.get(COL_TP_TAG).and_then(|s| s.chars().next()) {
                Some('P' | 'I') => num_aln_tpp += 1,
                Some('S' | 'i') => num_aln_tps += 1,
                _ => {}
            }
        }
    }

    // Count aligned rows tied with the best on ms, AS and MQ. The tp:A
    // tie-break is deliberately not a key here: a primary and a secondary with
    // equal scores still count as tied. n_alns_tied_top_score = 1 means the
    // representative won on score; > 1 means it was picked among equals
    // (primary first, then input order).
    let top_key = rows
        .iter()
        .filter(|r| r.is_aligned)
        .map(|r| (r.ms, r.aln_score, r.mapq))
        .max();
    let num_aln_maxscore: u64 = rows
        .iter()
        .filter(|r| r.is_aligned && Some((r.ms, r.aln_score, r.mapq)) == top_key)
        .count() as u64;

    // Build raw_fields for pass-through (all columns except query_name, query_len),
    // in READINFO_COLS order.
    let raw_fields: Vec<String> = vec![
        // Locus & span.
        target_chr.clone(),
        strand.to_string(),
        target_start.to_string(),
        target_end.to_string(),
        query_start.to_string(),
        query_end.to_string(),
        // Alignment selection & score.
        mapq_best.to_string(),
        tp_tag.to_string(),
        num_aln.to_string(),
        num_aln_tpp.to_string(),
        num_aln_tps.to_string(),
        num_aln_maxscore.to_string(),
        aln_score.to_string(),
        ms.to_string(),
        // Identity & coverage.
        fmt_float(seqid),
        fmt_float(query_aln_cov),
        query_aln_len.to_string(),
        // cs & strands.
        cs.clone(),
        query_transcript_strand.to_string(),
        genomic_intron_strand.to_string(),
        // cs-derived event counts.
        n_match_events.to_string(),
        n_match_bases.to_string(),
        n_sub_events.to_string(),
        n_sub_bases.to_string(),
        n_ins_events.to_string(),
        n_ins_bases.to_string(),
        n_del_events.to_string(),
        n_del_bases.to_string(),
        n_sc_events.to_string(),
        n_sc_start.to_string(),
        n_sc_end.to_string(),
        // Introns & junctions.
        n_splice_events.to_string(),
        n_splice_bases.to_string(),
        junctions.clone(),
        genomic_junctions.clone(),
    ];

    ReadInfoRow {
        read_name,
        read_len,
        target_chr,
        strand,
        query_transcript_strand,
        genomic_intron_strand,
        mapq_best,
        aln_score,
        ms,
        query_aln_cov,
        query_aln_len,
        seqid,
        junctions,
        num_aln,
        num_aln_maxscore,
        n_match_events,
        n_match_bases,
        n_sub_events,
        n_sub_bases,
        n_ins_events,
        n_ins_bases,
        n_del_events,
        n_del_bases,
        n_splice_events,
        n_splice_bases,
        n_sc_start,
        n_sc_end,
        n_sc_events,
        cs,
        genomic_junctions,
        query_start,
        query_end,
        target_start,
        target_end,
        tp_tag,
        num_aln_tpp,
        num_aln_tps,
        raw_fields,
    }
}

// ── CLI args ─────────────────────────────────────────────────────────────────
// Standalone alninfo → readinfo CLI; not registered in main.rs. Kept for possible reuse.

#[derive(clap::Args, Debug)]
#[allow(dead_code)]
pub struct ReadInfoArgs {
    /// Input alninfo TSV file (use '-' for stdin)
    #[arg(short = 'i', long = "input", value_name = "alninfo.tsv")]
    pub input: String,

    /// Output readinfo TSV file (default: stdout; '.gz' for gzip)
    #[arg(short = 'o', long = "output", value_name = "readinfo.tsv[.gz]")]
    pub output: Option<String>,
}

// ── Main entry point ─────────────────────────────────────────────────────────

#[allow(dead_code)]
pub fn run(args: &ReadInfoArgs) -> Result<()> {
    // Open input (transparently handles plain files, .gz, and stdin via '-')
    let reader = open_input(&args.input)?;

    // Open output
    let mut out = open_output(args.output.as_deref())?;
    writeln!(out, "{}", readinfo_header())?;

    let mut lines = reader.lines();

    // Skip header line
    match lines.next() {
        Some(Ok(h)) => {
            // Validate that it looks like a header
            if !h.starts_with(QUERY_NAME.name) {
                bail!("Expected header starting with '{}', got: {h}", QUERY_NAME.name);
            }
            // Build col index map for validation (optional)
            let _col_map: HashMap<String, usize> = h
                .split('\t')
                .enumerate()
                .map(|(i, s)| (s.to_owned(), i))
                .collect();
        }
        Some(Err(e)) => return Err(e.into()),
        None => return Ok(()), // empty file
    }

    let mut current_name: Option<String> = None;
    let mut group: Vec<AlnRow> = Vec::new();

    // One-shot lex-decrease guardrail. Constant memory (one prev-name).
    let mut warned_unsorted: bool = false;

    for line_res in lines {
        let line = line_res?;
        if line.is_empty() {
            continue;
        }

        let row = match parse_aln_row(&line) {
            Some(r) => r,
            None => continue,
        };

        let name = row.fields[COL_QUERY_NAME].clone();

        if let Some(ref cur) = current_name {
            if *cur != name {
                // Lex-decrease check at run boundaries (cheaper than per-row).
                if !warned_unsorted && name.as_bytes() < cur.as_bytes() {
                    eprintln!(
                        "WARNING: input alninfo is not byte-lex sorted by \
                         query_name (saw {cur:?} then {name:?}). `readinfo` \
                         groups by *contiguous* query_name runs and the \
                         downstream `compare` step requires byte-lex sort. \
                         If reads aren't contiguous, per-read summaries will \
                         be wrong. To pre-sort, either sort the alninfo:\n\
                         \n\
                         \t(head -1 in.alninfo.tsv; tail -n +2 in.alninfo.tsv \
                         | LC_ALL=C sort -t$'\\t' -k1,1) > sorted.alninfo.tsv\n\
                         \n\
                         or sort the PAF before toolkit paf2tables:\n\
                         \n\
                         \tLC_ALL=C sort -t$'\\t' -k1,1 in.paf | maligno \
                         toolkit paf2tables -i - --readinfo ...\n"
                    );
                    warned_unsorted = true;
                }
                // Flush group
                flush_group(&mut group, &mut *out)?;
                group.clear();
                current_name = Some(name);
            }
        } else {
            current_name = Some(name);
        }

        group.push(row);
    }

    // Flush last group
    if !group.is_empty() {
        flush_group(&mut group, &mut *out)?;
    }

    out.flush()?;
    Ok(())
}

#[allow(dead_code)]
fn flush_group<W: Write + ?Sized>(group: &mut Vec<AlnRow>, out: &mut W) -> Result<()> {
    if group.is_empty() {
        return Ok(());
    }
    collapse_group(group).write(out)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A full alninfo row for read `r1`: zeros everywhere except the fields
    /// `collapse_group` keys on. `target` "*" makes an unmapped placeholder; `tp`
    /// "" leaves the tp_tag cell empty.
    fn row(target: &str, ms: i64, aln_score: i64, mapq: u32, tp: &str) -> AlnRow {
        let mut f = vec!["0".to_string(); ALNINFO_COLS.len()];
        f[COL_QUERY_NAME] = "r1".into();
        f[COL_QUERY_LEN] = "100".into();
        f[COL_STRAND] = if target == "*" { "*" } else { "+" }.into();
        f[COL_QUERY_TRANSCRIPT_STRAND] = "*".into();
        f[COL_GENOMIC_INTRON_STRAND] = "*".into();
        f[COL_TARGET_NAME] = target.into();
        f[COL_MS] = ms.to_string();
        f[COL_AS] = aln_score.to_string();
        f[COL_MQ] = mapq.to_string();
        f[COL_CS] = String::new();
        f[COL_JUNCTIONS] = "()".into();
        f[COL_GENOMIC_JUNCTIONS] = "()".into();
        f[COL_TP_TAG] = tp.into();
        parse_aln_row(&f.join("\t")).expect("valid alninfo row")
    }

    /// Collapse `rows`, check the tpP + tpS <= n_alns invariant, and return
    /// (n_alns_primsupp, n_alns_secondary, n_alns).
    fn counts(mut rows: Vec<AlnRow>) -> (u64, u64, u64) {
        let ri = collapse_group(&mut rows);
        assert!(ri.num_aln_tpp + ri.num_aln_tps <= ri.num_aln);
        let at = |c: Col| ri.raw_fields[index_of(&READINFO_COLS, c) - 2].clone();
        assert_eq!(at(N_ALNS_PRIMSUPP), ri.num_aln_tpp.to_string());
        assert_eq!(at(N_ALNS_SECONDARY), ri.num_aln_tps.to_string());
        (ri.num_aln_tpp, ri.num_aln_tps, ri.num_aln)
    }

    #[test]
    fn header_has_tp_tag_after_mapq_and_matches_the_row_width() {
        let header = readinfo_header();
        let cols: Vec<&str> = header.split('\t').collect();
        assert_eq!(cols.len(), 37);
        assert_eq!(
            &cols[8..14],
            &["mapq", "tp_tag", "n_alns", "n_alns_primsupp", "n_alns_secondary", "n_alns_tied_top_score"]
        );
        let ri = collapse_group(&mut [row("chr1", 50, 50, 60, "P")]);
        assert_eq!(2 + ri.raw_fields.len(), cols.len());
    }

    #[test]
    fn single_primary() {
        assert_eq!(counts(vec![row("chr1", 50, 50, 60, "P")]), (1, 0, 1));
    }

    #[test]
    fn split_read_counts_every_primsupp_segment() {
        let rows = vec![row("chr1", 50, 50, 60, "P"), row("chr2", 30, 30, 60, "P")];
        assert_eq!(counts(rows), (2, 0, 2));
    }

    #[test]
    fn primary_plus_secondaries() {
        let rows = vec![
            row("chr1", 50, 50, 60, "P"),
            row("chr2", 50, 50, 0, "S"),
            row("chr3", 40, 40, 0, "S"),
        ];
        assert_eq!(counts(rows), (1, 2, 3));
    }

    #[test]
    fn inversions_fold_into_primsupp_and_secondary() {
        let rows = vec![
            row("chr1", 50, 50, 60, "P"),
            row("chr1", 20, 20, 60, "I"),
            row("chr2", 20, 20, 0, "i"),
        ];
        assert_eq!(counts(rows), (2, 1, 3));
    }

    #[test]
    fn mapped_rows_without_a_tp_tag_count_toward_neither() {
        let rows = vec![row("chr1", 50, 50, 60, "*"), row("chr2", 40, 40, 60, "")];
        assert_eq!(counts(rows), (0, 0, 2));
    }

    #[test]
    fn unmapped_read_is_zero_zero() {
        assert_eq!(counts(vec![row("*", 0, 0, 0, "*")]), (0, 0, 0));
    }

    #[test]
    fn representative_and_existing_aggregates_are_unchanged() {
        // A secondary ties the primary on (ms, AS) and wins on MQ: the representative
        // is still chosen by ms -> AS -> MQ, regardless of tp, and it won outright.
        let mut rows = vec![
            row("chr1", 50, 50, 10, "P"),
            row("chr2", 50, 50, 60, "S"),
            row("chr3", 20, 20, 60, "P"),
        ];
        let ri = collapse_group(&mut rows);
        assert_eq!(ri.target_chr, "chr2");
        assert_eq!(ri.tp_tag, 'S');
        assert_eq!((ri.num_aln, ri.num_aln_maxscore), (3, 1));
        assert_eq!((ri.num_aln_tpp, ri.num_aln_tps), (2, 1));
    }

    #[test]
    fn tied_count_includes_mapq_and_ignores_unmapped_rows() {
        // Two alignments tie on ms, AS and MQ; a third ties on ms and AS but has a
        // lower MQ, so it lost on a selection key and is not counted.
        let mut rows = vec![
            row("chr1", 50, 50, 60, "P"),
            row("chr2", 50, 50, 60, "S"),
            row("chr3", 50, 50, 10, "S"),
        ];
        let ri = collapse_group(&mut rows);
        assert_eq!((ri.num_aln, ri.num_aln_maxscore), (3, 2));
        let mut rows = vec![row("chr1", 50, 50, 60, "P"), row("chr2", 50, 50, 10, "S")];
        assert_eq!(collapse_group(&mut rows).num_aln_maxscore, 1, "MQ breaks the tie");
        assert_eq!(collapse_group(&mut [row("*", 0, 0, 0, "*")]).num_aln_maxscore, 0);
    }

    #[test]
    fn a_full_tie_prefers_the_primary_but_still_counts_as_tied() {
        // A secondary listed first ties the primary on ms, AS and MQ (a STAR
        // multimapper): the primary represents the read, and both count as tied.
        let mut rows = vec![row("chr2", 0, 50, 3, "S"), row("chr1", 0, 50, 3, "P")];
        let ri = collapse_group(&mut rows);
        assert_eq!((ri.target_chr.as_str(), ri.tp_tag), ("chr1", 'P'));
        assert_eq!((ri.num_aln, ri.num_aln_maxscore), (2, 2));
        // A higher-scoring secondary still wins: tp:A only breaks full ties.
        let mut rows = vec![row("chr1", 0, 50, 3, "P"), row("chr2", 0, 60, 3, "S")];
        assert_eq!(collapse_group(&mut rows).tp_tag, 'S');
        // An untagged mapped row ranks between primary and secondary, and any
        // mapped row beats the unmapped placeholder.
        let mut rows = vec![row("chr2", 0, 0, 0, "S"), row("chr3", 0, 0, 0, "*")];
        assert_eq!(collapse_group(&mut rows).target_chr, "chr3");
        let mut rows = vec![row("*", 0, 0, 0, "*"), row("chr2", 0, 0, 0, "S")];
        assert_eq!(collapse_group(&mut rows).target_chr, "chr2");
    }

    #[test]
    fn raw_fields_follow_the_readinfo_columns() {
        // One alignment whose every field holds a distinct value, so a field
        // emitted in the wrong position cannot match by accident.
        let vals: [(Col, &str); 36] = [
            (QUERY_NAME, "r9"), (QUERY_LEN, "101"), (QUERY_START, "2"), (QUERY_END, "99"),
            (STRAND, "-"), (TARGET_NAME, "chr9"), (TARGET_LEN, "5000"),
            (TARGET_START, "200"), (TARGET_END, "1300"), (N_RESIDUE_MATCHES, "90"),
            (ALN_BLOCK_LEN, "95"), (MAPQ, "42"), (TP_TAG, "S"), (AS_SCORE, "66"), (MS_SCORE, "77"),
            (ALN_IDENTITY, "0.5"), (QUERY_ALN_COV, "0.25"), (QUERY_ALN_LEN, "97"),
            (CS, ":10"), (QUERY_TRANSCRIPT_STRAND, "-"), (GENOMIC_INTRON_STRAND, "+"),
            (N_MATCH_EVENTS, "11"), (N_MATCH_BASES, "12"), (N_SUBSTITUTION_EVENTS, "13"),
            (N_SUBSTITUTION_BASES, "14"), (N_INSERTION_EVENTS, "15"), (N_INSERTION_BASES, "16"),
            (N_DELETION_EVENTS, "17"), (N_DELETION_BASES, "18"),
            (N_SOFTCLIP_EVENTS, "2"), (N_SOFTCLIP_BASES_START, "21"), (N_SOFTCLIP_BASES_END, "22"),
            (N_INTRON_EVENTS, "19"), (N_INTRON_BASES, "20"),
            (QUERY_JUNCTIONS, "(5, 6, 7)"), (GENOMIC_JUNCTIONS, "((300, 400),)"),
        ];
        let mut f = vec![String::new(); ALNINFO_COLS.len()];
        for (col, v) in vals {
            f[index_of(&ALNINFO_COLS, col)] = v.to_string();
        }
        let ri = collapse_group(&mut [parse_aln_row(&f.join("\t")).unwrap()]);

        let expected: [(Col, &str); 35] = [
            (TARGET_NAME, "chr9"), (STRAND, "-"), (TARGET_START, "200"), (TARGET_END, "1300"),
            (QUERY_START, "2"), (QUERY_END, "99"),
            (MAPQ, "42"), (TP_TAG, "S"), (N_ALNS, "1"), (N_ALNS_PRIMSUPP, "0"), (N_ALNS_SECONDARY, "1"),
            (N_ALNS_TIED_TOP_SCORE, "1"), (AS_SCORE, "66"), (MS_SCORE, "77"),
            (ALN_IDENTITY, "0.5"), (QUERY_ALN_COV, "0.25"), (QUERY_ALN_LEN, "97"),
            (CS, ":10"), (QUERY_TRANSCRIPT_STRAND, "-"), (GENOMIC_INTRON_STRAND, "+"),
            (N_MATCH_EVENTS, "11"), (N_MATCH_BASES, "12"), (N_SUBSTITUTION_EVENTS, "13"),
            (N_SUBSTITUTION_BASES, "14"), (N_INSERTION_EVENTS, "15"), (N_INSERTION_BASES, "16"),
            (N_DELETION_EVENTS, "17"), (N_DELETION_BASES, "18"),
            (N_SOFTCLIP_EVENTS, "2"), (N_SOFTCLIP_BASES_START, "21"), (N_SOFTCLIP_BASES_END, "22"),
            (N_INTRON_EVENTS, "19"), (N_INTRON_BASES, "20"),
            (QUERY_JUNCTIONS, "(5, 6, 7)"), (GENOMIC_JUNCTIONS, "((300, 400),)"),
        ];
        // Listed in readinfo order, so a reordered writer fails here.
        for (i, (col, _)) in expected.iter().enumerate() {
            assert_eq!(index_of(&READINFO_COLS, *col) - 2, i, "{:?} listed out of order", col.name);
        }
        assert_eq!(2 + ri.raw_fields.len(), READINFO_COLS.len());
        assert_eq!((ri.read_name.as_str(), ri.read_len), ("r9", 101));
        for (col, want) in expected {
            // raw_fields omits the two key columns that lead READINFO_COLS.
            let at = index_of(&READINFO_COLS, col) - 2;
            assert_eq!(ri.raw_fields[at], want, "column {:?}", col.name);
        }
    }

    #[test]
    fn scores_identity_and_coverage_come_from_the_representative() {
        // The representative wins on ms. A second alignment scores higher on AS,
        // identity, coverage and aligned length, and has the other transcript
        // strand; none of those may leak into the
        // readinfo row, which describes the representative alone.
        let aln = |target: &str, ms: i64, aln_score: i64, seqid: &str, len: &str, cov: &str| {
            let mut f = vec!["0".to_string(); ALNINFO_COLS.len()];
            f[COL_QUERY_NAME] = "r1".into();
            f[COL_QUERY_LEN] = "100".into();
            f[COL_STRAND] = "+".into();
            f[COL_QUERY_TRANSCRIPT_STRAND] = if target == "chr1" { "-" } else { "+" }.into();
            f[COL_GENOMIC_INTRON_STRAND] = if target == "chr1" { "-" } else { "+" }.into();
            f[COL_TARGET_NAME] = target.into();
            f[COL_MS] = ms.to_string();
            f[COL_AS] = aln_score.to_string();
            f[COL_ALN_IDENTITY] = seqid.into();
            f[COL_QUERY_ALN_LEN] = len.into();
            f[COL_QUERY_ALN_COV] = cov.into();
            f[COL_TP_TAG] = "P".into();
            parse_aln_row(&f.join("\t")).unwrap()
        };
        let ri = collapse_group(&mut [
            aln("chr2", 40, 90, "0.999", "95", "0.95"),
            aln("chr1", 50, 60, "0.9", "60", "0.6"),
        ]);
        assert_eq!(ri.target_chr, "chr1");
        assert_eq!((ri.ms, ri.aln_score, ri.query_aln_len), (50, 60, 60));
        assert_eq!((ri.seqid, ri.query_aln_cov), (0.9, 0.6));
        assert_eq!((ri.query_transcript_strand, ri.genomic_intron_strand), ('-', '-'));
        let at = |c: Col| ri.raw_fields[index_of(&READINFO_COLS, c) - 2].as_str();
        assert_eq!(
            [at(AS_SCORE), at(MS_SCORE), at(ALN_IDENTITY), at(QUERY_ALN_COV), at(QUERY_ALN_LEN),
             at(QUERY_TRANSCRIPT_STRAND), at(GENOMIC_INTRON_STRAND)],
            ["60", "50", "0.9", "0.6", "60", "-", "-"]
        );
        assert_eq!((ri.num_aln, ri.num_aln_maxscore), (2, 1));
    }
}
