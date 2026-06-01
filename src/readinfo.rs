use std::collections::HashMap;
use std::io::{BufRead, Write};

use anyhow::{bail, Result};

use crate::io_utils::{escape_tsv_field, fmt_float, open_input, open_output};
use crate::junction::junction_count_str;

// ── Column indices in alninfo TSV ────────────────────────────────────────────
const COL_QUERY_NAME: usize = 0;
const COL_QUERY_LEN: usize = 1;
const COL_QUERY_START: usize = 2;
const COL_QUERY_END: usize = 3;
const COL_STRAND: usize = 4;
const COL_TARGET_NAME: usize = 5;
const COL_TARGET_START: usize = 7;
const COL_TARGET_END: usize = 8;
const COL_MS: usize = 12;
const COL_AS: usize = 13;
const COL_CS: usize = 14;
const COL_N_MATCH_EVENTS: usize = 15;
const COL_N_MATCH_BASES: usize = 16;
const COL_N_SUB_EVENTS: usize = 17;
const COL_N_SUB_BASES: usize = 18;
const COL_N_INS_EVENTS: usize = 19;
const COL_N_INS_BASES: usize = 20;
const COL_N_DEL_EVENTS: usize = 21;
const COL_N_DEL_BASES: usize = 22;
const COL_N_SPLICE_EVENTS: usize = 23;
const COL_N_SPLICE_BASES: usize = 24;
const COL_N_SC_START: usize = 25;
const COL_N_SC_END: usize = 26;
const COL_N_SC_EVENTS: usize = 27;
const COL_JUNCTIONS: usize = 28;
const COL_SEQID: usize = 31;
const COL_QUERY_ALN_LEN: usize = 32;
const COL_QUERY_ALN_COV: usize = 33;
const COL_GENOMIC_JUNCTIONS: usize = 34;

pub const READINFO_HEADER: &str = "Read_Name\tRead_Len\t\
    TargetChr\tStrand\t\
    AS_Max\tms_Max\t\
    Query_Aln_Cov_Max\tQuery_Aln_Len_Max\t\
    seqid_Max\t\
    junctions\t\
    Num_Aln\tNum_Aln_MaxScore\tJuncCount\t\
    N_Match_Events\tN_Match_Bases\t\
    N_Substitution_Events\tN_Substitution_Bases\t\
    N_Insertion_Events\tN_Insertion_Bases\t\
    N_Deletion_Events\tN_Deletion_Bases\t\
    N_Splice_Junction_Events\tN_Splice_Junction_Bases\t\
    N_SoftClipped_Bases_Start\tN_SoftClipped_Bases_End\t\
    N_SoftClipped_Events\t\
    cs\t\
    genomic_junctions\t\
    Query_Start\tQuery_End\t\
    Target_Start\tTarget_End";

// ── ReadInfo row ─────────────────────────────────────────────────────────────

/// Per-read summary row, produced by the `readinfo` subcommand.
/// Fields are typed for use in the `compare` subcommand; `raw_fields` carries
/// all non-key columns pre-formatted for TSV pass-through.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct ReadInfoRow {
    pub read_name: String,
    pub read_len: u64,
    pub target_chr: String,
    pub strand: char,
    pub as_max: i64,
    pub ms_max: i64,
    pub query_aln_cov_max: f64,
    pub query_aln_len_max: u64,
    pub seqid_max: f64,
    pub junctions: String,
    pub num_aln: u64,
    pub num_aln_maxscore: u64,
    pub junc_count: usize,
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
    /// All fields except Read_Name and Read_Len, pre-formatted for pass-through.
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

// ── Helper: NaN-aware f64 max / min ─────────────────────────────────────────

fn f64_max_ignore_nan(a: f64, b: f64) -> f64 {
    if a.is_nan() {
        b
    } else if b.is_nan() {
        a
    } else if a >= b {
        a
    } else {
        b
    }
}

#[allow(dead_code)]
fn f64_min_ignore_nan(a: f64, b: f64) -> f64 {
    if a.is_nan() {
        b
    } else if b.is_nan() {
        a
    } else if a <= b {
        a
    } else {
        b
    }
}

// ── Alninfo row (lightweight, strings only) ───────────────────────────────────

struct AlnRow {
    fields: Vec<String>,
    ms: i64,
    aln_score: i64,
    is_aligned: bool, // Target_Name != "*"
}

fn parse_aln_row(line: &str) -> Option<AlnRow> {
    let fields: Vec<String> = line.split('\t').map(str::to_owned).collect();
    if fields.len() <= COL_GENOMIC_JUNCTIONS {
        return None;
    }
    let ms: i64 = fields[COL_MS].parse().unwrap_or(0);
    let aln_score: i64 = fields[COL_AS].parse().unwrap_or(0);
    let is_aligned = fields[COL_TARGET_NAME] != "*";
    Some(AlnRow {
        fields,
        ms,
        aln_score,
        is_aligned,
    })
}

fn collapse_group(rows: &mut [AlnRow]) -> ReadInfoRow {
    // Sort: ms desc, AS desc
    rows.sort_by(|a, b| {
        b.ms
            .cmp(&a.ms)
            .then_with(|| b.aln_score.cmp(&a.aln_score))
    });

    let best = &rows[0];
    let bf = &best.fields;

    let read_name = bf[COL_QUERY_NAME].clone();
    let read_len: u64 = bf[COL_QUERY_LEN].parse().unwrap_or(0);
    let target_chr = bf[COL_TARGET_NAME].clone();
    let strand: char = bf[COL_STRAND].chars().next().unwrap_or('*');
    // Best alignment's query span on the read and target span on the reference
    // (both 0-based half-open per PAF).
    let query_start: u64  = bf[COL_QUERY_START].parse().unwrap_or(0);
    let query_end: u64    = bf[COL_QUERY_END].parse().unwrap_or(0);
    let target_start: u64 = bf[COL_TARGET_START].parse().unwrap_or(0);
    let target_end: u64   = bf[COL_TARGET_END].parse().unwrap_or(0);

    // Best-alignment stats
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
    let junc_count = junction_count_str(&junctions);

    // Aggregates over all rows
    let mut as_max = i64::MIN;
    let mut ms_max = i64::MIN;
    let mut query_aln_cov_max: f64 = f64::NAN;
    let mut query_aln_len_max: u64 = 0;
    let mut seqid_max: f64 = f64::NAN;
    let mut num_aln: u64 = 0;

    for row in rows.iter() {
        let rf = &row.fields;
        let row_as: i64 = rf[COL_AS].parse().unwrap_or(0);
        let row_ms: i64 = rf[COL_MS].parse().unwrap_or(0);
        let row_cov: f64 = rf[COL_QUERY_ALN_COV]
            .parse()
            .unwrap_or(f64::NAN);
        let row_aln_len: u64 = rf[COL_QUERY_ALN_LEN].parse().unwrap_or(0);
        let row_seqid: f64 = rf[COL_SEQID].parse().unwrap_or(f64::NAN);

        if row_as > as_max {
            as_max = row_as;
        }
        if row_ms > ms_max {
            ms_max = row_ms;
        }
        query_aln_cov_max = f64_max_ignore_nan(query_aln_cov_max, row_cov);
        if row_aln_len > query_aln_len_max {
            query_aln_len_max = row_aln_len;
        }
        seqid_max = f64_max_ignore_nan(seqid_max, row_seqid);

        if row.is_aligned {
            num_aln += 1;
        }
    }

    // Handle single-row edge cases for min/max
    if as_max == i64::MIN {
        as_max = 0;
    }
    if ms_max == i64::MIN {
        ms_max = 0;
    }

    // Count of alignments tied at the chosen-best sort key — i.e., tied at
    // BOTH ms_max AND the highest AS among ms-tied rows. This matches the
    // full (ms desc, AS desc) selection rule used to pick the best alignment.
    // Num_Aln_MaxScore = 1 ⇒ a single unambiguous winner; > 1 ⇒ file order
    // ultimately broke the tie. Only aligned rows contribute, so unaligned
    // reads (Num_Aln = 0) get Num_Aln_MaxScore = 0.
    //
    // Practical note: STAR-style data writes ms=0 for every alignment, so
    // ms_max=0 and the AS tie-break does the real work. Counting at (ms, AS)
    // makes Num_Aln_MaxScore meaningful in that case.
    let as_at_ms_max: i64 = rows
        .iter()
        .filter(|r| r.is_aligned && r.ms == ms_max)
        .map(|r| r.aln_score)
        .max()
        .unwrap_or(0);
    let num_aln_maxscore: u64 = rows
        .iter()
        .filter(|r| r.is_aligned && r.ms == ms_max && r.aln_score == as_at_ms_max)
        .count() as u64;

    // Build raw_fields for pass-through (all columns except Read_Name, Read_Len)
    let raw_fields: Vec<String> = vec![
        target_chr.clone(),
        strand.to_string(),
        as_max.to_string(),
        ms_max.to_string(),
        fmt_float(query_aln_cov_max),
        query_aln_len_max.to_string(),
        fmt_float(seqid_max),
        junctions.clone(),
        num_aln.to_string(),
        num_aln_maxscore.to_string(),
        junc_count.to_string(),
        n_match_events.to_string(),
        n_match_bases.to_string(),
        n_sub_events.to_string(),
        n_sub_bases.to_string(),
        n_ins_events.to_string(),
        n_ins_bases.to_string(),
        n_del_events.to_string(),
        n_del_bases.to_string(),
        n_splice_events.to_string(),
        n_splice_bases.to_string(),
        n_sc_start.to_string(),
        n_sc_end.to_string(),
        n_sc_events.to_string(),
        cs.clone(),
        genomic_junctions.clone(),
        query_start.to_string(),
        query_end.to_string(),
        target_start.to_string(),
        target_end.to_string(),
    ];

    ReadInfoRow {
        read_name,
        read_len,
        target_chr,
        strand,
        as_max,
        ms_max,
        query_aln_cov_max,
        query_aln_len_max,
        seqid_max,
        junctions,
        num_aln,
        num_aln_maxscore,
        junc_count,
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
        raw_fields,
    }
}

// ── CLI args ─────────────────────────────────────────────────────────────────

#[derive(clap::Args, Debug)]
pub struct ReadInfoArgs {
    /// Input alninfo TSV file (use '-' for stdin)
    #[arg(short = 'i', long = "input", value_name = "alninfo.tsv")]
    pub input: String,

    /// Output readinfo TSV file (default: stdout; '.gz' for gzip)
    #[arg(short = 'o', long = "output", value_name = "readinfo.tsv[.gz]")]
    pub output: Option<String>,

    /// Skip sorting alignments within each read group by ms/AS (use when input
    /// is already ms-sorted per group)
    #[arg(long = "no-sort")]
    pub no_sort: bool,
}

// ── Main entry point ─────────────────────────────────────────────────────────

pub fn run(args: &ReadInfoArgs) -> Result<()> {
    // Open input (transparently handles plain files, .gz, and stdin via '-')
    let reader = open_input(&args.input)?;

    // Open output
    let mut out = open_output(args.output.as_deref())?;
    writeln!(out, "{READINFO_HEADER}")?;

    let mut lines = reader.lines();

    // Skip header line
    match lines.next() {
        Some(Ok(h)) => {
            // Validate that it looks like a header
            if !h.starts_with("Query_Name") {
                bail!("Expected header starting with 'Query_Name', got: {h}");
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

    let no_sort = args.no_sort;
    let mut current_name: Option<String> = None;
    let mut group: Vec<AlnRow> = Vec::new();

    // One-shot lex-decrease guardrail. Constant memory (one prev-name).
    // Fires once on the first decrease and stays quiet thereafter so it
    // doesn't spam the log on heavily-shuffled inputs.
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
                         Query_Name (saw {cur:?} then {name:?}). `readinfo` \
                         groups by *contiguous* Query_Name runs and the \
                         downstream `compare` step requires byte-lex sort. \
                         If reads aren't contiguous, per-read summaries will \
                         be wrong. To pre-sort, either sort the alninfo:\n\
                         \n\
                         \t(head -1 in.alninfo.tsv; tail -n +2 in.alninfo.tsv \
                         | LC_ALL=C sort -t$'\\t' -k1,1) > sorted.alninfo.tsv\n\
                         \n\
                         or sort the PAF before paf2alninfo:\n\
                         \n\
                         \tLC_ALL=C sort -t$'\\t' -k1,1 in.paf | maligno \
                         paf2alninfo -i - -o ...\n"
                    );
                    warned_unsorted = true;
                }
                // Flush group
                flush_group(&mut group, no_sort, &mut *out)?;
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
        flush_group(&mut group, no_sort, &mut *out)?;
    }

    out.flush()?;
    Ok(())
}

fn flush_group<W: Write + ?Sized>(
    group: &mut Vec<AlnRow>,
    no_sort: bool,
    out: &mut W,
) -> Result<()> {
    if group.is_empty() {
        return Ok(());
    }
    if !no_sort {
        // Sort handled inside collapse_group
    }
    let ri = if no_sort {
        // When no_sort: first row is considered best (already ms-sorted)
        collapse_group_nosort(group)
    } else {
        collapse_group(group)
    };
    ri.write(out)?;
    Ok(())
}

fn collapse_group_nosort(rows: &[AlnRow]) -> ReadInfoRow {
    // Same as collapse_group but don't sort; row[0] is best
    let best = &rows[0];
    let bf = &best.fields;

    let read_name = bf[COL_QUERY_NAME].clone();
    let read_len: u64 = bf[COL_QUERY_LEN].parse().unwrap_or(0);
    let target_chr = bf[COL_TARGET_NAME].clone();
    let strand: char = bf[COL_STRAND].chars().next().unwrap_or('*');
    // Best alignment's query span on the read and target span on the reference
    // (both 0-based half-open per PAF).
    let query_start: u64  = bf[COL_QUERY_START].parse().unwrap_or(0);
    let query_end: u64    = bf[COL_QUERY_END].parse().unwrap_or(0);
    let target_start: u64 = bf[COL_TARGET_START].parse().unwrap_or(0);
    let target_end: u64   = bf[COL_TARGET_END].parse().unwrap_or(0);

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
    let junc_count = junction_count_str(&junctions);

    let mut as_max = i64::MIN;
    let mut ms_max = i64::MIN;
    let mut query_aln_cov_max: f64 = f64::NAN;
    let mut query_aln_len_max: u64 = 0;
    let mut seqid_max: f64 = f64::NAN;
    let mut num_aln: u64 = 0;

    for row in rows.iter() {
        let rf = &row.fields;
        let row_as: i64 = rf[COL_AS].parse().unwrap_or(0);
        let row_ms: i64 = rf[COL_MS].parse().unwrap_or(0);
        let row_cov: f64 = rf[COL_QUERY_ALN_COV].parse().unwrap_or(f64::NAN);
        let row_aln_len: u64 = rf[COL_QUERY_ALN_LEN].parse().unwrap_or(0);
        let row_seqid: f64 = rf[COL_SEQID].parse().unwrap_or(f64::NAN);

        if row_as > as_max { as_max = row_as; }
        if row_ms > ms_max { ms_max = row_ms; }
        query_aln_cov_max = f64_max_ignore_nan(query_aln_cov_max, row_cov);
        if row_aln_len > query_aln_len_max { query_aln_len_max = row_aln_len; }
        seqid_max = f64_max_ignore_nan(seqid_max, row_seqid);
        if row.is_aligned { num_aln += 1; }
    }

    if as_max == i64::MIN { as_max = 0; }
    if ms_max == i64::MIN { ms_max = 0; }

    // Count of alignments tied at the chosen-best sort key — i.e., tied at
    // BOTH ms_max AND the highest AS among ms-tied rows. This matches the
    // full (ms desc, AS desc) selection rule used to pick the best alignment.
    // Num_Aln_MaxScore = 1 ⇒ a single unambiguous winner; > 1 ⇒ file order
    // ultimately broke the tie. Only aligned rows contribute, so unaligned
    // reads (Num_Aln = 0) get Num_Aln_MaxScore = 0.
    //
    // Practical note: STAR-style data writes ms=0 for every alignment, so
    // ms_max=0 and the AS tie-break does the real work. Counting at (ms, AS)
    // makes Num_Aln_MaxScore meaningful in that case.
    let as_at_ms_max: i64 = rows
        .iter()
        .filter(|r| r.is_aligned && r.ms == ms_max)
        .map(|r| r.aln_score)
        .max()
        .unwrap_or(0);
    let num_aln_maxscore: u64 = rows
        .iter()
        .filter(|r| r.is_aligned && r.ms == ms_max && r.aln_score == as_at_ms_max)
        .count() as u64;

    let raw_fields: Vec<String> = vec![
        target_chr.clone(),
        strand.to_string(),
        as_max.to_string(),
        ms_max.to_string(),
        fmt_float(query_aln_cov_max),
        query_aln_len_max.to_string(),
        fmt_float(seqid_max),
        junctions.clone(),
        num_aln.to_string(),
        num_aln_maxscore.to_string(),
        junc_count.to_string(),
        n_match_events.to_string(),
        n_match_bases.to_string(),
        n_sub_events.to_string(),
        n_sub_bases.to_string(),
        n_ins_events.to_string(),
        n_ins_bases.to_string(),
        n_del_events.to_string(),
        n_del_bases.to_string(),
        n_splice_events.to_string(),
        n_splice_bases.to_string(),
        n_sc_start.to_string(),
        n_sc_end.to_string(),
        n_sc_events.to_string(),
        cs.clone(),
        genomic_junctions.clone(),
        query_start.to_string(),
        query_end.to_string(),
        target_start.to_string(),
        target_end.to_string(),
    ];

    ReadInfoRow {
        read_name,
        read_len,
        target_chr,
        strand,
        as_max,
        ms_max,
        query_aln_cov_max,
        query_aln_len_max,
        seqid_max,
        junctions,
        num_aln,
        num_aln_maxscore,
        junc_count,
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
        raw_fields,
    }
}
