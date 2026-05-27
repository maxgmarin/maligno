use std::io::Write;

use crate::cs_parser::{parse_cs, CsStats};
use crate::paf::PafRecord;

/// Escape TSV field value: replace newlines, carriage returns, and tabs
/// to prevent breaking TSV format when fields contain special characters.
fn escape_tsv_field(field: &str) -> String {
    field
        .replace('\\', "\\\\")  // Backslash first to avoid double-escaping
        .replace('\n', "\\n")   // Newline
        .replace('\r', "\\r")   // Carriage return
        .replace('\t', "\\t")   // Tab
}

/// Fully computed per-alignment info row, matching `paf_to_df()` output.
#[derive(Debug)]
pub struct AlnInfo {
    // ── 12 standard PAF fields ────────────────────────────────────────────
    pub query_name:          String,
    pub query_len:           u64,
    pub query_start:         u64,
    pub query_end:           u64,
    pub strand:              char,
    pub target_name:         String,
    pub target_len:          u64,
    pub target_start:        u64,
    pub target_end:          u64,
    pub num_residue_matches: u64,
    pub aln_block_len:       u64,
    pub mapq:                u32,
    // ── Optional tags ─────────────────────────────────────────────────────
    pub ms:                  i64,
    pub aln_score:           i64,
    pub cs:                  String,
    // ── cs-derived stats ──────────────────────────────────────────────────
    pub n_match_events:           u32,
    pub n_match_bases:            u64,
    pub n_substitution_events:    u32,
    pub n_substitution_bases:     u64,
    pub n_insertion_events:       u32,
    pub n_insertion_bases:        u64,
    pub n_deletion_events:        u32,
    pub n_deletion_bases:         u64,
    pub n_splice_junction_events: u32,
    pub n_splice_junction_bases:  u64,
    pub n_softclipped_bases_start: u64,
    pub n_softclipped_bases_end:   u64,
    pub n_softclipped_events:      u32,
    /// Backwards-compat alias for N_Insertion_Bases.
    pub num_bp_inserted:          u64,
    /// Splice junction query-coordinate positions, plus-strand adjusted, sorted ascending.
    pub junctions:                Vec<i64>,
    pub splice_junction_count:    usize,
    // ── Derived scalars ───────────────────────────────────────────────────
    pub target_start_1based:  u64,
    pub seqid:                f64,
    pub query_aln_len:        u64,
    pub query_aln_cov:        f64,
}

impl AlnInfo {
    /// Build an `AlnInfo` from a parsed PAF record.
    pub fn from_paf(rec: &PafRecord<'_>) -> Self {
        // ── cs stats ──────────────────────────────────────────────────────
        let cs_stats: CsStats = if rec.cs.is_empty() {
            CsStats::default()
        } else {
            parse_cs(rec.cs)
        };

        // ── Soft-clip lengths (unaligned bases at read ends) ──────────────
        // Python's get_softclip_lengths():
        //   +strand: start = query_start,          end = query_len - query_end
        //   -strand: start = query_len - query_end, end = query_start
        let (sc_start, sc_end) = softclip_lengths(
            rec.query_len, rec.query_start, rec.query_end, rec.strand,
        );
        let n_softclipped_events = u32::from(sc_start > 0) + u32::from(sc_end > 0);

        // ── Junction coordinates ──────────────────────────────────────────
        // offset: same formula as calculate_offset() in junction_utils.py
        let offset: u64 = if rec.strand == '+' {
            rec.query_start
        } else {
            // '-' or unmapped
            rec.query_len.saturating_sub(rec.query_end)
        };

        // raw junctions: offset + q_pos for each ~ op
        // plus-strand adjustment for - strand: query_len - j
        let junctions: Vec<i64> = {
            let mut v: Vec<i64> = cs_stats.raw_junctions.iter()
                .map(|&q| (offset + q) as i64)
                .collect();

            if rec.strand == '-' {
                // adjust_junctions_negative_to_plus_strand():
                //   j → query_len - j, then sort ascending
                let ql = rec.query_len as i64;
                v = v.into_iter().map(|j| ql - j).collect();
                v.sort_unstable();
            }
            v
        };

        let splice_junction_count = junctions.len();

        // ── Derived scalars ───────────────────────────────────────────────
        let target_start_1based = rec.target_start + 1;
        let seqid = if rec.aln_block_len > 0 {
            rec.num_residue_matches as f64 / rec.aln_block_len as f64
        } else {
            f64::NAN
        };
        let query_aln_len = rec.query_end.saturating_sub(rec.query_start);
        let query_aln_cov = if rec.query_len > 0 {
            query_aln_len as f64 / rec.query_len as f64
        } else {
            f64::NAN
        };

        AlnInfo {
            query_name:          rec.query_name.to_owned(),
            query_len:           rec.query_len,
            query_start:         rec.query_start,
            query_end:           rec.query_end,
            strand:              rec.strand,
            target_name:         rec.target_name.to_owned(),
            target_len:          rec.target_len,
            target_start:        rec.target_start,
            target_end:          rec.target_end,
            num_residue_matches: rec.num_residue_matches,
            aln_block_len:       rec.aln_block_len,
            mapq:                rec.mapq,
            ms:                  rec.ms,
            aln_score:           rec.aln_score,
            cs:                  rec.cs.to_owned(),
            n_match_events:           cs_stats.n_match_events,
            n_match_bases:            cs_stats.n_match_bases,
            n_substitution_events:    cs_stats.n_substitution_events,
            n_substitution_bases:     cs_stats.n_substitution_bases,
            n_insertion_events:       cs_stats.n_insertion_events,
            n_insertion_bases:        cs_stats.n_insertion_bases,
            n_deletion_events:        cs_stats.n_deletion_events,
            n_deletion_bases:         cs_stats.n_deletion_bases,
            n_splice_junction_events: cs_stats.n_splice_junction_events,
            n_splice_junction_bases:  cs_stats.n_splice_junction_bases,
            n_softclipped_bases_start: sc_start,
            n_softclipped_bases_end:   sc_end,
            n_softclipped_events,
            num_bp_inserted:          cs_stats.n_insertion_bases,
            junctions,
            splice_junction_count,
            target_start_1based,
            seqid,
            query_aln_len,
            query_aln_cov,
        }
    }

    /// Write the TSV header row (column names).
    pub fn write_header<W: Write>(w: &mut W) -> std::io::Result<()> {
        writeln!(w,
            "Query_Name\tQuery_Len\tQuery_Start\tQuery_End\tStrand\t\
             Target_Name\tTarget_Len\tTarget_Start\tTarget_End\t\
             Num_Residue_Matches\tAln_Block_Len\tMQ\t\
             ms\tAS\tcs\t\
             N_Match_Events\tN_Match_Bases\t\
             N_Substitution_Events\tN_Substitution_Bases\t\
             N_Insertion_Events\tN_Insertion_Bases\t\
             N_Deletion_Events\tN_Deletion_Bases\t\
             N_Splice_Junction_Events\tN_Splice_Junction_Bases\t\
             N_SoftClipped_Bases_Start\tN_SoftClipped_Bases_End\t\
             N_SoftClipped_Events\tnum_bp_inserted\t\
             junctions\tsplice_junction_count\t\
             Target_Start_1based\tseqid\tQuery_Aln_Len\tQuery_Aln_Cov"
        )
    }

    /// Write one TSV data row.
    pub fn write_row<W: Write>(&self, w: &mut W) -> std::io::Result<()> {
        let escaped_cs = escape_tsv_field(&self.cs);
        write!(w,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t\
             {}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t",
            self.query_name,
            self.query_len, self.query_start, self.query_end, self.strand,
            self.target_name,
            self.target_len, self.target_start, self.target_end,
            self.num_residue_matches, self.aln_block_len, self.mapq,
            self.ms, self.aln_score, escaped_cs,
            self.n_match_events,        self.n_match_bases,
            self.n_substitution_events, self.n_substitution_bases,
            self.n_insertion_events,    self.n_insertion_bases,
            self.n_deletion_events,     self.n_deletion_bases,
            self.n_splice_junction_events, self.n_splice_junction_bases,
            self.n_softclipped_bases_start, self.n_softclipped_bases_end,
            self.n_softclipped_events,
            self.num_bp_inserted,
        )?;

        // junctions: Python's str(tuple) format
        write_junction_tuple(w, &self.junctions)?;

        writeln!(w, "\t{}\t{}\t{}\t{}\t{}",
            self.splice_junction_count,
            self.target_start_1based,
            fmt_float(self.seqid),
            self.query_aln_len,
            fmt_float(self.query_aln_cov),
        )
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Helpers
// ─────────────────────────────────────────────────────────────────────────────

/// Compute soft-clip lengths, matching Python's `get_softclip_lengths()`.
fn softclip_lengths(query_len: u64, query_start: u64, query_end: u64, strand: char)
    -> (u64, u64)
{
    match strand {
        '+' => (query_start, query_len.saturating_sub(query_end)),
        _   => (query_len.saturating_sub(query_end), query_start),
    }
}

/// Format a float to match Python's `repr(float)` exactly.
///
/// Python's repr() always includes a decimal point for finite non-integer
/// floats.  Rust's Display omits the `.0` for integer-valued floats (e.g.
/// `1.0` → `"1"`).  We restore the `.0` to match Python.
fn fmt_float(v: f64) -> String {
    if v.is_nan()      { return "NaN".to_owned(); }
    if v.is_infinite() { return if v > 0.0 { "inf".to_owned() } else { "-inf".to_owned() }; }

    let s = format!("{v}");
    // If Rust's display produced no decimal point and no exponent marker,
    // the value is an integer-valued float: append ".0".
    if s.bytes().all(|b| b == b'-' || b.is_ascii_digit()) {
        format!("{s}.0")
    } else {
        s
    }
}

/// Write junctions in Python's `str(tuple)` format.
///
/// | len  | Python repr          | example               |
/// |------|----------------------|-----------------------|
/// | 0    | `()`                 | `()`                  |
/// | 1    | `(n,)`               | `(108,)`              |
/// | 2+   | `(a, b, ...)`        | `(108, 168, 359)`     |
fn write_junction_tuple<W: Write>(w: &mut W, junctions: &[i64]) -> std::io::Result<()> {
    match junctions.len() {
        0 => write!(w, "()"),
        1 => write!(w, "({},)", junctions[0]),
        _ => {
            write!(w, "(")?;
            for (i, j) in junctions.iter().enumerate() {
                if i > 0 { write!(w, ", ")?; }
                write!(w, "{j}")?;
            }
            write!(w, ")")
        }
    }
}
