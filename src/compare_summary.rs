//! Aggregate summary statistics over a per-read comparison.
//!
//! Two delivery paths share the SAME classification logic (`classify`) and the
//! SAME accumulator (`CompareSummary`):
//!
//!   1. Built into `compare`: each matched read is `observe`d as its row streams
//!      out (O(1) memory — only counters are kept), and the summary is written as
//!      `{prefix}.compare.summary.tsv` plus an stderr block.
//!   2. The standalone `compare-summary` command: streams an existing comparison
//!      TSV (`compare` / `compare-readinfo` output) row-by-row and emits the same
//!      summary. Serves the manual `paf2tables` → `compare-readinfo` workflow.
//!
//! `classify` reads per-side values by **unsuffixed** readinfo column name through
//! two accessor closures, resolving them via a name→index map, so it is insensitive
//! to column order. The columns it needs are `cs`, `Strand`, `Query_Start`,
//! `Query_End`, `TargetChr`, `Target_Start`, `Target_End`.
//!
//! cs-tag equality is **motif-blind**: intron donor/acceptor letters are blanked
//! out (via `cs_strip_splice_motifs`) before comparing, so an intron with the same
//! length/position but a differently-reported motif (e.g. minimap2's `ct..ac` vs.
//! STAR's `nn..nn` placeholder) does not by itself make two alignments "different".
//! Every other structural detail — matches, substitutions, indels, intron
//! length/position — is still compared exactly.

use std::collections::HashMap;
use std::io::Write;

use anyhow::{bail, Context, Result};

use crate::cs_parser::{cs_revcomp, cs_strip_splice_motifs};
use crate::io_utils::open_output;
use crate::table_input::{open_table, InputFormat};

/// Whether each side's representative alignment is mapped.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum MapStatus {
    BothMapped,
    OnlyAMapped,
    OnlyBMapped,
    NeitherMapped,
}

/// The classification of one matched read (one comparison row).
#[derive(Debug, Clone, Copy)]
pub(crate) struct ReadClass {
    pub map_status: MapStatus,
    /// Same alignment relative to the query (same span + identical cs, possibly
    /// via reverse-complement on opposite strands).
    pub query_identical: bool,
    /// `query_identical` was reached via the reverse-complement (opposite-strand)
    /// branch — i.e. an inverted placement.
    pub query_identical_rc: bool,
    /// Same alignment relative to the reference (same-strand query identity plus
    /// identical TargetChr and target start/end).
    pub reference_identical: bool,
}

/// A `TargetChr` value indicates an unmapped read when it is empty or `"*"`.
#[inline]
fn is_mapped(target_chr: &str) -> bool {
    !(target_chr.is_empty() || target_chr == "*")
}

/// Classify one matched read from the two per-side column accessors. Pure and
/// mode-independent — the single source of truth for both call sites.
pub(crate) fn classify<'a, FA, FB>(get_a: &FA, get_b: &FB) -> ReadClass
where
    FA: Fn(&str) -> &'a str,
    FB: Fn(&str) -> &'a str,
{
    let mapped_a = is_mapped(get_a("TargetChr"));
    let mapped_b = is_mapped(get_b("TargetChr"));

    let map_status = match (mapped_a, mapped_b) {
        (true, true) => MapStatus::BothMapped,
        (true, false) => MapStatus::OnlyAMapped,
        (false, true) => MapStatus::OnlyBMapped,
        (false, false) => MapStatus::NeitherMapped,
    };

    let mut query_identical = false;
    let mut query_identical_rc = false;
    let mut reference_identical = false;

    if mapped_a && mapped_b {
        // Query span is in forward-read coordinates in PAF (strand-independent),
        // so it must match in both the same-strand and reverse-complement cases.
        let same_span = get_a("Query_Start") == get_b("Query_Start")
            && get_a("Query_End") == get_b("Query_End");
        if same_span {
            // Compare cs tags with intron donor/acceptor motif letters blanked out
            // (`~ct..ac` and `~nn..nn` compare equal if the intron length/position
            // match). Some aligners (e.g. STAR) report `nn` placeholders instead of
            // the true motif bases that minimap2 reports for the same splice site —
            // that's a limitation of the aligner's own output, not a real alignment
            // difference, so it must not by itself make two alignments "different".
            let cs_a = cs_strip_splice_motifs(get_a("cs"));
            if get_a("Strand") == get_b("Strand") {
                let cs_b = cs_strip_splice_motifs(get_b("cs"));
                if cs_a == cs_b {
                    query_identical = true;
                    // Reference identity only applies to the same-strand case.
                    reference_identical = get_a("TargetChr") == get_b("TargetChr")
                        && get_a("Target_Start") == get_b("Target_Start")
                        && get_a("Target_End") == get_b("Target_End");
                }
            } else if cs_a == cs_strip_splice_motifs(&cs_revcomp(get_b("cs"))) {
                // Opposite strands but the alignment is an exact reverse-complement
                // (e.g. an inverted locus between two assemblies).
                query_identical = true;
                query_identical_rc = true;
            }
        }
    }

    ReadClass {
        map_status,
        query_identical,
        query_identical_rc,
        reference_identical,
    }
}

/// Streaming accumulator of comparison summary statistics (all O(1) counters).
#[derive(Default, Debug)]
pub(crate) struct CompareSummary {
    pub reads_compared: u64,
    pub aligned_both: u64,
    pub aligned_only_a: u64,
    pub aligned_only_b: u64,
    pub aligned_neither: u64,
    pub query_identical: u64,
    pub query_identical_same_strand: u64,
    pub query_identical_rc: u64,
    pub reference_identical: u64,
    pub a_only_by_id: u64,
    pub b_only_by_id: u64,
}

impl CompareSummary {
    /// Tally one matched read.
    pub fn observe(&mut self, c: &ReadClass) {
        self.reads_compared += 1;
        match c.map_status {
            MapStatus::BothMapped => self.aligned_both += 1,
            MapStatus::OnlyAMapped => self.aligned_only_a += 1,
            MapStatus::OnlyBMapped => self.aligned_only_b += 1,
            MapStatus::NeitherMapped => self.aligned_neither += 1,
        }
        if c.query_identical {
            self.query_identical += 1;
            if c.query_identical_rc {
                self.query_identical_rc += 1;
            } else {
                self.query_identical_same_strand += 1;
            }
        }
        if c.reference_identical {
            self.reference_identical += 1;
        }
    }

    /// Record a read present only in set A by read-ID (not in B's PAF at all).
    pub fn note_a_only_id(&mut self) {
        self.a_only_by_id += 1;
    }

    /// Record a read present only in set B by read-ID.
    pub fn note_b_only_id(&mut self) {
        self.b_only_by_id += 1;
    }

    /// query-different among both-mapped reads.
    fn query_not_identical(&self) -> u64 {
        self.aligned_both - self.query_identical
    }

    /// The ordered (category, count) rows — the single layout used by both the
    /// TSV writer and the stderr renderer.
    ///
    /// Category names use the fixed `A` / `B` side identifiers, so summary keys
    /// are stable across datasets and machine-parseable without knowing the
    /// labels. The human-readable labels are emitted separately as `label_A` /
    /// `label_B` provenance rows (see `write_tsv` / `render_stderr`).
    fn rows(&self) -> Vec<(String, u64)> {
        vec![
            ("reads_compared".to_string(), self.reads_compared),
            ("aligned_both".to_string(), self.aligned_both),
            ("aligned_only_A".to_string(), self.aligned_only_a),
            ("aligned_only_B".to_string(), self.aligned_only_b),
            ("aligned_neither".to_string(), self.aligned_neither),
            ("query_identical".to_string(), self.query_identical),
            ("query_identical_same_strand".to_string(), self.query_identical_same_strand),
            ("query_identical_revcomp".to_string(), self.query_identical_rc),
            ("query_not_identical".to_string(), self.query_not_identical()),
            ("reference_identical".to_string(), self.reference_identical),
            ("present_only_in_A_by_id".to_string(), self.a_only_by_id),
            ("present_only_in_B_by_id".to_string(), self.b_only_by_id),
        ]
    }

    /// Write the summary as a 2-column TSV (`Category<TAB>Count`), preceded by
    /// the two label provenance rows so the file is self-describing.
    pub fn write_tsv(&self, path: &str, label_a: &str, label_b: &str) -> Result<()> {
        let mut w = open_output(Some(path))?;
        writeln!(w, "Category\tCount")?;
        writeln!(w, "label_A\t{label_a}")?;
        writeln!(w, "label_B\t{label_b}")?;
        for (k, v) in self.rows() {
            writeln!(w, "{k}\t{v}")?;
        }
        w.flush()?;
        Ok(())
    }

    /// Print a human-readable block to stderr.
    pub fn render_stderr(&self, label_a: &str, label_b: &str) {
        eprintln!("Comparison summary:");
        eprintln!("  {:<34} {label_a}", "label_A");
        eprintln!("  {:<34} {label_b}", "label_B");
        for (k, v) in self.rows() {
            eprintln!("  {k:<34} {v}");
        }
    }
}

// ── Standalone `compare-summary` command ────────────────────────────────────────

#[derive(clap::Args, Debug)]
pub struct CompareSummaryArgs {
    /// Comparison table from `compare` / `compare-readinfo`: TSV (`.gz` ok;
    /// `-` = stdin) or Parquet (`--format parquet` / `-o x.parquet`). See
    /// `--input-format`.
    #[arg(short = 'i', long = "input", value_name = "compare.tsv|compare.parquet")]
    input: String,

    /// Input serialization. `auto` (default) selects Parquet for a
    /// `.parquet`-named `--input` and TSV otherwise. Parquet requires a real
    /// file path — it cannot be read from stdin (`-`).
    #[arg(long = "input-format", value_enum, default_value = "auto")]
    input_format: InputFormat,

    /// Write the summary TSV here (default: stderr only). `.gz` ok; `-` for stdout.
    #[arg(short = 'o', long = "output", value_name = "summary.tsv")]
    output: Option<String>,
}

/// Confirm a compare-table header uses the fixed `_A` / `_B` side suffixes
/// introduced in v0.13.0. Shared by `compare-summary` and `find-query-diff`.
///
/// Pre-v0.13.0 tables suffixed per-side columns with the dataset *label*
/// (`TargetChr_Splice`), which made column names dataset-specific and ambiguous
/// whenever a label itself contained an underscore. Those tables are not
/// readable by this version — the error tells the user to regenerate.
pub(crate) fn require_ab_schema(cols: &[&str]) -> Result<()> {
    let has_a = cols.contains(&"TargetChr_A");
    let has_b = cols.contains(&"TargetChr_B");
    if has_a && has_b {
        return Ok(());
    }
    // A legacy table is recognizable by label-suffixed `TargetChr_*` columns;
    // name them in the error so the cause is obvious.
    let legacy: Vec<&str> = cols
        .iter()
        .copied()
        .filter(|c| c.starts_with("TargetChr_"))
        .collect();
    if !legacy.is_empty() {
        bail!(
            "this comparison table uses the pre-v0.13.0 label-suffixed schema ({}) \
             — regenerate it with maligno v0.13+ (`compare` / `compare-readinfo`), \
             which writes fixed `TargetChr_A` / `TargetChr_B` columns plus \
             `Label_A` / `Label_B`",
            legacy.join(", ")
        );
    }
    bail!(
        "comparison table is missing the `TargetChr_A` / `TargetChr_B` columns \
         — is this a maligno compare table?"
    );
}

/// Recover the human-readable set labels from a compare table's first data row
/// (the `Label_A` / `Label_B` columns). Falls back to `("A", "B")` when those
/// columns are absent or the table has no data rows, so label reporting is
/// best-effort and never blocks the actual comparison.
pub(crate) fn labels_from_row(col_index: &HashMap<&str, usize>, fields: &[&str]) -> (String, String) {
    let pick = |name: &str, fallback: &str| -> String {
        col_index
            .get(name)
            .and_then(|&i| fields.get(i))
            .copied()
            .filter(|s| !s.is_empty())
            .unwrap_or(fallback)
            .to_string()
    };
    (pick("Label_A", "A"), pick("Label_B", "B"))
}

pub fn run(args: &CompareSummaryArgs) -> Result<()> {
    // The unsuffixed columns `classify` reads; resolve each side's index up front.
    const NEEDED: [&str; 7] = [
        "TargetChr", "Strand", "cs", "Query_Start", "Query_End", "Target_Start", "Target_End",
    ];

    // For a Parquet input, project down to just the 16 columns this command
    // reads (of the comparison table's 96) — Parquet skips decoding the rest,
    // which is where its per-column storage actually pays off.
    let wanted: Vec<String> = ["Label_A", "Label_B"]
        .into_iter()
        .map(String::from)
        .chain(["A", "B"].iter().flat_map(|side| NEEDED.iter().map(move |b| format!("{b}_{side}"))))
        .collect();
    let wanted_refs: Vec<&str> = wanted.iter().map(String::as_str).collect();
    let mut source = open_table(&args.input, args.input_format, Some(&wanted_refs))?;

    // Header → column index map.
    let cols_owned = source
        .header()
        .with_context(|| format!("reading comparison table '{}'", args.input))?;
    let cols: Vec<&str> = cols_owned.iter().map(String::as_str).collect();
    let col_index: HashMap<&str, usize> =
        cols.iter().copied().enumerate().map(|(i, c)| (c, i)).collect();

    // Require the v0.13+ fixed `_A` / `_B` side suffixes.
    require_ab_schema(&cols)?;

    let resolve = |side: &str| -> Result<HashMap<&'static str, usize>> {
        let mut m = HashMap::new();
        for base in NEEDED {
            let name = format!("{base}_{side}");
            let idx = *col_index
                .get(name.as_str())
                .with_context(|| format!("comparison table is missing column '{name}'"))?;
            m.insert(base, idx);
        }
        Ok(m)
    };
    let idx_a = resolve("A")?;
    let idx_b = resolve("B")?;

    // The human-readable labels live in each row's `Label_A` / `Label_B`, so they
    // are picked up from the first data row and used for reporting only.
    let mut label_a = "A".to_string();
    let mut label_b = "B".to_string();
    let mut seen_row = false;

    let mut summary = CompareSummary::default();
    while let Some(fields_owned) = source.next_row()? {
        let fields: Vec<&str> = fields_owned.iter().map(String::as_str).collect();
        if !seen_row {
            (label_a, label_b) = labels_from_row(&col_index, &fields);
            seen_row = true;
        }
        let get_a = |c: &str| -> &str {
            idx_a.get(c).and_then(|&i| fields.get(i)).copied().unwrap_or("")
        };
        let get_b = |c: &str| -> &str {
            idx_b.get(c).and_then(|&i| fields.get(i)).copied().unwrap_or("")
        };
        summary.observe(&classify(&get_a, &get_b));
    }

    if let Some(path) = &args.output {
        summary.write_tsv(path, &label_a, &label_b)?;
        eprintln!("Wrote summary: {path}");
    }
    summary.render_stderr(&label_a, &label_b);
    eprintln!(
        "  (note: a comparison table holds only reads present in both sets, so \
         present_only_in_* by-ID counts are 0 here.)"
    );
    Ok(())
}

// ───────────────────────────────────────────────────────────────────────────
// Tests
// ───────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    fn index<'a>(cols: &[&'a str]) -> HashMap<&'a str, usize> {
        cols.iter().copied().enumerate().map(|(i, c)| (c, i)).collect()
    }

    // ─── require_ab_schema ─────────────────────────────────────────────────

    #[test]
    fn ab_schema_accepted() {
        let cols = ["Read_Name", "Label_A", "Label_B", "TargetChr_A", "TargetChr_B"];
        assert!(require_ab_schema(&cols).is_ok());
    }

    #[test]
    fn legacy_label_suffixed_schema_is_rejected_with_regenerate_hint() {
        // A pre-v0.13.0 table: per-side columns suffixed with the dataset label.
        let cols = ["Read_Name", "TargetChr_Splice", "TargetChr_SpliceHQ"];
        let err = require_ab_schema(&cols).unwrap_err().to_string();
        assert!(err.contains("pre-v0.13.0"), "unexpected error: {err}");
        assert!(err.contains("regenerate"), "unexpected error: {err}");
        // Names the offending columns so the cause is obvious.
        assert!(err.contains("TargetChr_Splice"), "unexpected error: {err}");
    }

    #[test]
    fn one_sided_ab_schema_is_rejected() {
        // `TargetChr_B` missing → not a usable comparison table.
        let cols = ["Read_Name", "TargetChr_A"];
        assert!(require_ab_schema(&cols).is_err());
    }

    #[test]
    fn non_compare_table_is_rejected() {
        // No `TargetChr_*` columns at all (e.g. an alninfo/readinfo table).
        let cols = ["Query_Name", "Query_Len", "TargetChr"];
        let err = require_ab_schema(&cols).unwrap_err().to_string();
        assert!(err.contains("missing the `TargetChr_A`"), "unexpected error: {err}");
    }

    // ─── labels_from_row ───────────────────────────────────────────────────

    #[test]
    fn labels_read_from_first_data_row() {
        let cols = ["Read_Name", "Read_Len", "Label_A", "Label_B", "TargetChr_A"];
        let fields = ["r1", "100", "Splice", "SpliceHQ", "chr1"];
        assert_eq!(
            labels_from_row(&index(&cols), &fields),
            ("Splice".to_string(), "SpliceHQ".to_string())
        );
    }

    #[test]
    fn labels_fall_back_when_columns_absent() {
        // No Label_A/Label_B columns → report the bare side identifiers.
        let cols = ["Read_Name", "TargetChr_A", "TargetChr_B"];
        let fields = ["r1", "chr1", "chr1"];
        assert_eq!(
            labels_from_row(&index(&cols), &fields),
            ("A".to_string(), "B".to_string())
        );
    }

    #[test]
    fn labels_fall_back_on_empty_values_and_short_rows() {
        let cols = ["Read_Name", "Label_A", "Label_B"];
        // Empty label value → fall back rather than reporting "".
        assert_eq!(
            labels_from_row(&index(&cols), &["r1", "", "SpliceHQ"]),
            ("A".to_string(), "SpliceHQ".to_string())
        );
        // Truncated row (fewer fields than the header) must not panic.
        assert_eq!(
            labels_from_row(&index(&cols), &["r1"]),
            ("A".to_string(), "B".to_string())
        );
    }

    // ─── label validation ──────────────────────────────────────────────────

    #[test]
    fn set_labels_reject_corrupting_characters() {
        use crate::compare_streaming::validate_set_label;
        assert!(validate_set_label("").is_err());
        assert!(validate_set_label("a\tb").is_err());
        assert!(validate_set_label("a\nb").is_err());
        assert!(validate_set_label("a/b").is_err());
        assert!(validate_set_label("a\\b").is_err());
    }

    #[test]
    fn set_labels_allow_underscores() {
        use crate::compare_streaming::validate_set_label;
        // The whole point of fixed `_A`/`_B` suffixes: an underscore in a label
        // is no longer ambiguous, so it must be accepted.
        assert_eq!(validate_set_label("my_run_1").unwrap(), "my_run_1");
        assert_eq!(validate_set_label("Splice").unwrap(), "Splice");
    }
}
