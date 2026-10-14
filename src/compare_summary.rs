//! Aggregate summary statistics over a per-read comparison.
//!
//! Two delivery paths share the SAME classification logic (`classify`) and the
//! SAME accumulator (`CompareSummary`):
//!
//!   1. Built into `compare`: each matched read is `observe`d as its row streams
//!      out (O(1) memory — only counters are kept), and the summary is written as
//!      `{prefix}.compare.summary.tsv` plus an stderr block.
//!   2. The `toolkit summary` command: streams an existing comparison
//!      table (`compare` / `toolkit merge-readinfo` output) row-by-row
//!      and emits the same summary. Serves the manual `toolkit
//!      paf2tables` → `toolkit merge-readinfo` workflow.
//!
//! `classify` reads per-side values by **unsuffixed** readinfo column name through
//! two accessor closures, resolving them via a name→index map, so it is insensitive
//! to column order. The columns it needs are `cs`, `strand`, `query_start`,
//! `query_end`, `target_name`, `target_start`, `target_end`.
//!
//! cs-tag equality is **motif-blind**: intron donor/acceptor letters are blanked
//! out (via `cs_strip_splice_motifs`) before comparing, so an intron with the same
//! length/position but a differently-reported motif (e.g. minimap2's `ct..ac` vs.
//! STAR's `nn..nn` placeholder) does not by itself make two alignments "different".
//! Every other structural detail — matches, substitutions, indels, intron
//! length/position — is still compared exactly.
//!
//! Reference-space classification (`RefClass`) is a second, independent axis for
//! both-mapped reads: whether the two sides land at the same genomic **position**
//! (`target_name` + `strand` + `target_start` — the cs tag's own operations
//! determine the alignment's length, so a matching start plus matching cs implies
//! a matching end too) and/or report the same **alignment** (`cs`, motif-blind).
//! This is a strict, literal comparison: unlike `query_identical`, there is no
//! reverse-complement accommodation — a real strand difference always means a
//! different position.
//!
//! `query_identical` also counts `aligned_neither` reads (neither side mapped
//! at all): both aligners agreeing a read doesn't map is itself a form of
//! agreement, not a disagreement to lump in with genuine `aligned_both`
//! mismatches. `RefClass`/junction-set identity stay a strictly "both mapped"
//! concept (`None` for `aligned_neither`, as for `aligned_only_A`/`_B`).

use std::collections::HashMap;
use std::io::Write;

use anyhow::{bail, Context, Result};

use crate::columns::{
    sided, Col, CLASSIFY_SIDE_COLS, CS, GENOMIC_JUNCTIONS, QUERY_JUNCTIONS, LABEL_A, LABEL_B,
    QUERY_END, QUERY_IDENTICAL_REVCOMP, QUERY_IDENTICAL_SAME_STRAND, QUERY_JUNCTIONS_IDENTICAL,
    QUERY_START, REF_DIFF_POSITION_DIFF_ALN, REF_DIFF_POSITION_SAME_ALN,
    REF_SAME_POSITION_DIFF_ALN, REF_SAME_POSITION_SAME_ALN, REF_SAME_POSITION_SAME_JUNCTIONS,
    SIDES, SIDE_A, SIDE_B, STRAND, SUMMARY_HEADER, TARGET_NAME, TARGET_START,
};
use crate::cs_parser::{cs_revcomp, cs_strip_splice_motifs};
use crate::io_utils::open_output;
use crate::junction::{
    genomic_junction_set_stats, junction_set_stats, parse_genomic_junction_str, parse_junction_str,
};
use crate::table_input::{open_table, InputFormat};

/// Whether each side's representative alignment is mapped.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum MapStatus {
    BothMapped,
    OnlyAMapped,
    OnlyBMapped,
    NeitherMapped,
}

/// Reference-space classification of one both-mapped read: whether the two
/// sides share the same genomic position (`target_name` + `strand` +
/// `target_start`) and/or report the same alignment (`cs`, motif-blind). See
/// the module doc comment for why this is a strict, literal comparison with
/// no reverse-complement accommodation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum RefClass {
    /// Same position, same alignment — reference-identical.
    SamePositionSameAln,
    /// Same position, different alignment (e.g. a different indel placement
    /// at the same site).
    SamePositionDiffAln,
    /// Same alignment, different position — the alignment was relocated.
    DiffPositionSameAln,
    /// Both position and alignment differ.
    DiffPositionDiffAln,
}

impl RefClass {
    /// Whether this classification's position axis matched.
    pub(crate) fn same_position(&self) -> bool {
        matches!(self, Self::SamePositionSameAln | Self::SamePositionDiffAln)
    }

    /// Whether this classification's alignment (`cs`) axis matched.
    pub(crate) fn same_aln(&self) -> bool {
        matches!(self, Self::SamePositionSameAln | Self::DiffPositionSameAln)
    }
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
    /// Reference-space classification. `None` unless both sides are mapped.
    pub ref_class: Option<RefClass>,
    /// Query-space splice-junction *set* identity (deduplicated per side).
    /// `None` unless both sides are mapped — same convention as `ref_class`.
    pub query_junctions_identical: Option<bool>,
    /// Genomic-coordinate splice-junction *set* identity, gated on
    /// `ref_class`'s position axis (`same_position()`) — mirrors
    /// `find-aln-diff`'s `ref_same_position_same_junctions` semantics.
    /// `None` unless both sides are mapped.
    pub ref_same_position_same_junctions: Option<bool>,
}

/// A `target_name` value indicates an unmapped read when it is empty or `"*"`.
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
    let mapped_a = is_mapped(get_a(TARGET_NAME.name));
    let mapped_b = is_mapped(get_b(TARGET_NAME.name));

    let map_status = match (mapped_a, mapped_b) {
        (true, true) => MapStatus::BothMapped,
        (true, false) => MapStatus::OnlyAMapped,
        (false, true) => MapStatus::OnlyBMapped,
        (false, false) => MapStatus::NeitherMapped,
    };

    let mut query_identical = false;
    let mut query_identical_rc = false;
    let mut ref_class = None;
    let mut query_junctions_identical = None;
    let mut ref_same_position_same_junctions = None;

    if mapped_a && mapped_b {
        // Compare cs tags with intron donor/acceptor motif letters blanked out
        // (`~ct..ac` and `~nn..nn` compare equal if the intron length/position
        // match). Some aligners (e.g. STAR) report `nn` placeholders instead of
        // the true motif bases that minimap2 reports for the same splice site —
        // that's a limitation of the aligner's own output, not a real alignment
        // difference, so it must not by itself make two alignments "different".
        let cs_a = cs_strip_splice_motifs(get_a(CS.name));
        let cs_b = cs_strip_splice_motifs(get_b(CS.name));

        // Query span is in forward-read coordinates in PAF (strand-independent),
        // so it must match in both the same-strand and reverse-complement cases.
        let same_span = get_a(QUERY_START.name) == get_b(QUERY_START.name)
            && get_a(QUERY_END.name) == get_b(QUERY_END.name);
        if same_span {
            if get_a(STRAND.name) == get_b(STRAND.name) {
                if cs_a == cs_b {
                    query_identical = true;
                }
            } else if cs_a == cs_strip_splice_motifs(&cs_revcomp(get_b(CS.name))) {
                // Opposite strands but the alignment is an exact reverse-complement
                // (e.g. an inverted locus between two assemblies).
                query_identical = true;
                query_identical_rc = true;
            }
        }

        // Reference-space classification: independent of query span, and a
        // strict literal comparison (no reverse-complement accommodation — a
        // real strand difference always means a different position). A
        // matching start plus matching cs implies a matching end too, since
        // the cs tag's own operations determine the alignment's length.
        let same_position = get_a(TARGET_NAME.name) == get_b(TARGET_NAME.name)
            && get_a(STRAND.name) == get_b(STRAND.name)
            && get_a(TARGET_START.name) == get_b(TARGET_START.name);
        let same_aln = cs_a == cs_b;
        ref_class = Some(match (same_position, same_aln) {
            (true, true) => RefClass::SamePositionSameAln,
            (true, false) => RefClass::SamePositionDiffAln,
            (false, true) => RefClass::DiffPositionSameAln,
            (false, false) => RefClass::DiffPositionDiffAln,
        });

        // Query-space splice-junction *set* identity — strand-agnostic (query
        // junctions are stored in plus-strand read coordinates), so no span
        // gate and no reverse-complement handling are needed here.
        let ja = parse_junction_str(get_a(QUERY_JUNCTIONS.name));
        let jb = parse_junction_str(get_b(QUERY_JUNCTIONS.name));
        let (_, only_a, only_b) = junction_set_stats(&ja, &jb);
        query_junctions_identical = Some(only_a == 0 && only_b == 0);

        // Genomic-coordinate junction-set identity, gated on the position axis
        // above. `genomic_junctions` cells carry only `(start, end)` pairs, so
        // each side's `target_name` is reattached before comparing.
        let chrom_a = get_a(TARGET_NAME.name);
        let chrom_b = get_b(TARGET_NAME.name);
        let ga: Vec<(String, u64, u64)> = parse_genomic_junction_str(get_a(GENOMIC_JUNCTIONS.name))
            .into_iter()
            .map(|(s, e)| (chrom_a.to_string(), s, e))
            .collect();
        let gb: Vec<(String, u64, u64)> = parse_genomic_junction_str(get_b(GENOMIC_JUNCTIONS.name))
            .into_iter()
            .map(|(s, e)| (chrom_b.to_string(), s, e))
            .collect();
        let (_, only_ga, only_gb) = genomic_junction_set_stats(&ga, &gb);
        ref_same_position_same_junctions = Some(same_position && only_ga == 0 && only_gb == 0);
    } else if map_status == MapStatus::NeitherMapped {
        // Both aligners agreeing a read doesn't map at all is itself a form
        // of agreement — count it as query_identical (same-strand branch, not
        // reverse-complement). `ref_class`/junction-identity stay `None`:
        // those remain a "both mapped" concept, untouched by this branch.
        query_identical = true;
    }

    ReadClass {
        map_status,
        query_identical,
        query_identical_rc,
        query_junctions_identical,
        ref_same_position_same_junctions,
        ref_class,
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
    pub ref_same_position_same_aln: u64,
    pub ref_same_position_diff_aln: u64,
    pub ref_diff_position_same_aln: u64,
    pub ref_diff_position_diff_aln: u64,
    pub query_junctions_identical: u64,
    pub ref_same_position_same_junctions: u64,
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
            // same_strand/rc are a strictly "both mapped" sub-split — a
            // NeitherMapped read is query_identical too (see `classify()`),
            // but must not inflate either sub-bucket, since
            // `query_not_identical()` depends on them staying an exact
            // partition of `aligned_both`.
            if c.map_status == MapStatus::BothMapped {
                if c.query_identical_rc {
                    self.query_identical_rc += 1;
                } else {
                    self.query_identical_same_strand += 1;
                }
            }
        }
        match c.ref_class {
            Some(RefClass::SamePositionSameAln) => self.ref_same_position_same_aln += 1,
            Some(RefClass::SamePositionDiffAln) => self.ref_same_position_diff_aln += 1,
            Some(RefClass::DiffPositionSameAln) => self.ref_diff_position_same_aln += 1,
            Some(RefClass::DiffPositionDiffAln) => self.ref_diff_position_diff_aln += 1,
            None => {}
        }
        if let Some(true) = c.query_junctions_identical {
            self.query_junctions_identical += 1;
        }
        if let Some(true) = c.ref_same_position_same_junctions {
            self.ref_same_position_same_junctions += 1;
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

    /// query-different among both-mapped reads. Computed from the two
    /// "both mapped" identity sub-counters (not `self.query_identical`
    /// directly), since `query_identical` also now includes `aligned_neither`
    /// reads — subtracting the grown total from `aligned_both` would
    /// underflow. This yields the exact same value as before for any dataset:
    /// `query_identical_same_strand` and `query_identical_rc` are still only
    /// ever incremented for `BothMapped` reads.
    fn query_not_identical(&self) -> u64 {
        self.aligned_both - self.query_identical_same_strand - self.query_identical_rc
    }

    /// query-junctions-different among both-mapped reads.
    fn query_junctions_not_identical(&self) -> u64 {
        self.aligned_both - self.query_junctions_identical
    }

    /// Same-position (reference) reads whose junction sets differ. Gated on
    /// `same_position` regardless of `same_aln` — matching how
    /// `ref_same_position_same_junctions` itself is gated — not on
    /// `aligned_both`, since a different-position read has no meaningful
    /// "junctions at the same position" comparison at all.
    fn ref_same_position_diff_junctions(&self) -> u64 {
        (self.ref_same_position_same_aln + self.ref_same_position_diff_aln)
            - self.ref_same_position_same_junctions
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
            (QUERY_IDENTICAL_SAME_STRAND.to_string(), self.query_identical_same_strand),
            (QUERY_IDENTICAL_REVCOMP.to_string(), self.query_identical_rc),
            ("query_not_identical".to_string(), self.query_not_identical()),
            (QUERY_JUNCTIONS_IDENTICAL.to_string(), self.query_junctions_identical),
            ("query_junctions_not_identical".to_string(), self.query_junctions_not_identical()),
            (REF_SAME_POSITION_SAME_ALN.to_string(), self.ref_same_position_same_aln),
            (REF_SAME_POSITION_DIFF_ALN.to_string(), self.ref_same_position_diff_aln),
            (REF_DIFF_POSITION_SAME_ALN.to_string(), self.ref_diff_position_same_aln),
            (REF_DIFF_POSITION_DIFF_ALN.to_string(), self.ref_diff_position_diff_aln),
            (REF_SAME_POSITION_SAME_JUNCTIONS.to_string(), self.ref_same_position_same_junctions),
            ("ref_same_position_diff_junctions".to_string(), self.ref_same_position_diff_junctions()),
            ("present_only_in_A_by_id".to_string(), self.a_only_by_id),
            ("present_only_in_B_by_id".to_string(), self.b_only_by_id),
        ]
    }

    /// Write the summary as a 2-column TSV (`category<TAB>count`), preceded by
    /// the two label provenance rows and any caller-supplied `extra`
    /// provenance rows (e.g. `find-aln-diff`'s `space`/`compare_by`), so the
    /// file is self-describing. This is the single column layout shared by
    /// `compare`, `toolkit summary`, and `find-aln-diff` — pass `&[]`
    /// for no extra rows.
    pub fn write_tsv(
        &self,
        path: &str,
        label_a: &str,
        label_b: &str,
        extra: &[(&str, &str)],
    ) -> Result<()> {
        let mut w = open_output(Some(path))?;
        writeln!(w, "{}", SUMMARY_HEADER.join("\t"))?;
        writeln!(w, "label_A\t{label_a}")?;
        writeln!(w, "label_B\t{label_b}")?;
        for (k, v) in extra {
            writeln!(w, "{k}\t{v}")?;
        }
        for (k, v) in self.rows() {
            writeln!(w, "{k}\t{v}")?;
        }
        w.flush()?;
        Ok(())
    }

    /// Print a human-readable block to stderr, in the same layout as
    /// `write_tsv` (see its doc comment).
    pub fn render_stderr(&self, label_a: &str, label_b: &str, extra: &[(&str, &str)]) {
        eprintln!("Comparison summary:");
        eprintln!("  {:<34} {label_a}", "label_A");
        eprintln!("  {:<34} {label_b}", "label_B");
        for (k, v) in extra {
            eprintln!("  {k:<34} {v}");
        }
        for (k, v) in self.rows() {
            eprintln!("  {k:<34} {v}");
        }
    }

    /// Print the abbreviated stderr block `compare` uses: just the alignment-
    /// status breakdown and the two headline query-identity stats. The TSV
    /// (`write_tsv`) and every other caller of `render_stderr` keep the full
    /// category set — this is only for `compare`'s own terminal output.
    pub fn render_stderr_brief(&self, label_a: &str, label_b: &str) {
        eprintln!("Comparison summary:");
        eprintln!("  {:<34} {label_a}", "label_A");
        eprintln!("  {:<34} {label_b}", "label_B");
        eprintln!("  {:<34} {}", "reads_compared", self.reads_compared);
        eprintln!("  {:<34} {}", "aligned_both", self.aligned_both);
        eprintln!("  {:<34} {}", "aligned_only_A", self.aligned_only_a);
        eprintln!("  {:<34} {}", "aligned_only_B", self.aligned_only_b);
        eprintln!("  {:<34} {}", "aligned_neither", self.aligned_neither);
        eprintln!("  {:<34} {}", "query_identical", self.query_identical);
        eprintln!("  {:<34} {}", "query_not_identical", self.query_not_identical());
    }
}

// ── `toolkit summary` command ──────────────────────────────────────────

#[derive(clap::Args, Debug)]
pub struct CompareSummaryArgs {
    /// Alignment comparison table (.parquet or .tsv[.gz]; `-` reads TSV from stdin).
    /// See `--input-format`.
    #[arg(short = 'i', long = "input", value_name = "compare.parquet|compare.tsv[.gz]")]
    input: String,

    /// Input serialization. `auto` (default) selects Parquet for a
    /// `.parquet`-named `--input` and TSV otherwise. Parquet requires a real
    /// file path — it cannot be read from stdin (`-`).
    #[arg(long = "input-format", value_enum, default_value = "auto")]
    input_format: InputFormat,

    /// Write the summary TSV here (default: stderr only). `.gz` ok; `-` for stdout.
    #[arg(short = 'o', long = "output", value_name = "summary.tsv[.gz]")]
    output: Option<String>,
}

/// Confirm a compare-table header uses the fixed `_A` / `_B` side suffixes.
/// Shared by `toolkit summary`, `find-aln-diff` and `query-junction-diff`.
pub(crate) fn require_ab_schema(cols: &[&str]) -> Result<()> {
    let (a, b) = (sided(TARGET_NAME, SIDE_A), sided(TARGET_NAME, SIDE_B));
    if cols.contains(&a.as_str()) && cols.contains(&b.as_str()) {
        return Ok(());
    }
    bail!(
        "comparison table is missing the `{a}` / `{b}` columns \
         — is this a maligno compare table?"
    );
}

/// Index of a (non-sided) column in a comparison table's header, or the
/// standard missing-column error.
pub(crate) fn require_col(col_index: &HashMap<&str, usize>, name: &str) -> Result<usize> {
    col_index
        .get(name)
        .copied()
        .with_context(|| format!("comparison table is missing column '{name}'"))
}

/// The columns to project from a comparison table: `keys` as named, then
/// `side_cols` suffixed for side A and for side B.
pub(crate) fn wanted_columns(keys: &[Col], side_cols: &[Col]) -> Vec<String> {
    keys.iter()
        .map(|c| c.name.to_string())
        .chain(SIDES.iter().flat_map(|&side| side_cols.iter().map(move |&c| sided(c, side))))
        .collect()
}

/// One side's column indices in a comparison table, keyed by **unsuffixed**
/// column name, so per-side code (`classify`, the interval and junction
/// builders) reads either side through the same names.
pub(crate) struct SideIndex(HashMap<&'static str, usize>);

impl SideIndex {
    /// Resolve every column in `cols` for `side`. A missing column is an error
    /// here, before any row is read.
    pub(crate) fn resolve(col_index: &HashMap<&str, usize>, cols: &[Col], side: &str) -> Result<Self> {
        let mut m = HashMap::new();
        for &c in cols {
            m.insert(c.name, require_col(col_index, &sided(c, side))?);
        }
        Ok(Self(m))
    }

    /// The cell for unsuffixed column `name` in `fields` (`""` if the row is
    /// short). Panics if `name` was not resolved: that is a bug in the caller's
    /// column list, not bad input, and must not read silently as an empty cell.
    pub(crate) fn get<'f>(&self, fields: &[&'f str], name: &str) -> &'f str {
        let i = *self
            .0
            .get(name)
            .unwrap_or_else(|| panic!("column {name:?} was read but never resolved"));
        fields.get(i).copied().unwrap_or("")
    }
}

/// Recover the human-readable set labels from a compare table's first data row
/// (the `label_A` / `label_B` columns). Falls back to `("A", "B")` when those
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
    (pick(LABEL_A.name, SIDE_A), pick(LABEL_B.name, SIDE_B))
}

pub fn run(args: &CompareSummaryArgs) -> Result<()> {
    // For a Parquet input, project down to just the 20 columns this command
    // reads (of the comparison table's 99) — Parquet skips decoding the rest,
    // which is where its per-column storage actually pays off.
    let wanted = wanted_columns(&[LABEL_A, LABEL_B], &CLASSIFY_SIDE_COLS);
    let wanted_refs: Vec<&str> = wanted.iter().map(String::as_str).collect();
    let mut source = open_table(&args.input, args.input_format, Some(&wanted_refs))?;

    // Header → column index map.
    let cols_owned = source
        .header()
        .with_context(|| format!("reading comparison table '{}'", args.input))?;
    let cols: Vec<&str> = cols_owned.iter().map(String::as_str).collect();
    let col_index: HashMap<&str, usize> =
        cols.iter().copied().enumerate().map(|(i, c)| (c, i)).collect();

    // Require the fixed `_A` / `_B` side suffixes.
    require_ab_schema(&cols)?;

    // The unsuffixed columns `classify` reads, resolved per side up front.
    let idx_a = SideIndex::resolve(&col_index, &CLASSIFY_SIDE_COLS, SIDE_A)?;
    let idx_b = SideIndex::resolve(&col_index, &CLASSIFY_SIDE_COLS, SIDE_B)?;

    // The human-readable labels live in each row's `label_A` / `label_B`, so they
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
        let get_a = |c: &str| idx_a.get(&fields, c);
        let get_b = |c: &str| idx_b.get(&fields, c);
        summary.observe(&classify(&get_a, &get_b));
    }

    if let Some(path) = &args.output {
        summary.write_tsv(path, &label_a, &label_b, &[])?;
        eprintln!("Wrote summary: {path}");
    }
    summary.render_stderr(&label_a, &label_b, &[]);
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
        let cols = ["query_name", "label_A", "label_B", "target_name_A", "target_name_B"];
        assert!(require_ab_schema(&cols).is_ok());
    }

    #[test]
    fn one_sided_ab_schema_is_rejected() {
        // `target_name_B` missing → not a usable comparison table.
        let cols = ["query_name", "target_name_A"];
        assert!(require_ab_schema(&cols).is_err());
    }

    #[test]
    fn non_compare_table_is_rejected() {
        // No `target_name_*` columns at all (e.g. an alninfo/readinfo table).
        let cols = ["query_name", "query_len", "target_name"];
        let err = require_ab_schema(&cols).unwrap_err().to_string();
        assert!(err.contains("missing the `target_name_A`"), "unexpected error: {err}");
    }

    // ─── labels_from_row ───────────────────────────────────────────────────

    #[test]
    fn labels_read_from_first_data_row() {
        let cols = ["query_name", "query_len", "label_A", "label_B", "target_name_A"];
        let fields = ["r1", "100", "Splice", "SpliceHQ", "chr1"];
        assert_eq!(
            labels_from_row(&index(&cols), &fields),
            ("Splice".to_string(), "SpliceHQ".to_string())
        );
    }

    #[test]
    fn labels_fall_back_when_columns_absent() {
        // No label_A/label_B columns → report the bare side identifiers.
        let cols = ["query_name", "target_name_A", "target_name_B"];
        let fields = ["r1", "chr1", "chr1"];
        assert_eq!(
            labels_from_row(&index(&cols), &fields),
            ("A".to_string(), "B".to_string())
        );
    }

    #[test]
    fn labels_fall_back_on_empty_values_and_short_rows() {
        let cols = ["query_name", "label_A", "label_B"];
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
        // Side identity comes from the fixed `_A`/`_B` suffixes, not the label,
        // so underscores in a label are fine and must be accepted.
        assert_eq!(validate_set_label("my_run_1").unwrap(), "my_run_1");
        assert_eq!(validate_set_label("Splice").unwrap(), "Splice");
    }

    // ─── classify(): junction-identity fields ──────────────────────────────
    // Tests for `classify()`'s junction-identity fields, driven through its full
    // accessor interface.

    /// Build a `get_a`/`get_b`-shaped accessor from a fixed field list. Both
    /// sides default to "mapped, same reference position, same cs" so a test
    /// can vary only the field(s) it cares about.
    fn getter(fields: &'static [(&'static str, &'static str)]) -> impl Fn(&str) -> &'static str {
        move |c: &str| fields.iter().find(|(k, _)| *k == c).map(|(_, v)| *v).unwrap_or("")
    }

    const BOTH_MAPPED_SAME_POSITION: &[(&str, &str)] = &[
        ("target_name", "chr1"),
        ("strand", "+"),
        ("cs", ":10"),
        ("query_start", "0"),
        ("query_end", "10"),
        ("target_start", "100"),
    ];

    /// `extra` entries take precedence over `BOTH_MAPPED_SAME_POSITION`'s
    /// defaults (the getter returns the first match for a given key).
    fn mapped_with(extra: &'static [(&'static str, &'static str)]) -> Vec<(&'static str, &'static str)> {
        let mut v = extra.to_vec();
        v.extend_from_slice(BOTH_MAPPED_SAME_POSITION);
        v
    }

    #[test]
    fn query_junctions_identical_empty_equals_empty() {
        let a = mapped_with(&[("query_junctions", "()")]);
        let b = mapped_with(&[("query_junctions", "()")]);
        let class = classify(&getter(Box::leak(a.into_boxed_slice())), &getter(Box::leak(b.into_boxed_slice())));
        assert_eq!(class.query_junctions_identical, Some(true));
    }

    #[test]
    fn query_junctions_identical_same_set() {
        let a = mapped_with(&[("query_junctions", "(10, 20)")]);
        let b = mapped_with(&[("query_junctions", "(10, 20)")]);
        let class = classify(&getter(Box::leak(a.into_boxed_slice())), &getter(Box::leak(b.into_boxed_slice())));
        assert_eq!(class.query_junctions_identical, Some(true));
    }

    #[test]
    fn query_junctions_identical_ignores_order() {
        let a = mapped_with(&[("query_junctions", "(10, 20)")]);
        let b = mapped_with(&[("query_junctions", "(20, 10)")]);
        let class = classify(&getter(Box::leak(a.into_boxed_slice())), &getter(Box::leak(b.into_boxed_slice())));
        assert_eq!(class.query_junctions_identical, Some(true));
    }

    #[test]
    fn query_junctions_identical_disjoint_is_different() {
        let a = mapped_with(&[("query_junctions", "(10,)")]);
        let b = mapped_with(&[("query_junctions", "(20,)")]);
        let class = classify(&getter(Box::leak(a.into_boxed_slice())), &getter(Box::leak(b.into_boxed_slice())));
        assert_eq!(class.query_junctions_identical, Some(false));
    }

    #[test]
    fn query_junctions_identical_one_side_empty_is_different() {
        let a = mapped_with(&[("query_junctions", "(10,)")]);
        let b = mapped_with(&[("query_junctions", "()")]);
        let class = classify(&getter(Box::leak(a.into_boxed_slice())), &getter(Box::leak(b.into_boxed_slice())));
        assert_eq!(class.query_junctions_identical, Some(false));
    }

    #[test]
    fn query_junctions_identical_subset_is_different() {
        let a = mapped_with(&[("query_junctions", "(10, 20)")]);
        let b = mapped_with(&[("query_junctions", "(10,)")]);
        let class = classify(&getter(Box::leak(a.into_boxed_slice())), &getter(Box::leak(b.into_boxed_slice())));
        assert_eq!(class.query_junctions_identical, Some(false));
    }

    #[test]
    fn ref_same_position_same_junctions_true_when_position_and_set_match() {
        let a = mapped_with(&[("genomic_junctions", "((100, 250),)")]);
        let b = mapped_with(&[("genomic_junctions", "((100, 250),)")]);
        let class = classify(&getter(Box::leak(a.into_boxed_slice())), &getter(Box::leak(b.into_boxed_slice())));
        assert_eq!(class.ref_same_position_same_junctions, Some(true));
    }

    #[test]
    fn ref_same_position_same_junctions_false_when_set_differs() {
        let a = mapped_with(&[("genomic_junctions", "((100, 250),)")]);
        let b = mapped_with(&[("genomic_junctions", "((400, 800),)")]);
        let class = classify(&getter(Box::leak(a.into_boxed_slice())), &getter(Box::leak(b.into_boxed_slice())));
        assert_eq!(class.ref_same_position_same_junctions, Some(false));
    }

    #[test]
    fn ref_same_position_same_junctions_false_when_position_differs() {
        // Same genomic-junction set, but a different target_start (position).
        let a = mapped_with(&[("genomic_junctions", "((100, 250),)"), ("target_start", "100")]);
        let b = mapped_with(&[("genomic_junctions", "((100, 250),)"), ("target_start", "200")]);
        let class = classify(&getter(Box::leak(a.into_boxed_slice())), &getter(Box::leak(b.into_boxed_slice())));
        assert_eq!(class.ref_same_position_same_junctions, Some(false));
    }

    #[test]
    fn junction_identity_fields_are_none_unless_both_mapped() {
        let a = mapped_with(&[]);
        let unmapped: Vec<(&str, &str)> = vec![("target_name", "*")];
        let class = classify(&getter(Box::leak(a.into_boxed_slice())), &getter(Box::leak(unmapped.into_boxed_slice())));
        assert_eq!(class.query_junctions_identical, None);
        assert_eq!(class.ref_same_position_same_junctions, None);
    }

    #[test]
    fn neither_mapped_is_query_identical_but_not_via_revcomp() {
        // Both aligners agreeing a read is unmapped is itself agreement.
        let unmapped_a: Vec<(&str, &str)> = vec![("target_name", "*")];
        let unmapped_b: Vec<(&str, &str)> = vec![("target_name", "")];
        let class = classify(
            &getter(Box::leak(unmapped_a.into_boxed_slice())),
            &getter(Box::leak(unmapped_b.into_boxed_slice())),
        );
        assert_eq!(class.map_status, MapStatus::NeitherMapped);
        assert!(class.query_identical);
        assert!(!class.query_identical_rc);
        // Reference-space and junction-identity axes stay a strictly
        // "both mapped" concept — unaffected by this branch.
        assert_eq!(class.ref_class, None);
        assert_eq!(class.query_junctions_identical, None);
        assert_eq!(class.ref_same_position_same_junctions, None);
    }

    #[test]
    fn one_sided_mapped_reads_are_not_query_identical() {
        // Sanity check the NeitherMapped branch didn't leak into
        // OnlyAMapped/OnlyBMapped — those must stay `false`, unchanged.
        let mapped = mapped_with(&[]);
        let unmapped: Vec<(&str, &str)> = vec![("target_name", "*")];
        let class = classify(
            &getter(Box::leak(mapped.into_boxed_slice())),
            &getter(Box::leak(unmapped.into_boxed_slice())),
        );
        assert_eq!(class.map_status, MapStatus::OnlyAMapped);
        assert!(!class.query_identical);
    }

    // ─── CompareSummary: new counters + rows() ─────────────────────────────

    #[test]
    fn aligned_neither_counts_as_query_identical_without_underflow() {
        let mut s = CompareSummary::default();
        // One both-mapped, same-strand-identical read.
        s.observe(&ReadClass {
            map_status: MapStatus::BothMapped,
            query_identical: true,
            query_identical_rc: false,
            ref_class: Some(RefClass::SamePositionSameAln),
            query_junctions_identical: Some(true),
            ref_same_position_same_junctions: Some(true),
        });
        // One both-mapped read whose alignment genuinely differs.
        s.observe(&ReadClass {
            map_status: MapStatus::BothMapped,
            query_identical: false,
            query_identical_rc: false,
            ref_class: Some(RefClass::SamePositionDiffAln),
            query_junctions_identical: Some(false),
            ref_same_position_same_junctions: Some(false),
        });
        // Two neither-mapped reads — query_identical per classify()'s new
        // branch, but not via the reverse-complement sub-bucket.
        for _ in 0..2 {
            s.observe(&ReadClass {
                map_status: MapStatus::NeitherMapped,
                query_identical: true,
                query_identical_rc: false,
                ref_class: None,
                query_junctions_identical: None,
                ref_same_position_same_junctions: None,
            });
        }

        let rows: HashMap<String, u64> = s.rows().into_iter().collect();
        assert_eq!(rows["aligned_both"], 2);
        assert_eq!(rows["aligned_neither"], 2);
        // query_identical now includes the 2 neither-mapped reads.
        assert_eq!(rows["query_identical"], 3);
        // ...but query_not_identical must still reflect only the genuine
        // both-mapped mismatch (1), not underflow/panic from the subtraction.
        assert_eq!(rows["query_not_identical"], 1);
    }

    #[test]
    fn compare_summary_rows_include_junction_counters() {
        let mut s = CompareSummary::default();
        // Read 1: both mapped, query-junctions identical, ref same-position-same-junctions.
        s.observe(&ReadClass {
            map_status: MapStatus::BothMapped,
            query_identical: true,
            query_identical_rc: false,
            ref_class: Some(RefClass::SamePositionSameAln),
            query_junctions_identical: Some(true),
            ref_same_position_same_junctions: Some(true),
        });
        // Read 2: both mapped, junctions differ on both axes.
        s.observe(&ReadClass {
            map_status: MapStatus::BothMapped,
            query_identical: false,
            query_identical_rc: false,
            ref_class: Some(RefClass::SamePositionDiffAln),
            query_junctions_identical: Some(false),
            ref_same_position_same_junctions: Some(false),
        });

        let rows: HashMap<String, u64> = s.rows().into_iter().collect();
        assert_eq!(rows["aligned_both"], 2);
        assert_eq!(rows["query_junctions_identical"], 1);
        assert_eq!(rows["query_junctions_not_identical"], 1);
        assert_eq!(rows["ref_same_position_same_junctions"], 1);
        // ref_same_position_diff_junctions denominator is same_position (both
        // reads here), not aligned_both — both are same_position, one same_junctions.
        assert_eq!(rows["ref_same_position_diff_junctions"], 1);
    }

    #[test]
    fn write_tsv_extra_rows_appear_between_labels_and_counters() {
        let s = CompareSummary::default();
        let path = std::env::temp_dir().join("maligno_test_write_tsv_extra.tsv");
        let path_str = path.to_str().unwrap();
        s.write_tsv(path_str, "Splice", "SpliceHQ", &[("space", "query"), ("compare_by", "all")])
            .unwrap();
        let content = std::fs::read_to_string(&path).unwrap();
        let lines: Vec<&str> = content.lines().collect();
        let _ = std::fs::remove_file(&path);
        assert_eq!(lines[0], "category\tcount");
        assert_eq!(lines[1], "label_A\tSplice");
        assert_eq!(lines[2], "label_B\tSpliceHQ");
        assert_eq!(lines[3], "space\tquery");
        assert_eq!(lines[4], "compare_by\tall");
        assert_eq!(lines[5], "reads_compared\t0");
    }
}
