//! `find-aln-diff` — from a `compare` / `compare-toolkit merge-readinfo`
//! table, find every read whose alignment differs between sets A and B, in
//! query space or reference space, and the genomic regions where those reads
//! cluster.
//!
//! `--space` selects the coordinate space that defines a "difference":
//!   - `query` (default) — compares each side's alignment **relative to the
//!     read** (via the shared `classify`): a read is different iff it is in
//!     exactly one of `diff_aln_to_both` (mapped both sides, NOT
//!     query-identical), `diff_aln_only_A`, `diff_aln_only_B`.
//!     Reverse-complement matches count as **identical** (not a difference).
//!   - `reference` — compares each side's placement **on the reference
//!     genome**: same `TargetChr`/`Strand`/`Target_Start` and same alignment
//!     (`cs`, or the genomic-coordinate junction set under `--compare-by
//!     junctions`). Strict and literal — no reverse-complement accommodation,
//!     so a real strand difference always means a different position. A
//!     difference is `reference_diff` (mapped both sides, not
//!     reference-identical), `diff_aln_only_A`, or `diff_aln_only_B`.
//! In both spaces, reads unmapped on both sides are excluded.
//!
//! `--compare-by` selects what "identical" means for both-mapped reads, within
//! whichever space is active:
//!   - `all` (default) — the full cs tag must match, motif-blind (intron
//!     donor/acceptor motif letters are ignored, so a differently-reported
//!     motif at the same intron position/length is not by itself a
//!     difference); any other mismatch/indel/soft-clip/junction-position
//!     difference counts.
//!   - `junctions` — only the splice-junction set must match (query-space
//!     under `--space query`, genomic-coordinate under `--space reference`);
//!     reads with identical junctions but differing mismatches/indels/
//!     soft-clips count as the same. Reads aligned in only one set are still
//!     reported as differences (they have no junctions to compare on the
//!     missing side).
//!
//! Outputs (to `--outdir`, `--prefix`-named; `query_diff`/`query_identical` or
//! `reference_diff`/`reference_identical` depending on `--space`; gzipped by
//! default, `--no-gzip` to opt out):
//!   1. `{prefix}.{stem}_reads.tsv[.gz]`      — one row per differing read:
//!      `Read_Name`, `outcome`, plus 8 classification booleans (`1`/`0`)
//!      computed the same way regardless of `--space`/`--compare-by` —
//!      `query_identical_same_strand`, `query_identical_revcomp`,
//!      `query_junctions_identical`, `ref_same_position_same_aln`,
//!      `ref_same_position_diff_aln`, `ref_diff_position_same_aln`,
//!      `ref_diff_position_diff_aln`, `ref_same_position_same_junctions`.
//!      `{prefix}.{identical_stem}_reads.tsv[.gz]` (`--emit-identical-reads`)
//!      carries the same 8 columns for the complementary read set.
//!   2. `{prefix}.{stem}_regions.A.bed[.gz]`  — merged A-coordinate loci
//!   3. `{prefix}.{stem}_regions.B.bed[.gz]`  — merged B-coordinate loci
//!   4. `{prefix}.{stem}_summary.tsv`         — the same category-tally schema
//!      `compare`/`compare-toolkit summary` write (`CompareSummary::rows()`),
//!      with `space`/`compare_by` provenance rows prepended, so all three
//!      commands share one column layout. Its counters are mode-independent
//!      (always the raw `classify()` output) — mode-specific behavior only
//!      selects which reads land in `{stem}_reads`/region-BED outputs.
//!
//! The per-row core of this command (bool-cols + differing-read selection +
//! genomic interval building) lives in `AlnDiffAccumulator` below so that
//! `compare` can drive the identical logic inline, in its own single merge
//! pass, for its default (`--space query --compare-by all`) fused output.

use std::fs;
use std::io::Write;
use std::path::Path;

use anyhow::{Context, Result};

use crate::compare_summary::{classify, labels_from_row, require_ab_schema, CompareSummary, MapStatus, ReadClass, RefClass};
use crate::interval_merge::{merge_and_count, Ivl, Locus};
use crate::io_utils::open_output;
use crate::table_input::{open_table, InputFormat};
use std::collections::HashMap;

/// What aspect of the alignment defines a difference between A and B, within
/// whichever `--space` is active.
#[derive(Clone, Copy, Debug, PartialEq, Eq, clap::ValueEnum)]
pub enum CompareBy {
    /// Flag any difference across the whole alignment (compares the full cs tag,
    /// motif-blind: intron donor/acceptor letters are ignored).
    All,
    /// Flag only reads whose splice-junction set differs — the query-space set
    /// under `--space query`, the genomic-coordinate set under `--space reference`.
    Junctions,
}

/// Which coordinate space defines a "difference" between A and B.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, clap::ValueEnum)]
pub enum DiffSpace {
    /// Compare each side's alignment relative to the read.
    #[default]
    Query,
    /// Compare each side's placement on the reference genome: same
    /// `TargetChr`/`Strand`/`Target_Start` and same alignment. Strict and
    /// literal — no reverse-complement accommodation, so a real strand
    /// difference always means a different position.
    Reference,
}

/// Find reads whose alignment differs between A and B, in query or reference
/// space, and the genomic regions where they cluster.
#[derive(clap::Args, Debug)]
pub struct FindAlnDiffArgs {
    /// Comparison table from `compare` / `compare-toolkit merge-readinfo`:
    /// TSV (`.gz` ok; `-` = stdin) or Parquet (`--format parquet` / `-o
    /// x.parquet`). See `--input-format`.
    #[arg(short = 'i', long = "input", value_name = "compare.tsv|compare.parquet")]
    input: String,

    /// Input serialization. `auto` (default) selects Parquet for a
    /// `.parquet`-named `--input` and TSV otherwise. Parquet requires a real
    /// file path — it cannot be read from stdin (`-`).
    #[arg(long = "input-format", value_enum, default_value = "auto")]
    input_format: InputFormat,

    /// Which coordinate space defines a "difference". `query` (default)
    /// compares each side's alignment relative to the read; `reference`
    /// compares each side's placement on the reference genome (strict,
    /// literal — no reverse-complement accommodation).
    #[arg(long = "space", value_enum, default_value = "query")]
    space: DiffSpace,

    /// Output directory (created if it does not exist).
    #[arg(long = "outdir", value_name = "DIR")]
    outdir: String,

    /// Filename prefix for all outputs.
    #[arg(long = "prefix", value_name = "STR")]
    prefix: String,

    /// Do not gzip the read and region output tables (gzipped by default).
    #[arg(long = "no-gzip")]
    no_gzip: bool,

    /// What defines a difference: `all` (default) compares the full cs tag,
    /// motif-blind (intron donor/acceptor letters ignored) — any other
    /// mismatch/indel/soft-clip/junction-position difference counts; `junctions`
    /// compares only the splice-junction set (query-space or genomic-coordinate,
    /// depending on `--space`; reads with identical junctions but differing
    /// mismatches/indels/soft-clips count as the same). In `junctions` mode the
    /// outputs gain a `.junctions` filename segment.
    #[arg(long = "compare-by", value_enum, default_value = "all")]
    compare_by: CompareBy,

    /// Also write a TSV of Read_Names considered identical between A and B
    /// under the active `--space`. Off by default.
    #[arg(long = "emit-identical-reads")]
    emit_identical_reads: bool,
}

/// Whether a differing read is mapped on both sides or only this coordinate
/// space's side (rendered as `n_only_A` / `n_only_B` depending on the table).
#[derive(Clone, Copy)]
enum Outcome {
    Both,
    OnlySide,
}

/// Per-interval payload folded into each merged locus.
struct DiffMeta {
    strand: char,
    outcome: Outcome,
}

/// Per-locus accumulator.
#[derive(Default)]
struct DiffAcc {
    n_reads: u64,
    n_both: u64,
    n_only: u64,
    n_plus: u64,
    n_minus: u64,
}

fn fold_diff(a: &mut DiffAcc, m: &DiffMeta) {
    a.n_reads += 1;
    match m.outcome {
        Outcome::Both => a.n_both += 1,
        Outcome::OnlySide => a.n_only += 1,
    }
    match m.strand {
        '+' => a.n_plus += 1,
        '-' => a.n_minus += 1,
        _ => {}
    }
}

/// Build a genomic interval from one side's accessor. Returns `None` (caller
/// counts it as a bad interval) if the coord is unmapped, unparseable, or `end <= start`.
fn build_ivl<'a>(get: impl Fn(&str) -> &'a str, outcome: Outcome) -> Option<Ivl<DiffMeta>> {
    let chrom = get("TargetChr");
    if chrom.is_empty() || chrom == "*" {
        return None;
    }
    let start: u64 = get("Target_Start").parse().ok()?;
    let end: u64 = get("Target_End").parse().ok()?;
    if end <= start {
        return None;
    }
    let strand = get("Strand").chars().next().unwrap_or('.');
    Some(Ivl {
        chrom: chrom.to_string(),
        start,
        end,
        meta: DiffMeta { strand, outcome },
    })
}

/// The 8 classification-boolean column names, computed the same way
/// regardless of `--space`/`--compare-by` (see the module doc comment).
const BOOL_COLS: &str = "query_identical_same_strand\tquery_identical_revcomp\t\
    query_junctions_identical\tref_same_position_same_aln\tref_same_position_diff_aln\t\
    ref_diff_position_same_aln\tref_diff_position_diff_aln\tref_same_position_same_junctions";

/// Result of `AlnDiffAccumulator::finish`.
pub(crate) struct AlnDiffStats {
    pub n_diff_rows: u64,
    pub n_bad_interval: u64,
}

/// The per-row core of `find-aln-diff`: given one matched read's raw
/// `classify()` output (`base`), decides whether it's a "difference" under
/// this accumulator's fixed `--space`/`--compare-by`, writes it to the reads
/// table if so, and buffers its genomic interval(s) for the region-BED pass.
///
/// Shared by the standalone `find-aln-diff` command (`run`, below, built with
/// its CLI's `--space`/`--compare-by`) and by `compare`'s fused default output
/// (built fixed at `DiffSpace::Query`/`CompareBy::All` — the passthrough case
/// in `observe_row`'s mode-match — so `compare` can drive it inline, in its
/// own single merge pass, with no re-read of its own output table).
pub(crate) struct AlnDiffAccumulator {
    reads_w: Box<dyn Write>,
    vec_a: Vec<Ivl<DiffMeta>>,
    vec_b: Vec<Ivl<DiffMeta>>,
    space: DiffSpace,
    compare_by: CompareBy,
    n_diff_rows: u64,
    n_bad_interval: u64,
}

impl AlnDiffAccumulator {
    pub(crate) fn new(reads_out: &str, space: DiffSpace, compare_by: CompareBy) -> Result<Self> {
        let mut reads_w = open_output(Some(reads_out))?;
        writeln!(reads_w, "Read_Name\toutcome\t{BOOL_COLS}")?;
        Ok(Self {
            reads_w,
            vec_a: Vec::new(),
            vec_b: Vec::new(),
            space,
            compare_by,
            n_diff_rows: 0,
            n_bad_interval: 0,
        })
    }

    /// Classify one matched read, write it to the reads table if it's a
    /// difference under this accumulator's active mode, and buffer its
    /// interval(s). Returns the *effective* (mode-derived) classification —
    /// `query_identical`/`query_identical_rc` reflect the active
    /// `--space`/`--compare-by`, everything else is passed through from
    /// `base` — for the caller's own use (e.g. `--emit-identical-reads`).
    pub(crate) fn observe_row<'a>(
        &mut self,
        read_name: &str,
        get_a: impl Fn(&str) -> &'a str,
        get_b: impl Fn(&str) -> &'a str,
        base: &ReadClass,
    ) -> Result<ReadClass> {
        // The 8 classification booleans, computed unconditionally (independent
        // of `--space`/`--compare-by`) so every emitted row carries the full
        // picture regardless of which mode selected it.
        let bool_cols = format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            (base.query_identical && !base.query_identical_rc) as u8,
            (base.query_identical && base.query_identical_rc) as u8,
            base.query_junctions_identical.unwrap_or(false) as u8,
            matches!(base.ref_class, Some(RefClass::SamePositionSameAln)) as u8,
            matches!(base.ref_class, Some(RefClass::SamePositionDiffAln)) as u8,
            matches!(base.ref_class, Some(RefClass::DiffPositionSameAln)) as u8,
            matches!(base.ref_class, Some(RefClass::DiffPositionDiffAln)) as u8,
            base.ref_same_position_same_junctions.unwrap_or(false) as u8,
        );

        let (identical, identical_rc) = match (self.space, self.compare_by) {
            (DiffSpace::Query, CompareBy::All) => (base.query_identical, base.query_identical_rc),
            (DiffSpace::Query, CompareBy::Junctions) => {
                (base.query_junctions_identical.unwrap_or(false), false)
            }
            (DiffSpace::Reference, compare_by) => {
                let same_position = base.ref_class.map(|rc| rc.same_position()).unwrap_or(false);
                let same_aln = match compare_by {
                    CompareBy::All => base.ref_class.map(|rc| rc.same_aln()).unwrap_or(false),
                    CompareBy::Junctions => base.ref_same_position_same_junctions.unwrap_or(false),
                };
                (same_position && same_aln, false)
            }
        };
        let class = ReadClass {
            map_status: base.map_status,
            query_identical: identical,
            query_identical_rc: identical_rc,
            ref_class: base.ref_class,
            query_junctions_identical: base.query_junctions_identical,
            ref_same_position_same_junctions: base.ref_same_position_same_junctions,
        };

        // Select differing reads and their placement(s).
        let both_mapped_diff_category = match self.space {
            DiffSpace::Query => "diff_aln_to_both",
            DiffSpace::Reference => "reference_diff",
        };
        let (category, in_a, in_b) = match class.map_status {
            MapStatus::BothMapped if !class.query_identical => {
                (both_mapped_diff_category, true, true)
            }
            MapStatus::OnlyAMapped => ("diff_aln_only_A", true, false),
            MapStatus::OnlyBMapped => ("diff_aln_only_B", false, true),
            // query-identical (incl. reverse-complement) or unmapped-both → not a difference
            _ => return Ok(class),
        };

        writeln!(self.reads_w, "{read_name}\t{category}\t{bool_cols}")?;
        self.n_diff_rows += 1;

        let outcome = if in_a && in_b { Outcome::Both } else { Outcome::OnlySide };
        if in_a {
            match build_ivl(&get_a, outcome) {
                Some(iv) => self.vec_a.push(iv),
                None => self.n_bad_interval += 1,
            }
        }
        if in_b {
            match build_ivl(&get_b, outcome) {
                Some(iv) => self.vec_b.push(iv),
                None => self.n_bad_interval += 1,
            }
        }
        Ok(class)
    }

    /// Flush the reads table, merge the buffered intervals into region-BED
    /// tables, and return the run's stats. `source_desc` is a human-readable
    /// description of where the rows came from, used only in the region
    /// tables' `[INFO]` log lines.
    pub(crate) fn finish(
        mut self,
        regions_a_out: &str,
        regions_b_out: &str,
        label_a: &str,
        label_b: &str,
        source_desc: &str,
    ) -> Result<AlnDiffStats> {
        self.reads_w.flush()?;
        let loci_a = merge_and_count(self.vec_a, DiffAcc::default, fold_diff);
        write_region_table(regions_a_out, &loci_a, "A", label_a, "n_only_A", source_desc)?;
        let loci_b = merge_and_count(self.vec_b, DiffAcc::default, fold_diff);
        write_region_table(regions_b_out, &loci_b, "B", label_b, "n_only_B", source_desc)?;
        Ok(AlnDiffStats {
            n_diff_rows: self.n_diff_rows,
            n_bad_interval: self.n_bad_interval,
        })
    }
}

pub fn run(args: &FindAlnDiffArgs) -> Result<()> {
    let outdir = Path::new(&args.outdir);
    fs::create_dir_all(outdir)
        .with_context(|| format!("cannot create --outdir '{}'", args.outdir))?;

    let ext = if args.no_gzip { "" } else { ".gz" };
    // Non-default mode gets a `.junctions` filename segment so its outputs never
    // clobber the default (`all`) run at the same --outdir/--prefix, and so the
    // default run stays byte-for-byte backward-compatible.
    let tag = match args.compare_by {
        CompareBy::All => "",
        CompareBy::Junctions => ".junctions",
    };
    // `--space` picks the filename stem, so a query-space and reference-space
    // run at the same --outdir/--prefix never clobber each other, and the
    // default (`query`) run stays byte-for-byte backward-compatible.
    let stem = match args.space {
        DiffSpace::Query => "query_diff",
        DiffSpace::Reference => "reference_diff",
    };
    let identical_stem = match args.space {
        DiffSpace::Query => "query_identical",
        DiffSpace::Reference => "reference_identical",
    };
    let path = |name: String| outdir.join(name).to_string_lossy().into_owned();
    let reads_out = path(format!("{}.{}_reads{}.tsv{}", args.prefix, stem, tag, ext));
    let regions_a_out = path(format!("{}.{}_regions.A{}.bed{}", args.prefix, stem, tag, ext));
    let regions_b_out = path(format!("{}.{}_regions.B{}.bed{}", args.prefix, stem, tag, ext));
    let summary_out = path(format!("{}.{}_summary{}.tsv", args.prefix, stem, tag));
    let identical_out = path(format!("{}.{}_reads{}.tsv{}", args.prefix, identical_stem, tag, ext));

    // `junctions`/`genomic_junctions` are only read by `--compare-by junctions`,
    // but required unconditionally: they're always present in the comparison
    // table, so demanding them up front turns a mid-stream surprise into an
    // early, clear error.
    const NEEDED: [&str; 9] = [
        "TargetChr", "Strand", "cs", "Query_Start", "Query_End", "Target_Start", "Target_End",
        "junctions", "genomic_junctions",
    ];

    // ── Input format: resolve, then open ───────────────────────────────────────
    // For a Parquet input, project down to just the columns this command reads
    // (19 of the comparison table's 96) — Parquet skips decoding the rest, which
    // is where its per-column storage actually pays off.
    let wanted: Vec<String> = ["Read_Name", "Label_A", "Label_B"]
        .into_iter()
        .map(String::from)
        .chain(["A", "B"].iter().flat_map(|side| NEEDED.iter().map(move |b| format!("{b}_{side}"))))
        .collect();
    let wanted_refs: Vec<&str> = wanted.iter().map(String::as_str).collect();
    let mut source = open_table(&args.input, args.input_format, Some(&wanted_refs))?;

    // ── Header: column index, labels, per-side indices ────────────────────────
    let cols_owned = source
        .header()
        .with_context(|| format!("reading comparison table '{}'", args.input))?;
    let cols: Vec<&str> = cols_owned.iter().map(String::as_str).collect();
    let col_index: HashMap<&str, usize> =
        cols.iter().copied().enumerate().map(|(i, c)| (c, i)).collect();
    require_ab_schema(&cols)?;
    let read_name_idx = *col_index
        .get("Read_Name")
        .context("comparison table is missing column 'Read_Name'")?;

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

    // Human-readable labels come from each row's `Label_A` / `Label_B` columns
    // (picked up from the first data row); used for reporting only.
    let mut label_a = "A".to_string();
    let mut label_b = "B".to_string();
    let mut seen_row = false;

    let mut acc = AlnDiffAccumulator::new(&reads_out, args.space, args.compare_by)?;

    let mut identical_w: Option<Box<dyn Write>> = if args.emit_identical_reads {
        let mut w = open_output(Some(&identical_out))?;
        writeln!(w, "Read_Name\tcategory\t{BOOL_COLS}")?;
        Some(w)
    } else {
        None
    };

    // The one, mode-independent summary — same schema `compare`/
    // `compare-toolkit summary` write (see the module doc comment).
    let mut summary = CompareSummary::default();

    // ── Pass 1: stream rows → read TSV + differing-interval vectors ────────────
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

        // `base` is the raw, mode-independent classification — observed as-is
        // into `summary` so its meaning matches `compare`'s and
        // `compare-toolkit summary`'s (always "query+all" semantics).
        let base = classify(&get_a, &get_b);
        summary.observe(&base);

        let read_name = fields.get(read_name_idx).copied().unwrap_or("");
        let class = acc.observe_row(read_name, get_a, get_b, &base)?;

        if class.query_identical {
            if let Some(w) = identical_w.as_mut() {
                let bool_cols = format!(
                    "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                    (base.query_identical && !base.query_identical_rc) as u8,
                    (base.query_identical && base.query_identical_rc) as u8,
                    base.query_junctions_identical.unwrap_or(false) as u8,
                    matches!(base.ref_class, Some(RefClass::SamePositionSameAln)) as u8,
                    matches!(base.ref_class, Some(RefClass::SamePositionDiffAln)) as u8,
                    matches!(base.ref_class, Some(RefClass::DiffPositionSameAln)) as u8,
                    matches!(base.ref_class, Some(RefClass::DiffPositionDiffAln)) as u8,
                    base.ref_same_position_same_junctions.unwrap_or(false) as u8,
                );
                let cat = match args.space {
                    DiffSpace::Reference => "reference_identical",
                    DiffSpace::Query => match args.compare_by {
                        CompareBy::Junctions => "query_identical_junctions",
                        CompareBy::All if class.query_identical_rc => "query_identical_revcomp",
                        CompareBy::All => "query_identical_same_strand",
                    },
                };
                writeln!(w, "{read_name}\t{cat}\t{bool_cols}")?;
            }
        }
    }
    if let Some(w) = identical_w.as_mut() {
        w.flush()?;
    }

    // ── Pass 2: merge each coordinate space → region tables ────────────────────
    let stats = acc.finish(&regions_a_out, &regions_b_out, &label_a, &label_b, &args.input)?;

    // ── Summary (TSV + stderr), shared schema with `compare` ───────────────────
    let space_str = match args.space {
        DiffSpace::Query => "query",
        DiffSpace::Reference => "reference",
    };
    let compare_by_str = match args.compare_by {
        CompareBy::All => "all",
        CompareBy::Junctions => "junctions",
    };
    let extra = [("space", space_str), ("compare_by", compare_by_str)];
    summary.write_tsv(&summary_out, &label_a, &label_b, &extra)?;

    eprintln!(
        "Diff summary (space={space_str}, compare-by={compare_by_str}, A={label_a}, B={label_b}): \
         {} differing reads written to {reads_out}",
        stats.n_diff_rows
    );
    summary.render_stderr(&label_a, &label_b, &extra);
    if stats.n_bad_interval > 0 {
        eprintln!("  ({} intervals skipped: unparseable or degenerate coordinates)", stats.n_bad_interval);
    }
    eprintln!("Outputs in {}:", args.outdir);
    eprintln!("  {reads_out}");
    eprintln!("  {regions_a_out}");
    eprintln!("  {regions_b_out}");
    eprintln!("  {summary_out}");
    if args.emit_identical_reads {
        eprintln!("  {identical_out}");
    }
    Ok(())
}

fn write_region_table(
    path: &str,
    loci: &[Locus<DiffAcc>],
    coord: &str,
    label: &str,
    only_col: &str,
    input: &str,
) -> Result<()> {
    eprintln!(
        "[INFO] find-aln-diff: coord_space={coord}(label={label})  input={input}  -> {path}"
    );
    let mut w = open_output(Some(path))?;
    writeln!(w, "#chrom\tstart\tend\tn_reads\tn_both\t{only_col}\tn_plus\tn_minus")?;
    for l in loci {
        writeln!(
            w,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            l.chrom, l.start, l.end, l.acc.n_reads, l.acc.n_both, l.acc.n_only, l.acc.n_plus, l.acc.n_minus
        )?;
    }
    w.flush()?;
    Ok(())
}
