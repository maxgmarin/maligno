//! `find-aln-diff` — from a `compare` / `toolkit compare-readinfo` table
//! (Parquet), find every read whose core alignment differs between sets A and
//! B, and the genomic regions where those reads cluster.
//!
//! `--mode` defines when a read's two alignments are the same. MAPQ and scores
//! are never compared.
//!   - `any-diff-query` (default) — same query span and same `cs` (intron
//!     donor/acceptor motif letters ignored). Opposite strands count as the
//!     same when B's `cs` is the exact reverse complement of A's.
//!   - `any-diff-reference` — the same as `any-diff-query` on the same strand,
//!     and also at the same `target_name` / `target_start`. For A and B aligned
//!     to the same reference.
//!
//! In both modes a read mapped on one side only differs, and a read mapped on
//! neither side is identical. `--aligned-both-only` compares only the reads
//! mapped on both sides; the others are counted in the summary as excluded.
//!
//! Outputs (to `--outdir`, `--prefix`-named; `{stem}` is `query_aln_diff` or
//! `reference_aln_diff`; gzipped unless `--no-gzip`):
//!   1. `{prefix}.{stem}_reads.tsv[.gz]` — one row per differing read:
//!      `query_name`, `outcome` (`diff_aln_both`, `diff_aln_only_A`,
//!      `diff_aln_only_B`), then the 7 classification columns.
//!      `--emit-identical-reads` adds `{prefix}.{query|reference}_aln_identical_reads.tsv[.gz]`
//!      for the identical reads (`identical_same_strand`, `identical_revcomp`,
//!      `neither_mapped`), in the same layout.
//!   2. `{prefix}.{stem}_regions.{A,B}.bed[.gz]` — each differing read's
//!      representative span on that side, merged into regions.
//!   3. `{prefix}.{stem}_summary.tsv` — the mode and parameters, then this
//!      run's selection counts.
//!
//! The per-read decision ([`decide`]) and the per-read core
//! ([`AlnDiffAccumulator`]) are independent of where rows come from: `run`
//! feeds them rows of the Parquet table.

use std::fs;
use std::io::Write;
use std::path::Path;

use anyhow::{bail, Context, Result};

use crate::columns::{
    CLASSIFICATION_COLS, CLASSIFY_SIDE_COLS, LABEL_A, LABEL_B, QUERY_NAME, SIDE_A, SIDE_B, STRAND,
    SUMMARY_HEADER, TARGET_NAME, TARGET_END, TARGET_START,
};
use crate::compare_summary::{
    classify, labels_from_row, require_ab_schema, require_col, wanted_columns, MapStatus,
    ReadClass, RefClass, SideIndex,
};
use crate::interval_merge::{merge_and_count, Ivl, Locus};
use crate::io_utils::open_output;
use crate::table_input::{open_table, InputFormat};
use std::collections::HashMap;

/// When a read's A and B alignments count as the same.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, clap::ValueEnum)]
pub enum DiffMode {
    /// Same query span and same `cs` (motif-blind); an exact reverse complement
    /// on opposite strands also counts as the same.
    #[default]
    AnyDiffQuery,
    /// Same as `any-diff-query` on the same strand, and at the same
    /// `target_name` / `target_start`.
    AnyDiffReference,
}

impl DiffMode {
    fn name(self) -> &'static str {
        match self {
            DiffMode::AnyDiffQuery => "any-diff-query",
            DiffMode::AnyDiffReference => "any-diff-reference",
        }
    }

    /// Filename stems: (differing reads, identical reads).
    fn stems(self) -> (&'static str, &'static str) {
        match self {
            DiffMode::AnyDiffQuery => ("query_aln_diff", "query_aln_identical"),
            DiffMode::AnyDiffReference => ("reference_aln_diff", "reference_aln_identical"),
        }
    }
}

/// Find reads whose core alignment differs between A and B, and the genomic
/// regions where they cluster.
#[derive(clap::Args, Debug)]
pub struct FindAlnDiffArgs {
    /// Alignment comparison table in Parquet format.
    #[arg(short = 'i', long = "input", value_name = "compare.parquet")]
    input: String,

    /// Output directory (created if it does not exist).
    #[arg(short = 'o', long = "outdir", value_name = "DIR")]
    outdir: String,

    /// Filename prefix for all outputs.
    #[arg(short = 'p', long = "prefix", value_name = "STR")]
    prefix: String,

    /// When a read's A and B alignments count as the same.
    #[arg(long = "mode", value_enum, default_value = "any-diff-query")]
    mode: DiffMode,

    /// Compare only the reads mapped in both A and B. Reads mapped on one side
    /// or neither are left out of every output and counted in the summary.
    #[arg(long = "aligned-both-only")]
    aligned_both_only: bool,

    /// Also write a table of the reads whose alignments are the same.
    #[arg(long = "emit-identical-reads")]
    emit_identical_reads: bool,

    /// Do not gzip the read and region output tables (gzipped by default).
    #[arg(long = "no-gzip")]
    no_gzip: bool,
}

/// How a differing read differs.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Outcome {
    Both,
    OnlyA,
    OnlyB,
}

impl Outcome {
    fn label(self) -> &'static str {
        match self {
            Outcome::Both => "diff_aln_both",
            Outcome::OnlyA => "diff_aln_only_A",
            Outcome::OnlyB => "diff_aln_only_B",
        }
    }
}

/// Why an identical read is identical.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum IdenticalKind {
    SameStrand,
    Revcomp,
    NeitherMapped,
}

impl IdenticalKind {
    fn label(self) -> &'static str {
        match self {
            IdenticalKind::SameStrand => "identical_same_strand",
            IdenticalKind::Revcomp => "identical_revcomp",
            IdenticalKind::NeitherMapped => "neither_mapped",
        }
    }
}

/// The decision for one read.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Decision {
    /// Left out by `--aligned-both-only`.
    Excluded,
    Identical(IdenticalKind),
    Differ(Outcome),
}

/// Decide whether one read's alignments differ under `mode`, from its
/// `classify()` output.
pub(crate) fn decide(mode: DiffMode, aligned_both_only: bool, base: &ReadClass) -> Decision {
    match base.map_status {
        MapStatus::BothMapped => {}
        _ if aligned_both_only => return Decision::Excluded,
        MapStatus::OnlyAMapped => return Decision::Differ(Outcome::OnlyA),
        MapStatus::OnlyBMapped => return Decision::Differ(Outcome::OnlyB),
        MapStatus::NeitherMapped => return Decision::Identical(IdenticalKind::NeitherMapped),
    }
    let identical = match mode {
        DiffMode::AnyDiffQuery => base.query_identical,
        DiffMode::AnyDiffReference => {
            base.query_identical
                && !base.query_identical_rc
                && base.ref_class.is_some_and(|rc| rc.same_position())
        }
    };
    match (identical, base.query_identical_rc) {
        (false, _) => Decision::Differ(Outcome::Both),
        (true, false) => Decision::Identical(IdenticalKind::SameStrand),
        (true, true) => Decision::Identical(IdenticalKind::Revcomp),
    }
}

/// The 7 classification columns of a reads-table row, the same in every mode.
/// Same-strand / reverse-complement identity is only set for reads mapped on
/// both sides.
fn classification_cells(base: &ReadClass) -> String {
    let both_mapped = base.map_status == MapStatus::BothMapped;
    format!(
        "{}\t{}\t{}\t{}\t{}\t{}\t{}",
        (both_mapped && base.query_identical && !base.query_identical_rc) as u8,
        (both_mapped && base.query_identical && base.query_identical_rc) as u8,
        base.query_junctions_identical.unwrap_or(false) as u8,
        matches!(base.ref_class, Some(RefClass::SamePositionSameAln)) as u8,
        matches!(base.ref_class, Some(RefClass::SamePositionDiffAln)) as u8,
        matches!(base.ref_class, Some(RefClass::DiffPositionSameAln)) as u8,
        matches!(base.ref_class, Some(RefClass::DiffPositionDiffAln)) as u8,
    )
}

/// Header of the reads tables: the read name, `outcome`, then the 7
/// classification columns.
fn reads_header() -> String {
    format!("{}\toutcome\t{}", QUERY_NAME.name, CLASSIFICATION_COLS.join("\t"))
}

/// Per-interval payload folded into each merged locus.
struct DiffMeta {
    strand: char,
    both: bool,
}

/// Per-locus accumulator.
#[derive(Default)]
struct DiffAcc {
    n_diff_aln_total: u64,
    n_diff_aln_both: u64,
    n_diff_aln_only: u64,
    n_diff_aln_plus_strand: u64,
    n_diff_aln_minus_strand: u64,
}

fn fold_diff(a: &mut DiffAcc, m: &DiffMeta) {
    a.n_diff_aln_total += 1;
    if m.both {
        a.n_diff_aln_both += 1;
    } else {
        a.n_diff_aln_only += 1;
    }
    match m.strand {
        '+' => a.n_diff_aln_plus_strand += 1,
        '-' => a.n_diff_aln_minus_strand += 1,
        _ => {}
    }
}

/// Build a genomic interval from one side's accessor. Returns `None` (caller
/// counts it as a bad interval) if the coord is unmapped, unparseable, or `end <= start`.
fn build_ivl<'a>(get: impl Fn(&str) -> &'a str, both: bool) -> Option<Ivl<DiffMeta>> {
    let chrom = get(TARGET_NAME.name);
    if chrom.is_empty() || chrom == "*" {
        return None;
    }
    let start: u64 = get(TARGET_START.name).parse().ok()?;
    let end: u64 = get(TARGET_END.name).parse().ok()?;
    if end <= start {
        return None;
    }
    let strand = get(STRAND.name).chars().next().unwrap_or('.');
    Some(Ivl { chrom: chrom.to_string(), start, end, meta: DiffMeta { strand, both } })
}

/// Per-decision read counts for the summary.
#[derive(Default, Debug, PartialEq, Eq)]
pub(crate) struct DiffCounts {
    pub reads_in_table: u64,
    pub excluded: u64,
    pub identical_same_strand: u64,
    pub identical_revcomp: u64,
    pub neither_mapped: u64,
    pub diff_aln_both: u64,
    pub diff_aln_only_a: u64,
    pub diff_aln_only_b: u64,
}

impl DiffCounts {
    fn note(&mut self, d: Decision) {
        self.reads_in_table += 1;
        match d {
            Decision::Excluded => self.excluded += 1,
            Decision::Identical(IdenticalKind::SameStrand) => self.identical_same_strand += 1,
            Decision::Identical(IdenticalKind::Revcomp) => self.identical_revcomp += 1,
            Decision::Identical(IdenticalKind::NeitherMapped) => self.neither_mapped += 1,
            Decision::Differ(Outcome::Both) => self.diff_aln_both += 1,
            Decision::Differ(Outcome::OnlyA) => self.diff_aln_only_a += 1,
            Decision::Differ(Outcome::OnlyB) => self.diff_aln_only_b += 1,
        }
    }

    fn identical(&self) -> u64 {
        self.identical_same_strand + self.identical_revcomp + self.neither_mapped
    }

    fn differ(&self) -> u64 {
        self.diff_aln_both + self.diff_aln_only_a + self.diff_aln_only_b
    }
}

/// Result of `AlnDiffAccumulator::finish`.
pub(crate) struct AlnDiffStats {
    pub counts: DiffCounts,
    pub n_regions_a: u64,
    pub n_regions_b: u64,
    pub n_bad_interval: u64,
}

/// The per-read core of `find-aln-diff`: decides each read with [`decide`],
/// writes it to the differing (or, when requested, identical) reads table, and
/// buffers its genomic interval(s) for the region tables. Independent of where
/// rows come from.
pub(crate) struct AlnDiffAccumulator {
    mode: DiffMode,
    aligned_both_only: bool,
    reads_w: Box<dyn Write>,
    identical_w: Option<Box<dyn Write>>,
    vec_a: Vec<Ivl<DiffMeta>>,
    vec_b: Vec<Ivl<DiffMeta>>,
    counts: DiffCounts,
    n_bad_interval: u64,
}

impl AlnDiffAccumulator {
    pub(crate) fn new(
        mode: DiffMode,
        aligned_both_only: bool,
        reads_out: &str,
        identical_out: Option<&str>,
    ) -> Result<Self> {
        let mut reads_w = open_output(Some(reads_out))?;
        writeln!(reads_w, "{}", reads_header())?;
        let identical_w = match identical_out {
            Some(p) => {
                let mut w = open_output(Some(p))?;
                writeln!(w, "{}", reads_header())?;
                Some(w)
            }
            None => None,
        };
        Ok(Self {
            mode,
            aligned_both_only,
            reads_w,
            identical_w,
            vec_a: Vec::new(),
            vec_b: Vec::new(),
            counts: DiffCounts::default(),
            n_bad_interval: 0,
        })
    }

    /// Decide one read from its `classify()` output, write its reads-table row,
    /// and buffer its interval(s) when it differs.
    pub(crate) fn observe_row<'a>(
        &mut self,
        read_name: &str,
        get_a: impl Fn(&str) -> &'a str,
        get_b: impl Fn(&str) -> &'a str,
        base: &ReadClass,
    ) -> Result<Decision> {
        let decision = decide(self.mode, self.aligned_both_only, base);
        self.counts.note(decision);
        match decision {
            Decision::Excluded => {}
            Decision::Identical(kind) => {
                if let Some(w) = self.identical_w.as_mut() {
                    writeln!(w, "{read_name}\t{}\t{}", kind.label(), classification_cells(base))?;
                }
            }
            Decision::Differ(outcome) => {
                writeln!(self.reads_w, "{read_name}\t{}\t{}", outcome.label(), classification_cells(base))?;
                let both = outcome == Outcome::Both;
                if outcome != Outcome::OnlyB {
                    match build_ivl(&get_a, both) {
                        Some(iv) => self.vec_a.push(iv),
                        None => self.n_bad_interval += 1,
                    }
                }
                if outcome != Outcome::OnlyA {
                    match build_ivl(&get_b, both) {
                        Some(iv) => self.vec_b.push(iv),
                        None => self.n_bad_interval += 1,
                    }
                }
            }
        }
        Ok(decision)
    }

    /// Flush the reads tables, merge the buffered intervals into the region
    /// tables, and return the run's stats.
    pub(crate) fn finish(mut self, regions_a_out: &str, regions_b_out: &str) -> Result<AlnDiffStats> {
        self.reads_w.flush()?;
        if let Some(w) = self.identical_w.as_mut() {
            w.flush()?;
        }
        let loci_a = merge_and_count(self.vec_a, DiffAcc::default, fold_diff);
        write_region_table(regions_a_out, &loci_a, "n_diff_aln_only_A")?;
        let loci_b = merge_and_count(self.vec_b, DiffAcc::default, fold_diff);
        write_region_table(regions_b_out, &loci_b, "n_diff_aln_only_B")?;
        Ok(AlnDiffStats {
            counts: self.counts,
            n_regions_a: loci_a.len() as u64,
            n_regions_b: loci_b.len() as u64,
            n_bad_interval: self.n_bad_interval,
        })
    }
}

/// The summary rows: the mode and parameters, then the selection counts.
fn summary_rows(
    label_a: &str,
    label_b: &str,
    mode: DiffMode,
    aligned_both_only: bool,
    s: &AlnDiffStats,
) -> Vec<(&'static str, String)> {
    let c = &s.counts;
    vec![
        ("label_A", label_a.to_string()),
        ("label_B", label_b.to_string()),
        ("mode", mode.name().to_string()),
        ("aligned_both_only", aligned_both_only.to_string()),
        ("reads_in_table", c.reads_in_table.to_string()),
        ("reads_not_aligned_both_excluded", c.excluded.to_string()),
        ("reads_compared", (c.reads_in_table - c.excluded).to_string()),
        ("reads_identical", c.identical().to_string()),
        ("identical_same_strand", c.identical_same_strand.to_string()),
        ("identical_revcomp", c.identical_revcomp.to_string()),
        ("neither_mapped", c.neither_mapped.to_string()),
        ("reads_differ", c.differ().to_string()),
        ("diff_aln_both", c.diff_aln_both.to_string()),
        ("diff_aln_only_A", c.diff_aln_only_a.to_string()),
        ("diff_aln_only_B", c.diff_aln_only_b.to_string()),
        ("n_regions_A", s.n_regions_a.to_string()),
        ("n_regions_B", s.n_regions_b.to_string()),
    ]
}

pub fn run(args: &FindAlnDiffArgs) -> Result<()> {
    if !args.input.ends_with(".parquet") {
        bail!(
            "--input must be a Parquet file ('.parquet') — got '{}'. \
             find-aln-diff doesn't support TSV input.",
            args.input
        );
    }
    let outdir = Path::new(&args.outdir);
    fs::create_dir_all(outdir)
        .with_context(|| format!("cannot create --outdir '{}'", args.outdir))?;

    let ext = if args.no_gzip { "" } else { ".gz" };
    let (stem, identical_stem) = args.mode.stems();
    let path = |name: String| outdir.join(name).to_string_lossy().into_owned();
    let reads_out = path(format!("{}.{stem}_reads.tsv{ext}", args.prefix));
    let regions_a_out = path(format!("{}.{stem}_regions.A.bed{ext}", args.prefix));
    let regions_b_out = path(format!("{}.{stem}_regions.B.bed{ext}", args.prefix));
    let summary_out = path(format!("{}.{stem}_summary.tsv", args.prefix));
    let identical_out = path(format!("{}.{identical_stem}_reads.tsv{ext}", args.prefix));

    // Read only the columns this command needs.
    let wanted = wanted_columns(&[QUERY_NAME, LABEL_A, LABEL_B], &CLASSIFY_SIDE_COLS);
    let wanted_refs: Vec<&str> = wanted.iter().map(String::as_str).collect();
    let mut source = open_table(&args.input, InputFormat::Parquet, Some(&wanted_refs))?;

    let cols_owned = source
        .header()
        .with_context(|| format!("reading comparison table '{}'", args.input))?;
    let cols: Vec<&str> = cols_owned.iter().map(String::as_str).collect();
    let col_index: HashMap<&str, usize> =
        cols.iter().copied().enumerate().map(|(i, c)| (c, i)).collect();
    require_ab_schema(&cols)?;
    let read_name_idx = require_col(&col_index, QUERY_NAME.name)?;
    let idx_a = SideIndex::resolve(&col_index, &CLASSIFY_SIDE_COLS, SIDE_A)?;
    let idx_b = SideIndex::resolve(&col_index, &CLASSIFY_SIDE_COLS, SIDE_B)?;

    // Labels come from the first row's `label_A` / `label_B`.
    let mut label_a = "A".to_string();
    let mut label_b = "B".to_string();
    let mut seen_row = false;

    let mut acc = AlnDiffAccumulator::new(
        args.mode,
        args.aligned_both_only,
        &reads_out,
        args.emit_identical_reads.then_some(identical_out.as_str()),
    )?;

    while let Some(fields_owned) = source.next_row()? {
        let fields: Vec<&str> = fields_owned.iter().map(String::as_str).collect();
        if !seen_row {
            (label_a, label_b) = labels_from_row(&col_index, &fields);
            seen_row = true;
        }
        let get_a = |c: &str| idx_a.get(&fields, c);
        let get_b = |c: &str| idx_b.get(&fields, c);
        let base = classify(&get_a, &get_b);
        let read_name = fields.get(read_name_idx).copied().unwrap_or("");
        acc.observe_row(read_name, get_a, get_b, &base)?;
    }

    let stats = acc.finish(&regions_a_out, &regions_b_out)?;
    let rows = summary_rows(&label_a, &label_b, args.mode, args.aligned_both_only, &stats);

    let mut w = open_output(Some(&summary_out))?;
    writeln!(w, "{}", SUMMARY_HEADER.join("\t"))?;
    for (k, v) in &rows {
        writeln!(w, "{k}\t{v}")?;
    }
    w.flush()?;

    eprintln!("find-aln-diff summary:");
    for (k, v) in &rows {
        eprintln!("  {k:<34} {v}");
    }
    if stats.n_bad_interval > 0 {
        eprintln!("  ({} intervals skipped: unparseable or degenerate coordinates)", stats.n_bad_interval);
    }
    eprintln!("Outputs in {}:", args.outdir);
    eprintln!("  {reads_out}");
    if args.emit_identical_reads {
        eprintln!("  {identical_out}");
    }
    eprintln!("  {regions_a_out}");
    eprintln!("  {regions_b_out}");
    eprintln!("  {summary_out}");
    Ok(())
}

fn write_region_table(path: &str, loci: &[Locus<DiffAcc>], only_col: &str) -> Result<()> {
    let mut w = open_output(Some(path))?;
    writeln!(
        w,
        "#chrom\tstart\tend\tn_diff_aln_total\tn_diff_aln_both\t{only_col}\t\
         n_diff_aln_plus_strand\tn_diff_aln_minus_strand"
    )?;
    for l in loci {
        writeln!(
            w,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            l.chrom,
            l.start,
            l.end,
            l.acc.n_diff_aln_total,
            l.acc.n_diff_aln_both,
            l.acc.n_diff_aln_only,
            l.acc.n_diff_aln_plus_strand,
            l.acc.n_diff_aln_minus_strand
        )?;
    }
    w.flush()?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn class(map_status: MapStatus, query_identical: bool, rc: bool, rc_class: Option<RefClass>) -> ReadClass {
        ReadClass {
            map_status,
            query_identical,
            query_identical_rc: rc,
            ref_class: rc_class,
            query_junctions_identical: None,
        }
    }

    const Q: DiffMode = DiffMode::AnyDiffQuery;
    const R: DiffMode = DiffMode::AnyDiffReference;

    #[test]
    fn one_sided_reads_differ_and_unmapped_reads_are_identical() {
        for mode in [Q, R] {
            let only_a = class(MapStatus::OnlyAMapped, false, false, None);
            let only_b = class(MapStatus::OnlyBMapped, false, false, None);
            let neither = class(MapStatus::NeitherMapped, true, false, None);
            assert_eq!(decide(mode, false, &only_a), Decision::Differ(Outcome::OnlyA));
            assert_eq!(decide(mode, false, &only_b), Decision::Differ(Outcome::OnlyB));
            assert_eq!(decide(mode, false, &neither), Decision::Identical(IdenticalKind::NeitherMapped));
            for c in [only_a, only_b, neither] {
                assert_eq!(decide(mode, true, &c), Decision::Excluded);
            }
        }
    }

    #[test]
    fn identical_alignment_at_the_same_place_is_identical_in_both_modes() {
        let c = class(MapStatus::BothMapped, true, false, Some(RefClass::SamePositionSameAln));
        for mode in [Q, R] {
            for both_only in [false, true] {
                assert_eq!(decide(mode, both_only, &c), Decision::Identical(IdenticalKind::SameStrand));
            }
        }
    }

    #[test]
    fn reverse_complement_is_identical_only_in_query_mode() {
        let c = class(MapStatus::BothMapped, true, true, Some(RefClass::DiffPositionDiffAln));
        assert_eq!(decide(Q, false, &c), Decision::Identical(IdenticalKind::Revcomp));
        assert_eq!(decide(R, false, &c), Decision::Differ(Outcome::Both));
    }

    #[test]
    fn same_alignment_elsewhere_is_identical_only_in_query_mode() {
        let c = class(MapStatus::BothMapped, true, false, Some(RefClass::DiffPositionSameAln));
        assert_eq!(decide(Q, true, &c), Decision::Identical(IdenticalKind::SameStrand));
        assert_eq!(decide(R, true, &c), Decision::Differ(Outcome::Both));
    }

    #[test]
    fn same_position_needs_query_identity_in_reference_mode() {
        // Same position and cs, but a different query span: not query-identical.
        let c = class(MapStatus::BothMapped, false, false, Some(RefClass::SamePositionSameAln));
        assert_eq!(decide(Q, false, &c), Decision::Differ(Outcome::Both));
        assert_eq!(decide(R, false, &c), Decision::Differ(Outcome::Both));
    }

    #[test]
    fn region_table_merges_spans_and_counts_outcomes() {
        let iv = |start, end, strand, both| Ivl {
            chrom: "chr1".to_string(),
            start,
            end,
            meta: DiffMeta { strand, both },
        };
        let loci = merge_and_count(
            vec![iv(100, 200, '+', true), iv(150, 300, '-', false), iv(500, 600, '+', true)],
            DiffAcc::default,
            fold_diff,
        );
        let path = std::env::temp_dir().join(format!("maligno_test_regions_{}.bed", std::process::id()));
        let p = path.to_str().unwrap();
        write_region_table(p, &loci, "n_diff_aln_only_A").unwrap();
        let content = std::fs::read_to_string(&path).unwrap();
        let _ = std::fs::remove_file(&path);
        let lines: Vec<&str> = content.lines().collect();
        assert_eq!(
            lines[0],
            "#chrom\tstart\tend\tn_diff_aln_total\tn_diff_aln_both\tn_diff_aln_only_A\tn_diff_aln_plus_strand\tn_diff_aln_minus_strand"
        );
        assert_eq!(lines[1], "chr1\t100\t300\t2\t1\t1\t1\t1");
        assert_eq!(lines[2], "chr1\t500\t600\t1\t1\t0\t1\t0");
        assert_eq!(lines.len(), 3);
    }
}
