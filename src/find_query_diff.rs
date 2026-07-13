//! `find-query-diff` — from a `compare` / `compare-readinfo` TSV, find
//! every read whose **query-space alignment differs** between sets A and B, and
//! the genomic regions where those reads cluster.
//!
//! "Query-different" is defined exactly as in `compare-summary` (via the shared
//! `classify`): a read is different iff it is in exactly one of
//!   - `diff_aln_to_both`  (mapped both sides, NOT query-identical),
//!   - `diff_aln_only_A`   (mapped only in A),
//!   - `diff_aln_only_B`   (mapped only in B).
//! Reverse-complement matches count as **identical** (not a difference), and
//! reads unmapped on both sides are excluded.
//!
//! `--compare-by` selects what "identical" means for both-mapped reads:
//!   - `all` (default) — the full cs tag must match, motif-blind (the `classify`
//!     definition above: intron donor/acceptor motif letters are ignored, so a
//!     differently-reported motif at the same intron position/length is not by
//!     itself a difference); any other mismatch/indel/soft-clip/junction-position
//!     difference counts.
//!   - `junctions` — only the **query-space splice-junction set** must match; reads
//!     with identical junctions but differing mismatches/indels/soft-clips count as
//!     the same. Reads aligned in only one set are still reported as differences
//!     (they have no junctions to compare on the missing side).
//!
//! Outputs (to `--outdir`, `--prefix`-named):
//!   1. `{prefix}.query_diff_reads.tsv[.gz]`      — one row per differing read + category
//!   2. `{prefix}.query_diff_regions.A.bed[.gz]`  — merged A-coordinate loci
//!   3. `{prefix}.query_diff_regions.B.bed[.gz]`  — merged B-coordinate loci
//!   4. `{prefix}.query_diff_summary.tsv`         — category tally (+ stderr)

use std::collections::HashMap;
use std::fs;
use std::io::{BufRead, Write};
use std::path::Path;

use anyhow::{bail, Context, Result};

use crate::compare_summary::{classify, detect_labels, CompareSummary, MapStatus, ReadClass};
use crate::interval_merge::{merge_and_count, Ivl, Locus};
use crate::io_utils::{open_input, open_output};
use crate::junction::{junction_set_stats, parse_junction_str};

#[derive(Clone, Copy, Debug, clap::ValueEnum)]
pub enum CoordSide {
    A,
    B,
    Both,
}

/// What aspect of the alignment defines a difference between A and B.
#[derive(Clone, Copy, Debug, PartialEq, Eq, clap::ValueEnum)]
pub enum CompareBy {
    /// Flag any difference across the whole alignment (compares the full cs tag,
    /// motif-blind: intron donor/acceptor letters are ignored).
    All,
    /// Flag only reads whose query-space splice-junction set differs.
    Junctions,
}

/// Find reads whose query-space alignment differs between A and B, and the
/// genomic regions where they cluster.
#[derive(clap::Args, Debug)]
pub struct FindQueryDiffArgs {
    /// Comparison TSV from `compare` / `compare-readinfo` (`.gz` ok; `-` = stdin).
    #[arg(short = 'i', long = "input", value_name = "compare.tsv")]
    input: String,

    /// Output directory (created if it does not exist).
    #[arg(long = "outdir", value_name = "DIR")]
    outdir: String,

    /// Filename prefix for all outputs.
    #[arg(long = "prefix", value_name = "STR")]
    prefix: String,

    /// Which coordinate-space region table(s) to emit. `both` (default) writes
    /// both A and B tables. The read TSV and summary are always written in full.
    #[arg(long = "coord-side", value_enum, default_value = "both")]
    coord_side: CoordSide,

    /// Gzip the read TSV and region tables (appends `.gz`).
    #[arg(long = "gzip")]
    gzip: bool,

    /// What defines a difference: `all` (default) compares the full cs tag,
    /// motif-blind (intron donor/acceptor letters ignored) — any other
    /// mismatch/indel/soft-clip/junction-position difference counts; `junctions`
    /// compares only the query-space splice-junction set (reads with identical
    /// junctions but differing mismatches/indels/soft-clips count as the same). In
    /// `junctions` mode the outputs gain a `.junctions` filename segment.
    #[arg(long = "compare-by", value_enum, default_value = "all")]
    compare_by: CompareBy,

    /// Also write a TSV of Read_Names whose query alignment is identical between
    /// A and B (same-strand or reverse-complement). Off by default.
    #[arg(long = "emit-identical-reads")]
    emit_identical_reads: bool,
}

impl FindQueryDiffArgs {
    /// Construct args for `compare`'s built-in invocation of this command: always
    /// both coordinate sides, always gzip'd output — `compare` does not expose
    /// either as its own option, so those defaults live here in exactly one place.
    /// `emit_identical_reads` stays off and is not exposed via `compare`.
    /// `compare_by` is pinned to `All` so the compare-fused run is byte-identical
    /// to the historical behavior; the junctions mode is not exposed via `compare`.
    pub(crate) fn for_compare(input: String, outdir: String, prefix: String) -> Self {
        Self {
            input,
            outdir,
            prefix,
            coord_side: CoordSide::Both,
            gzip: true,
            compare_by: CompareBy::All,
            emit_identical_reads: false,
        }
    }
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

/// Junction-space query identity for one both-mapped read: `true` iff the two
/// query-space junction *sets* are equal. Strand-agnostic — query junctions are
/// stored in plus-strand read coordinates (see `record.rs`), so no span gate and
/// no reverse-complement handling are needed (unlike the cs-tag comparison). Two
/// reads with no junctions on either side compare equal.
fn junctions_identical<'a>(
    get_a: &impl Fn(&str) -> &'a str,
    get_b: &impl Fn(&str) -> &'a str,
) -> bool {
    let ja = parse_junction_str(get_a("junctions"));
    let jb = parse_junction_str(get_b("junctions"));
    let (_overlap, only_a, only_b) = junction_set_stats(&ja, &jb);
    only_a == 0 && only_b == 0
}

pub fn run(args: &FindQueryDiffArgs) -> Result<()> {
    let outdir = Path::new(&args.outdir);
    fs::create_dir_all(outdir)
        .with_context(|| format!("cannot create --outdir '{}'", args.outdir))?;

    let want_a = matches!(args.coord_side, CoordSide::A | CoordSide::Both);
    let want_b = matches!(args.coord_side, CoordSide::B | CoordSide::Both);
    let ext = if args.gzip { ".gz" } else { "" };
    // Non-default mode gets a `.junctions` filename segment so its outputs never
    // clobber the default (`all`) run at the same --outdir/--prefix, and so the
    // default run stays byte-for-byte backward-compatible.
    let tag = match args.compare_by {
        CompareBy::All => "",
        CompareBy::Junctions => ".junctions",
    };
    let path = |name: String| outdir.join(name).to_string_lossy().into_owned();
    let reads_out = path(format!("{}.query_diff_reads{}.tsv{}", args.prefix, tag, ext));
    let regions_a_out = path(format!("{}.query_diff_regions.A{}.bed{}", args.prefix, tag, ext));
    let regions_b_out = path(format!("{}.query_diff_regions.B{}.bed{}", args.prefix, tag, ext));
    let summary_out = path(format!("{}.query_diff_summary{}.tsv", args.prefix, tag));
    let identical_out = path(format!("{}.query_identical_reads{}.tsv{}", args.prefix, tag, ext));

    // ── Header: column index, labels, per-side indices ────────────────────────
    let mut reader = open_input(&args.input)
        .with_context(|| format!("opening comparison table '{}'", args.input))?;
    let mut header = String::new();
    if reader.read_line(&mut header)? == 0 {
        bail!("comparison table '{}' is empty", args.input);
    }
    let header = header.trim_end_matches(['\n', '\r']);
    let cols: Vec<&str> = header.split('\t').collect();
    let col_index: HashMap<&str, usize> =
        cols.iter().copied().enumerate().map(|(i, c)| (c, i)).collect();
    let (label_a, label_b) = detect_labels(&cols)?;
    let read_name_idx = *col_index
        .get("Read_Name")
        .context("comparison table is missing column 'Read_Name'")?;

    // `junctions` (query-space set) is needed by `--compare-by junctions`; it is
    // present in both the full (94-col) and junctions (47-col) compare tables, so
    // requiring it unconditionally never breaks either input.
    const NEEDED: [&str; 8] = [
        "TargetChr", "Strand", "cs", "Query_Start", "Query_End", "Target_Start", "Target_End",
        "junctions",
    ];
    let resolve = |label: &str| -> Result<HashMap<&'static str, usize>> {
        let mut m = HashMap::new();
        for base in NEEDED {
            let name = format!("{base}_{label}");
            let idx = *col_index
                .get(name.as_str())
                .with_context(|| format!("comparison table is missing column '{name}'"))?;
            m.insert(base, idx);
        }
        Ok(m)
    };
    let idx_a = resolve(&label_a)?;
    let idx_b = resolve(&label_b)?;

    // ── Pass 1: stream rows → read TSV + differing-interval vectors ────────────
    let mut reads_w = open_output(Some(&reads_out))?;
    writeln!(reads_w, "Read_Name\toutcome")?;

    let mut identical_w: Option<Box<dyn Write>> = if args.emit_identical_reads {
        let mut w = open_output(Some(&identical_out))?;
        writeln!(w, "Read_Name\tcategory")?;
        Some(w)
    } else {
        None
    };

    let mut summary = CompareSummary::default();
    let mut vec_a: Vec<Ivl<DiffMeta>> = Vec::new();
    let mut vec_b: Vec<Ivl<DiffMeta>> = Vec::new();
    let mut n_bad_interval: u64 = 0;

    for line in reader.lines() {
        let line = line?;
        if line.is_empty() {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        let get_a = |c: &str| -> &str {
            idx_a.get(c).and_then(|&i| fields.get(i)).copied().unwrap_or("")
        };
        let get_b = |c: &str| -> &str {
            idx_b.get(c).and_then(|&i| fields.get(i)).copied().unwrap_or("")
        };

        // `classify` gives map_status (mode-independent) + the cs-based identity.
        // In `junctions` mode we keep the map_status but redefine "identical" as
        // "query-space junction sets are equal"; the rest of the pipeline (read
        // selection, summary, intervals) is driven by this single `class`, so the
        // outputs and summary stay consistent by construction.
        let base = classify(&get_a, &get_b);
        let class = match args.compare_by {
            CompareBy::All => base,
            CompareBy::Junctions => ReadClass {
                map_status: base.map_status,
                query_identical: matches!(base.map_status, MapStatus::BothMapped)
                    && junctions_identical(&get_a, &get_b),
                query_identical_rc: false,
                reference_identical: false,
            },
        };
        summary.observe(&class);
        let read_name = fields.get(read_name_idx).copied().unwrap_or("");

        if class.query_identical {
            if let Some(w) = identical_w.as_mut() {
                let cat = match args.compare_by {
                    CompareBy::Junctions => "query_identical_junctions",
                    CompareBy::All if class.query_identical_rc => "query_identical_revcomp",
                    CompareBy::All => "query_identical_same_strand",
                };
                writeln!(w, "{read_name}\t{cat}")?;
            }
        }

        // Select query-different reads and their placement(s).
        let (category, in_a, in_b) = match class.map_status {
            MapStatus::BothMapped if !class.query_identical => ("diff_aln_to_both", true, true),
            MapStatus::OnlyAMapped => ("diff_aln_only_A", true, false),
            MapStatus::OnlyBMapped => ("diff_aln_only_B", false, true),
            // query-identical (incl. reverse-complement) or unmapped-both → not a difference
            _ => continue,
        };

        writeln!(reads_w, "{read_name}\t{category}")?;

        let outcome = if in_a && in_b { Outcome::Both } else { Outcome::OnlySide };
        if in_a && want_a {
            match build_ivl(&get_a, outcome) {
                Some(iv) => vec_a.push(iv),
                None => n_bad_interval += 1,
            }
        }
        if in_b && want_b {
            match build_ivl(&get_b, outcome) {
                Some(iv) => vec_b.push(iv),
                None => n_bad_interval += 1,
            }
        }
    }
    reads_w.flush()?;
    if let Some(w) = identical_w.as_mut() {
        w.flush()?;
    }

    // ── Pass 2: merge each selected coordinate space → region tables ──────────
    if want_a {
        let loci = merge_and_count(vec_a, DiffAcc::default, fold_diff);
        write_region_table(&regions_a_out, &loci, "A", &label_a, "n_only_A", &args.input)?;
    }
    if want_b {
        let loci = merge_and_count(vec_b, DiffAcc::default, fold_diff);
        write_region_table(&regions_b_out, &loci, "B", &label_b, "n_only_B", &args.input)?;
    }

    // ── Summary (TSV + stderr) ────────────────────────────────────────────────
    let compare_by_str = match args.compare_by {
        CompareBy::All => "all",
        CompareBy::Junctions => "junctions",
    };
    let rows = summary_rows(&summary);
    let mut sw = open_output(Some(&summary_out))?;
    writeln!(sw, "Category\tCount")?;
    // Record the active mode in every summary (both `all` and `junctions`) so the
    // file is self-describing regardless of how it was produced.
    writeln!(sw, "compare_by\t{compare_by_str}")?;
    for (k, v) in &rows {
        writeln!(sw, "{k}\t{v}")?;
    }
    sw.flush()?;

    eprintln!("Query-diff summary (compare-by={compare_by_str}):");
    for (k, v) in &rows {
        eprintln!("  {k:<24} {v}");
    }
    if n_bad_interval > 0 {
        eprintln!("  ({n_bad_interval} intervals skipped: unparseable or degenerate coordinates)");
    }
    eprintln!("Outputs in {}:", args.outdir);
    eprintln!("  {reads_out}");
    if want_a {
        eprintln!("  {regions_a_out}");
    }
    if want_b {
        eprintln!("  {regions_b_out}");
    }
    eprintln!("  {summary_out}");
    if args.emit_identical_reads {
        eprintln!("  {identical_out}");
    }
    Ok(())
}

/// Derive this command's category rows from the shared `CompareSummary` counters.
/// `diff_aln_to_both == aligned_both - query_identical` (== `query_not_identical`).
fn summary_rows(s: &CompareSummary) -> Vec<(String, u64)> {
    let diff_both = s.aligned_both - s.query_identical;
    let diff_total = diff_both + s.aligned_only_a + s.aligned_only_b;
    vec![
        ("reads_compared".to_string(), s.reads_compared),
        ("query_different_total".to_string(), diff_total),
        ("diff_aln_to_both".to_string(), diff_both),
        ("diff_aln_only_A".to_string(), s.aligned_only_a),
        ("diff_aln_only_B".to_string(), s.aligned_only_b),
        ("query_identical_total".to_string(), s.query_identical),
        ("aligned_neither".to_string(), s.aligned_neither),
    ]
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
        "[INFO] find-query-diff: coord_space={coord}(label={label})  input={input}  -> {path}"
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

#[cfg(test)]
mod tests {
    use super::*;

    /// Build a column accessor that returns `juncs` for the `junctions` column and
    /// "" otherwise — enough to exercise `junctions_identical`.
    fn getter(juncs: &'static str) -> impl Fn(&str) -> &'static str {
        move |c| if c == "junctions" { juncs } else { "" }
    }

    #[test]
    fn junctions_identical_empty_equals_empty() {
        // Both reads unspliced → junction sets both empty → identical.
        assert!(junctions_identical(&getter("()"), &getter("()")));
    }

    #[test]
    fn junctions_identical_same_set() {
        assert!(junctions_identical(&getter("(10, 20)"), &getter("(10, 20)")));
    }

    #[test]
    fn junctions_identical_ignores_order() {
        // Set semantics: order does not matter.
        assert!(junctions_identical(&getter("(10, 20)"), &getter("(20, 10)")));
    }

    #[test]
    fn junctions_identical_disjoint_is_different() {
        assert!(!junctions_identical(&getter("(10,)"), &getter("(20,)")));
    }

    #[test]
    fn junctions_identical_one_side_empty_is_different() {
        assert!(!junctions_identical(&getter("(10,)"), &getter("()")));
    }

    #[test]
    fn junctions_identical_subset_is_different() {
        assert!(!junctions_identical(&getter("(10, 20)"), &getter("(10,)")));
    }
}
