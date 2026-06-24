//! Aggregate summary statistics over a per-read comparison.
//!
//! Two delivery paths share the SAME classification logic (`classify`) and the
//! SAME accumulator (`CompareSummary`):
//!
//!   1. Built into `compare`: each matched read is `observe`d as its row streams
//!      out (O(1) memory — only counters are kept), and the summary is written as
//!      `{prefix}.compare[.junctions].summary.tsv` plus an stderr block.
//!   2. The standalone `compare-summary` command: streams an existing comparison
//!      TSV (`compare` / `compare-readinfo` output) row-by-row and emits the same
//!      summary. Serves the manual `paf2tables` → `compare-readinfo` workflow.
//!
//! `classify` reads per-side values by **unsuffixed** readinfo column name through
//! two accessor closures, so it is independent of the compare `--mode` (the columns
//! it needs — `cs`, `Strand`, `Query_Start`, `Query_End`, `TargetChr`,
//! `Target_Start`, `Target_End` — are present in both the 94-col `full` and the
//! 47-col `junctions` outputs).

use std::collections::HashMap;
use std::io::{BufRead, Write};

use anyhow::{bail, Context, Result};

use crate::cs_parser::cs_revcomp;
use crate::io_utils::{open_input, open_output};

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
            let cs_a = get_a("cs");
            let cs_b = get_b("cs");
            if get_a("Strand") == get_b("Strand") {
                if cs_a == cs_b {
                    query_identical = true;
                    // Reference identity only applies to the same-strand case.
                    reference_identical = get_a("TargetChr") == get_b("TargetChr")
                        && get_a("Target_Start") == get_b("Target_Start")
                        && get_a("Target_End") == get_b("Target_End");
                }
            } else if cs_a == cs_revcomp(cs_b) {
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
    fn rows(&self, label_a: &str, label_b: &str) -> Vec<(String, u64)> {
        vec![
            ("reads_compared".to_string(), self.reads_compared),
            ("aligned_both".to_string(), self.aligned_both),
            (format!("aligned_only_{label_a}"), self.aligned_only_a),
            (format!("aligned_only_{label_b}"), self.aligned_only_b),
            ("aligned_neither".to_string(), self.aligned_neither),
            ("query_identical".to_string(), self.query_identical),
            ("query_identical_same_strand".to_string(), self.query_identical_same_strand),
            ("query_identical_revcomp".to_string(), self.query_identical_rc),
            ("query_not_identical".to_string(), self.query_not_identical()),
            ("reference_identical".to_string(), self.reference_identical),
            (format!("present_only_in_{label_a}_by_id"), self.a_only_by_id),
            (format!("present_only_in_{label_b}_by_id"), self.b_only_by_id),
        ]
    }

    /// Write the summary as a 2-column TSV (`Category<TAB>Count`).
    pub fn write_tsv(&self, path: &str, label_a: &str, label_b: &str) -> Result<()> {
        let mut w = open_output(Some(path))?;
        writeln!(w, "Category\tCount")?;
        for (k, v) in self.rows(label_a, label_b) {
            writeln!(w, "{k}\t{v}")?;
        }
        w.flush()?;
        Ok(())
    }

    /// Print a human-readable block to stderr.
    pub fn render_stderr(&self, label_a: &str, label_b: &str) {
        eprintln!("Comparison summary:");
        for (k, v) in self.rows(label_a, label_b) {
            eprintln!("  {k:<34} {v}");
        }
    }
}

// ── Standalone `compare-summary` command ────────────────────────────────────────

#[derive(clap::Args, Debug)]
pub struct CompareSummaryArgs {
    /// Comparison TSV from `compare` / `compare-readinfo` (`.gz` ok; `-` for stdin).
    #[arg(short = 'i', long = "input", value_name = "compare.tsv")]
    input: String,

    /// Write the summary TSV here (default: stderr only). `.gz` ok; `-` for stdout.
    #[arg(short = 'o', long = "output", value_name = "summary.tsv")]
    output: Option<String>,
}

pub fn run(args: &CompareSummaryArgs) -> Result<()> {
    let mut reader = open_input(&args.input)
        .with_context(|| format!("opening comparison table '{}'", args.input))?;

    // Header → column index map.
    let mut header = String::new();
    if reader.read_line(&mut header)? == 0 {
        bail!("comparison table '{}' is empty", args.input);
    }
    let header = header.trim_end_matches(['\n', '\r']);
    let cols: Vec<&str> = header.split('\t').collect();
    let col_index: HashMap<&str, usize> =
        cols.iter().copied().enumerate().map(|(i, c)| (c, i)).collect();

    // Detect label_a / label_b from the two `TargetChr_<label>` columns, in order.
    let labels: Vec<String> = cols
        .iter()
        .filter_map(|c| c.strip_prefix("TargetChr_").map(|s| s.to_string()))
        .collect();
    if labels.len() != 2 {
        bail!(
            "expected exactly two `TargetChr_<label>` columns in the header, found {} \
             — is '{}' a maligno compare table?",
            labels.len(),
            args.input
        );
    }
    let (label_a, label_b) = (labels[0].clone(), labels[1].clone());

    // The unsuffixed columns `classify` reads; resolve each side's index up front.
    const NEEDED: [&str; 7] = [
        "TargetChr", "Strand", "cs", "Query_Start", "Query_End", "Target_Start", "Target_End",
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

    let mut summary = CompareSummary::default();
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
