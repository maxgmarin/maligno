//! Subcommand `paf2alninfo`: PAF → per-alignment info TSV.
//!
//! Reads a PAF file (optionally gzipped, or stdin via `-`), parses each
//! alignment record, computes cs-tag–derived statistics and derived scalars,
//! and writes a 35-column TSV (one row per alignment).
//!
//! Unaligned reads (PAF `Target_Name == "*"` or `Strand == "*"`) are **kept**:
//! they produce a full row with zeroed alignment statistics. This is what lets
//! unaligned reads flow downstream into readinfo and the comparison.
//!
//! Output is sorted by (Query_Name, Query_Start, Query_End) unless `--no-sort`,
//! matching the Python `paf_to_df()` ordering.

use anyhow::Result;
use std::io::{BufRead, Write};

use crate::io_utils::{open_input, open_output};
use crate::paf::parse_line;
use crate::record::AlnInfo;

#[derive(clap::Args, Debug)]
pub struct Paf2AlnInfoArgs {
    /// Input PAF file. Use '-' for stdin; '.gz' is auto-decompressed.
    #[arg(short = 'i', long = "input", value_name = "in.paf")]
    input: String,

    /// Output TSV file. Append '.gz' for gzip output. Defaults to stdout.
    #[arg(short = 'o', long = "output", value_name = "out.tsv[.gz]")]
    output: Option<String>,

    /// Disable the default sort by (Query_Name, Query_Start, Query_End).
    #[arg(long = "no-sort")]
    no_sort: bool,
}

pub fn run(args: &Paf2AlnInfoArgs) -> Result<()> {
    // ── Read + parse every PAF record ─────────────────────────────────────
    let reader = open_input(&args.input)?;
    let mut records = collect_records(reader)?;

    // ── Optional sort (default: on) ───────────────────────────────────────
    if !args.no_sort {
        records.sort_by(|a, b| {
            a.query_name
                .cmp(&b.query_name)
                .then_with(|| a.query_start.cmp(&b.query_start))
                .then_with(|| a.query_end.cmp(&b.query_end))
        });
    }

    // ── Write header + rows ───────────────────────────────────────────────
    let mut w = open_output(args.output.as_deref())?;
    AlnInfo::write_header(&mut w)?;
    for rec in &records {
        rec.write_row(&mut w)?;
    }
    w.flush()?;
    Ok(())
}

/// Read a PAF stream line by line; parse and compute each record into AlnInfo.
fn collect_records<R: BufRead>(mut reader: R) -> Result<Vec<AlnInfo>> {
    let mut records: Vec<AlnInfo> = Vec::new();
    let mut line = String::with_capacity(4096);
    let mut lineno = 0u64;

    while reader.read_line(&mut line)? > 0 {
        lineno += 1;
        let trimmed = line.trim_end_matches(|c| c == '\n' || c == '\r');

        // Skip blank/header lines (rare in PAF, but defensive).
        if trimmed.is_empty() || trimmed.starts_with('#') {
            line.clear();
            continue;
        }

        match parse_line(trimmed, lineno) {
            Ok(rec) => records.push(AlnInfo::from_paf(&rec)),
            Err(e) => eprintln!("WARNING: {e}"),
        }

        line.clear();
    }

    Ok(records)
}
