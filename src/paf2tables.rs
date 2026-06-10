//! Subcommand `paf2tables`: the primary PAF entry point.
//!
//! Takes a PAF and writes the per-alignment **alninfo** table (35 cols), the
//! per-read **readinfo** summary (33 cols), or **both in a single pass**,
//! depending on which output paths are supplied:
//!
//! ```text
//!   --alninfo  only   → stream every alignment row (no grouping required)
//!   --readinfo only   → group by Query_Name, collapse to the best alignment
//!   --alninfo+--readinfo → do both in one pass (tee alninfo while collapsing)
//! ```
//!
//! This module owns the shared streaming logic; the legacy `paf2alninfo` and
//! `paf2readinfo` subcommands are thin deprecated aliases that delegate here, so
//! all three produce byte-identical output.
//!
//! ## Grouping requirement
//!
//! readinfo correctness depends on every read's alignments being **contiguous**
//! in the input (the same precondition `readinfo` has always had); alninfo never
//! cares about order. By default the reader groups contiguous `Query_Name` runs
//! and emits a one-shot warning if the input isn't byte-lex sorted.
//! `--strict-grouping` adds a rigorous guard that errors if a `Query_Name`
//! reappears non-contiguously (at the cost of O(#distinct reads) memory).

use std::collections::HashSet;
use std::io::{BufRead, Write};

use anyhow::{bail, Result};

use crate::io_utils::{open_input, open_output};
use crate::paf::parse_line;
use crate::paf_groups::PafGroups;
use crate::readinfo::{collapse_group, READINFO_HEADER};
use crate::record::AlnInfo;

#[derive(clap::Args, Debug)]
#[command(group(
    clap::ArgGroup::new("outputs")
        .required(true)
        .multiple(true)
        .args(["alninfo", "readinfo"])
))]
pub struct Paf2TablesArgs {
    /// Input PAF file. Use '-' for stdin; '.gz' is auto-decompressed.
    #[arg(short = 'i', long = "input", value_name = "in.paf")]
    input: String,

    /// Write the per-alignment info table (35 cols) here. '.gz' for gzip.
    /// Omit to skip the alninfo output.
    #[arg(short = 'a', long = "alninfo", value_name = "alninfo.tsv[.gz]")]
    alninfo: Option<String>,

    /// Write the per-read best-alignment summary (33 cols) here. '.gz' for gzip.
    /// Omit to skip the readinfo output. Requires the PAF be grouped by Query_Name.
    #[arg(short = 'r', long = "readinfo", value_name = "readinfo.tsv[.gz]")]
    readinfo: Option<String>,

    /// Error if a Query_Name reappears non-contiguously (rigorous grouping check;
    /// uses memory proportional to the number of distinct reads). Only affects
    /// the readinfo output.
    #[arg(long = "strict-grouping")]
    strict_grouping: bool,
}

pub fn run(args: &Paf2TablesArgs) -> Result<()> {
    match (args.alninfo.as_deref(), args.readinfo.as_deref()) {
        // Unreachable: clap's `outputs` ArgGroup (required = true) already rejects
        // an invocation with neither --alninfo nor --readinfo before run() is called.
        (None, None) => {
            bail!(
                "nothing to do: specify --alninfo <path> and/or --readinfo <path> \
                 to choose which table(s) to write."
            );
        }
        // alninfo-only: pure streaming, no grouping required (works on unsorted input).
        (Some(alninfo_path), None) => {
            let mut reader = open_input(&args.input)?;
            let mut w = open_output(Some(alninfo_path))?;
            stream_alninfo(&mut reader, &mut w)?;
            w.flush()?;
            Ok(())
        }
        // readinfo-only or both: group by Query_Name; tee alninfo when requested.
        (alninfo_path, Some(readinfo_path)) => {
            let reader = open_input(&args.input)?;
            let mut readinfo_out = open_output(Some(readinfo_path))?;
            match alninfo_path {
                Some(p) => {
                    let mut alninfo_out = open_output(Some(p))?;
                    stream_readinfo(
                        reader,
                        &mut readinfo_out,
                        Some(&mut *alninfo_out),
                        args.strict_grouping,
                    )?;
                    alninfo_out.flush()?;
                }
                None => {
                    stream_readinfo(reader, &mut readinfo_out, None, args.strict_grouping)?;
                }
            }
            readinfo_out.flush()?;
            Ok(())
        }
    }
}

/// Stream a PAF into the alninfo table, one row per alignment, in input order.
/// No grouping is performed — this is the single source of truth for the
/// `paf2alninfo` behavior. Malformed lines emit a WARNING and are skipped.
pub(crate) fn stream_alninfo<R: BufRead, W: Write>(reader: &mut R, w: &mut W) -> Result<()> {
    AlnInfo::write_header(w)?;

    let mut line = String::with_capacity(4096);
    let mut lineno: u64 = 0;
    loop {
        line.clear();
        let n = reader.read_line(&mut line)?;
        if n == 0 {
            break;
        }
        lineno += 1;
        let trimmed = line.trim_end_matches(|c| c == '\n' || c == '\r');
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        match parse_line(trimmed, lineno) {
            Ok(rec) => AlnInfo::from_paf(&rec).write_row(w)?,
            Err(e) => eprintln!("WARNING: {e}"),
        }
    }
    Ok(())
}

/// Stream a PAF into the readinfo table by grouping contiguous `Query_Name`
/// runs and collapsing each to its best alignment. If `alninfo_out` is `Some`,
/// every alignment's alninfo row is tee'd to it in the same pass (byte-identical
/// to `paf2alninfo`). With `strict_grouping`, errors if a `Query_Name` reappears
/// non-contiguously.
pub(crate) fn stream_readinfo<W: Write + ?Sized>(
    reader: Box<dyn BufRead>,
    readinfo_out: &mut W,
    mut alninfo_out: Option<&mut dyn Write>,
    strict_grouping: bool,
) -> Result<()> {
    if let Some(w) = alninfo_out.as_deref_mut() {
        AlnInfo::write_header(w)?;
    }
    writeln!(readinfo_out, "{READINFO_HEADER}")?;

    let mut groups = PafGroups::new(reader, /* warn_unsorted = */ true);
    let mut seen: Option<HashSet<String>> = if strict_grouping {
        Some(HashSet::new())
    } else {
        None
    };

    while let Some(mut group) = groups.next_group_tee(&mut alninfo_out)? {
        if let Some(seen) = seen.as_mut() {
            // group is non-empty by construction; first row carries the name.
            let name = group[0].query_name();
            if !seen.insert(name.to_owned()) {
                bail!(
                    "read {name:?} appears non-contiguously; input PAF is not \
                     grouped by Query_Name. Pre-sort with \
                     `LC_ALL=C sort -t$'\\t' -k1,1`, or drop --strict-grouping to \
                     group only contiguous runs.",
                );
            }
        }
        collapse_group(&mut group).write(readinfo_out)?;
    }

    Ok(())
}
