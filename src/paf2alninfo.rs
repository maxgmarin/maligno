//! Subcommand `paf2alninfo`: PAF → per-alignment info TSV (streaming).
//!
//! Reads a PAF file (optionally gzipped, or stdin via `-`), parses each
//! alignment record, computes cs-tag–derived statistics and derived scalars,
//! and writes a 35-column TSV (one row per alignment).
//!
//! Unaligned reads (PAF `Target_Name == "*"` or `Strand == "*"`) are **kept**:
//! they produce a full row with zeroed alignment statistics, allowing them to
//! flow downstream into readinfo and the comparison.
//!
//! **Pure streaming, constant memory.** As of v0.2.7 paf2alninfo no longer
//! collects the entire input in memory and no longer sorts. Each PAF line is
//! parsed and emitted independently in input order. For the standard pipeline
//! (`paf2alninfo` → `readinfo` → `compare`), pre-sort the PAF by Query_Name
//! once upstream:
//!
//!     LC_ALL=C sort -t$'\t' -k1,1 in.paf > sorted.paf
//!
//! Unix `sort` does external-sort with bounded memory and handles files
//! larger than RAM. The pre-sort satisfies both `readinfo`'s contiguity
//! requirement and `compare`'s byte-lex sort requirement in one pass.

use anyhow::Result;

use crate::io_utils::{open_input, open_output};
use crate::paf2tables::stream_alninfo;

#[derive(clap::Args, Debug)]
pub struct Paf2AlnInfoArgs {
    /// Input PAF file. Use '-' for stdin; '.gz' is auto-decompressed.
    #[arg(short = 'i', long = "input", value_name = "in.paf")]
    input: String,

    /// Output TSV file. Append '.gz' for gzip output. Defaults to stdout.
    #[arg(short = 'o', long = "output", value_name = "out.tsv[.gz]")]
    output: Option<String>,
}

pub fn run(args: &Paf2AlnInfoArgs) -> Result<()> {
    eprintln!(
        "NOTE: 'paf2alninfo' is deprecated; use 'paf2tables --alninfo <out>' instead."
    );
    let mut reader = open_input(&args.input)?;
    let mut w = open_output(args.output.as_deref())?;
    stream_alninfo(&mut reader, &mut w)?;
    w.flush()?;
    Ok(())
}
