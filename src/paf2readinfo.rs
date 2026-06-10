//! Subcommand `paf2readinfo` *(deprecated alias)*: fused `PAF → readinfo TSV`.
//!
//! Superseded by `paf2tables --readinfo <out>`, which can also emit the alninfo
//! table in the same pass. This thin wrapper delegates to the shared
//! `paf2tables::stream_readinfo` logic so the output is identical; it is kept
//! only so existing scripts keep working.

use anyhow::Result;

use crate::io_utils::{open_input, open_output};
use crate::paf2tables::stream_readinfo;

#[derive(clap::Args, Debug)]
pub struct Paf2ReadInfoArgs {
    /// Input PAF file. Use '-' for stdin; '.gz' is auto-decompressed.
    /// Must be grouped by Query_Name (each read's alignments contiguous).
    #[arg(short = 'i', long = "input", value_name = "in.paf")]
    input: String,

    /// Output readinfo TSV file. Append '.gz' for gzip. Defaults to stdout.
    #[arg(short = 'o', long = "output", value_name = "readinfo.tsv[.gz]")]
    output: Option<String>,
}

pub fn run(args: &Paf2ReadInfoArgs) -> Result<()> {
    eprintln!(
        "NOTE: 'paf2readinfo' is deprecated; use 'paf2tables --readinfo <out>' instead."
    );
    let reader = open_input(&args.input)?;
    let mut out = open_output(args.output.as_deref())?;
    stream_readinfo(reader, &mut out, None, /* strict_grouping = */ false)?;
    out.flush()?;
    Ok(())
}
