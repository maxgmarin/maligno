//! Subcommand `paf2readinfo`: fused `PAF → readinfo TSV` in one streaming pass.
//!
//! Equivalent to `paf2alninfo | readinfo`, but never materializes the
//! intermediate alninfo file. The output is **byte-identical** to running the
//! two steps separately, because every PAF record is routed through the exact
//! same conversion + serialization path (`AlnInfo::from_paf` → `write_row`)
//! before being collapsed by the unchanged `readinfo` logic.
//!
//! Requires the input PAF to be **grouped** by `Query_Name` (each read's
//! alignments contiguous) — the same precondition `readinfo` already has.
//! Constant memory: only one read's group is held at a time.

use std::io::Write;

use anyhow::Result;

use crate::io_utils::{open_input, open_output};
use crate::paf_groups::PafGroups;
use crate::readinfo::{collapse_group, READINFO_HEADER};

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
    let mut out = open_output(args.output.as_deref())?;
    writeln!(out, "{READINFO_HEADER}")?;

    let reader = open_input(&args.input)?;
    let mut groups = PafGroups::new(reader, /* warn_unsorted = */ true);

    while let Some(mut group) = groups.next_group()? {
        collapse_group(&mut group).write(&mut out)?;
    }

    out.flush()?;
    Ok(())
}
