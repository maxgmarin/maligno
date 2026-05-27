//! Subcommand `sam2paf`: SAM → PAF converter.
//!
//! A high-performance Rust port of the `sam2paf` sub-command from
//! paftools.MGM.js. Reads a SAM file (or stdin) and writes PAF records to
//! stdout. Output is byte-for-byte compatible with paftools.js sam2paf.
//!
//! This is a utility subcommand intended for use *before* the main pipeline:
//!   SAM ──sam2paf──▶ PAF ──paf2alninfo──▶ alninfo ──readinfo──▶ readinfo ──compare──▶ compare

mod cigar;
mod convert;
mod cs_generator;
mod md;

use std::fs::File;
use std::io::{self, BufReader, BufWriter, Write};

use anyhow::Result;

use convert::Options;

#[derive(clap::Args, Debug)]
pub struct Sam2pafArgs {
    /// Input SAM file; use '-' to read from stdin.
    #[arg(value_name = "in.sam")]
    pub input: String,

    /// Convert primary and supplementary alignments only
    /// (skip secondary, FLAG 0x100).
    #[arg(short = 'p', long = "primary-or-supp")]
    pub pri_only: bool,

    /// Convert primary alignments only
    /// (skip secondary FLAG 0x100 and supplementary FLAG 0x800).
    /// Implies -p.
    #[arg(short = 'P', long = "primary-only")]
    pub pri_pri_only: bool,

    /// Output the cs tag in long form (`=ACGT`).
    /// By default the short form (`:N`) is used.
    #[arg(short = 'L', long = "long-cs")]
    pub long_cs: bool,

    /// Emit placeholder PAF records for unmapped reads.
    /// By default unmapped records are silently discarded.
    #[arg(short = 'U', long = "unaligned")]
    pub convert_unaligned: bool,
}

pub fn run(args: &Sam2pafArgs) -> Result<()> {
    let opts = Options {
        pri_only:          args.pri_only || args.pri_pri_only,
        pri_pri_only:      args.pri_pri_only,
        long_cs:           args.long_cs,
        convert_unaligned: args.convert_unaligned,
    };

    // Large write buffer on stdout: flush in bulk rather than per-line.
    let stdout = io::stdout();
    let mut writer = BufWriter::with_capacity(1 << 20, stdout.lock());

    if args.input == "-" {
        let reader = BufReader::with_capacity(1 << 20, io::stdin().lock());
        convert::convert(reader, &mut writer, &opts)?;
    } else {
        let file = File::open(&args.input)
            .map_err(|e| anyhow::anyhow!("cannot open '{}': {}", args.input, e))?;
        let reader = BufReader::with_capacity(1 << 20, file);
        convert::convert(reader, &mut writer, &opts)?;
    }

    writer.flush()?;
    Ok(())
}
