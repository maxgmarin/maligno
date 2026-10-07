//! Subcommand `sam2paf`: SAM → PAF converter.
//!
//! A high-performance Rust port of the `sam2paf` sub-command from
//! paftools.MGM.js. Reads a SAM or BAM file (or SAM on stdin) and writes PAF records to
//! stdout. Output is byte-for-byte compatible with paftools.js sam2paf.
//!
//! A standalone utility; `compare` also runs this converter internally on SAM/BAM
//! inputs:
//!   SAM/BAM ──sam2paf──▶ PAF ──compare──▶ comparison table
//!   SAM/BAM ──sam2paf──▶ PAF ──toolkit paf2tables──▶ alninfo / readinfo

mod cigar;
pub(crate) mod convert;
mod cs_generator;
mod md;

use std::io::{self, BufReader, BufWriter, Write};

use anyhow::{bail, Result};

use convert::Options;

use crate::aln_input;

#[derive(clap::Args, Debug)]
pub struct Sam2pafArgs {
    /// Input SAM or BAM file (auto-detected); use '-' to read SAM from stdin.
    #[arg(value_name = "in.sam|in.bam")]
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
        // SAM or BAM (sniffed); a BAM is read as the SAM text `samtools view -h` prints.
        let fmt = aln_input::detect(&args.input)?;
        if fmt.is_paf() {
            bail!("'{}' is not a SAM or BAM file", args.input);
        }
        convert::convert(aln_input::open_sam_text(&args.input, fmt)?, &mut writer, &opts)?;
    }

    writer.flush()?;
    Ok(())
}
