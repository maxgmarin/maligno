//! maligno — unified PAF alignment-comparison toolkit.
//!
//! **Start here:** `paf2tables` is the primary PAF entry point — it turns a PAF
//! into the per-alignment info table, the per-read summary table, or both in a
//! single pass:
//!
//!   0. `paf2tables`         PAF → alninfo TSV (35 cols) and/or readinfo TSV (33 cols)
//!
//! Core pipeline subcommands:
//!
//!   1. `readinfo`           alninfo TSV        → per-read summary TSV   (33 cols)
//!   2. `compare`            two readinfo TSVs  → per-read comparison TSV (88 cols by default;
//!                                                94 with `--compare-genomic-junctions`)
//!   3. `compare-junctions`  two readinfo TSVs  → streamlined splice-focused comparison (47 cols)
//!   4. `sam2paf`            SAM                → PAF  (utility; use before paf2tables)
//!   5. `pafcompare`         two PAFs           → comparison TSV (one pass; requires
//!                                                identical read-name order)
//!
//! Deprecated aliases (delegate to `paf2tables`; kept for backward compatibility):
//!
//!   - `paf2alninfo`         PAF → alninfo TSV   (= `paf2tables --alninfo`)
//!   - `paf2readinfo`        PAF → readinfo TSV  (= `paf2tables --readinfo`)
//!
//! Full pipeline (including optional SAM conversion):
//!
//!     SAM ──sam2paf──▶ PAF ──paf2alninfo──▶ alninfo.tsv ──readinfo──▶ readinfo.tsv ─┐
//!                                                                                    ├─compare─▶ compare.tsv
//!     SAM ──sam2paf──▶ PAF ──paf2alninfo──▶ alninfo.tsv ──readinfo──▶ readinfo.tsv ─┘
//!
//! The comparison is a streaming merge-join (constant memory) and an INNER
//! join: only reads present in BOTH readinfo files (matched on Read_Name +
//! Read_Len) produce an output row. Reads present in only one file are dropped
//! but counted in the end-of-run summary. Note that unaligned reads still
//! appear in each readinfo file (with Num_Aln=0), so they DO get compared.

// ── Pipeline modules ──────────────────────────────────────────────────────────
mod cigar_junctions;    // CIGAR-based intron extractor (utility; not yet wired in)
mod compare_junctions;  // streamlined splice-focused comparison
mod compare_streaming;
mod cs_parser;          // cs-tag parser  (PAF → alninfo path; also extracts genomic junctions)
mod io_utils;
mod junction;
mod paf;
mod paf2alninfo;        // deprecated alias → paf2tables --alninfo
mod paf2readinfo;       // deprecated alias → paf2tables --readinfo
mod paf2tables;         // primary PAF entry point (alninfo and/or readinfo, one pass)
mod paf_groups;         // shared PAF → per-read group reader (paf2tables / pafcompare)
mod pafcompare;         // fused paired-PAF → comparison (one pass)
mod readinfo;
mod record;

// ── sam2paf utility submodule ─────────────────────────────────────────────────
mod sam2paf; // SAM → PAF converter (self-contained; owns cigar/convert/cs_generator/md)

use anyhow::Result;
use clap::{Parser, Subcommand};

use compare_junctions::CompareJunctionsArgs;
use compare_streaming::CompareStreamingArgs;
use paf2alninfo::Paf2AlnInfoArgs;
use paf2readinfo::Paf2ReadInfoArgs;
use paf2tables::Paf2TablesArgs;
use pafcompare::PafCompareArgs;
use readinfo::ReadInfoArgs;
use sam2paf::Sam2pafArgs;

/// Unified PAF alignment-comparison toolkit.
#[derive(Parser, Debug)]
#[command(name = "maligno", version, about)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand, Debug)]
enum Commands {
    /// PAF → alninfo (35 cols) and/or readinfo (33 cols). The primary PAF entry point.
    Paf2tables(Paf2TablesArgs),
    /// alninfo TSV → per-read summary TSV (best alignment by ms, then AS, then MQ).
    Readinfo(ReadInfoArgs),
    /// Streaming comparison of two readinfo TSVs (strict order by default; constant memory).
    Compare(CompareStreamingArgs),
    /// Streamlined splice-junction-focused comparison (47 cols).
    CompareJunctions(CompareJunctionsArgs),
    /// SAM → PAF converter (utility; use before paf2tables to start the pipeline).
    Sam2paf(Sam2pafArgs),
    /// Fused paired-PAF → comparison in one pass (requires identical read-name order).
    Pafcompare(PafCompareArgs),
    /// [DEPRECATED: use `paf2tables --alninfo`] PAF → per-alignment info TSV.
    Paf2alninfo(Paf2AlnInfoArgs),
    /// [DEPRECATED: use `paf2tables --readinfo`] PAF → per-read summary TSV.
    Paf2readinfo(Paf2ReadInfoArgs),
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    match &cli.command {
        Commands::Paf2tables(args)       => paf2tables::run(args),
        Commands::Readinfo(args)         => readinfo::run(args),
        Commands::Compare(args)          => compare_streaming::run(args),
        Commands::CompareJunctions(args) => compare_junctions::run(args),
        Commands::Sam2paf(args)          => sam2paf::run(args),
        Commands::Pafcompare(args)       => pafcompare::run(args),
        Commands::Paf2alninfo(args)      => paf2alninfo::run(args),
        Commands::Paf2readinfo(args)     => paf2readinfo::run(args),
    }
}
