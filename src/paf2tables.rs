//! Subcommand `paf2tables`: one PAF → the per-alignment **alninfo** table (36
//! cols), the per-read **readinfo** table (37 cols), or both in a single pass,
//! depending on which output paths are supplied.
//!
//! The input is first sorted by `query_name` in byte order, with the same
//! external sort `compare` uses (`--sort-mem` in RAM, sorted runs spilled to
//! `--sort-tmp-dir`). The outputs are therefore byte-identical to `compare
//! --emit-alninfo` / `--emit-readinfo`, and the readinfo is ready for `toolkit
//! compare-readinfo`. For SAM/BAM, convert first with `maligno sam2paf -p -U`.
//!
//! `--presorted` skips the sort for input that is already sorted by
//! `query_name` in byte order (`LC_ALL=C sort -t$'\t' -k1,1`), and errors as
//! soon as a `query_name` is lower than the previous one.

use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};

use anyhow::{bail, Context, Result};

use crate::external_sort::{parse_mem, sort_lines_to_file};
use crate::io_utils::{open_input, open_output};
use crate::paf_groups::PafGroups;
use crate::readinfo::{collapse_group, readinfo_header};
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
    #[arg(short = 'i', long = "input", value_name = "in.paf[.gz]")]
    input: String,

    /// Write the per-alignment info table (36 cols) here. '.gz' for gzip.
    /// Omit to skip the alninfo output.
    #[arg(long = "alninfo", value_name = "alninfo.tsv[.gz]")]
    alninfo: Option<String>,

    /// Write the per-read best-alignment summary (37 cols) here. '.gz' for gzip.
    /// Omit to skip the readinfo output.
    #[arg(short = 'r', long = "readinfo", value_name = "readinfo.tsv[.gz]")]
    readinfo: Option<String>,

    /// Skip the sort: the input is already sorted by query_name in byte order
    /// (`LC_ALL=C sort -t$'\t' -k1,1`). Errors if it is not.
    #[arg(long = "presorted", conflicts_with = "keep_sorted_paf")]
    presorted: bool,

    /// In-memory sort buffer (K/M/G suffix, or plain bytes).
    #[arg(long = "sort-mem", value_name = "SIZE", default_value = "1G")]
    sort_mem: String,

    /// Directory for the sort's temporary spill files (default: the folder of
    /// --readinfo, or of --alninfo when --readinfo is not given).
    #[arg(long = "sort-tmp-dir", value_name = "DIR")]
    tmp_dir: Option<String>,

    /// Number of sort threads.
    #[arg(long = "sort-threads", value_name = "N", default_value_t = 1)]
    sort_threads: usize,

    /// Keep the sorted PAF as `<input stem>.sorted.paf.gz` in the folder of
    /// --readinfo (or --alninfo) instead of deleting it at the end. Needs a file
    /// --input, not stdin.
    #[arg(long = "keep-sorted-paf")]
    keep_sorted_paf: bool,
}

pub fn run(args: &Paf2TablesArgs) -> Result<()> {
    if args.keep_sorted_paf && args.input == "-" {
        bail!("--keep-sorted-paf names the sorted PAF after --input, so it needs a file, not stdin ('-')");
    }

    if args.presorted {
        return write_tables(&args.input, args, true);
    }

    // The output folder: where the kept (or temporary) sorted PAF goes, and the
    // default --sort-tmp-dir.
    let first_output = args.readinfo.as_deref().or(args.alninfo.as_deref()).unwrap_or("-");
    let out_dir = output_folder(first_output);
    let tmp_dir = args.tmp_dir.as_deref().map(PathBuf::from).unwrap_or_else(|| out_dir.clone());
    fs::create_dir_all(&tmp_dir)
        .with_context(|| format!("cannot create --sort-tmp-dir '{}'", tmp_dir.display()))?;
    let mem = parse_mem(&args.sort_mem)?;

    let sorted = if args.keep_sorted_paf {
        out_dir.join(format!("{}.sorted.paf.gz", input_stem(&args.input)))
    } else {
        out_dir.join(format!(".paf2tables.{}.sorted.paf", std::process::id()))
    };
    let sorted = sorted.to_string_lossy().into_owned();

    eprintln!("[INFO] sorting the input by query_name...");
    let result = open_input(&args.input)
        .and_then(|lines| sort_lines_to_file(lines, &sorted, mem, &tmp_dir, Some(args.sort_threads)))
        .with_context(|| format!("sorting input ({})", args.input))
        .and_then(|()| write_tables(&sorted, args, false));
    if !args.keep_sorted_paf || result.is_err() {
        let _ = fs::remove_file(&sorted);
    }
    result
}

/// Write the requested tables from a PAF sorted by `query_name`, in one grouped
/// pass. Removes any partial output on failure.
fn write_tables(paf: &str, args: &Paf2TablesArgs, require_sorted: bool) -> Result<()> {
    let result = (|| -> Result<()> {
        let mut alninfo_out = args.alninfo.as_deref().map(|p| open_output(Some(p))).transpose()?;
        let mut readinfo_out = args.readinfo.as_deref().map(|p| open_output(Some(p))).transpose()?;
        stream_tables(
            open_input(paf)?,
            alninfo_out.as_deref_mut().map(|w| w as &mut dyn Write),
            readinfo_out.as_deref_mut().map(|w| w as &mut dyn Write),
            require_sorted,
        )
    })();
    if result.is_err() {
        for p in [args.alninfo.as_deref(), args.readinfo.as_deref()].into_iter().flatten() {
            if p != "-" {
                let _ = fs::remove_file(p);
            }
        }
    }
    result
}

/// Stream a PAF grouped by `query_name` into the alninfo table (every alignment,
/// in input order) and/or the readinfo table (each group collapsed to its
/// representative alignment). Malformed lines emit a WARNING and are skipped.
pub(crate) fn stream_tables(
    reader: Box<dyn std::io::BufRead>,
    mut alninfo_out: Option<&mut dyn Write>,
    mut readinfo_out: Option<&mut dyn Write>,
    require_sorted: bool,
) -> Result<()> {
    if let Some(w) = alninfo_out.as_deref_mut() {
        AlnInfo::write_header(w)?;
    }
    if let Some(w) = readinfo_out.as_deref_mut() {
        writeln!(w, "{}", readinfo_header())?;
    }

    let mut groups = PafGroups::new(reader, require_sorted);
    while let Some(mut group) = groups.next_group_tee(&mut alninfo_out)? {
        if let Some(w) = readinfo_out.as_deref_mut() {
            collapse_group(&mut group).write(w)?;
        }
    }
    if let Some(w) = alninfo_out {
        w.flush()?;
    }
    if let Some(w) = readinfo_out {
        w.flush()?;
    }
    Ok(())
}

/// The folder an output path is written into (`.` for stdout or a bare name).
fn output_folder(path: &str) -> PathBuf {
    match Path::new(path).parent() {
        Some(p) if path != "-" && !p.as_os_str().is_empty() => p.to_path_buf(),
        _ => PathBuf::from("."),
    }
}

/// The input file name without its folder and its `.gz` and `.paf` extensions.
fn input_stem(path: &str) -> String {
    let name = Path::new(path).file_name().map(|n| n.to_string_lossy().into_owned()).unwrap_or_default();
    let name = name.strip_suffix(".gz").unwrap_or(&name);
    name.strip_suffix(".paf").unwrap_or(name).to_string()
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::{BufReader, Cursor};

    fn paf_line(name: &str, tstart: u64) -> String {
        format!("{name}\t100\t0\t100\t+\tchr1\t1000\t{tstart}\t{}\t100\t100\t60\ttp:A:P\tcs:Z::100\n", tstart + 100)
    }

    fn run_stream(paf: &str, require_sorted: bool) -> Result<(String, String)> {
        let reader: Box<dyn std::io::BufRead> = Box::new(BufReader::new(Cursor::new(paf.to_owned())));
        let (mut aln, mut ri) = (Vec::new(), Vec::new());
        stream_tables(reader, Some(&mut aln), Some(&mut ri), require_sorted)?;
        Ok((String::from_utf8(aln)?, String::from_utf8(ri)?))
    }

    #[test]
    fn presorted_input_must_be_in_byte_order() {
        let sorted = [paf_line("r1", 0), paf_line("r1", 500), paf_line("r2", 0)].concat();
        let (aln, ri) = run_stream(&sorted, true).unwrap();
        assert_eq!(aln.lines().count(), 4);
        assert_eq!(ri.lines().count(), 3);

        // Grouped but not in byte order.
        let unsorted = [paf_line("r2", 0), paf_line("r1", 0)].concat();
        let err = run_stream(&unsorted, true).unwrap_err().to_string();
        assert!(err.contains("not sorted by query_name"), "{err}");

        // A read whose rows are split up.
        let split = [paf_line("r1", 0), paf_line("r2", 0), paf_line("r1", 500)].concat();
        assert!(run_stream(&split, true).is_err());
    }

    #[test]
    fn alninfo_keeps_every_row_in_input_order() {
        let paf = [paf_line("r1", 0), paf_line("r1", 500), paf_line("r2", 0)].concat();
        let reader: Box<dyn std::io::BufRead> = Box::new(BufReader::new(Cursor::new(paf)));
        let mut aln = Vec::new();
        stream_tables(reader, Some(&mut aln), None, true).unwrap();
        let starts: Vec<String> = String::from_utf8(aln)
            .unwrap()
            .lines()
            .skip(1)
            .map(|l| l.split('\t').nth(7).unwrap().to_string())
            .collect();
        assert_eq!(starts, ["0", "500", "0"]);
    }

    #[test]
    fn input_stem_drops_folder_and_extensions() {
        assert_eq!(input_stem("/a/b/reads.paf.gz"), "reads");
        assert_eq!(input_stem("reads.paf"), "reads");
        assert_eq!(input_stem("x.sorted.paf.gz"), "x.sorted");
        assert_eq!(input_stem("reads"), "reads");
    }

    #[test]
    fn output_folder_defaults_to_cwd() {
        assert_eq!(output_folder("-"), PathBuf::from("."));
        assert_eq!(output_folder("ri.tsv.gz"), PathBuf::from("."));
        assert_eq!(output_folder("out/ri.tsv.gz"), PathBuf::from("out"));
    }
}
