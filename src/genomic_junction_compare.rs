//! `toolkit genomic-junction-compare` — compare A's and B's splice junctions in
//! **genomic (reference) coordinates**, from an alignment comparison table.
//!
//! A genomic junction is `(contig, intron_start, intron_end, genome_strand)`
//! ([`GenomicJunction`]): 0-based half-open intron coordinates on the reference,
//! plus the strand of the transcript on the reference (`paf::genome_strand`,
//! from each side's `strand` + `transcript_strand`). The same coordinates on `+`
//! and `-` are different junctions; an unknown strand (`*`) matches only `*`.
//!
//! Each side's junctions are rebuilt from its representative alignment's `cs`
//! tag (`target_start` + the `cs` walk, the same arithmetic that produces the
//! table's `genomic_junctions` columns), which also yields each intron's
//! splice-site motif letters. `nn` letters (STAR, or `cs` rebuilt from MD) are
//! kept as `nn-nn`.
//!
//! Outputs (to `--outdir`, `--prefix`-named):
//!   1. `{prefix}.genomic_junction_compare.per_read.tsv.gz` — one row per read
//!      with at least one genomic junction on either side. Junction lists hold
//!      complete `(contig, start, end, strand)` tuples.
//!   2. `{prefix}.genomic_junction_compare.per_junction.tsv.gz` — one row per
//!      distinct junction, with read counts and motifs per side.
//!   3. `{prefix}.genomic_junction_compare.summary.tsv` — `category`/`count`.
//!
//! Coordinates are only comparable when A and B were aligned to the same
//! reference; with different references, don't run this command.

use std::collections::HashMap;
use std::fs;
use std::io::Write;
use std::path::Path;

use anyhow::{Context, Result};

use crate::columns::{
    gjc_per_read_header, header_line, CS, GJC_INPUT_SIDE_COLS, GJC_PER_JUNCTION_COLS, LABEL_A,
    LABEL_B, QUERY_NAME, SIDE_A, SIDE_B, STRAND, SUMMARY_HEADER, TARGET_NAME, TARGET_START,
    TRANSCRIPT_STRAND,
};
use crate::compare_summary::{labels_from_row, require_ab_schema, require_col, wanted_columns, SideIndex};
use crate::cs_parser::{oriented_motif, parse_cs};
use crate::io_utils::open_output;
use crate::junction::{
    format_genomic_junction_tuple, genomic_junction_set_diffs, genomic_junction_set_stats,
    GenomicJunction,
};
use crate::paf::genome_strand;
use crate::table_input::{open_table, InputFormat};

/// Compare A's and B's splice junctions in genomic coordinates, with strand and
/// splice-site motifs, from an alignment comparison table.
#[derive(clap::Args, Debug)]
pub struct GenomicJunctionCompareArgs {
    /// Alignment comparison table (.parquet or .tsv[.gz]; `-` reads TSV from stdin).
    /// See `--input-format`.
    #[arg(short = 'i', long = "input", value_name = "compare.parquet|compare.tsv[.gz]")]
    input: String,

    /// Input serialization. `auto` (default) selects Parquet for a
    /// `.parquet`-named `--input` and TSV otherwise. Parquet requires a real
    /// file path — it cannot be read from stdin (`-`).
    #[arg(long = "input-format", value_enum, default_value = "auto")]
    input_format: InputFormat,

    /// Output directory (created if it does not exist).
    #[arg(long = "outdir", value_name = "DIR")]
    outdir: String,

    /// Filename prefix for all outputs.
    #[arg(long = "prefix", value_name = "STR")]
    prefix: String,
}

/// A `target_name` value indicates an unmapped read when it is empty or `"*"`
/// (same convention as `compare_summary` and `query_junction_diff`).
#[inline]
fn is_mapped(target_name: &str) -> bool {
    !(target_name.is_empty() || target_name == "*")
}

/// One side's genomic junctions, rebuilt from its `cs`.
#[derive(Debug, Default, PartialEq)]
struct SideJunctions {
    genome_strand: char,
    /// In `cs` order, i.e. reference-ascending. An alignment's introns never
    /// repeat (its reference position only moves forward), so these are distinct.
    junctions: Vec<GenomicJunction>,
    /// Reference-forward motif letters, parallel to `junctions`.
    motifs: Vec<[u8; 4]>,
}

fn side_junctions(
    target_name: &str,
    target_start: &str,
    strand: &str,
    transcript_strand: &str,
    cs: &str,
) -> SideJunctions {
    if !is_mapped(target_name) {
        return SideJunctions { genome_strand: '*', ..Default::default() };
    }
    let first = |s: &str| s.chars().next().unwrap_or('*');
    let gs = genome_strand(first(strand), first(transcript_strand));
    let tstart: u64 = target_start.parse().unwrap_or(0);
    let stats = parse_cs(cs);
    let junctions = stats
        .raw_genomic_junctions
        .iter()
        .map(|&(s, e)| (target_name.to_string(), tstart + s, tstart + e, gs))
        .collect();
    SideJunctions { genome_strand: gs, junctions, motifs: stats.raw_junction_motifs }
}

/// The motifs as a Python tuple of strings, in the junction's orientation:
/// `()`, `('gt-ag',)`, `('gt-ag', 'gc-ag')`.
fn format_motif_list(motifs: &[[u8; 4]], genome_strand: char) -> String {
    let items: Vec<String> =
        motifs.iter().map(|&m| format!("'{}'", oriented_motif(m, genome_strand))).collect();
    match items.len() {
        0 => "()".to_string(),
        1 => format!("({},)", items[0]),
        _ => format!("({})", items.join(", ")),
    }
}

/// Per-junction accumulator: reads supporting it, and the motif seen, per side.
#[derive(Debug, Default)]
struct JunctionAcc {
    n_reads_a: u64,
    n_reads_b: u64,
    motif_a: Option<[u8; 4]>,
    motif_b: Option<[u8; 4]>,
}

/// Keep the first motif seen, but prefer real letters over `nnnn`.
fn note_motif(slot: &mut Option<[u8; 4]>, m: [u8; 4]) {
    match slot {
        None => *slot = Some(m),
        Some(cur) if cur == b"nnnn" && &m != b"nnnn" => *slot = Some(m),
        _ => {}
    }
}

/// The `category`/`count` summary.
#[derive(Debug, Default)]
struct Summary {
    reads_compared: u64,
    reads_with_junctions_both: u64,
    reads_with_junctions_only_a: u64,
    reads_with_junctions_only_b: u64,
    reads_with_junctions_neither: u64,
    reads_genomic_junctions_identical: u64,
    reads_genomic_junctions_differ: u64,
    reads_genome_strand_differ: u64,
    n_distinct_junctions_a: u64,
    n_distinct_junctions_b: u64,
    n_distinct_junctions_shared: u64,
    n_distinct_junctions_only_a: u64,
    n_distinct_junctions_only_b: u64,
    n_distinct_junctions_unknown_strand_a: u64,
    n_distinct_junctions_unknown_strand_b: u64,
}

impl Summary {
    fn rows(&self) -> Vec<(&'static str, u64)> {
        vec![
            ("reads_compared", self.reads_compared),
            ("reads_with_junctions_both", self.reads_with_junctions_both),
            ("reads_with_junctions_only_A", self.reads_with_junctions_only_a),
            ("reads_with_junctions_only_B", self.reads_with_junctions_only_b),
            ("reads_with_junctions_neither", self.reads_with_junctions_neither),
            ("reads_genomic_junctions_identical", self.reads_genomic_junctions_identical),
            ("reads_genomic_junctions_differ", self.reads_genomic_junctions_differ),
            ("reads_genome_strand_differ", self.reads_genome_strand_differ),
            ("n_distinct_junctions_A", self.n_distinct_junctions_a),
            ("n_distinct_junctions_B", self.n_distinct_junctions_b),
            ("n_distinct_junctions_shared", self.n_distinct_junctions_shared),
            ("n_distinct_junctions_only_A", self.n_distinct_junctions_only_a),
            ("n_distinct_junctions_only_B", self.n_distinct_junctions_only_b),
            ("n_distinct_junctions_unknown_strand_A", self.n_distinct_junctions_unknown_strand_a),
            ("n_distinct_junctions_unknown_strand_B", self.n_distinct_junctions_unknown_strand_b),
        ]
    }

    fn write_tsv(&self, path: &str, label_a: &str, label_b: &str) -> Result<()> {
        let mut w = open_output(Some(path))?;
        writeln!(w, "{}", SUMMARY_HEADER.join("\t"))?;
        writeln!(w, "label_A\t{label_a}")?;
        writeln!(w, "label_B\t{label_b}")?;
        for (k, v) in self.rows() {
            writeln!(w, "{k}\t{v}")?;
        }
        w.flush()?;
        Ok(())
    }
}

/// The per-read comparison of two sides' junction sets.
#[derive(Debug, PartialEq)]
struct ReadComparison {
    shared: u64,
    only_a: u64,
    only_b: u64,
    only_a_list: Vec<GenomicJunction>,
    only_b_list: Vec<GenomicJunction>,
}

fn compare_sides(a: &SideJunctions, b: &SideJunctions) -> ReadComparison {
    let (shared, only_a, only_b) = genomic_junction_set_stats(&a.junctions, &b.junctions);
    let (only_a_list, only_b_list) = genomic_junction_set_diffs(&a.junctions, &b.junctions);
    ReadComparison { shared, only_a, only_b, only_a_list, only_b_list }
}

pub fn run(args: &GenomicJunctionCompareArgs) -> Result<()> {
    let outdir = Path::new(&args.outdir);
    fs::create_dir_all(outdir).with_context(|| format!("cannot create --outdir '{}'", args.outdir))?;
    let path = |name: &str| {
        outdir
            .join(format!("{}.genomic_junction_compare.{name}", args.prefix))
            .to_string_lossy()
            .into_owned()
    };
    let per_read_out = path("per_read.tsv.gz");
    let per_junction_out = path("per_junction.tsv.gz");
    let summary_out = path("summary.tsv");

    let wanted = wanted_columns(&[QUERY_NAME, LABEL_A, LABEL_B], &GJC_INPUT_SIDE_COLS);
    let wanted_refs: Vec<&str> = wanted.iter().map(String::as_str).collect();
    let mut source = open_table(&args.input, args.input_format, Some(&wanted_refs))?;

    let cols_owned = source
        .header()
        .with_context(|| format!("reading comparison table '{}'", args.input))?;
    let cols: Vec<&str> = cols_owned.iter().map(String::as_str).collect();
    let col_index: HashMap<&str, usize> = cols.iter().copied().enumerate().map(|(i, c)| (c, i)).collect();
    require_ab_schema(&cols)?;
    let query_name_idx = require_col(&col_index, QUERY_NAME.name)?;
    let idx_a = SideIndex::resolve(&col_index, &GJC_INPUT_SIDE_COLS, SIDE_A)?;
    let idx_b = SideIndex::resolve(&col_index, &GJC_INPUT_SIDE_COLS, SIDE_B)?;

    let mut per_read_w = open_output(Some(&per_read_out))?;
    writeln!(per_read_w, "{}", gjc_per_read_header().join("\t"))?;

    let mut label_a = SIDE_A.to_string();
    let mut label_b = SIDE_B.to_string();
    let mut seen_row = false;
    let mut summary = Summary::default();
    let mut by_junction: HashMap<GenomicJunction, JunctionAcc> = HashMap::new();

    while let Some(fields_owned) = source.next_row()? {
        let fields: Vec<&str> = fields_owned.iter().map(String::as_str).collect();
        if !seen_row {
            (label_a, label_b) = labels_from_row(&col_index, &fields);
            seen_row = true;
        }
        let side = |idx: &SideIndex| {
            let get = |c: &str| idx.get(&fields, c);
            let s = side_junctions(
                get(TARGET_NAME.name),
                get(TARGET_START.name),
                get(STRAND.name),
                get(TRANSCRIPT_STRAND.name),
                get(CS.name),
            );
            (s, get(TARGET_NAME.name), get(STRAND.name), get(CS.name))
        };
        let (sa, target_a, strand_a, cs_a) = side(&idx_a);
        let (sb, target_b, strand_b, cs_b) = side(&idx_b);

        summary.reads_compared += 1;
        let (has_a, has_b) = (!sa.junctions.is_empty(), !sb.junctions.is_empty());
        match (has_a, has_b) {
            (true, true) => summary.reads_with_junctions_both += 1,
            (true, false) => summary.reads_with_junctions_only_a += 1,
            (false, true) => summary.reads_with_junctions_only_b += 1,
            (false, false) => {
                summary.reads_with_junctions_neither += 1;
                continue;
            }
        }

        let cmp = compare_sides(&sa, &sb);
        if has_a && has_b {
            if cmp.only_a == 0 && cmp.only_b == 0 {
                summary.reads_genomic_junctions_identical += 1;
            } else {
                summary.reads_genomic_junctions_differ += 1;
            }
            if sa.genome_strand != '*' && sb.genome_strand != '*' && sa.genome_strand != sb.genome_strand {
                summary.reads_genome_strand_differ += 1;
            }
        }

        let query_name = fields.get(query_name_idx).copied().unwrap_or("");
        writeln!(
            per_read_w,
            "{query_name}\t{target_a}\t{target_b}\t{strand_a}\t{strand_b}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{cs_a}\t{cs_b}",
            sa.genome_strand,
            sb.genome_strand,
            sa.junctions.len(),
            sb.junctions.len(),
            cmp.shared,
            cmp.only_a + cmp.only_b,
            cmp.only_a,
            cmp.only_b,
            format_genomic_junction_tuple(&sa.junctions),
            format_genomic_junction_tuple(&sb.junctions),
            format_motif_list(&sa.motifs, sa.genome_strand),
            format_motif_list(&sb.motifs, sb.genome_strand),
            format_genomic_junction_tuple(&cmp.only_a_list),
            format_genomic_junction_tuple(&cmp.only_b_list),
        )?;

        for (j, &m) in sa.junctions.into_iter().zip(sa.motifs.iter()) {
            let acc = by_junction.entry(j).or_default();
            acc.n_reads_a += 1;
            note_motif(&mut acc.motif_a, m);
        }
        for (j, &m) in sb.junctions.into_iter().zip(sb.motifs.iter()) {
            let acc = by_junction.entry(j).or_default();
            acc.n_reads_b += 1;
            note_motif(&mut acc.motif_b, m);
        }
    }
    per_read_w.flush()?;

    // ── Per-junction table + distinct-junction counts ───────────────────────
    let mut junctions: Vec<(GenomicJunction, JunctionAcc)> = by_junction.into_iter().collect();
    junctions.sort_by(|x, y| x.0.cmp(&y.0));
    let mut per_junction_w = open_output(Some(&per_junction_out))?;
    writeln!(per_junction_w, "{}", header_line(&GJC_PER_JUNCTION_COLS))?;
    for ((contig, start, end, strand), acc) in &junctions {
        let motif = |m: Option<[u8; 4]>| m.map_or(".".to_string(), |m| oriented_motif(m, *strand));
        writeln!(
            per_junction_w,
            "{contig}\t{start}\t{end}\t{strand}\t{}\t{}\t{}\t{}",
            acc.n_reads_a,
            acc.n_reads_b,
            motif(acc.motif_a),
            motif(acc.motif_b),
        )?;
        let (in_a, in_b) = (acc.n_reads_a > 0, acc.n_reads_b > 0);
        summary.n_distinct_junctions_a += in_a as u64;
        summary.n_distinct_junctions_b += in_b as u64;
        summary.n_distinct_junctions_shared += (in_a && in_b) as u64;
        summary.n_distinct_junctions_only_a += (in_a && !in_b) as u64;
        summary.n_distinct_junctions_only_b += (in_b && !in_a) as u64;
        if *strand == '*' {
            summary.n_distinct_junctions_unknown_strand_a += in_a as u64;
            summary.n_distinct_junctions_unknown_strand_b += in_b as u64;
        }
    }
    per_junction_w.flush()?;

    summary.write_tsv(&summary_out, &label_a, &label_b)?;

    eprintln!(
        "Genomic junction compare (A={label_a}, B={label_b}): {} of {} reads have junctions on both sides; \
         {} identical, {} differ",
        summary.reads_with_junctions_both,
        summary.reads_compared,
        summary.reads_genomic_junctions_identical,
        summary.reads_genomic_junctions_differ,
    );
    eprintln!(
        "  distinct junctions: A={}, B={} ({} shared, {} only in A, {} only in B)",
        summary.n_distinct_junctions_a,
        summary.n_distinct_junctions_b,
        summary.n_distinct_junctions_shared,
        summary.n_distinct_junctions_only_a,
        summary.n_distinct_junctions_only_b,
    );
    eprintln!("Outputs in {}:", args.outdir);
    eprintln!("  {per_read_out}");
    eprintln!("  {per_junction_out}");
    eprintln!("  {summary_out}");
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn junctions_are_rebuilt_from_cs_with_strand_and_motif() {
        // '-' alignment, read antisense ('-'): the transcript is on the reference '+'.
        let s = side_junctions("chr1", "1000", "-", "-", ":10~gt100ag:5~nn50nn:5");
        assert_eq!(s.genome_strand, '+');
        assert_eq!(
            s.junctions,
            vec![("chr1".to_string(), 1010, 1110, '+'), ("chr1".to_string(), 1115, 1165, '+')]
        );
        assert_eq!(s.motifs, vec![*b"gtag", *b"nnnn"]);
        assert_eq!(
            format_genomic_junction_tuple(&s.junctions),
            "(('chr1', 1010, 1110, '+'), ('chr1', 1115, 1165, '+'))"
        );
        assert_eq!(format_motif_list(&s.motifs, s.genome_strand), "('gt-ag', 'nn-nn')");
    }

    #[test]
    fn unmapped_or_unknown_strand_sides() {
        let unmapped = side_junctions("*", "0", "*", "*", "");
        assert_eq!(unmapped.genome_strand, '*');
        assert!(unmapped.junctions.is_empty());
        let unknown = side_junctions("chr1", "0", "+", "*", ":5~ct20ac:5");
        assert_eq!(unknown.genome_strand, '*');
        assert_eq!(format_motif_list(&unknown.motifs, '*'), "('ct-ac',)", "kept reference-forward");
    }

    #[test]
    fn strand_disagreement_makes_every_junction_differ() {
        let a = side_junctions("chr1", "0", "+", "+", ":5~gt20ag:5");
        let b = side_junctions("chr1", "0", "+", "-", ":5~gt20ag:5");
        let c = compare_sides(&a, &b);
        assert_eq!((c.shared, c.only_a, c.only_b), (0, 1, 1));
        assert_eq!(format_genomic_junction_tuple(&c.only_a_list), "(('chr1', 5, 25, '+'),)");
        assert_eq!(format_genomic_junction_tuple(&c.only_b_list), "(('chr1', 5, 25, '-'),)");
        let same = compare_sides(&a, &a);
        assert_eq!((same.shared, same.only_a, same.only_b), (1, 0, 0));
    }

    #[test]
    fn motif_prefers_real_letters_over_nn() {
        let mut slot = None;
        note_motif(&mut slot, *b"nnnn");
        note_motif(&mut slot, *b"gtag");
        note_motif(&mut slot, *b"gcag");
        assert_eq!(slot, Some(*b"gtag"));
    }

    #[test]
    fn end_to_end_on_a_small_table() {
        let dir = std::env::temp_dir().join(format!("maligno_gjc_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let input = dir.join("t.tsv");
        let mut header = vec!["query_name".to_string(), "label_A".into(), "label_B".into()];
        for side in [SIDE_A, SIDE_B] {
            header.extend(GJC_INPUT_SIDE_COLS.iter().map(|&c| crate::columns::sided(c, side)));
        }
        // Columns: query_name, labels, then target_name, target_start, strand,
        // transcript_strand, cs for A and for B.
        let rows = [
            "r1\tX\tY\tchr1\t100\t+\t+\t:5~gt20ag:5\tchr1\t100\t+\t+\t:5~nn20nn:5",
            "r2\tX\tY\tchr1\t100\t+\t+\t:5~gt20ag:5\t*\t0\t*\t*\t",
            "r3\tX\tY\tchr1\t100\t+\t*\t:10\tchr1\t100\t+\t*\t:10",
        ];
        let body = format!("{}\n{}\n", header.join("\t"), rows.join("\n"));
        std::fs::write(&input, body).unwrap();
        let args = GenomicJunctionCompareArgs {
            input: input.to_string_lossy().into_owned(),
            input_format: InputFormat::Auto,
            outdir: dir.to_string_lossy().into_owned(),
            prefix: "t".into(),
        };
        run(&args).unwrap();

        let summary = std::fs::read_to_string(dir.join("t.genomic_junction_compare.summary.tsv")).unwrap();
        let get = |k: &str| {
            summary.lines().find_map(|l| l.strip_prefix(&format!("{k}\t"))).unwrap().to_string()
        };
        assert_eq!(get("label_A"), "X");
        assert_eq!(get("reads_compared"), "3");
        assert_eq!(get("reads_with_junctions_both"), "1");
        assert_eq!(get("reads_with_junctions_only_A"), "1");
        assert_eq!(get("reads_with_junctions_neither"), "1");
        assert_eq!(get("reads_genomic_junctions_identical"), "1");
        assert_eq!(get("n_distinct_junctions_shared"), "1");

        let gz = |name: &str| {
            let f = std::fs::File::open(dir.join(name)).unwrap();
            let mut s = String::new();
            std::io::Read::read_to_string(&mut flate2::read::MultiGzDecoder::new(f), &mut s).unwrap();
            s
        };
        let per_junction = gz("t.genomic_junction_compare.per_junction.tsv.gz");
        assert_eq!(per_junction.lines().nth(1).unwrap(), "chr1\t105\t125\t+\t2\t1\tgt-ag\tnn-nn");
        let per_read = gz("t.genomic_junction_compare.per_read.tsv.gz");
        assert_eq!(per_read.lines().count(), 3, "header + the two reads with junctions");
        let header_len = gjc_per_read_header().len();
        assert!(per_read.lines().all(|l| l.split('\t').count() == header_len));
        let _ = std::fs::remove_dir_all(&dir);
    }
}
