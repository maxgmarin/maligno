<picture>
  <source media="(prefers-color-scheme: dark)" srcset="docs/maligno-dark.png">
  <source media="(prefers-color-scheme: light)" srcset="docs/maligno-light.png">
  <img alt="maligno" src="docs/maligno-light.png" width="300">
</picture>

# maligno

**maligno** is a toolkit for systematically comparing two sets of alignments read-by-read.
The primary command, **`maligno compare`**, takes two alignment files (e.g. PAF) for the same set of reads
and produces a detailed per-read alignment comparison table. This table enables easy comparison of alignment stats (alignment scores, coverage, indels, mismatches,
soft-clipping, and splice-junction agreement) with a difference column for each metric.

This makes it easy to compare alignment results across conditions such as:
1) Different aligners (`minimap2` vs `bwa mem`)
2) Different parameters of the same aligner (`splice` vs `splice:hq` for transcript alignment with `minimap2`)
3) Different reference genomes (`GRCh38` vs `CHM13`)

With a single comparison table, you can quickly ask:
- "How many reads aligned exactly the same?"
- "Which read IDs differ in their alignment between set A and B, and how?"
- "How many reads improved their alignment score by X points?"
- "How do splice junctions differ, in read or genomic coordinate space?"

In short, `maligno` gives you an efficient framework for detailed alignment comparisons
across different parameters, aligners, and references.

---

## Install

Requires a [Rust toolchain](https://rustup.rs).

```bash
cargo build --release          # binary at target/release/maligno
cargo install --path .         # or install to ~/.cargo/bin (usually on PATH)
```

A static Linux (musl) build for HPC is described in the
[reference](docs/REFERENCE.md#build).

---

## Quick start

Compare the two bundled test PAFs (the same `GRCh38-Gencode-Chr22` transcripts aligned with
differing minimap2 paramters. (`--x splice` vs `--x splice:hq`). This test dataset includes all GENCODE reference transcripts from human chromosome 22 aligned with different `minimap2` alignment parameters. The set of sequenced aligned (ReadIDs) are identical between the two PAF files.

```bash
maligno compare \
  -a test_data/Splice.AlnToHG38.PriAln.paf.gz   --label-a Splice \
  -b test_data/SpliceHQ.AlnToHG38.PriAln.paf.gz --label-b SpliceHQ \
  -o results/ --prefix Splice_vs_SpliceHQ
```

This writes a results directory with the per-set tables and the comparison table:

```
results/Splice_vs_SpliceHQ.Splice.alninfo.tsv.gz     results/Splice_vs_SpliceHQ.SpliceHQ.alninfo.tsv.gz
results/Splice_vs_SpliceHQ.Splice.readinfo.tsv.gz    results/Splice_vs_SpliceHQ.SpliceHQ.readinfo.tsv.gz
results/Splice_vs_SpliceHQ.compare.tsv.gz
```

For the splice-focused view (fewer columns, junction-centric):

```bash
maligno compare --mode junctions \
  -a test_data/Splice.AlnToHG38.PriAln.paf.gz   --label-a Splice \
  -b test_data/SpliceHQ.AlnToHG38.PriAln.paf.gz --label-b SpliceHQ \
  -o results/ --prefix Splice_vs_SpliceHQ
# → results/Splice_vs_SpliceHQ.compare.junctions.tsv.gz
```

All inputs/outputs transparently support gzip (`.gz`) and stdin/stdout (`-`).

---

## The `compare` command

`compare` runs the whole pipeline in three main steps:

1. **Sort** both PAFs by `Query_Name` 
2. **Verify** both PAFs carry the **same read-ID set**. By default it **errors**
   if they differ, reporting how many IDs are shared / only in A / only in B.
3. **Compare** the representative alignment for each readID between sets A and B.
For each read to be compared, the following is done.
  - Select representative alignment for each read (in cases where a read has multiple alignments reported)
    - The best alignment per read is selected by the following alignment scores: `ms` tag, `AS` tag, alignment `MQ`), 
  - The selected alignments for each read are then systematically compared and results are written to a final "compare.tsv" that keep track of each alignments info and differences between them.



| Flag | Purpose |
|------|---------|
| `-a`, `-b` | input PAF for set A / B (`.gz` and `-` ok) |
| `--label-a`, `--label-b` | names used in filenames and per-side column suffixes |
| `-o`, `--outdir` | output directory |
| `-p`, `--prefix` | filename prefix for all outputs |
| `--mode` | `full` (default, 94 cols) or `junctions` (47-col splice view) |
| `--sort-mem` | in-RAM sort buffer per file (default `1G`; `K`/`M`/`G`) |
| `--sort-threads` | sort threads (default `1`) |
| `--presorted` | skip the sort — inputs already hold the same reads in the same order |
| `--allow-id-mismatch` | compare the shared intersection instead of erroring when read-ID sets differ |
| `--no-alninfo`, `--no-readinfo` | skip writing those per-set tables |
| `--skip-find-query-diff` | skip the automatic `find-query-diff` step run after the comparison table is written |
| `--keep-sorted-paf` | keep the intermediate sorted PAFs |


---

## Outputs

A `compare` run writes, under `--outdir`, files prefixed with `--prefix`:

| File | Cols | Contents |
|------|------|----------|
| `{prefix}.{label}.alninfo.tsv.gz` | 35 | **per-alignment** table — one row per PAF alignment (every alignment, per set) |
| `{prefix}.{label}.readinfo.tsv.gz` | 33 | **per-read** table — the chosen best alignment for each read (per set) |
| `{prefix}.compare.tsv.gz` | 94 | the **comparison** table (`full` mode) |
| `{prefix}.compare.junctions.tsv.gz` | 47 | the **comparison** table (`junctions` mode) |
| `{prefix}.compare[.junctions].summary.tsv` | 2 | **aggregate summary statistics** (see below) |
| `{prefix}.query_diff_reads.tsv.gz`, `{prefix}.query_diff_regions.{A,B}.bed.gz`, `{prefix}.query_diff_summary.tsv` | — | **query-different reads + regions** (see below); skip with `--skip-find-query-diff` |

### The comparison table

One row per read, organized in column groups (left to right):

| Group | `full` (94) | `junctions` (47) | What it holds |
|-------|:-----------:|:----------------:|---------------|
| **Join keys** | 1–2 | 1–2 | `Read_Name`, `Read_Len` |
| **Per-side data — A** | 3–33 | 3–17 | the best alignment's stats for set A, each column suffixed `_<label-a>` |
| **Per-side data — B** | 34–64 | 18–32 | the same columns for set B, suffixed `_<label-b>` |
| **Comparison metrics** | 65–90 | 33–43 | A-vs-B differences/ratios: `Strand_Match`, score diffs (`AS_Diff`, `ms_Diff`, …), `seqid_Diff`, coverage/indel/soft-clip diffs (full only), and junction-set counts |
| **Non-overlap objects** | 91–94 | 44–47 | the actual junctions that failed to overlap: `Junctions_OnlyA/B` and `Genomic_Junctions_OnlyA/B` |

**Junctions are compared as sets** of coordinates, deduplicated per side, in two
coordinate systems:

- **Query coordinates** — `N_Matched_Junctions`, `N_Unmatched_Junctions`,
  `N_Junctions_OnlyA`, `N_Junctions_OnlyB`.
- **Genomic coordinates** — the parallel `Genomic_N_*` columns (always emitted;
  computed only when both sides map to the same reference). `chrom` is tracked via
  the per-side `TargetChr_<label>` column, so junctions on different contigs never
  falsely match.

The `junctions` / `genomic_junctions` columns (and the trailing `*_OnlyA/B`
object lists) use a Python tuple-of-tuples format — parse them in Python with
`ast.literal_eval`.

For the **exhaustive column-by-column dictionary**, the genomic-junction format
details, and schema-migration notes, see the [reference](docs/REFERENCE.md#compare-readinfo-and-the-comparison-core).

### Summary statistics

Alongside the comparison table, `compare` writes a small `…summary.tsv` with
predefined aggregate counts (tallied as rows stream, so memory stays constant).
The headline is the **per-read alignment status**, followed by **identity** stats:

| Category | Meaning |
|----------|---------|
| `aligned_both` / `aligned_only_<label>` / `aligned_neither` | how the read's representative alignment maps in each set (an unmapped side is `TargetChr == "*"`) |
| `query_identical` | both sides mapped over the same query span with the **same alignment relative to the read** (identical `cs`). Split into `…_same_strand` and `…_revcomp` (an inverted, opposite-strand match — `cs_A == reverse_complement(cs_B)`) |
| `reference_identical` | query-identical **and** same `TargetChr` + target start/end (same placement on the reference) |
| `present_only_in_<label>_by_id` | reads found in only one set's PAF (0 unless `--allow-id-mismatch`) |

To get the same summary from an existing comparison table (e.g. from the manual
`paf2tables` → `compare-readinfo` workflow), use **`compare-summary`**:

```bash
maligno compare-summary -i AvsB.compare.tsv.gz -o AvsB.compare.summary.tsv
```

It works on either `--mode` (the identity check uses columns present in both).
Full definitions are in the [reference](docs/REFERENCE.md#compare-readinfo-and-the-comparison-core).

### Query-different reads & regions

After writing the comparison table, `compare` also runs **`find-query-diff`** —
it finds every read whose alignment is **not** `query_identical` (using the exact
same definition as above), and reports where those reads cluster on the genome.
This runs by default (always both coordinate sides, gzip'd — not configurable);
skip it with `--skip-find-query-diff`.

Outputs: `{prefix}.query_diff_reads.tsv.gz` (one row per differing read + its
category — `diff_aln_to_both` / `diff_aln_only_<label>`), a merged, `bedtools
merge`-style region table per side (`{prefix}.query_diff_regions.{A,B}.bed.gz` —
`chrom, start, end, n_reads, n_both, n_only_<label>, n_plus, n_minus`), and a
category-tally `{prefix}.query_diff_summary.tsv`.

To run it on an existing comparison table (e.g. from the manual workflow, or with
different `--coord-side`/`--gzip` choices), use the standalone command:

```bash
maligno find-query-diff -i AvsB.compare.tsv.gz --outdir results/ --prefix AvsB
```

Add **`--emit-identical-reads`** to also write `{prefix}.query_identical_reads.tsv`
(`Read_Name` + `query_identical_same_strand`/`query_identical_revcomp`) — the
complement of the diff-reads file. Off by default; not used by `compare`'s
built-in invocation.

---

## Test data

`test_data/` holds two Chr22-scale PAFs (~0.5 MB each) — the same 11,578
transcripts aligned with minimap2 `--x splice` vs `--x splice:hq`.

Run from the repo root; outputs go to `test_data/test_results/` (gitignored):

```bash
# Full comparison (sort → verify read-IDs → per-set tables + comparison table).
maligno compare \
  -a test_data/Splice.AlnToHG38.PriAln.paf.gz   --label-a Splice \
  -b test_data/SpliceHQ.AlnToHG38.PriAln.paf.gz --label-b SpliceHQ \
  -o test_data/test_results/ --prefix Splice_vs_SpliceHQ

# Splice-focused (47-col) view.
maligno compare --mode junctions \
  -a test_data/Splice.AlnToHG38.PriAln.paf.gz   --label-a Splice \
  -b test_data/SpliceHQ.AlnToHG38.PriAln.paf.gz --label-b SpliceHQ \
  -o test_data/test_results/ --prefix Splice_vs_SpliceHQ

# Inspect a comparison header (column number → name).
zcat < test_data/test_results/Splice_vs_SpliceHQ.compare.tsv.gz | head -1 | tr '\t' '\n' | nl

# Sanity-check column counts (expect 35, 33, 94, 47).
for f in test_data/test_results/Splice_vs_SpliceHQ.*.tsv.gz; do
  printf '%s\t' "$f"; zcat < "$f" | awk -F'\t' '{print NF}' | sort -u | paste -sd, -
done

# Count reads whose query-coordinate junctions don't all agree between the two sets.
zcat < test_data/test_results/Splice_vs_SpliceHQ.compare.tsv.gz \
  | awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)c[$i]=i;next} $c["N_Unmatched_Junctions"]>0' | wc -l
```

---

## Extended documentation

The full manual lives in **[`docs/REFERENCE.md`](docs/REFERENCE.md)**:

- Manual building blocks — `paf2tables` (PAF → alninfo/readinfo) and
  `compare-readinfo` (the comparison engine).
- Utilities — `sam2paf` (SAM → PAF) and `utils-readinfo`.
- The complete column dictionary for every table, genomic-junction format,
  and schema-migration notes.
- Sort-order troubleshooting and the overlap diagnostic script.
- Static HPC build and the source layout.

---

## License

See [LICENSE](LICENSE).
