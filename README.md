<picture>
  <source media="(prefers-color-scheme: dark)" srcset="docs/maligno-dark.png">
  <source media="(prefers-color-scheme: light)" srcset="docs/maligno-light.png">
  <img alt="maligno" src="docs/maligno-light.png" width="300">
</picture>

# maligno

**maligno** is a toolkit for systematically comparing two sets of alignments read-by-read.
The primary command, **`maligno compare`**, takes two alignment files ([PAF](https://github.com/lh3/miniasm/blob/master/PAF.md)) for the same set of reads
and produces a detailed per-read alignment comparison table. This table enables easy comparison of alignment stats (alignment scores, coverage, indels, mismatches,
soft-clipping, and splice-junction agreement) across the analyzed reads.

This makes it easy to compare alignment results across conditions such as:
1) Different aligners (`minimap2` vs `bwa mem`)
2) Different parameters of the same aligner (`splice` vs `splice:hq` for transcript alignment with `minimap2`)
3) Different reference genomes (`GRCh38` vs `CHM13`)

With a single comparison table, you can quickly ask:
- "How many reads aligned identically?"
- "Which read IDs differ in their alignment between set A and B, and how do they differ?"
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

### Step 1: Use `maligno compare` to generate detailed comparisons of each sequence's alignment across the two input alignment files

```bash
maligno compare \
  -a test_data/Splice.AlnToHG38.PriAln.paf.gz   --label-a Splice \
  -b test_data/SpliceHQ.AlnToHG38.PriAln.paf.gz --label-b SpliceHQ \
  -o results/ --prefix Splice_vs_SpliceHQ
```

This writes a results directory with the comparison table and — by default —
the differing reads and the genomic regions where they cluster
(`find-aln-diff`'s default-mode output, computed in the same pass; see
Step 2):

```
results/Splice_vs_SpliceHQ.compare.tsv.gz
results/Splice_vs_SpliceHQ.compare.parquet
results/Splice_vs_SpliceHQ.compare.summary.tsv
results/Splice_vs_SpliceHQ.query_diff_reads.tsv.gz
results/Splice_vs_SpliceHQ.query_diff_regions.A.bed.gz
results/Splice_vs_SpliceHQ.query_diff_regions.B.bed.gz
```

The per-set `alninfo`/`readinfo` tables are opt-in (`--emit-alninfo`,
`--emit-readinfo`) — pass them if you want the per-alignment or per-read
detail tables alongside the comparison output; they roughly double the run's
disk footprint and add meaningfully to its runtime, so they're off by default.

### Step 2: Use `maligno find-aln-diff` for reference-space or junctions-only differences

Step 1 already found the reads that align differently in query space (the
default, and most common, comparison mode) — no separate command needed for
that. Run **`find-aln-diff`** standalone only when you need a different
`--space`/`--compare-by` combination (e.g. reference-space differences, or
differences restricted to the splice-junction set), or want to regenerate the
diff outputs from an existing comparison table without re-running `compare`:

```bash
maligno find-aln-diff \
  -i results/Splice_vs_SpliceHQ.compare.tsv.gz \
  --space reference \
  --outdir results/ --prefix Splice_vs_SpliceHQ
```

---

## The `compare` command

`compare` runs the whole pipeline in four main steps:

1. **Sort** both PAFs by `Query_Name`.
2. **Verify** both PAFs carry the **same read-ID set**. By default it **errors**
   if they differ, reporting how many IDs are shared / only in A / only in B.
3. **Select representative alignment for each readID within each read set (A and B)** —
   in cases where a read has multiple alignments reported, the best alignment is
   chosen by alignment score: `ms` tag, then `AS` tag, then alignment `MQ`.
4. **Compare the representative alignment across sets A and B** — the selected
   alignments for each read are systematically compared and the results are
   written to a final `compare.tsv` that keeps track of each alignment's info
   and the differences between them.

### Options

| Flag | Purpose |
|------|---------|
| `-a`, `-b` | input PAF for set A / B (`.gz` and `-` ok) |
| `--label-a`, `--label-b` | names for each set, used in filenames and recorded in the comparison table's `Label_A` / `Label_B` columns |
| `-o`, `--outdir` | output directory |
| `-p`, `--prefix` | filename prefix for all outputs |
| `--sort-mem` | in-RAM sort buffer per file (default `1G`; `K`/`M`/`G`) |
| `--sort-threads` | sort threads (default `1`) |
| `--presorted` | skip the sort — inputs already hold the same reads in the same order |
| `--allow-id-mismatch` | compare the shared intersection instead of erroring when read-ID sets differ |
| `--format` | which serialization(s) of the comparison table: `tsv`, `parquet`, or `both` (default) |
| `--emit-alninfo`, `--emit-readinfo` | write those per-set tables (off by default) |
| `--keep-sorted-paf` | keep the intermediate sorted PAFs |
| `--skip-find-aln-diff` | don't also emit `find-aln-diff`'s default-mode output (differing reads + region tables) |


---

## Outputs

A `compare` run writes, under `--outdir`, files prefixed with `--prefix`:

| File | Cols | Contents |
|------|------|----------|
| `{prefix}.{label}.alninfo.tsv.gz` | 35 | **per-alignment** table — one row per PAF alignment (every alignment, per set); opt-in via `--emit-alninfo` |
| `{prefix}.{label}.readinfo.tsv.gz` | 33 | **per-read** table — the chosen best alignment for each read (per set); opt-in via `--emit-readinfo` |
| `{prefix}.compare.tsv.gz` | 96 | the **comparison** table (unless `--format parquet`) |
| `{prefix}.compare.parquet` | 96 | the same table as Parquet (unless `--format tsv`) |
| `{prefix}.compare.summary.tsv` | 2 | **aggregate summary statistics** (see below) |
| `{prefix}.query_diff_reads.tsv.gz` | 10 | differing reads (default `find-aln-diff` mode — see below; unless `--skip-find-aln-diff`) |
| `{prefix}.query_diff_regions.{A,B}.bed.gz` | 8 | genomic regions where differing reads cluster, per set (unless `--skip-find-aln-diff`) |

### The comparison table

One row per read, 96 columns: the representative alignment's stats for set A and
for set B (suffixed `_A` / `_B`), plus a block of columns comparing how those two
alignments differ (score, coverage, indels, soft-clipping, and junction agreement
in both query and genomic space).

For the full column-by-column layout and dictionary, junction/format details, and
schema-migration notes, see [`docs/COMPARE_TABLE.md`](docs/COMPARE_TABLE.md).

### Summary statistics

Alongside the comparison table, `compare` writes a small `…summary.tsv` with
predefined aggregate counts (tallied as rows stream, so memory stays constant).
This is the **same schema** written by `compare-pipeline summary` and by
standalone `find-aln-diff`'s own summary output (see below) — one set of
category names shared across all three. The headline is the **per-read
alignment status**, followed by **identity** stats:

| Category | Meaning |
|----------|---------|
| `label_A` / `label_B` | which dataset each side is (the `--label-a` / `--label-b` values), so the summary is self-describing |
| `aligned_both` / `aligned_only_A` / `aligned_only_B` / `aligned_neither` | how the read's representative alignment maps in each set (an unmapped side is `TargetChr == "*"`) |
| `query_identical` / `query_not_identical` | both sides mapped over the same query span with the **same alignment relative to the read** (identical `cs` tag operations), and its complement among `aligned_both` reads |
| `query_junctions_identical` / `query_junctions_not_identical` | same, but comparing only the **query-space splice-junction set** — a looser criterion than `query_identical` (ignores mismatches/indels/soft-clips) |
| `ref_same_position_same_aln` | reference-identical: same `TargetChr` + `Strand` + `Target_Start` (same genomic position) **and** same `cs` |
| `ref_same_position_diff_aln` | same reference position, different alignment (e.g. a different indel placement at the same site) |
| `ref_diff_position_same_aln` | same alignment, different reference position |
| `ref_diff_position_diff_aln` | both reference position and alignment differ |
| `ref_same_position_same_junctions` / `ref_same_position_diff_junctions` | among same-position reads, whether the **genomic-coordinate splice-junction set** also matches |
| `present_only_in_A_by_id` / `present_only_in_B_by_id` | sequences found in only one set's PAF (will be 0 unless `--allow-id-mismatch` is used) |

To get the same summary from an existing comparison table, use
**`compare-pipeline summary`**:

```bash
maligno compare-pipeline summary -i AvsB.compare.tsv.gz -o AvsB.compare.summary.tsv
```

Full definitions are in the [reference](docs/REFERENCE.md#compare-pipeline-merge-readinfo-and-the-comparison-core).

### Finding reads with differing alignments (`maligno find-aln-diff`)

**`find-aln-diff`** reads a comparison table and finds every read whose
alignment differs between A and B, then reports where those reads cluster on
the genome.

`compare` already runs this **by default**, at its default settings
(`--space query --compare-by all`), in the same pass that builds the
comparison table — see [Step 2](#step-2-use-maligno-find-aln-diff-for-reference-space-or-junctions-only-differences)
above (`--skip-find-aln-diff` opts out). Run `find-aln-diff` standalone when
you need a different `--space`/`--compare-by` combination, or want to
regenerate these outputs from an existing comparison table without
re-running `compare`.

`--space` selects the coordinate space: `query` (default) compares each side's
alignment relative to the read; `reference` compares each side's placement on
the reference genome instead — only meaningful when both input alignment sets
were aligned to the same reference genome.

```bash
maligno find-aln-diff -i results/AvsB.compare.tsv.gz \
  --outdir results/ --prefix AvsB
```

#### Outputs of `maligno find-aln-diff`:

| File | Contents |
|------|----------|
| `{prefix}.query_diff_reads.tsv.gz` | one row per differing read: `Read_Name`, `outcome` (category — `diff_aln_to_both` / `diff_aln_only_A` / `diff_aln_only_B`), plus 8 classification booleans (`1`/`0`) |
| `{prefix}.query_diff_regions.{A,B}.bed.gz` | Tables (BED format) of all regions with differing alignments. A table is generate for each alignment set (A and B). The columns are `chrom, start, end, n_reads, n_both, n_only_A/n_only_B, n_plus, n_minus` |
| `{prefix}.query_diff_summary.tsv` | the **same summary schema** as `compare.summary.tsv` (see [Summary statistics](#summary-statistics) above), with `space`/`compare_by` provenance rows prepended — its counters are mode-independent; only which reads land in the other two files above is mode-dependent |

The 8 classification booleans in `query_diff_reads.tsv.gz` are: `query_identical_same_strand`, `query_identical_revcomp`, `query_junctions_identical`, `ref_same_position_same_aln`, `ref_same_position_diff_aln`, `ref_diff_position_same_aln`, `ref_diff_position_diff_aln`, `ref_same_position_same_junctions` — so a single run can show, e.g., a read that's query-different but reference-identical, without a second run in the other `--space`.



**`--compare-by`** selects what counts as a difference in alignment between the two sets:

- `all` (default) — compare the full `cs` tag operations. 
- `junctions` — compare only the **query-space splice-junction set**; reads with
  identical junctions but differing mismatches/indels/soft-clips count as the
  **same**. Reads aligned in only one set are still reported (they have no
  junctions to compare on the missing side). In this mode the outputs gain a
  `.junctions` filename segment (e.g. `{prefix}.query_diff_reads.junctions.tsv.gz`),
  so a `junctions` run never clobbers an `all` run at the same prefix. The
  junctions-different read set is always a subset of the `all`-different set.
  `--compare-by` is a standalone-only option — `compare` always uses `all`.

---

## Test data

`test_data/` holds two Chr22-scale PAFs (~0.5 MB each) — the same 11,578
transcripts aligned with minimap2 `--x splice` vs `--x splice:hq`.

If run from the repo root, the outputs will go to `test_data/test_results/`:

```bash
# Full comparison (sort → verify read-IDs → comparison table), plus the
# opt-in per-set alninfo/readinfo tables (--emit-alninfo/--emit-readinfo) so
# this walkthrough can sanity-check all three table shapes below.
maligno compare \
  -a test_data/Splice.AlnToHG38.PriAln.paf.gz   --label-a Splice \
  -b test_data/SpliceHQ.AlnToHG38.PriAln.paf.gz --label-b SpliceHQ \
  -o test_data/test_results/ --prefix Splice_vs_SpliceHQ \
  --emit-alninfo --emit-readinfo

# Inspect a comparison header (column number → name).
zcat < test_data/test_results/Splice_vs_SpliceHQ.compare.tsv.gz | head -1 | tr '\t' '\n' | nl

# Sanity-check column counts (expect 35, 33, 96, plus 10 for the default
# query_diff_reads.tsv.gz — unless --skip-find-aln-diff was passed).
for f in test_data/test_results/Splice_vs_SpliceHQ.*.tsv.gz; do
  printf '%s\t' "$f"; zcat < "$f" | awk -F'\t' '{print NF}' | sort -u | paste -sd, -
done

# Count reads whose query-coordinate junctions don't exactly agree between the two sets.
zcat < test_data/test_results/Splice_vs_SpliceHQ.compare.tsv.gz \
  | awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)c[$i]=i;next} $c["N_Unmatched_Junctions"]>0' | wc -l
```

---

## Extended documentation

The comparison table's column-by-column format lives in
**[`docs/COMPARE_TABLE.md`](docs/COMPARE_TABLE.md)**.

The full manual lives in **[`docs/REFERENCE.md`](docs/REFERENCE.md)**:

- `compare-pipeline`'s manual building blocks — `paf2tables` (PAF → alninfo/readinfo),
  `merge-readinfo` (the comparison engine), and `summary`.
- Utilities — `sam2paf` (SAM → PAF).
- The complete column dictionary for every table, genomic-junction format,
  and schema-migration notes.
- Sort-order troubleshooting and the overlap diagnostic script.
- Static HPC build and the source layout.

---

## License

See [LICENSE](LICENSE).
