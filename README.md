<picture>
  <source media="(prefers-color-scheme: dark)" srcset="docs/maligno-dark.png">
  <source media="(prefers-color-scheme: light)" srcset="docs/maligno-light.png">
  <img alt="maligno" src="docs/maligno-light.png" width="300">
</picture>

# maligno

**maligno** is a toolkit for systematically comparing two sets of alignments read-by-read.
The primary command, **`maligno compare`**, takes two alignment files ([PAF](https://github.com/lh3/miniasm/blob/master/PAF.md), [SAM/BAM](https://en.wikipedia.org/wiki/SAM_(file_format))) for the same set of reads
and produces a detailed per-read alignment comparison table. This table enables easy comparison of alignment stats (alignment scores, coverage, indels, mismatches,
soft-clipping, and splice-junction agreement) across the aligned sequences.

This means we can use `maligno` to compare alignment results between: <br>
1) Different aligners (`minimap2` vs `bwa mem`) <br>
2) Different parameters of the same aligner (`splice` vs `splice:hq` for transcript alignment with `minimap2`) <br>
3) Different reference genomes (`GRCh38` vs `CHM13`) <br>

With a single alignment comparison table, you can quickly ask:
- "How many reads aligned identically?"
- "Which read IDs differ in their alignment between sets A and B, and how do they specifically differ?"
- "How many reads improved their alignment score by X points?"
- "How do splice junctions differ between sets A and B? For which specific reads?"

In short, `maligno` gives you an framework for detailed alignment comparisons
across different parameters, aligners, and references.

---

## Install

### Installing prebuilt binaries

Prebuilt binaries for linux and mac OS are attached to each
[GitHub release](https://github.com/maxgmarin/maligno/releases).
Each tarball holds the `maligno` binary and this README.

```bash
# Go to the directory where you want to place the maligno binary
cd <directory_on_PATH>

# Download the latest release (Linux x86_64)
wget https://github.com/maxgmarin/maligno/releases/download/v0.42.1/maligno-v0.42.1-x86_64-unknown-linux-musl.tar.gz

# Uncompress the .tar.gz to get the binary
tar -xzf maligno-v0.42.1-x86_64-unknown-linux-musl.tar.gz

# Test that the binary runs
./maligno --help
./maligno --version

# Then add the directory with the binary to your $PATH variable
```

### Building from source

Requires a [Rust toolchain](https://rustup.rs).

```bash
cargo build --release          # binary at target/release/maligno
cargo install --path .         # or install to ~/.cargo/bin (usually on PATH)
```

A static Linux (musl) build for HPC is described in the
[reference](docs/REFERENCE.md#build).

---

## Quick start

Compare the two bundled test PAFs: all GENCODE reference transcripts from human chromosome 22
(`GRCh38-Gencode-Chr22`), aligned with two different minimap2 parameter sets (`-x splice` vs
`-x splice:hq`). Both PAF files contain the same set of aligned sequences (read IDs).

####  Use `maligno compare` to generate detailed comparisons of each sequence's alignment across the two input alignment files

```bash
maligno compare \
  -a test_data/Splice.AlnToHG38.PriAln.paf.gz   --label-a Splice \
  -b test_data/SpliceHQ.AlnToHG38.PriAln.paf.gz --label-b SpliceHQ \
  -o test_data/test_results/ --prefix Splice_vs_SpliceHQ
```

This will write a results directory with the alignment comparison table and
its summary statistics:

```
test_data/test_results/Splice_vs_SpliceHQ.compare.parquet
test_data/test_results/Splice_vs_SpliceHQ.compare.summary.tsv
```

#### Use `maligno toolkit find-aln-diff` to list the reads whose alignments differ, and where they cluster in the coordinate systems of A and B

```bash
maligno toolkit find-aln-diff \
  -i test_data/test_results/Splice_vs_SpliceHQ.compare.parquet \
  -o test_data/test_results/ -p Splice_vs_SpliceHQ
```

```
test_data/test_results/Splice_vs_SpliceHQ.query_aln_diff_reads.tsv.gz
test_data/test_results/Splice_vs_SpliceHQ.query_aln_diff_regions.A.bed.gz
test_data/test_results/Splice_vs_SpliceHQ.query_aln_diff_regions.B.bed.gz
test_data/test_results/Splice_vs_SpliceHQ.query_aln_diff_summary.tsv
```

By default (`--mode any-diff-query`) a read differs when its alignment relative
to the read differs: a different query span or a different `cs` (an exact
reverse complement on the opposite strand counts as the same). With
`--mode any-diff-reference` it also differs when A and B place it at a
different position or strand on the reference.

---

## Usage

### `maligno compare` command

`compare` takes two alignment files, each **PAF, SAM or BAM** (format auto-detected
per file, so mixing e.g. a BAM with a PAF is fine; CRAM is not supported). SAM/BAM
inputs are converted internally with the same code as `maligno sam2paf -p -U`
(primary + supplementary records, unmapped reads kept), so results are identical to
converting them yourself first. Any BAM sort order works.

```bash
maligno compare \
  -a sample.bwamem.sorted.bam  --label-a bwamem \
  -b sample.minibwa.sorted.bam --label-b minibwa \
  -o results/ --prefix bwamem_vs_minibwa
```

It then runs the whole pipeline in four main steps:

1. **Sort** both inputs by `query_name`.
2. **Verify** both inputs carry the **same read-ID set**. By default it **errors**
   if they differ, reporting how many IDs are shared / only in A / only in B.
3. **Select representative alignment for each readID within each read set (A and B)**.
   In cases where a read has multiple alignments reported, the best alignment is
   chosen by alignment score: `ms` tag, then `AS` tag, then mapping quality (`mapq`);
   remaining ties prefer the aligner's primary alignment (`tp:A`), then input order.
   Every value reported for the read describes that one alignment.
4. **Compare the representative alignment across sets A and B**. The selected
   alignments for each read are systematically compared and the results are
   written to a final comparison table ([Parquet](https://www.databricks.com/blog/what-is-parquet) by default) that keeps track of each alignment's info
   and the differences between them.

#### Options

| Flag | Purpose |
|------|---------|
| `-a`/`--aln-a`, `-b`/`--aln-b` | input alignments for set A / B: PAF (`.paf[.gz]`), SAM or BAM |
| `--sam-records-to-convert` | SAM/BAM inputs only: `primary-supp` (default, = `sam2paf -p`), `primary` (= `-P`), or `all`; unmapped reads are always kept. This also decides which alignment types the `n_alns_primsupp`/`n_alns_secondary` counts can see (the default drops secondaries, so `n_alns_secondary` is `0`) |
| `--label-a`, `--label-b` | names for each set, used in filenames and recorded in the comparison table's `label_A` / `label_B` columns |
| `-o`, `--outdir` | output directory |
| `-p`, `--prefix` | filename prefix for all outputs |
| `--sort-mem` | in-RAM sort buffer per file (default `1G`; `K`/`M`/`G`) |
| `--sort-threads` | sort threads (default `1`) |
| `--sort-tmp-dir` | directory for the sort's temporary files (default: `--outdir`) |
| `--presorted` | skip the sorting step; assumes the inputs already hold the same reads in the same order (for SAM/BAM: name-sort both with `samtools sort -n`; coordinate-sorted BAMs are rejected). Sorting differently than maligno can cause small differences in the selected representative alignment |
| `--allow-id-mismatch` | compare the shared intersection instead of erroring when read-ID sets differ |
| `--format` | which serialization(s) of the comparison table: `parquet` (default), `tsv` (gzipped), or `both` |
| `--emit-alninfo`, `--emit-readinfo` | write those per-set tables (off by default) |
| `--keep-sorted-paf` | keep the intermediate sorted PAFs |


---

### Outputs

A `compare` run can write the following files in the user defined output directory:

| File | Cols | Contents |
|------|------|----------|
| `{prefix}.{label}.alninfo.tsv.gz` | 36 | **per-alignment** table: one row per PAF alignment (every alignment, per set); opt-in via `--emit-alninfo` |
| `{prefix}.{label}.readinfo.tsv.gz` | 37 | **per-read** table: the chosen best alignment for each read (per set); opt-in via `--emit-readinfo` |
| `{prefix}.compare.parquet` | 93 | the **comparison** table (unless `--format tsv`) |
| `{prefix}.compare.tsv.gz` | 93 | the same table as gzipped TSV (only with `--format tsv` or `both`) |
| `{prefix}.compare.summary.tsv` | 2 | **aggregate summary statistics** |

Column-by-column definitions of the comparison table are in
[`docs/output-tables/compare.md`](docs/output-tables/compare.md), and the summary
categories (alignment status, query and reference identity) in
[`docs/output-tables/compare-summary.md`](docs/output-tables/compare-summary.md).
`toolkit compare-summary` recomputes the summary from an existing comparison table.

To compare A's and B's splice junctions in **genomic (reference)
coordinates** (each junction with its strand and splice-site motif, per read
and per distinct junction), run **`toolkit genomic-junction-compare`** on the
comparison table ([Parquet](https://www.databricks.com/blog/what-is-parquet) or TSV). It only makes sense when A and B were
aligned to the same reference:

```bash
maligno toolkit genomic-junction-compare \
  -i test_data/test_results/Splice_vs_SpliceHQ.compare.parquet \
  --outdir test_data/test_results/ --prefix Splice_vs_SpliceHQ
```

```
test_data/test_results/Splice_vs_SpliceHQ.genomic_junction_compare.per_read.tsv.gz
test_data/test_results/Splice_vs_SpliceHQ.genomic_junction_compare.per_junction.tsv.gz
test_data/test_results/Splice_vs_SpliceHQ.genomic_junction_compare.summary.tsv
```


## Included test dataset (Annotated Gencode v49 Human Transcripts from Chr22)

`test_data/` holds two Chr22-scale PAFs (~0.5 MB each) — the same 11,578
transcripts aligned with minimap2 `-x splice` vs `-x splice:hq`.

If run from the repo root, the outputs will go to `test_data/test_results/`:

```bash
# Full comparison: the comparison table and its summary, plus the opt-in
# per-set alninfo/readinfo tables.
maligno compare \
  -a test_data/Splice.AlnToHG38.PriAln.paf.gz   --label-a Splice \
  -b test_data/SpliceHQ.AlnToHG38.PriAln.paf.gz --label-b SpliceHQ \
  -o test_data/test_results/ --prefix Splice_vs_SpliceHQ \
  --emit-alninfo --emit-readinfo

# Reconstruct the splice junctions of the reads whose query-space junctions differ.
maligno toolkit query-junction-compare \
  -i test_data/test_results/Splice_vs_SpliceHQ.compare.parquet \
  -o test_data/test_results/ -p Splice_vs_SpliceHQ
```

---

## Extended documentation

A per-column spec for every output file `compare`, `toolkit find-aln-diff`, `toolkit
query-junction-compare` and `toolkit genomic-junction-compare` can write — alninfo,
readinfo, the comparison table, the summary TSV, the find-aln-diff tables,
the query-junction-compare tables and the genomic-junction-compare tables — lives
in **[`docs/output-tables/`](docs/output-tables/)**.

The full manual lives in **[`docs/REFERENCE.md`](docs/REFERENCE.md)**:

- Every command, including `sam2paf` and `toolkit`'s individual building blocks.
- How the representative alignment is chosen and how reads are classified.
- Junction comparison, TSV vs [Parquet](https://www.databricks.com/blog/what-is-parquet), and troubleshooting.
- The static HPC build.

---

## FAQ

### Preprocessing BAMs produced by the STAR aligner (adding the `MD` tag)

**STAR does not emit a `cs` tag or `MD` on its own, so its BAMs need preprocessing with `samtools calmd` before conversion from SAM to PAF:**


When `maligno` is parsing SAM/BAM files it uses an aligner-supplied `cs:Z:`
SAM tag as-is when one is present. When there's no `cs` tag, `maligno` falls back to deriving one from **`CIGAR` +
`MD` + `SEQ`**. A mapped record with neither a `cs` tag nor a usable `MD` tag
(or with `MD` but `SEQ` = `*`, common for secondary alignments) is an **error**. If the aligner (like STAR) does NOT emit an `MD` tag it will need to be generated with `samtools calmd`.

```bash
samtools faidx reference.fasta   # only if reference.fasta.fai doesn't already exist
samtools calmd -b star_output.bam reference.fasta > star_output.calmd.bam
```

- `reference.fasta` must be the **exact same reference** STAR aligned
  against (same contig names and sequence) — `calmd` looks up the reference
  base at each alignment position to recompute `MD`, and a mismatched
  reference silently produces wrong tags rather than an error.
- Avoid `-e` — it rewrites reference-matching bases in `SEQ` as `=`, which
  would corrupt the base calls `sam2paf`'s `cs` derivation needs to read back
  out.
- `calmd` doesn't require any particular sort order (it recomputes each
  record independently), but isn't multithreaded either — budget time for it
  on large BAMs.
- Any aligner that doesn't emit a `cs` tag hits this same requirement for
  `MD` — e.g. `bwa mem` always emits `MD` by default, but `minibwa map` needs
  an explicit `-b MD` flag to emit it at all. (or set `minibwa map` to emit the cs tag with `-b cs`)
- Check that your aligner's BAM actually carries `cs` or `MD` before running `sam2paf`
  or `compare` on it.
- `compare` compares `cs` strings as text, so both sides should use the same `cs`
  form. SAM/BAM conversion produces the short form (minimap2's default); don't pair
  it with a PAF made with `minimap2 --cs=long`.

---

