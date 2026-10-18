<picture>
  <source media="(prefers-color-scheme: dark)" srcset="docs/maligno-dark.png">
  <source media="(prefers-color-scheme: light)" srcset="docs/maligno-light.png">
  <img alt="maligno" src="docs/maligno-light.png" width="300">
</picture>

# maligno

**maligno** is a toolkit for systematically comparing two sets of alignments read-by-read.
The primary command, **`maligno compare`**, takes two alignment files ([PAF](https://github.com/lh3/miniasm/blob/master/PAF.md), SAM or BAM) for the same set of reads
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
the differing reads and the genomic regions where they cluster :

```
test_data/test_results/Splice_vs_SpliceHQ.compare.parquet
test_data/test_results/Splice_vs_SpliceHQ.compare.summary.tsv
test_data/test_results/Splice_vs_SpliceHQ.query_diff_reads.tsv.gz
test_data/test_results/Splice_vs_SpliceHQ.query_diff_regions.A.bed.gz
test_data/test_results/Splice_vs_SpliceHQ.query_diff_regions.B.bed.gz
```

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
   written to a final comparison table (Parquet by default) that keeps track of each alignment's info
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
| `--skip-find-aln-diff` | don't also emit `find-aln-diff`'s default-mode output (differing reads + region tables) |


---

### Outputs

A `compare` run can write the following files in the user defined output directory:

| File | Cols | Contents |
|------|------|----------|
| `{prefix}.{label}.alninfo.tsv.gz` | 36 | **per-alignment** table: one row per PAF alignment (every alignment, per set); opt-in via `--emit-alninfo` |
| `{prefix}.{label}.readinfo.tsv.gz` | 37 | **per-read** table: the chosen best alignment for each read (per set); opt-in via `--emit-readinfo` |
| `{prefix}.compare.parquet` | 93 | the **comparison** table (unless `--format tsv`) |
| `{prefix}.compare.tsv.gz` | 93 | the same table as gzipped TSV (only with `--format tsv` or `both`) |
| `{prefix}.compare.summary.tsv` | 2 | **aggregate summary statistics** (see below) |
| `{prefix}.query_diff_reads.tsv.gz` | 9 | table of all reads with difference in alignment between set A and B  |
| `{prefix}.query_diff_regions.{A,B}.bed.gz` | 10 | genomic regions where differing reads cluster, per set |

Separately, **`toolkit query-junction-compare`** takes an existing
`compare.parquet` (Parquet input only), selects reads whose
**query-space** splice junctions differ, and reconstructs those junctions
per side, paired in both query and genomic coordinate space, plus a rollup
of which specific junctions are unsupported by the other side and how many
reads *total* (across the whole table) carry each one:

```bash
maligno toolkit query-junction-compare \
  -i test_data/test_results/Splice_vs_SpliceHQ.compare.parquet \
  --outdir test_data/test_results/ --prefix Splice_vs_SpliceHQ
```

```
test_data/test_results/Splice_vs_SpliceHQ.query_junction_compare.summary.tsv
test_data/test_results/Splice_vs_SpliceHQ.query_junction_compare.per_read_per_junc_info.tsv.gz
test_data/test_results/Splice_vs_SpliceHQ.query_junction_compare.unmatched_junctions.A.tsv.gz
test_data/test_results/Splice_vs_SpliceHQ.query_junction_compare.unmatched_junctions.B.tsv.gz
```

To compare A's and B's splice junctions in **genomic (reference)
coordinates** (each junction with its strand and splice-site motif, per read
and per distinct junction), run **`toolkit genomic-junction-compare`** on the
comparison table (Parquet or TSV). It only makes sense when A and B were
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
# Full comparison (sort → verify read-IDs → comparison table), plus the
# opt-in per-set alninfo/readinfo tables (--emit-alninfo/--emit-readinfo) so
# this walkthrough can sanity-check all three table shapes below. The
# comparison table is Parquet-only by default; --format both also writes the
# gzipped TSV that the zcat/awk checks below read.
maligno compare \
  -a test_data/Splice.AlnToHG38.PriAln.paf.gz   --label-a Splice \
  -b test_data/SpliceHQ.AlnToHG38.PriAln.paf.gz --label-b SpliceHQ \
  -o test_data/test_results/ --prefix Splice_vs_SpliceHQ \
  --emit-alninfo --emit-readinfo --format both

# Inspect a comparison header (column number → name).
zcat < test_data/test_results/Splice_vs_SpliceHQ.compare.tsv.gz | head -1 | tr '\t' '\n' | nl

# Sanity-check column counts (expect 36 for alninfo, 37 for readinfo, 93 for compare, 9 for query_diff_reads).
for f in test_data/test_results/Splice_vs_SpliceHQ.*.tsv.gz; do
  printf '%s\t' "$f"; zcat < "$f" | awk -F'\t' '{print NF}' | sort -u | paste -sd, -
done

# Count reads whose query-coordinate junctions don't exactly agree between the two sets.
zcat < test_data/test_results/Splice_vs_SpliceHQ.compare.tsv.gz \
  | awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)c[$i]=i;next} $c["n_query_junctions_not_shared"]>0' | wc -l
```

---

## Extended documentation

A per-column spec for every output file `compare`, `toolkit
query-junction-compare` and `toolkit genomic-junction-compare` can write — alninfo,
readinfo, the comparison table, the summary TSV, the diff-reads/regions tables,
the query-junction-compare tables and the genomic-junction-compare tables — lives
in **[`docs/output-tables/`](docs/output-tables/)**.

The full manual lives in **[`docs/REFERENCE.md`](docs/REFERENCE.md)**:

- Every command, including `sam2paf` and `toolkit`'s individual building blocks.
- How the representative alignment is chosen and how reads are classified.
- Junction comparison, TSV vs Parquet, and troubleshooting.
- The static HPC build.

---

## Extra details for the output tables


#### The alignment comparison table (`{prefix}.compare.parquet` or `{prefix}.compare.tsv.gz`)

One row per read, 93 columns: the representative alignment's stats for set A and
for set B (suffixed `_A` / `_B`), plus a block of columns comparing how those two
alignments differ (score, coverage, indels, soft-clipping, and junction agreement
in query space, i.e. positions in the read). Each side also carries the
transcript strand from minimap2's `ts:A` or STAR's `XS:A`, twice: relative to the
read (`query_transcript_strand`) and on the reference (`genomic_intron_strand`)
(see [REFERENCE.md](docs/REFERENCE.md#transcript-strand)).

For the full column-by-column layout and dictionary, see
[`docs/output-tables/compare.md`](docs/output-tables/compare.md).

---

#### Summary statistics in `{prefix}.compare.summary.tsv`

Alongside the comparison table, `compare` writes a small `…summary.tsv` with
predefined aggregate counts. 

The headline is the **per-read alignment status**, followed by **identity** stats:

| Category | Meaning |
|----------|---------|
| `label_A` / `label_B` | which dataset each side is (the `--label-a` / `--label-b` values), so the summary is self-describing |
| `reads_compared` | reads written to the comparison table |
| `aligned_both` / `aligned_only_A` / `aligned_only_B` / `aligned_neither` | how the read's representative alignment maps in each set (an unmapped side is `target_name == "*"`) |
| `query_identical` / `query_not_identical` | `query_identical`: either both sides mapped over the same query span with the **same alignment relative to the read** (the same `cs` on the same strand, or reverse-complement `cs` on opposite strands; intron motif letters ignored), or neither side mapped at all (`aligned_neither`: both aligners agreeing a read doesn't map is agreement, not disagreement); `query_not_identical` is its complement among `aligned_both` reads only |
| `query_identical_same_strand` / `query_identical_revcomp` | the both-mapped `query_identical` reads, split by same-strand vs reverse-complement match |
| `query_junctions_identical` / `query_junctions_not_identical` | same, but comparing only the **query-space splice-junction set** (ignores mismatches/indels/soft-clips) |
| `ref_same_position_same_aln` | reference-identical: same `target_name` + `strand` + `target_start` (same genomic position) **and** same `cs` |
| `ref_same_position_diff_aln` | same reference position, different alignment (e.g. a different indel placement at the same site) |
| `ref_diff_position_same_aln` | same alignment, different reference position |
| `ref_diff_position_diff_aln` | both reference position and alignment differ |
| `present_only_in_A_by_id` / `present_only_in_B_by_id` | sequences found in only one set's PAF (will be 0 unless `--allow-id-mismatch` is used) |

To get the same summary from an existing comparison table, use
**`toolkit compare-summary`**:

```bash
maligno toolkit compare-summary -i AvsB.compare.parquet -o AvsB.compare.summary.tsv
```

Full definitions are in the [reference](docs/REFERENCE.md#classification).

---

## Preprocessing STAR BAMs for SAM/BAM input (adding the `MD` tag)

`maligno sam2paf` (and `compare` on SAM/BAM input) uses an aligner-supplied `cs:Z:`
SAM tag as-is when one is present (e.g. minimap2 emits it natively).
When there's no `cs` tag, it falls back to deriving one from **`CIGAR` +
`MD` + `SEQ`**. A mapped record with neither a `cs` tag nor a usable `MD` tag
(or with `MD` but `SEQ` = `*`, common for secondary alignments) is an **error**.

If the aligner does NOT emit an `MD` tag it will need to be generated with `samtools calmd`.

STAR does not emit a `cs` tag or `MD` on its own, so its BAMs need one preprocessing pass with `samtools calmd` before conversion from SAM to PAF:

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

