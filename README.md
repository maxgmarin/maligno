<picture>
  <source media="(prefers-color-scheme: dark)" srcset="docs/maligno-dark.png">
  <source media="(prefers-color-scheme: light)" srcset="docs/maligno-light.png">
  <img alt="maligno" src="docs/maligno-light.png" width="300">
</picture>

# maligno

**maligno** is a toolkit for comparing two sets of alignments read-by-read. The
headline command, **`maligno compare`**, takes two PAF files for the same reads
(e.g. the same transcripts aligned with two parameter sets, or against two
references) and produces a per-read comparison table — scores, coverage, indels,
and splice-junction agreement — in a single command.

It owns its preconditions: it sorts the inputs, verifies they describe the same
reads, and does the comparison in one pass, so you can't accidentally mis-compare
unsorted or mismatched files.

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

Compare the two bundled test PAFs (the same Chr22 transcripts aligned with
minimap2 `--x splice` vs `--x splice:hq`):

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

`compare` runs the whole pipeline in three steps:

1. **Sort** both PAFs by `Query_Name` with an in-process external sort (bounded
   memory — buffers up to `--sort-mem`, spilling to temp files), guaranteeing the
   grouping and matching order rather than trusting the input.
2. **Verify** both PAFs carry the **same read-ID set**. By default it **errors**
   if they differ, reporting how many IDs are shared / only in A / only in B.
3. **Compare** in a single pass: it picks the best alignment per read (by `ms`,
   then `AS`, then `MQ`), then merge-joins the two sides to write the comparison
   table, teeing out the per-set tables as it goes.

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
| `--keep-sorted-paf` | keep the intermediate sorted PAFs |

> Prefer to drive each step yourself (`paf2tables` → `compare-readinfo`)? The
> manual building blocks produce the identical comparison table and are
> documented in the [reference](docs/REFERENCE.md).

---

## Outputs

A `compare` run writes, under `--outdir`, files prefixed with `--prefix`:

| File | Cols | Contents |
|------|------|----------|
| `{prefix}.{label}.alninfo.tsv.gz` | 35 | **per-alignment** table — one row per PAF alignment (every alignment, per set) |
| `{prefix}.{label}.readinfo.tsv.gz` | 33 | **per-read** table — the chosen best alignment for each read (per set) |
| `{prefix}.compare.tsv.gz` | 94 | the **comparison** table (`full` mode) |
| `{prefix}.compare.junctions.tsv.gz` | 47 | the **comparison** table (`junctions` mode) |

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

---

## Test data

`test_data/` holds two Chr22-scale PAFs (~0.5 MB each) — the same 11,578
transcripts aligned with minimap2 `--x splice` vs `--x splice:hq`.

```bash
# Full comparison (sort → verify read-IDs → per-set tables + comparison table).
maligno compare \
  -a test_data/Splice.AlnToHG38.PriAln.paf.gz   --label-a Splice \
  -b test_data/SpliceHQ.AlnToHG38.PriAln.paf.gz --label-b SpliceHQ \
  -o results/ --prefix Splice_vs_SpliceHQ

# Splice-focused (47-col) view.
maligno compare --mode junctions \
  -a test_data/Splice.AlnToHG38.PriAln.paf.gz   --label-a Splice \
  -b test_data/SpliceHQ.AlnToHG38.PriAln.paf.gz --label-b SpliceHQ \
  -o results/ --prefix Splice_vs_SpliceHQ

# Inspect a comparison header (column number → name).
zcat < results/Splice_vs_SpliceHQ.compare.tsv.gz | head -1 | tr '\t' '\n' | nl

# Sanity-check column counts (expect 35, 33, 94, 47).
for f in results/Splice_vs_SpliceHQ.*.tsv.gz; do
  printf '%s\t' "$f"; zcat < "$f" | awk -F'\t' '{print NF}' | sort -u | paste -sd, -
done

# Count reads whose query-coordinate junctions don't all agree between the two sets.
zcat < results/Splice_vs_SpliceHQ.compare.tsv.gz \
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
