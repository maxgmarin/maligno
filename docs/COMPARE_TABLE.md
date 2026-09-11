# The comparison table

`compare` (and `compare-toolkit merge-readinfo`) write one row per read, 96
columns, organized in column groups (left to right):

| Group | Cols | What it holds |
|-------|:----:|---------------|
| **Join keys** | 1–2 | `Read_Name`, `Read_Len` |
| **Set labels** | 3–4 | `Label_A`, `Label_B` — the `--label-a` / `--label-b` values, repeated on every row |
| **Per-side data — A** | 5–35 | the best alignment's stats for set A, each column suffixed `_A` |
| **Per-side data — B** | 36–66 | the same columns for set B, suffixed `_B` |
| **Comparison metrics** | 67–92 | A-vs-B differences/ratios: `Strand_Match`, `seqid_Diff`, coverage/length diffs, score diffs (`AS_Diff`, `ms_Diff`, …), indel/soft-clip diffs, and junction-set counts in both query and genomic space |
| **Non-overlap objects** | 93–96 | the actual junctions that failed to overlap: `Junctions_OnlyA/B` and `Genomic_Junctions_OnlyA/B` |

Within each per-side block the columns are grouped by topic — locus and span,
alignment selection and score, identity and coverage, junction counts, cs-derived
event counts, then the three long strings (`junctions`, `genomic_junctions`, `cs`)
last. Inspect the exact layout of any table with
`gzip -dc … | head -1 | tr '\t' '\n' | nl`.

Per-side columns always use the **fixed** `_A` / `_B` suffixes — never the
dataset label — so column names are identical for every comparison and stay
unambiguous even when a label itself contains an underscore. Which dataset each
side *is* comes from the `Label_A` / `Label_B` columns, carried on every row so
that any row subset of the table remains self-describing.

## Junction comparison

**Junctions are compared as sets** of coordinates, deduplicated per side, in two
coordinate systems:

- **Query coordinates** — `N_Matched_Junctions`, `N_Unmatched_Junctions`,
  `N_Junctions_OnlyA`, `N_Junctions_OnlyB`.
- **Genomic coordinates** — the parallel `Genomic_N_*` columns (always emitted;
  computed only when both sides map to the same reference). `chrom` is tracked via
  the per-side `TargetChr_A` / `TargetChr_B` column, so junctions on different
  contigs never falsely match.

The `junctions` / `genomic_junctions` columns (and the trailing `*_OnlyA/B`
object lists) use a Python tuple-of-tuples format — parse them in Python with
`ast.literal_eval`.

## Nulls

Nulls mean **undefined**, not zero: an unmapped side has a null `cs` and a null
`seqid_Max`, and a ratio over a zero denominator is null. Values that mean something
are kept — `TargetChr` stays `*` for unmapped, an empty junction set stays `()`, and
a real zero stays `0`.

For the exhaustive column-by-column dictionary, the genomic-junction format
details, and schema-migration notes, see the
[reference](REFERENCE.md#compare-toolkit-merge-readinfo-and-the-comparison-core).
