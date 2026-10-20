# `query_junction_compare.all_junctions.{A,B}.tsv`

`{prefix}.query_junction_compare.all_junctions.A.tsv.gz` /
`{prefix}.query_junction_compare.all_junctions.B.tsv.gz` — every distinct
genomic splice junction on that side, across **every** read of the comparison
table, with how many reads carry it. Written by `toolkit query-junction-compare`
with `--emit-all-junctions`, one file per side.

Per-side tables don't need A and B to share a reference: with GRCh38 vs a
personal assembly, for example, each file still lists that side's junctions and
counts. (For junctions shared between A and B on the same reference, see
`genomic-junction-compare`'s [per-junction table](genomic-junction-compare-per-junction.md).)

A junction is `(target_name, intron_start, intron_end, intron_strand)`, and
rows are sorted by it. The columns are the same as the
[unmatched tables](query-junction-compare-unmatched.md).

7 columns.

| # | Column | Description |
|---|---|---|
| 1 | `target_name` | Genomic contig |
| 2 | `intron_start` | Intron start (0-based half-open) |
| 3 | `intron_end` | Intron end |
| 4 | `intron_strand` | The intron's strand on the reference (this side's `genomic_intron_strand`; `*` when unknown) |
| 5 | `motif` | Splice-site motif, donor-acceptor in the junction's orientation (e.g. `gt-ag`); see Notes for which occurrence it comes from |
| 6 | `n_reads_query_junctions_different` | Differing reads where this junction is unmatched on the other side; `0` for a junction that never is |
| 7 | `n_reads_with_junction_total` | Reads, across every row of the comparison table, with this junction on this side |

**Invariants**
- The unmatched table for the same side is exactly the rows of this table
  with `n_reads_query_junctions_different > 0`.
- The number of rows equals the summary's `n_unique_junctions_A` /
  `n_unique_junctions_B`.

**Notes**
- `motif`: for a junction with unmatched occurrences, the same value as in the
  unmatched table (the first among the unmatched reads, preferring real
  letters over `nn-nn`); otherwise the first among all reads with the
  junction, preferring real letters.
- Counts are reads, using each read's representative alignment only, as in
  the rest of the comparison table.
