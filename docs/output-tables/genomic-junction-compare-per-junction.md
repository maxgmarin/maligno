# `genomic_junction_compare.per_junction.tsv.gz`

`{prefix}.genomic_junction_compare.per_junction.tsv.gz` — written by
[`toolkit genomic-junction-compare`](../REFERENCE.md#toolkit-genomic-junction-compare).
One row per **distinct genomic junction** seen on either side, across the
whole comparison table, sorted by `target_name`, `intron_start`,
`intron_end`, `intron_strand`.

11 columns.

| # | Column | Description |
|---|---|---|
| 1 | `target_name` | Reference/target name |
| 2 | `intron_start` | Intron start, 0-based |
| 3 | `intron_end` | Intron end, half-open |
| 4 | `intron_strand` | Strand of the junction on the reference (`+`, `-`, or `*` unknown) |
| 5 | `n_reads_A` | Number of reads whose A alignment has this junction |
| 6 | `n_reads_B` | Number of reads whose B alignment has this junction |
| 7 | `n_reads_both` | Number of reads whose A **and** B alignments both have this junction: the two alignments of the read agree on it |
| 8 | `n_reads_only_A` | Reads with the junction in A's alignment but not B's, `n_reads_A − n_reads_both` |
| 9 | `n_reads_only_B` | Reads with the junction in B's alignment but not A's, `n_reads_B − n_reads_both` |
| 10 | `motif_A` | Splice-site motif reported by A, in the junction's orientation (e.g. `gt-ag`); `nn-nn` if A didn't report the bases; `.` if no A read has the junction |
| 11 | `motif_B` | The same for B |

**Notes**
- A junction is `(target_name, intron_start, intron_end, intron_strand)`: the
  same coordinates on `+` and `-` (or `*`) are separate rows.
- `n_reads_A > 0 && n_reads_B > 0` is a shared junction; one of them `0`
  means it was seen on one side only.
- `n_reads_both` counts per read. A junction can be shared overall (seen by
  both A and B somewhere) yet have reads where only one alignment calls it:
  `n_reads_only_A` / `n_reads_only_B` > 0. Summed over all junctions,
  `n_reads_both` equals the per-read table's total `n_junctions_shared`.
- When reads on one side report different motif letters for the same
  junction, the first one with real letters is kept (real letters are
  preferred over `nn-nn`). On the same reference, real letters from A and B
  must agree.
- `motif_*` follows the junction's orientation: reverse-complemented on `-`
  junctions, reference-forward on `*` junctions. See
  [genomic-junction-compare-per-read.md](genomic-junction-compare-per-read.md).
