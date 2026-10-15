# `genomic_junction_compare.per_junction.tsv.gz`

`{prefix}.genomic_junction_compare.per_junction.tsv.gz` — written by
[`toolkit genomic-junction-compare`](../REFERENCE.md#toolkit-genomic-junction-compare).
One row per **distinct genomic junction** seen on either side, across the
whole comparison table, sorted by `target_name`, `intron_start`,
`intron_end`, `intron_strand`.

8 columns.

| # | Column | Description |
|---|---|---|
| 1 | `target_name` | Reference/target name |
| 2 | `intron_start` | Intron start, 0-based |
| 3 | `intron_end` | Intron end, half-open |
| 4 | `intron_strand` | Strand of the junction on the reference (`+`, `-`, or `*` unknown) |
| 5 | `n_reads_A` | Number of reads whose A alignment has this junction |
| 6 | `n_reads_B` | Number of reads whose B alignment has this junction |
| 7 | `motif_A` | Splice-site motif reported by A, in the junction's orientation (e.g. `gt-ag`); `nn-nn` if A didn't report the bases; `.` if no A read has the junction |
| 8 | `motif_B` | The same for B |

**Notes**
- A junction is `(target_name, intron_start, intron_end, intron_strand)`: the
  same coordinates on `+` and `-` (or `*`) are separate rows.
- `n_reads_A > 0 && n_reads_B > 0` is a shared junction; one of them `0`
  means it was seen on one side only.
- When reads on one side report different motif letters for the same
  junction, the first one with real letters is kept (real letters are
  preferred over `nn-nn`). On the same reference, real letters from A and B
  must agree.
- `motif_*` follows the junction's orientation: reverse-complemented on `-`
  junctions, reference-forward on `*` junctions. See
  [genomic-junction-compare-per-read.md](genomic-junction-compare-per-read.md).
