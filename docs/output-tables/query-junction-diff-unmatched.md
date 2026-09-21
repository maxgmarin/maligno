# `query_junction_diff_unmatched.{A,B}.tsv`

`{prefix}.query_junction_diff_unmatched.A.tsv.gz` /
`{prefix}.query_junction_diff_unmatched.B.tsv.gz` — distinct genomic splice
junctions, per side, that were never matched on the other side, with how many
reads support each. Written by `compare-toolkit query-junction-diff`, one
file per side (mirrors the `query_diff_regions.{A,B}.bed.gz` per-side-file
convention).

Each file is a rollup over
[`per_read_query_junction_diff.summary.tsv.gz`](per-read-query-junction-diff.md)'s
`matched_in_genomic == false` rows for that side, grouped by `(chrom,
genomic_start, genomic_end, strand)`. The A file covers junctions called in A
unsupported by B; the B file, junctions called in B unsupported by A.

5 columns.

| # | Column | Description |
|---|---|---|
| 1 | `chrom` | Genomic contig |
| 2 | `genomic_start` | Junction start (0-based half-open) |
| 3 | `genomic_end` | Junction end |
| 4 | `strand` | This side's alignment strand at this junction |
| 5 | `n_reads_query_junctions_different` | Count of **distinct reads** on this side reporting this exact junction, unmatched on the other side |

**Notes**
- A read that reports the same unmatched junction more than once (from a
  single alignment) is only counted once toward
  `n_reads_query_junctions_different` — this is a distinct-reads count, not a
  raw occurrence tally.
- A junction with a high count is a systematic disagreement between A and B
  at that locus (recurring across many reads); a count of `1` is more likely
  per-read alignment noise.
- Side is fixed by the filename, not a column — "A unsupported by B" and "B
  unsupported by A" at the same coordinates are different, non-interchangeable
  events, so they are never merged into one row.
