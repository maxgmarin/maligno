# `query_junction_diff.unmatched_junctions.{A,B}.tsv`

`{prefix}.query_junction_diff.unmatched_junctions.A.tsv.gz` /
`{prefix}.query_junction_diff.unmatched_junctions.B.tsv.gz` — distinct genomic splice
junctions, per side, that were never matched on the other side, with how many
reads support each — plus how many reads *total*, across the whole
comparison table, carry that same junction on this side. Written by
`toolkit query-junction-diff`, one file per side (mirrors the
`query_diff_regions.{A,B}.bed.gz` per-side-file convention).

Columns 1–5 roll up the `matched_in_genomic == false` rows of
[`query_junction_diff.per_read_per_junc_info.tsv.gz`](per-read-query-junction-diff.md)
for that side, grouped by `(chrom, genomic_start, genomic_end, strand)`. The A
file covers junctions called in A but unsupported by B, and the B file the
reverse. Column 6 counts how many reads in the whole comparison table carry the
same junction on this side.

6 columns.

| # | Column | Description |
|---|---|---|
| 1 | `chrom` | Genomic contig |
| 2 | `genomic_start` | Junction start (0-based half-open) |
| 3 | `genomic_end` | Junction end |
| 4 | `strand` | This side's alignment strand at this junction |
| 5 | `n_reads_query_junctions_different` | Count of **distinct reads** on this side reporting this exact junction, unmatched on the other side |
| 6 | `n_reads_with_junction_total` | Count of reads, across **every** row of the comparison table, with this exact junction on this side — differing or not |

**Notes**
- A read that reports the same unmatched junction more than once (from a
  single alignment) is only counted once toward
  `n_reads_query_junctions_different` — this is a distinct-reads count, not a
  raw occurrence tally.
- A junction with a high count is a systematic disagreement between A and B
  at that locus (recurring across many reads); a count of `1` is more likely
  per-read alignment noise.
- `n_reads_with_junction_total >= n_reads_query_junctions_different` always
  holds — every unmatched-and-differing occurrence of a junction is
  necessarily also counted in its total. The difference between the two is
  how many reads carried this junction and it *wasn't* a problem — a junction
  whose total is much larger than its unmatched count is common and
  well-supported elsewhere; a junction whose total equals its unmatched count
  never shows up anywhere else in the dataset.
- Side is fixed by the filename, not a column — "A unsupported by B" and "B
  unsupported by A" at the same coordinates are different, non-interchangeable
  events, so they are never merged into one row.
