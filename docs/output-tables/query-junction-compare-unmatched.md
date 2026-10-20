# `query_junction_compare.unmatched_junctions.{A,B}.tsv`

`{prefix}.query_junction_compare.unmatched_junctions.A.tsv.gz` /
`{prefix}.query_junction_compare.unmatched_junctions.B.tsv.gz` — distinct genomic splice
junctions, per side, that were never matched on the other side, with how many
reads support each — plus how many reads *total*, across the whole
comparison table, carry that same junction on this side. Written by
`toolkit query-junction-compare`, one file per side (like `find-aln-diff`'s
per-side region tables).

Columns 1–6 roll up the `matched_in_genomic == false` rows of
[`query_junction_compare.per_read_per_junc_info.tsv.gz`](per-read-query-junction-compare.md)
for that side, grouped by junction: `(target_name, intron_start, intron_end,
intron_strand)`, the same identity and column names as `genomic-junction-compare`'s
[per-junction table](genomic-junction-compare-per-junction.md). The A file
covers junctions called in A but unsupported by B, and the B file the reverse.
Column 7 counts how many reads in the whole comparison table carry the same
junction on this side. With `--emit-all-junctions`, every junction of each side
is also written, in the same columns; this table is its subset with
`n_reads_query_junctions_different > 0` (see
[query-junction-compare-all-junctions.md](query-junction-compare-all-junctions.md)).

7 columns.

| # | Column | Description |
|---|---|---|
| 1 | `target_name` | Genomic contig |
| 2 | `intron_start` | Intron start (0-based half-open) |
| 3 | `intron_end` | Intron end |
| 4 | `intron_strand` | The intron's strand on the reference (this side's `genomic_intron_strand`; `*` when unknown) |
| 5 | `motif` | Splice-site motif, donor-acceptor in the junction's orientation (e.g. `gt-ag`): the first seen among this side's unmatched reads, preferring real letters over `nn-nn` |
| 6 | `n_reads_query_junctions_different` | Count of **distinct reads** on this side reporting this exact junction, unmatched on the other side |
| 7 | `n_reads_with_junction_total` | Count of reads, across **every** row of the comparison table, with this exact junction on this side — differing or not |

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
- For scale, the number of rows in each file and the total number of distinct
  junctions per side are in
  [`query_junction_compare.summary.tsv`](query-junction-compare-summary.md)
  (`n_unmatched_junctions_{A,B}` and `n_unique_junctions_{A,B}`).
- Side is fixed by the filename, not a column — "A unsupported by B" and "B
  unsupported by A" at the same coordinates are different, non-interchangeable
  events, so they are never merged into one row.
