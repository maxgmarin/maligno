# `per_read_query_junction_diff.summary.tsv`

`{prefix}.per_read_query_junction_diff.summary.tsv.gz` — one row per
reconstructed splice junction, per side (A/B), per read selected as
"differing" (see [query-junction-diff-summary.md](query-junction-diff-summary.md)
for the selection formula). Written by `compare-toolkit query-junction-diff`.

11 columns.

| # | Column | Description |
|---|---|---|
| 1 | `Read_Name` | Read ID |
| 2 | `side` | Which set this junction was reconstructed from: `A` or `B` |
| 3 | `junction_index` | 1-based ordinal position of this junction along the read, 5'→3' in original-read orientation (see Notes) |
| 4 | `strand` | This alignment's strand (`+`/`-`) |
| 5 | `query_pos` | Query-space junction coordinate (plus-strand-normalized, same convention as [compare.md](compare.md)'s `junctions` column) |
| 6 | `chrom` | Genomic contig of this junction |
| 7 | `genomic_start` | Junction's absolute genomic interval start (0-based half-open) |
| 8 | `genomic_end` | Junction's absolute genomic interval end |
| 9 | `matched_in_query` | `1` if `query_pos` also appears in the *other* side's query-space junction set for this read |
| 10 | `matched_in_genomic` | `1` if the genomic interval also appears in the *other* side's genomic junction set for this read |
| 11 | `other_side_aligned` | `1` if the other side has any alignment at all for this read; `0` means this read only mapped on the `side` shown, so columns 9/10 are trivially `0` (nothing to compare against) |

**Notes**
- **Why this table re-derives junctions from `cs` instead of reusing
  `compare.tsv`'s stored `junctions`/`genomic_junctions` columns**: for
  `-`-strand alignments, `record.rs` sorts the stored `junctions` column into
  ascending query order but never reorders the paired `genomic_junctions`
  column to match, so the two stored columns are not reliably index-paired
  per junction (every individual value in each is still correct — only the
  positional correspondence breaks). This has never affected any existing
  output, since nothing else in the crate reads the two columns index-paired.
  This table is the first consumer that needs that pairing, so it re-parses
  each side's `cs` tag fresh and sorts the query and genomic coordinates
  together, with one shared permutation, guaranteeing `query_pos`/`(chrom,
  genomic_start, genomic_end)` on the same row describe the same intron.
- `junction_index` is computed independently per side — it is a display/sort
  key, not a matching key. A's junction #3 and B's junction #3 are not
  guaranteed to be the same biological junction (an extra/missing junction
  partway through the read shifts every later index on that side).
- For a `+`-strand alignment, `genomic_start` increases monotonically with
  `junction_index`; for a `-`-strand alignment it decreases monotonically
  (query 5'→3' order runs opposite to the reference's forward direction on
  the minus strand).
- Every distinct `matched_in_genomic == false` junction, aggregated across
  reads, is rolled up in
  [`query_junction_diff_unmatched.A/B.tsv.gz`](query-junction-diff-unmatched.md).
