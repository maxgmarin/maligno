# `query_junction_compare.summary.tsv`

`{prefix}.query_junction_compare.summary.tsv` — read-selection funnel counts, then
distinct-junction counts per side, written by `toolkit query-junction-compare`
over an existing `compare.parquet`. Always uncompressed. 2 columns:
`category`, `count`.

Every `category` below is one **row** in the file, in this order:

| Category | Meaning |
|---|---|
| `label_A` | Which dataset side A is (from the input table's `label_A` column). String, not a count — provenance row |
| `label_B` | Which dataset side B is |
| `n_total` | Total rows read from the input comparison table |
| `n_neither_aligned` | Excluded: unmapped on both sides |
| `n_both_aligned_query_junctions_identical` | Excluded: mapped both sides, query-space junction sets match (`n_query_junctions_only_A == 0 && n_query_junctions_only_B == 0`) |
| `n_aligned_only_no_junctions` | Excluded: mapped on exactly one side, but that side has zero splice junctions — nothing to report |
| `n_query_junctions_different` | Selected — these reads got full junction reconstruction, written to [`query_junction_compare.per_read_per_junc_info.tsv.gz`](per-read-query-junction-compare.md) |
| `n_query_junctions_different_aligned_both` | ...of those, mapped on both sides with differing query-space junction sets |
| `n_query_junctions_different_aligned_only_A` | ...of those, mapped only in A, with `n_intron_events_A > 0` |
| `n_query_junctions_different_aligned_only_B` | ...of those, mapped only in B, with `n_intron_events_B > 0` |
| `n_unique_junctions_A` | Distinct genomic junctions `(target_name, intron_start, intron_end, intron_strand)` across every row's side-A alignment |
| `n_unique_junctions_B` | The same, for side B |
| `n_unique_junctions_shared` | Distinct junctions present on both sides, anywhere in the table |
| `n_unique_junctions_only_A` | Distinct junctions on side A that appear in no side-B alignment |
| `n_unique_junctions_only_B` | Distinct junctions on side B that appear in no side-A alignment |
| `n_unmatched_junctions_A` | Rows in [`unmatched_junctions.A.tsv.gz`](query-junction-compare-unmatched.md): distinct side-A junctions unmatched within at least one differing read |
| `n_unmatched_junctions_B` | Rows in `unmatched_junctions.B.tsv.gz` |

**Invariants**
- `n_total = n_neither_aligned + n_both_aligned_query_junctions_identical + n_aligned_only_no_junctions + n_query_junctions_different`
- `n_query_junctions_different = n_query_junctions_different_aligned_both + n_query_junctions_different_aligned_only_A + n_query_junctions_different_aligned_only_B`
- `n_unique_junctions_A = n_unique_junctions_shared + n_unique_junctions_only_A` (symmetric for B)

**Notes**
- Every row of the input table is tallied into exactly one bucket, whether or
  not it passes the "differing" filter — so this file is a complete
  accounting of the input, not just of what was selected.
- Reconciles with `compare.summary.tsv`
  ([spec](compare-summary.md)): `n_neither_aligned` ==
  its `aligned_neither`; `n_both_aligned_query_junctions_identical +
  n_query_junctions_different_aligned_both` == its `aligned_both`;
  `n_query_junctions_different_aligned_only_A` +
  (only-A reads with zero junctions) == its `aligned_only_A` (symmetric for B).
- The `n_unique_junctions_*` rows count each read's representative alignment
  only (the one per read per side that the comparison table keeps). Secondary
  or supplementary alignments that were not chosen are not in the table and
  are not counted.
- The junction key matches the unmatched tables and `genomic-junction-compare`:
  `intron_strand` is the intron's strand on the reference, so the same intron
  counts once whichever strand its reads aligned to. The
  `n_unique_junctions_*` rows equal `genomic-junction-compare`'s
  `n_distinct_junctions_*` on the same table.
- `n_unique_junctions_only_A` and `n_unmatched_junctions_A` answer different
  questions. A junction is *unmatched* when a differing read's A alignment
  has it and that same read's B alignment doesn't; it may still appear in
  other reads' B alignments. A junction is *only in A* when no B alignment in
  the table has it, including junctions from reads that were not selected as
  differing (for example, reads whose query-space junctions agree but are
  placed at a different genomic locus).
