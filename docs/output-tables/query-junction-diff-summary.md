# `query_junction_diff.summary.tsv`

`{prefix}.query_junction_diff.summary.tsv` — parse-time funnel counts written
by `compare-toolkit query-junction-diff`, over an existing
`compare.tsv`/`.parquet`. Always uncompressed. 2 columns: `Category`, `Count`.

Every `Category` below is one **row** in the file, in this order:

| Category | Meaning |
|---|---|
| `label_A` | Which dataset side A is (from the input table's `Label_A` column). String, not a count — provenance row |
| `label_B` | Which dataset side B is |
| `n_total` | Total rows read from the input comparison table |
| `n_neither_aligned` | Excluded: unmapped on both sides |
| `n_both_aligned_query_junctions_identical` | Excluded: mapped both sides, query-space junction sets match (`N_Junctions_OnlyA == 0 && N_Junctions_OnlyB == 0`) |
| `n_aligned_only_no_junctions` | Excluded: mapped on exactly one side, but that side has zero splice junctions — nothing to report |
| `n_query_junctions_different` | Selected — these reads got full junction reconstruction, written to [`per_read_query_junction_diff.summary.tsv.gz`](per-read-query-junction-diff.md) |
| `n_query_junctions_different_aligned_both` | ...of those, mapped on both sides with differing query-space junction sets |
| `n_query_junctions_different_aligned_only_A` | ...of those, mapped only in A, with `JuncCount_A > 0` |
| `n_query_junctions_different_aligned_only_B` | ...of those, mapped only in B, with `JuncCount_B > 0` |

**Invariants**
- `n_total = n_neither_aligned + n_both_aligned_query_junctions_identical + n_aligned_only_no_junctions + n_query_junctions_different`
- `n_query_junctions_different = n_query_junctions_different_aligned_both + n_query_junctions_different_aligned_only_A + n_query_junctions_different_aligned_only_B`

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
