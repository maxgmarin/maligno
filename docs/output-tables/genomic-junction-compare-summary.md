# `genomic_junction_compare.summary.tsv`

`{prefix}.genomic_junction_compare.summary.tsv` — written by
[`toolkit genomic-junction-compare`](../REFERENCE.md#toolkit-genomic-junction-compare).
2 columns: `category`, `count` (the same layout as the other summaries).

Every `category` below is one **row**, in this order:

| Category | Meaning |
|---|---|
| `label_A` | Which dataset side A is (from the table's `label_A` column). String, not a count |
| `label_B` | Which dataset side B is |
| `reads_compared` | Rows of the comparison table |
| `reads_with_junctions_both` | Reads with at least one genomic junction on both sides |
| `reads_with_junctions_only_A` | Reads with junctions in A and none in B (B unmapped, or mapped without a junction) |
| `reads_with_junctions_only_B` | The reverse |
| `reads_with_junctions_neither` | Reads with no junction on either side |
| `reads_genomic_junctions_identical` | Of `reads_with_junctions_both`: the two junction sets are identical |
| `reads_genomic_junctions_differ` | Of `reads_with_junctions_both`: the sets differ |
| `reads_genomic_intron_strand_differ` | Of `reads_with_junctions_both`: both intron strands are known and different |
| `n_distinct_junctions_A` | Distinct junctions seen in A across all reads |
| `n_distinct_junctions_B` | Distinct junctions seen in B |
| `n_distinct_junctions_shared` | Distinct junctions seen in both |
| `n_distinct_junctions_only_A` | Distinct junctions seen only in A |
| `n_distinct_junctions_only_B` | Distinct junctions seen only in B |
| `n_distinct_junctions_unknown_strand_A` | Of `n_distinct_junctions_A`, those with intron strand `*` |
| `n_distinct_junctions_unknown_strand_B` | The same for B |

**Notes**
- The four `reads_with_junctions_*` rows add up to `reads_compared`, and
  `reads_genomic_junctions_identical + reads_genomic_junctions_differ =
  reads_with_junctions_both`.
- A junction is `(target_name, intron_start, intron_end, intron_strand)`, so a
  strand disagreement, or an unknown strand on one side only, makes the
  junctions differ even when the coordinates agree.
  `reads_genomic_intron_strand_differ` and the `unknown_strand` rows show how much of
  the difference comes from strand.
- "Identical" is set equality only; it doesn't require the alignments to
  start at the same position. Position is covered by the `ref_*` rows of
  [compare-summary.md](compare-summary.md).
- The distinct-junction counts differ from query-junction-diff's
  `n_unique_junctions_*`, which key junctions by alignment strand instead of
  intron strand.
