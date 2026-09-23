# Output-table specs

Column-by-column specifications for every file `maligno compare` (and its
`compare-toolkit` building blocks) can write. Each doc gives, for its file:
what produces it, when it's written, and a table with one row per column.

| File | Spec |
|---|---|
| `{prefix}.{label}.alninfo.tsv.gz` | [alninfo.md](alninfo.md) |
| `{prefix}.{label}.readinfo.tsv.gz` | [readinfo.md](readinfo.md) |
| `{prefix}.compare.tsv.gz` / `.compare.parquet` | [compare.md](compare.md) |
| `{prefix}.compare.summary.tsv` | [compare-summary.md](compare-summary.md) |
| `{prefix}.query_diff_reads.tsv.gz` | [query-diff-reads.md](query-diff-reads.md) |
| `{prefix}.query_diff_regions.{A,B}.bed.gz` | [query-diff-regions.md](query-diff-regions.md) |
| `{prefix}.query_junction_diff.summary.tsv` | [query-junction-diff-summary.md](query-junction-diff-summary.md) |
| `{prefix}.query_junction_diff.per_read_per_junc_info.tsv.gz` | [per-read-query-junction-diff.md](per-read-query-junction-diff.md) |
| `{prefix}.query_junction_diff.unmatched_junctions.{A,B}.tsv.gz` | [query-junction-diff-unmatched.md](query-junction-diff-unmatched.md) |

These are per-column references. For prose explanations, schema-migration
history, and the classification logic behind the derived columns, see the
main [REFERENCE.md](../REFERENCE.md) and [COMPARE_TABLE.md](../COMPARE_TABLE.md).
