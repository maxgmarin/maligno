# `query_diff_reads.tsv`

`{prefix}.query_diff_reads.tsv.gz` — one row per read whose alignment
**differs** between set A and set B, plus 8 classification booleans.
Written by `compare` by default (fused into its single merge pass, fixed at
`--space query --compare-by all`; `--skip-find-aln-diff` opts out) and by
standalone `compare-toolkit find-aln-diff` (any `--space`/`--compare-by`
combination — filenames gain a `.junctions`/`reference_diff` segment
accordingly; see [REFERENCE.md](../REFERENCE.md#--compare-by--what-defines-a-difference-within-the-active---space)
for the full mode matrix).

10 columns.

| # | Column | Description |
|---|---|---|
| 1 | `Read_Name` | Read ID |
| 2 | `outcome` | Category: `diff_aln_to_both` (mapped in both sets, not identical under the active mode; `reference_diff` under `--space reference`), `diff_aln_only_A` (mapped in A only), or `diff_aln_only_B` (mapped in B only) |
| 3 | `query_identical_same_strand` | `1` if `query_identical` reached via the same-strand branch |
| 4 | `query_identical_revcomp` | `1` if `query_identical` reached via the reverse-complement branch |
| 5 | `query_junctions_identical` | `1` if both mapped and the **query-space** splice-junction sets match — computed unconditionally, not just under `--compare-by junctions` |
| 6 | `ref_same_position_same_aln` | `1` if reference-identical: same `TargetChr`/`Strand`/`Target_Start` and same `cs` |
| 7 | `ref_same_position_diff_aln` | `1` if same position, different `cs` |
| 8 | `ref_diff_position_same_aln` | `1` if same `cs`, different position (relocated) |
| 9 | `ref_diff_position_diff_aln` | `1` if both position and `cs` differ |
| 10 | `ref_same_position_same_junctions` | `1` if both mapped, same position, and the **genomic-coordinate** junction sets also match — computed unconditionally |

**Notes**
- Columns 3–10 are computed **the same way regardless of the active
  `--space`/`--compare-by`** — so one run can show, e.g., a read that's
  query-different but reference-identical, without a second run. All 8 are
  `0` for a read mapped on only one side.
- For a both-mapped read, exactly one of columns 6–9 is `1` (they partition
  the four reference-space outcomes).
- A read is only written here if it's a *difference* under the accumulator's
  active mode — a read considered "identical" is excluded (see the
  complementary, opt-in `{prefix}.query_identical_reads.tsv.gz`, produced
  only by standalone `find-aln-diff --emit-identical-reads`, never by
  `compare`; it carries the same 10-column schema for the complementary read
  set).
- Every genomic locus this table's reads cluster into is summarized in
  [query-diff-regions.md](query-diff-regions.md).
