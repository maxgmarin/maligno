# Edge-case comparison fixture

A 16-read PAF pair that exercises every branch of `compare`'s per-read
classifier. It exists because the main `test_data/` chr22 pair, despite having
11,578 reads, reaches **only one** of the four mapping-status branches.

## Files

| File pair | How it was made | Reaches |
|---|---|---|
| `edge_cases.{Splice,SpliceHQ}.paf.gz` | 16 selected reads, extracted as-is | most branches |
| `edge_cases_idmismatch.{Splice,SpliceHQ}.paf.gz` | one *different* read deleted from each side | `present_only_in_{A,B}_by_id`, and the default read-ID-set error |

`aligned_only_A` needs no extra data: run `compare` with the two files swapped
(the fixture is deliberately asymmetric).

## Coverage

| Branch / hazard | chr22 pair | this fixture |
|---|:--:|:--:|
| `aligned_both` | ✅ | ✅ 13 |
| `aligned_only_B` | ❌ | ✅ 1 (`ENST00000578854.1`) |
| `aligned_only_A` | ❌ | ✅ 1 (swap the sides) |
| `aligned_neither` (counts as `query_identical`) | ❌ | ✅ 2 |
| `query_identical_same_strand` | ✅ | ✅ 6 |
| `query_identical_revcomp` | ❌ | ✅ 1 (`ENST00000619436.1`) |
| `query_not_identical` | ✅ | ✅ 6 |
| `ref_same_position_same_aln` | ✅ | ✅ 6 |
| `ref_same_position_diff_aln` | ✅ | ✅ 6 |
| `ref_diff_position_same_aln` (relocated) | ❌ | ❌ — not yet covered |
| `ref_diff_position_diff_aln` | ✅ | ✅ 1 (`ENST00000619436.1`) |
| junction-set difference | ✅ | ✅ 3 |
| cs difference with **identical** junctions | ✅ | ✅ 3 |
| `present_only_in_{A,B}_by_id` | ❌ | ✅ (idmismatch variant) |
| read-ID-set mismatch **error** | ❌ | ✅ (idmismatch variant, default flags) |
| empty-string cells (`cs` of an unmapped read) | ❌ | ✅ |
| `*` target, `()` junctions, `NaN` identity | ❌ | ✅ |
| integral floats rendered `N.0` | ✅ | ✅ |
| single-exon / multi-junction / minus-strand reads | ✅ | ✅ 2 each |

### Still not covered

- **Multiple alignments per read** (`Num_Aln > 1`). Both this fixture and the
  chr22 pair are primary-alignment-only, one row per read, so
  `readinfo.rs::collapse_group`'s tie-breaking (ms → AS → MQ) is never
  exercised. Would need a non-`PriAln` source.

## Use

Run `compare` on each scenario from the repo root:

```bash
EC=test_data/edge_cases

# Most branches
maligno compare -a $EC/edge_cases.Splice.paf.gz -b $EC/edge_cases.SpliceHQ.paf.gz \
  --label-a Splice --label-b SpliceHQ --outdir /tmp/edge --prefix edge

# aligned_only_A: the same pair with the sides swapped
maligno compare -a $EC/edge_cases.SpliceHQ.paf.gz -b $EC/edge_cases.Splice.paf.gz \
  --label-a SpliceHQ --label-b Splice --outdir /tmp/edge_swapped --prefix edge

# present_only_in_{A,B}_by_id (without --allow-id-mismatch this must error)
maligno compare --allow-id-mismatch \
  -a $EC/edge_cases_idmismatch.Splice.paf.gz -b $EC/edge_cases_idmismatch.SpliceHQ.paf.gz \
  --label-a Splice --label-b SpliceHQ --outdir /tmp/edge_idmismatch --prefix edge
```
