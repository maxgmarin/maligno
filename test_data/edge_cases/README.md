# Edge-case comparison fixture

A 16-read PAF pair that exercises every branch of `compare`'s per-read
classifier. It exists because the main `test_data/` chr22 pair, despite having
11,578 reads, reaches **only one** of the four mapping-status branches.

## Why this fixture is needed

Measured on `test_data/*.AlnToHG38.PriAln.paf.gz`:

```
aligned_both     11578
aligned_only_A       0
aligned_only_B       0
aligned_neither      0
```

Zero unmapped reads, zero one-sided reads, zero rows where `TargetChr_*` is
`*`, and **zero empty-string cells**. So the entire unmapped code path — empty
`cs`, `junctions` rendered as `()`, `NaN` identity, zeroed counters — was
untested. That matters most for the planned `ComparisonRow` refactor: today the
per-side columns are emitted as raw pass-through strings, and parsing them into
typed fields then re-serializing would silently normalize an empty cell to `0`.
Nothing in the repo would have caught it.

## Provenance — every record is real

No row is hand-written. All records were extracted verbatim from the
507,365-transcript GENCODE v49 PriChr comparison
(`minimap2 -x splice` vs `-x splice:hq`, aligned to GRCh38), which is the same
software and parameters that produced the chr22 pair. `read_categories.tsv`
records which read was chosen for which reason.

Two derived variants, both pure subsetting — no invented content:

| File pair | How it was made | Reaches |
|---|---|---|
| `edge_cases.{Splice,SpliceHQ}.paf.gz` | 16 selected reads, extracted as-is | most branches |
| `edge_cases_idmismatch.{Splice,SpliceHQ}.paf.gz` | one *different* read deleted from each side | `present_only_in_{A,B}_by_id`, and the default read-ID-set error |

`aligned_only_A` needs no extra data: run `compare` with the two files swapped
(the fixture is deliberately asymmetric). `scripts/schema-stability-manifest.sh`
does exactly that as its `edge_swapped` scenario.

## Coverage

| Branch / hazard | chr22 pair | this fixture |
|---|:--:|:--:|
| `aligned_both` | ✅ | ✅ 13 |
| `aligned_only_B` | ❌ | ✅ 1 (`ENST00000578854.1`) |
| `aligned_only_A` | ❌ | ✅ 1 (swap the sides) |
| `aligned_neither` | ❌ | ✅ 2 |
| `query_identical_same_strand` | ✅ | ✅ 6 |
| `query_identical_revcomp` | ❌ | ✅ 1 (`ENST00000619436.1`) |
| `query_not_identical` | ✅ | ✅ 6 |
| junction-set difference | ✅ | ✅ 3 |
| cs difference with **identical** junctions | ✅ | ✅ 3 |
| `present_only_in_{A,B}_by_id` | ❌ | ✅ (idmismatch variant) |
| read-ID-set mismatch **error** | ❌ | ✅ (idmismatch variant, default flags) |
| empty-string cells (`cs` of an unmapped read) | ❌ | ✅ |
| `*` target, `()` junctions, `NaN` identity | ❌ | ✅ |
| integral floats rendered `N.0` (see below) | ✅ | ✅ |
| single-exon / multi-junction / minus-strand reads | ✅ | ✅ 2 each |

`ENST00000619436.1` is the one revcomp case in all 507,365 transcripts: the same
transcript aligned to chrY on `+` by `splice` and on `−` by `splice:hq`, at
different loci, with cs tags that are exact reverse complements. It also
separates `query_identical` (7) from `reference_identical` (6).

### Still not covered

- **Multiple alignments per read** (`Num_Aln > 1`). Both this fixture and the
  chr22 pair are primary-alignment-only, one row per read, so
  `readinfo.rs::collapse_group`'s tie-breaking (ms → AS → MQ) is never
  exercised. Would need a non-`PriAln` source.
- **Float formatting** is *covered* but is worth calling out as the main
  refactor hazard: `io_utils::fmt_float` appends `.0` to integral values and
  Rust's default `{}` for `f64` does not. On the chr22 table `QueryAlnCov_Diff`
  is integral in 11,573 of 11,578 rows, so a generic `serde` serializer would
  rewrite ~99% of rows across eight columns unless `f64` goes through
  `fmt_float`.

## Use

```bash
./scripts/schema-stability-manifest.sh ./target/release/maligno /tmp/w > after.txt
diff test_data/schema_manifest.v0.15.0.txt after.txt && echo "OUTPUT UNCHANGED"
```

`test_data/schema_manifest.v0.15.0.txt` is the committed baseline: 46 SHA-256
fingerprints over the decompressed outputs of six scenarios. The filename carries
a version on purpose — an output change must rename it, which makes regenerating
the baseline a deliberate act rather than an invisible overwrite.

Its history is a good illustration of what the gate is for:

- **v0.14.0** (column regroup) moved exactly 4 of 45 fingerprints — the four
  `compare.tsv` files — and left every summary and region output identical, which
  is what a naming-and-ordering-only change should look like.
- **v0.14.1** (typed `ComparisonRow`) moved **nothing**: 46 of 46 identical, as a
  pure refactor must.
- **v0.15.0** (unmapped soft-clip fix) moved exactly the 15 `alninfo`/`readinfo`/
  `compare` fingerprints of the three edge scenarios and **not one chr22
  fingerprint** — because chr22 contains no unmapped reads, so only this fixture
  could see the change. Every summary and `find-query-diff` output stayed
  identical too, confirming the classifier never reads soft-clip.

## Regenerating

`make-edge-case-fixture.sh` rebuilds the pair from the full GENCODE comparison.
It needs that dataset present and is **not** part of the normal test loop — the
committed `.paf.gz` files are the artifact; the script only records how they
were derived.
