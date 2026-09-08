#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# Schema-stability manifest: fingerprint every `compare` output across all
# classification branches, so a refactor can be proven output-preserving.
#
# Usage:
#   ./scripts/schema-stability-manifest.sh <maligno-binary> <workdir> > manifest.txt
#
# Typical use (the Phase 1/2 gate for the ComparisonRow / serde refactor):
#   ./scripts/schema-stability-manifest.sh ./old-maligno /tmp/w1 > before.txt
#   ./scripts/schema-stability-manifest.sh ./target/release/maligno /tmp/w2 > after.txt
#   diff before.txt after.txt && echo "OUTPUT UNCHANGED"
#
# Hashes are taken over DECOMPRESSED bytes: gzip embeds an mtime, so hashing
# the .gz files directly would differ on every run and prove nothing.
# -----------------------------------------------------------------------------
set -euo pipefail

BIN="${1:?usage: $0 <maligno-binary> <workdir>}"
WORK="${2:?usage: $0 <maligno-binary> <workdir>}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TD="$HERE/test_data"
EC="$TD/edge_cases"

# Fail loudly on a missing/unusable binary. Without this the manifest comes out
# empty and the script exits 0 — a silent pass, which is worse than no check.
command -v "$BIN" >/dev/null 2>&1 || [ -x "$BIN" ] || {
  echo "error: '$BIN' is not an executable maligno binary" >&2; exit 2; }
"$BIN" --version >/dev/null 2>&1 || {
  echo "error: '$BIN --version' failed; is this a maligno binary?" >&2; exit 2; }
for f in "$TD/Splice.AlnToHG38.PriAln.paf.gz" "$EC/edge_cases.Splice.paf.gz"; do
  [ -r "$f" ] || { echo "error: missing test data: $f" >&2; exit 2; }
done
rm -rf "$WORK"; mkdir -p "$WORK"

sha () { shasum -a 256 | cut -d' ' -f1; }

# Emit "<hash>  <label>" for every file a run produced, decompressing .gz.
fingerprint () {
  local dir="$1" label="$2"
  [ -n "$(find "$dir" -type f -print -quit)" ] || {
    echo "error: nothing to fingerprint in $dir" >&2; exit 3; }
  find "$dir" -type f | LC_ALL=C sort | while read -r f; do
    local rel="${f#"$dir"/}"
    case "$f" in
      *.gz) printf '%s  %s/%s\n' "$(gzip -dc "$f" | sha)" "$label" "${rel%.gz}" ;;
      *)    printf '%s  %s/%s\n' "$(sha < "$f")"          "$label" "$rel" ;;
    esac
  done
}

run () {                      # run <label> <argv...>
  local label="$1"; shift
  local out="$WORK/$label"
  mkdir -p "$out"
  if ! "$@" --outdir "$out" --prefix run >"$WORK/$label.log" 2>&1; then
    echo "error: scenario '$label' failed; see $WORK/$label.log" >&2; exit 3
  fi
  [ -n "$(find "$out" -type f -print -quit)" ] || {
    echo "error: scenario '$label' produced no output files" >&2; exit 3; }
  fingerprint "$out" "$label"
}

# 1. chr22 bulk data — volume, real-world column values.
run chr22 "$BIN" compare \
  -a "$TD/Splice.AlnToHG38.PriAln.paf.gz"   --label-a Splice \
  -b "$TD/SpliceHQ.AlnToHG38.PriAln.paf.gz" --label-b SpliceHQ

# 2. Edge cases — unmapped, revcomp, junction/cs diffs. See edge_cases/README.md.
run edge "$BIN" compare \
  -a "$EC/edge_cases.Splice.paf.gz"   --label-a Splice \
  -b "$EC/edge_cases.SpliceHQ.paf.gz" --label-b SpliceHQ

# 3. Same data, sides swapped — the only way to reach aligned_only_A.
run edge_swapped "$BIN" compare \
  -a "$EC/edge_cases.SpliceHQ.paf.gz" --label-a SpliceHQ \
  -b "$EC/edge_cases.Splice.paf.gz"   --label-b Splice

# 4. Differing read-ID sets — reaches present_only_in_{A,B}_by_id.
run edge_idmismatch "$BIN" compare --allow-id-mismatch \
  -a "$EC/edge_cases_idmismatch.Splice.paf.gz"   --label-a Splice \
  -b "$EC/edge_cases_idmismatch.SpliceHQ.paf.gz" --label-b SpliceHQ

# 5. find-query-diff --compare-by junctions over the edge-case table.
"$BIN" find-query-diff -i "$WORK/edge/run.compare.tsv.gz" \
  --outdir "$WORK/fqd_junctions" --prefix run --compare-by junctions >/dev/null 2>&1
fingerprint "$WORK/fqd_junctions" fqd_junctions

# 6. Standalone compare-summary (the non-fused classifier path).
mkdir -p "$WORK/summary"
"$BIN" compare-summary -i "$WORK/edge/run.compare.tsv.gz" \
  -o "$WORK/summary/run.summary.tsv" >/dev/null 2>&1 \
  || { echo "error: compare-summary failed" >&2; exit 3; }
fingerprint "$WORK/summary" summary

# 7. The default read-ID-mismatch error must stay an error (text not hashed,
#    only the fact that it fails, so wording stays free to improve).
if "$BIN" compare \
     -a "$EC/edge_cases_idmismatch.Splice.paf.gz" \
     -b "$EC/edge_cases_idmismatch.SpliceHQ.paf.gz" \
     --outdir "$WORK/should_fail" --prefix x >/dev/null 2>&1; then
  echo "MISMATCH_GUARD  FAIL-did-not-error"
else
  echo "MISMATCH_GUARD  ok-errored"
fi
