#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# Rebuild the edge-case fixture from the full GENCODE Splice-vs-SpliceHQ dataset.
#
# NOT part of the normal test loop. The committed .paf.gz files are the
# artifact; this script exists so their derivation is reproducible and auditable.
# Every emitted record is extracted verbatim — nothing is synthesized.
#
# Usage:
#   ./make-edge-case-fixture.sh <gencode-results-dir> [outdir]
#
# <gencode-results-dir> must contain:
#   Splice.AlnToHG38.PriAln.paf.gz, SpliceHQ.AlnToHG38.PriAln.paf.gz
#   Splice_vs_SpliceHQ.query_diff_reads.tsv.gz
#   Splice_vs_SpliceHQ.query_diff_reads.junctions.tsv.gz
# -----------------------------------------------------------------------------
set -euo pipefail

SRC="${1:?usage: $0 <gencode-results-dir> [outdir]}"
OUT="${2:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
A="$SRC/Splice.AlnToHG38.PriAln.paf.gz"
B="$SRC/SpliceHQ.AlnToHG38.PriAln.paf.gz"
for f in "$A" "$B" \
         "$SRC/Splice_vs_SpliceHQ.query_diff_reads.tsv.gz" \
         "$SRC/Splice_vs_SpliceHQ.query_diff_reads.junctions.tsv.gz"; do
  [ -r "$f" ] || { echo "error: missing input: $f" >&2; exit 2; }
done

# The 16 chosen reads and why. Selected once by scanning the full comparison
# table for varied properties (junction count, strand, mapping status); pinned
# here by ID so the fixture is stable and reviewable rather than re-sampled.
cat > "$OUT/read_categories.tsv" <<'CATS'
identical_single_exon	ENST00000194152.4
identical_single_exon	ENST00000194155.7
identical_multi_junction	ENST00000001008.6
identical_multi_junction	ENST00000002125.9
identical_minus_strand	ENST00000000412.8
identical_minus_strand	ENST00000001146.7
juncdiff	ENST00000054668.5
juncdiff	ENST00000233836.5
juncdiff	ENST00000254043.8
csdiff_only	ENST00000229771.11
csdiff_only	ENST00000242066.10
csdiff_only	ENST00000242770.9
unmapped_both	ENST00000290239.7
unmapped_both	ENST00000362102.3
aligned_only_B	ENST00000578854.1
identical_revcomp	ENST00000619436.1
CATS

KEEP="$(mktemp)"; trap 'rm -f "$KEEP"' EXIT
cut -f2 "$OUT/read_categories.tsv" | LC_ALL=C sort > "$KEEP"

extract () {  # extract <src.paf.gz> <dest.paf.gz>
  gzip -dc "$1" \
    | awk -F'\t' 'NR==FNR{k[$1];next} $1 in k' "$KEEP" - \
    | LC_ALL=C sort -t$'\t' -k1,1 \
    | gzip -c > "$2"
}

extract "$A" "$OUT/edge_cases.Splice.paf.gz"
extract "$B" "$OUT/edge_cases.SpliceHQ.paf.gz"

# ID-mismatch variant: drop one *different* real read from each side, so the
# read-ID sets differ in both directions. Reaches present_only_in_{A,B}_by_id
# and, with default flags, the read-ID-set mismatch error.
gzip -dc "$OUT/edge_cases.Splice.paf.gz"   | grep -v $'^ENST00000194152.4\t' | gzip -c > "$OUT/edge_cases_idmismatch.Splice.paf.gz"
gzip -dc "$OUT/edge_cases.SpliceHQ.paf.gz" | grep -v $'^ENST00000001008.6\t' | gzip -c > "$OUT/edge_cases_idmismatch.SpliceHQ.paf.gz"

n_a=$(gzip -dc "$OUT/edge_cases.Splice.paf.gz" | wc -l | tr -d ' ')
n_b=$(gzip -dc "$OUT/edge_cases.SpliceHQ.paf.gz" | wc -l | tr -d ' ')
[ "$n_a" -eq 16 ] && [ "$n_b" -eq 16 ] || {
  echo "error: expected 16 reads per side, got A=$n_a B=$n_b" >&2; exit 3; }
echo "wrote 4 PAFs + read_categories.tsv to $OUT (16 reads per side)"
