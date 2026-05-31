#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# Sort-order / overlap diagnostic for two readinfo TSVs.
#
# The `compare` and `compare-junctions` subcommands use a streaming merge-join
# on the (Read_Name, Read_Len) key, which assumes both inputs are sorted in
# the same byte-lexicographic order. If they aren't, matches are silently
# missed and the `matched` count in the comparison summary will be lower than
# it should be. This script tells you what the matched count SHOULD be by
# computing the full-key set intersection — a much cheaper check than redoing
# the join.
#
# Usage:
#   ./check-readinfo-overlap.sh a.readinfo.tsv[.gz] b.readinfo.tsv[.gz]
#
# Exit codes:
#   0 — diagnostic ran successfully
#   1 — bad arguments
# -----------------------------------------------------------------------------
set -euo pipefail

if [ $# -ne 2 ]; then
  echo "usage: $0 a.readinfo.tsv[.gz] b.readinfo.tsv[.gz]" >&2
  exit 1
fi

a="$1"
b="$2"

decompress() {
  case "$1" in
    *.gz) zcat < "$1" ;;
    *)    cat    "$1" ;;
  esac
}

# Build unique (Read_Name, Read_Len) keys from each file.
a_keys=$(mktemp); b_keys=$(mktemp)
trap 'rm -f "$a_keys" "$b_keys"' EXIT

decompress "$a" | tail -n +2 | awk -F'\t' '{print $1"\t"$2}' \
  | LC_ALL=C sort -u > "$a_keys"
decompress "$b" | tail -n +2 | awk -F'\t' '{print $1"\t"$2}' \
  | LC_ALL=C sort -u > "$b_keys"

n_a=$(wc -l < "$a_keys" | tr -d ' ')
n_b=$(wc -l < "$b_keys" | tr -d ' ')
n_full=$(LC_ALL=C comm -12 "$a_keys" "$b_keys" | wc -l | tr -d ' ')

# Also compute intersection on Read_Name alone — if this is larger than the
# full-key intersection, some reads share names but differ on Read_Len.
n_name=$(LC_ALL=C comm -12 \
  <(awk -F'\t' '{print $1}' "$a_keys" | LC_ALL=C sort -u) \
  <(awk -F'\t' '{print $1}' "$b_keys" | LC_ALL=C sort -u) \
  | wc -l | tr -d ' ')

cat <<EOF
Sort/overlap diagnostic
  A: $a
  B: $b

  rows in A (unique full-key):     $n_a
  rows in B (unique full-key):     $n_b
  intersection by (Name, Len):     $n_full
  intersection by Name only:       $n_name
EOF

if [ "$n_name" -gt "$n_full" ]; then
  diff_count=$(( n_name - n_full ))
  cat <<EOF

  ⚠  $diff_count read names appear in both files but with different Read_Len
     values. Those rows will NOT match in compare. Usually indicates an
     upstream soft-clip / qlen difference between the two BAMs.
EOF
fi

cat <<EOF

  ⇒ \`compare\`'s 'matched' count should equal $n_full. If maligno's reported
    matched count is lower than that, the two readinfo files are not sorted
    in the same byte-lexicographic order. Re-sort each with:

        LC_ALL=C sort -t\$'\\t' -k1,1 -k2,2n input.readinfo.tsv > sorted.tsv
        # (keep the header separately)

    or rerun the maligno pipeline from PAF level — paf2alninfo + readinfo
    produce consistently-sorted output by construction.
EOF
