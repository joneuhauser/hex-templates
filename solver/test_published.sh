#!/usr/bin/env bash
set -euo pipefail
case_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
bash "$case_dir/run.sh" build
binary="$case_dir/build/mirror_search"
out=$(mktemp -d "$case_dir/published-tests-XXXXXX")
for n in 36 44; do
  plane=diagonal; interior=33; orbits=12; vertices=51
  if [[ $n == 44 ]]; then plane=axial; interior=43; orbits=14; vertices=61; fi
  for mode in reference overlap budget all; do
    flags=()
    case "$mode" in
      reference) flags=(--early-overlap 0 --component-bound 0) ;;
      overlap) flags=(--early-budget 0 --early-parity 0) ;;
      budget) flags=(--early-parity 0) ;;
    esac
    stem="$out/seed$n-$mode"
    "$binary" --input "$case_dir/input/pyramid.mesh" --cap "$n" --interior "$interior" --threads 1 --planes "$plane" \
      --replay "$case_dir/input/published$n-normalized.mesh" --replay-actions "$case_dir/input/published$n-$plane.actions" \
      --output "$stem.jsonl" "${flags[@]}" 2> "$stem.log"
    grep -q "REPLAY_OK cells=$n cell_orbits=$orbits vertices=$vertices" "$stem.log"
    grep -q '^SEARCH_FINISHED candidates=1$' "$stem.log"
    cmp "$out/seed$n-reference.jsonl" "$stem.jsonl"
  done
  # The exactly symmetric geometry must also pass automatic coordinate actions.
  stem="$out/seed$n-symmetric"
  "$binary" --input "$case_dir/input/pyramid.mesh" --cap "$n" --interior "$interior" --threads 1 --planes "$plane" \
    --replay "$case_dir/input/published$n-symmetric.mesh" --output "$stem.jsonl" 2> "$stem.log"
  grep -q '^SEARCH_FINISHED candidates=1$' "$stem.log"
  cmp "$out/seed$n-reference.jsonl" "$stem.jsonl"
  # A seed with more cells than its cap must still fail.
  if "$binary" --input "$case_dir/input/pyramid.mesh" --cap "$((n-1))" --interior "$interior" --threads 1 --planes "$plane" \
    --replay "$case_dir/input/published$n-normalized.mesh" --replay-actions "$case_dir/input/published$n-$plane.actions" \
    --output "$out/seed$n-too-small.jsonl" 2> "$out/seed$n-too-small.log"; then
    echo "ERROR: seed$n accepted at cap $((n-1))." >&2; exit 1
  else [[ $? == 3 ]]; fi
  echo "PASS: published $n cells, $plane planes, $orbits cell orbits, $vertices vertices; all pruning modes."
done
# Explicit actions are checked, not trusted: the wrong boundary action must fail.
if "$binary" --input "$case_dir/input/pyramid.mesh" --cap 44 --interior 43 --threads 1 --planes diagonal \
  --replay "$case_dir/input/published44-normalized.mesh" --replay-actions "$case_dir/input/published44-axial.actions" \
  --output "$out/wrong-planes.jsonl" 2> "$out/wrong-planes.log"; then
  echo 'ERROR: incorrect mirror action accepted.' >&2; exit 1
else [[ $? == 2 ]]; fi
grep -q 'Replay boundary action mismatch' "$out/wrong-planes.log"
echo "Published-seed replay tests passed. Logs: $out"
