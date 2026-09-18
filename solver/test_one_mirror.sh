#!/usr/bin/env bash
set -euo pipefail
case_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
bash "$case_dir/run.sh" build
binary="$case_dir/build/mirror_search"
out=$(mktemp -d "$case_dir/one-tests-XXXXXX")
for n in 36 44; do
  plane=diagonal; interior=33
  if [[ $n == 44 ]]; then plane=axial; interior=43; fi
  for mode in reference all; do
    flags=(); [[ $mode != reference ]] || flags=(--early-overlap 0 --component-bound 0)
    stem="$out/seed$n-$mode"
    "$binary" --boundary-canonical 0 --mirrors 1 --planes "$plane" --input "$case_dir/input/pyramid.mesh" --cap "$n" --interior "$interior" --threads 1 \
      --replay "$case_dir/input/published$n-normalized.mesh" --replay-actions "$case_dir/input/published$n-$plane-one.actions" \
      --output "$stem.jsonl" "${flags[@]}" 2> "$stem.log"
    grep -q '^SEARCH_FINISHED candidates=1$' "$stem.log"
  done
  cmp "$out/seed$n-reference.jsonl" "$out/seed$n-all.jsonl"
  "$binary" --boundary-canonical 0 --mirrors 1 --planes "$plane" --input "$case_dir/input/pyramid.mesh" --cap "$n" --interior "$interior" --threads 1 \
    --replay "$case_dir/input/published$n-symmetric.mesh" --output "$out/seed$n-symmetric.jsonl" 2> "$out/seed$n-symmetric.log"
  cmp "$out/seed$n-all.jsonl" "$out/seed$n-symmetric.jsonl"
  if "$binary" --boundary-canonical 0 --mirrors 1 --planes "$plane" --input "$case_dir/input/pyramid.mesh" --cap "$((n-1))" --interior "$interior" --threads 1 \
    --replay "$case_dir/input/published$n-symmetric.mesh" --output "$out/seed$n-small.jsonl" 2> "$out/seed$n-small.log"; then
    echo 'ERROR: oversized seed accepted.'; exit 1
  else [[ $? == 3 ]]; fi
  tail -n 2 "$out/seed$n-all.log"
done
for plane in axial diagonal; do
  "$binary" --boundary-canonical 0 --mirrors 1 --planes "$plane" --input "$case_dir/input/one-$plane-boundary.mesh" --cap 8 --threads 1 \
    --replay "$case_dir/input/one-$plane-seed.mesh" --output "$out/one-$plane.jsonl" 2> "$out/one-$plane.log"
  orbits=4; [[ $plane != diagonal ]] || orbits=6
  grep -q "REPLAY_OK cells=8 cell_orbits=$orbits" "$out/one-$plane.log"
  if "$binary" --boundary-canonical 0 --mirrors 2 --planes "$plane" --input "$case_dir/input/one-$plane-boundary.mesh" --cap 8 --threads 1 \
    --replay "$case_dir/input/one-$plane-seed.mesh" --output "$out/two-$plane.jsonl" 2> "$out/two-$plane.log"; then
    echo 'ERROR: asymmetric fixture accepted in two-plane mode.'; exit 1
  else [[ $? == 2 ]]; fi
  grep -q 'Input coordinates do not have the requested reflections' "$out/two-$plane.log"
  # Isolate the older overlap optimization: cover packing can reject doomed
  # prefixes before they reach the traced DFS entry.
  for fast in 0 1; do
    "$binary" --boundary-canonical 0 --mirrors 1 --planes "$plane" --input "$case_dir/input/pyramid.mesh" --cap 10 --threads 4 \
      --face-cover 0 --component-bound 0 --early-overlap "$fast" --state-trace "$out/H10-$plane-$fast.states" \
      --output "$out/H10-$plane-$fast.jsonl" 2> "$out/H10-$plane-$fast.log"
    LC_ALL=C sort "$out/H10-$plane-$fast.states" > "$out/H10-$plane-$fast.sorted"
  done
  cmp "$out/H10-$plane-0.sorted" "$out/H10-$plane-1.sorted"
  for bound in 0 1; do
    "$binary" --boundary-canonical 0 --mirrors 1 --planes "$plane" --input "$case_dir/input/cube.mesh" --cap 7 --threads 4 \
      --component-bound "$bound" --output "$out/cube7-$plane-$bound.jsonl" 2> "$out/cube7-$plane-$bound.log"
    LC_ALL=C sort "$out/cube7-$plane-$bound.jsonl" > "$out/cube7-$plane-$bound.sorted"
  done
  cmp "$out/cube7-$plane-0.sorted" "$out/cube7-$plane-1.sorted"
  grep -q '^SEARCH_FINISHED candidates=2$' "$out/cube7-$plane-1.log"
done
echo "One-mirror tests passed. Logs: $out"
