#!/usr/bin/env bash
set -euo pipefail
case_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
bash "$case_dir/run.sh" build
cmake --build "$case_dir/build" --target pattern_test parity_test completion_test geometry_test -- -j "${BUILD_JOBS:-8}"
"$case_dir/build/pattern_test"
"$case_dir/build/parity_test"
"$case_dir/build/completion_test" "$case_dir/input"
"$case_dir/build/geometry_test" "$case_dir/input"
out=$(mktemp -d "$case_dir/test-results-XXXXXX")
binary="$case_dir/build/mirror_search"
"$binary" --input "$case_dir/input/pyramid.mesh" --cap 36 --threads 1 \
  --replay "$case_dir/input/symmetric36.mesh" --output "$out/seed36.jsonl" 2> "$out/seed36.log"
grep -q 'REPLAY_OK cells=36 cell_orbits=12 vertices=51' "$out/seed36.log"
grep -q '^SEARCH_FINISHED candidates=1$' "$out/seed36.log"
if "$binary" --input "$case_dir/input/pyramid.mesh" --cap 35 --threads 1 \
  --replay "$case_dir/input/symmetric36.mesh" --output "$out/seed35.jsonl" 2> "$out/seed35.log"; then
  echo 'ERROR: cap 35 accepted the 36-cell seed.' >&2; exit 1
else [[ $? == 3 ]]; fi
for planes in diagonal axial; do
  "$binary" --input "$case_dir/input/cube.mesh" --cap 1 --threads 2 --planes "$planes" \
    --output "$out/cube-$planes.jsonl" 2> "$out/cube-$planes.log"
  grep -q '^SEARCH_FINISHED candidates=1$' "$out/cube-$planes.log"
  for n in 2 3; do
    "$binary" --input "$case_dir/input/grid$n-boundary.mesh" --cap "$((n*n*n))" --threads 1 --planes "$planes" \
      --replay "$case_dir/input/grid$n-seed.mesh" --output "$out/grid$n-$planes.jsonl" 2> "$out/grid$n-$planes.log"
    grep -q '^SEARCH_FINISHED candidates=1$' "$out/grid$n-$planes.log"
  done
done
# Positive exhaustive searches exercising cell-orbit sizes 1, 2 and 4.
for fixture in stack2 slab4; do
  cap=2; [[ $fixture != slab4 ]] || cap=4
  for planes in diagonal axial; do
    for fast in 0 1; do
      stem="$out/$fixture-$planes-fast$fast"
      "$binary" --input "$case_dir/input/$fixture-boundary.mesh" --cap "$cap" --threads 1 --planes "$planes" \
        --component-bound 0 --early-overlap "$fast" --state-trace "$stem.states" --output "$stem.jsonl" 2> "$stem.log"
      grep -q '^SEARCH_FINISHED candidates=1$' "$stem.log"
    done
    cmp "$out/$fixture-$planes-fast0.states" "$out/$fixture-$planes-fast1.states"
    cmp "$out/$fixture-$planes-fast0.jsonl" "$out/$fixture-$planes-fast1.jsonl"
  done
done
# Compare complete partial-state multisets, not only the final candidate count.
for planes in diagonal axial; do
  for fast in 0 1; do
    stem="$out/H12-$planes-fast$fast"
    "$binary" --input "$case_dir/input/pyramid.mesh" --cap 12 --threads 4 --planes "$planes" \
      --component-bound 0 --early-overlap "$fast" --seconds 120 --state-trace "$stem.states" --output "$stem.jsonl" 2> "$stem.log"
    grep -q '^SEARCH_FINISHED candidates=0$' "$stem.log"
    LC_ALL=C sort "$stem.states" > "$stem.sorted"
  done
  cmp "$out/H12-$planes-fast0.sorted" "$out/H12-$planes-fast1.sorted"
done
# Isolate the new parity optimization from the earlier overlap/budget checks.
for planes in diagonal axial; do
  for parity in 0 1; do
    stem="$out/H15-$planes-parity$parity"
    "$binary" --input "$case_dir/input/pyramid.mesh" --cap 15 --threads 4 --planes "$planes" \
      --component-bound 0 --early-parity "$parity" --seconds 120 --state-trace "$stem.states" --output "$stem.jsonl" 2> "$stem.log"
    grep -q '^SEARCH_FINISHED candidates=0$' "$stem.log"
    LC_ALL=C sort "$stem.states" > "$stem.sorted"
  done
  cmp "$out/H15-$planes-parity0.sorted" "$out/H15-$planes-parity1.sorted"
  grep -q ' reject_odd_cycle=0 ' "$out/H15-$planes-parity1.log"
done
# Completion precheck versus the full component calculation, with fixed ordering.
for early in 0 1; do
  stem="$out/H23-completion$early"
  "$binary" --input "$case_dir/input/pyramid.mesh" --cap 23 --threads 4 --early-completion "$early" \
    --seconds 120 --state-trace "$stem.states" --output "$stem.jsonl" 2> "$stem.log"
  grep -q '^SEARCH_FINISHED candidates=0$' "$stem.log"
  LC_ALL=C sort "$stem.states" > "$stem.sorted"
done
cmp "$out/H23-completion0.sorted" "$out/H23-completion1.sorted"
# Both known cube fillings must survive the new ordering and bounds.
for planes in diagonal axial; do
  for order in 0 3; do
    for bound in 0 1; do
      stem="$out/cube7-$planes-o$order-b$bound"
      "$binary" --input "$case_dir/input/cube.mesh" --cap 7 --threads 4 --planes "$planes" \
        --face-order "$order" --component-bound "$bound" --output "$stem.jsonl" 2> "$stem.log"
      grep -q '^SEARCH_FINISHED candidates=2$' "$stem.log"
      grep -q '^{"hexes":1,' "$stem.jsonl"
      grep -q '^{"hexes":7,' "$stem.jsonl"
    done
  done
done
for threads in 1 4; do
  "$binary" --input "$case_dir/input/pyramid.mesh" --cap 8 --threads "$threads" --seconds 120 \
    --output "$out/H8-t$threads.jsonl" 2> "$out/H8-t$threads.log"
  grep -q '^SEARCH_FINISHED candidates=0$' "$out/H8-t$threads.log"
done
for counter in branches orbit_trials states candidates; do
  one=$(sed -n "s/.* $counter=\([0-9]*\).*/\1/p" "$out/H8-t1.log" | tail -n 1)
  four=$(sed -n "s/.* $counter=\([0-9]*\).*/\1/p" "$out/H8-t4.log" | tail -n 1)
  [[ -n $one && $one == "$four" ]]
done
echo "Tests passed: seed replay, cap rejection, reflected grids, cube enumeration, thread consistency, overlap tables, rollback parity, completion bounds, and reference-state equivalence. Logs: $out"
