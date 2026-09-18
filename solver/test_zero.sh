#!/usr/bin/env bash
set -euo pipefail
case_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
bash "$case_dir/run.sh" build
cmake --build "$case_dir/build" --target zero_test -- -j "${BUILD_JOBS:-2}"
"$case_dir/build/zero_test" "$case_dir/input"
binary="$case_dir/build/mirror_search"
out=$(mktemp -d "$case_dir/zero-tests-XXXXXX")
for n in 36 44; do
  interior=33; [[ $n != 44 ]] || interior=43
  for order in 0 1 2 3 4; do
    for work in 1 1000; do
      stem="$out/replay$n-o$order-w$work"
      "$binary" --input "$case_dir/input/pyramid.mesh" --mirrors 0 --planes none --cap "$n" --interior "$interior" \
        --threads 1 --face-order "$order" --cover-work "$work" --boundary-canonical 0 --replay "$case_dir/input/published$n-normalized.mesh" \
        --output "$stem.jsonl" 2> "$stem.log"
      grep -q "REPLAY_OK cells=$n cell_orbits=$n" "$stem.log"
      grep -q '^SEARCH_FINISHED candidates=1$' "$stem.log"
    done
  done
done
# Disable the new packing warm start in these legacy exact-state checks:
# its stronger bound depends on ancestor reuse. test_mrv.sh checks complete
# candidate equivalence with it enabled and with reuse disabled.
# The orientation selector must have no effect when there is no reflection.
for plane in none diagonal axial; do
  stem="$out/H16-$plane"
  "$binary" --input "$case_dir/input/pyramid.mesh" --mirrors 0 --planes "$plane" --cap 16 --threads 4 \
    --quick-packing 0 --state-trace "$stem.states" --output "$stem.jsonl" 2> "$stem.log"
  grep -q '^SEARCH_FINISHED candidates=0$' "$stem.log"
  LC_ALL=C sort "$stem.states" > "$stem.sorted"
  cmp "$out/H16-none.sorted" "$stem.sorted"
done
for flag in edge-completion completion-domains cover-incremental cover-flat early-parity; do
  stem="$out/reference-$flag"
  "$binary" --input "$case_dir/input/pyramid.mesh" --mirrors 0 --cap 16 --threads 4 "--$flag" 0 \
    --quick-packing 0 --state-trace "$stem.states" --output "$stem.jsonl" 2> "$stem.log"
  grep -q '^SEARCH_FINISHED candidates=0$' "$stem.log"
  LC_ALL=C sort "$stem.states" > "$stem.sorted"
  cmp "$out/H16-none.sorted" "$stem.sorted"
done
for fixture in cube skew-cube-boundary; do
  for cover in 0 1; do
    stem="$out/$fixture-cover$cover"
    "$binary" --input "$case_dir/input/$fixture.mesh" --mirrors 0 --cap 7 --interior 8 --threads 4 --face-order 3 \
      --face-cover "$cover" --output "$stem.jsonl" 2> "$stem.log"
    grep -q '^SEARCH_FINISHED candidates=2$' "$stem.log"
    LC_ALL=C sort "$stem.jsonl" > "$stem.sorted"
  done
  cmp "$out/$fixture-cover0.sorted" "$out/$fixture-cover1.sorted"
  "$binary" --input "$case_dir/input/$fixture.mesh" --mirrors 0 --cap 7 --interior 8 --threads 4 --face-order 3 \
    --cover-fixed 0 --output "$out/$fixture-fixed0.jsonl" 2> "$out/$fixture-fixed0.log"
  grep -q '^SEARCH_FINISHED candidates=2$' "$out/$fixture-fixed0.log"
  LC_ALL=C sort "$out/$fixture-fixed0.jsonl" > "$out/$fixture-fixed0.sorted"
  cmp "$out/$fixture-cover1.sorted" "$out/$fixture-fixed0.sorted"
  "$binary" --input "$case_dir/input/$fixture.mesh" --mirrors 0 --cap 7 --interior 8 --threads 4 --root-face 2 \
    --output "$out/$fixture-root2.jsonl" 2> "$out/$fixture-root2.log"
  grep -q '^SEARCH_FINISHED candidates=2$' "$out/$fixture-root2.log"
  python3 "$case_dir/compare_meshes.py" 8 "$out/$fixture-cover1.jsonl" "$out/$fixture-root2.jsonl"
done
for mirrors in 1 2; do
  for plane in axial diagonal; do
    if "$binary" --input "$case_dir/input/skew-cube-boundary.mesh" --mirrors "$mirrors" --planes "$plane" --cap 7 \
      --threads 1 --output "$out/skew-rejected.jsonl" 2> "$out/skew-rejected.log"; then
      echo 'ERROR: a mirror was found in the asymmetric fixture.'; exit 1
    else [[ $? == 2 ]]; fi
    grep -q 'Input coordinates do not have the requested reflections' "$out/skew-rejected.log"
  done
done
echo "Zero-mirror tests passed: published seeds, asymmetric geometry, one-at-a-time allocation, completion prefixes, orientation independence, and exact state comparisons. Logs: $out"
