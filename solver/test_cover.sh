#!/usr/bin/env bash
set -euo pipefail
case_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
bash "$case_dir/run.sh" build
cmake --build "$case_dir/build" --target cover_test -- -j "${BUILD_JOBS:-2}"
"$case_dir/build/cover_test" "$case_dir/input"
binary="$case_dir/build/mirror_search"
out=$(mktemp -d "$case_dir/cover-tests-XXXXXX")
for mirrors in 1 2; do
  for plane in diagonal axial; do
    for cache in 0 1; do
      stem="$out/H18-m$mirrors-$plane-cache$cache"
      "$binary" --boundary-canonical 0 --input "$case_dir/input/pyramid.mesh" --cap 18 --mirrors "$mirrors" --planes "$plane" --threads 4 \
        --cover-cache "$cache" --state-trace "$stem.states" --output "$stem.jsonl" 2> "$stem.log"
      grep -q '^SEARCH_FINISHED candidates=0$' "$stem.log"
      LC_ALL=C sort "$stem.states" > "$stem.sorted"
    done
    cmp "$out/H18-m$mirrors-$plane-cache0.sorted" "$out/H18-m$mirrors-$plane-cache1.sorted"
    # Fix the traversal: order 4 uses the graph disabled by --face-cover 0.
    for cover in 0 1; do
      stem="$out/cube7-m$mirrors-$plane-cover$cover"
      "$binary" --boundary-canonical 0 --input "$case_dir/input/cube.mesh" --cap 7 --mirrors "$mirrors" --planes "$plane" --threads 4 \
        --face-order 3 --face-cover "$cover" --output "$stem.jsonl" 2> "$stem.log"
      grep -q '^SEARCH_FINISHED candidates=2$' "$stem.log"
      LC_ALL=C sort "$stem.jsonl" > "$stem.sorted"
    done
    cmp "$out/cube7-m$mirrors-$plane-cover0.sorted" "$out/cube7-m$mirrors-$plane-cover1.sorted"
  done
done
# A tiny work allowance must relax the cover check, preserving positive meshes.
for n in 36 44; do
  plane=diagonal; interior=33
  if [[ $n == 44 ]]; then plane=axial; interior=43; fi
  for mirrors in 1 2; do
    for work in 1 1000; do
      for order in 0 1 2 3; do
        stem="$out/replay$n-m$mirrors-w$work-o$order"
        "$binary" --boundary-canonical 0 --input "$case_dir/input/pyramid.mesh" --cap "$n" --interior "$interior" --mirrors "$mirrors" --planes "$plane" \
          --threads 1 --cover-work "$work" --face-order "$order" --replay "$case_dir/input/published$n-symmetric.mesh" \
          --output "$stem.jsonl" 2> "$stem.log"
        grep -q '^SEARCH_FINISHED candidates=1$' "$stem.log"
      done
    done
  done
done
echo "Cover tests passed: masked actions, known partial fillings, cache equivalence, positive candidates, and conservative work limits. Logs: $out"
