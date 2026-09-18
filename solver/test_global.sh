#!/usr/bin/env bash
set -euo pipefail
case_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
bash "$case_dir/run.sh" build
cmake --build "$case_dir/build" --target global_test -- -j "${BUILD_JOBS:-2}"
"$case_dir/build/global_test" "$case_dir/input"
binary="$case_dir/build/mirror_search"
out=$(mktemp -d "$case_dir/global-tests-XXXXXX")
for mirrors in 1 2; do
  for plane in diagonal axial; do
    for mode in reference symmetry incremental flat; do
      flags=(--cover-symmetry 0 --cover-incremental 0 --cover-cache 0)
      [[ $mode == reference ]] || flags=(--cover-symmetry 1 --cover-incremental 0 --cover-cache 0)
      [[ $mode != incremental ]] || flags=(--cover-symmetry 1 --cover-incremental 1 --cover-cache 0)
      # Force frequent replacement/collisions in a tiny exact-key cache.
      [[ $mode != flat ]] || flags=(--cover-symmetry 1 --cover-incremental 1 --cover-cache 1 --cover-flat 1 --cover-cache-bits 10)
      stem="$out/H18-m$mirrors-$plane-$mode"
      "$binary" --boundary-canonical 0 --input "$case_dir/input/pyramid.mesh" --cap 18 --mirrors "$mirrors" --planes "$plane" --threads 4 \
        --quick-packing 0 --face-order 3 --state-trace "$stem.states" --output "$stem.jsonl" "${flags[@]}" 2> "$stem.log"
      grep -q '^SEARCH_FINISHED candidates=0$' "$stem.log"
      LC_ALL=C sort "$stem.states" > "$stem.sorted"
      cmp "$out/H18-m$mirrors-$plane-reference.sorted" "$stem.sorted"
    done
    for order in 3 4; do
      for global in 0 1; do
        stem="$out/cube7-m$mirrors-$plane-o$order-g$global"
        "$binary" --boundary-canonical 0 --input "$case_dir/input/cube.mesh" --cap 7 --mirrors "$mirrors" --planes "$plane" --threads 4 \
          --face-order "$order" --cover-global "$global" --cover-seed "$global" --output "$stem.jsonl" 2> "$stem.log"
        grep -q '^SEARCH_FINISHED candidates=2$' "$stem.log"
        LC_ALL=C sort "$stem.jsonl" > "$stem.sorted"
      done
      cmp "$out/cube7-m$mirrors-$plane-o$order-g0.sorted" "$out/cube7-m$mirrors-$plane-o$order-g1.sorted"
    done
  done
done
for n in 36 44; do
  plane=diagonal; interior=33
  if [[ $n == 44 ]]; then plane=axial; interior=43; fi
  for mirrors in 1 2; do
    for work in 1 1000; do
      stem="$out/replay$n-m$mirrors-w$work-o4"
      "$binary" --boundary-canonical 0 --input "$case_dir/input/pyramid.mesh" --cap "$n" --interior "$interior" --mirrors "$mirrors" --planes "$plane" \
        --threads 1 --face-order 4 --cover-work "$work" --replay "$case_dir/input/published$n-symmetric.mesh" \
        --output "$stem.jsonl" 2> "$stem.log"
      grep -q '^SEARCH_FINISHED candidates=1$' "$stem.log"
    done
  done
done
echo "Global tests passed: independent intersections, incremental graph recomputation, exact state multisets, cache collisions, and positive meshes. Logs: $out"
