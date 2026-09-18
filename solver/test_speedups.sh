#!/usr/bin/env bash
set -euo pipefail
case_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
bash "$case_dir/run.sh" build
cmake --build "$case_dir/build" --target speedup_test -- -j "${BUILD_JOBS:-2}"
"$case_dir/build/speedup_test"
binary="$case_dir/build/mirror_search"
out=$(mktemp -d "$case_dir/speedup-tests-XXXXXX")
for mirrors in 1 2; do
  for plane in diagonal axial; do
    for fast in 0 1; do
      stem="$out/H16-m$mirrors-$plane-fast$fast"
      "$binary" --boundary-canonical 0 --input "$case_dir/input/pyramid.mesh" --cap 16 --mirrors "$mirrors" --planes "$plane" --threads 4 \
        --edge-completion "$fast" --completion-domains "$fast" --incremental-faces "$fast" --component-index "$fast" --overlap-domains "$fast" \
        --state-trace "$stem.states" --output "$stem.jsonl" 2> "$stem.log"
      grep -q '^SEARCH_FINISHED candidates=0$' "$stem.log"
      LC_ALL=C sort "$stem.states" > "$stem.sorted"
      stem="$out/cube7-m$mirrors-$plane-fast$fast"
      "$binary" --boundary-canonical 0 --input "$case_dir/input/cube.mesh" --cap 7 --mirrors "$mirrors" --planes "$plane" --threads 4 \
        --edge-completion "$fast" --completion-domains "$fast" --incremental-faces "$fast" --component-index "$fast" --overlap-domains "$fast" \
        --output "$stem.jsonl" 2> "$stem.log"
      grep -q '^SEARCH_FINISHED candidates=2$' "$stem.log"
      LC_ALL=C sort "$stem.jsonl" > "$stem.sorted"
    done
    cmp "$out/H16-m$mirrors-$plane-fast0.sorted" "$out/H16-m$mirrors-$plane-fast1.sorted"
    cmp "$out/cube7-m$mirrors-$plane-fast0.sorted" "$out/cube7-m$mirrors-$plane-fast1.sorted"
  done
done
# Published positive meshes exercise the bound under several insertion orders.
for n in 36 44; do
  plane=diagonal; interior=33
  if [[ $n == 44 ]]; then plane=axial; interior=43; fi
  for order in 0 1 2 3; do
    stem="$out/prefix-replay-$n-o$order"
    "$binary" --boundary-canonical 0 --input "$case_dir/input/pyramid.mesh" --cap "$n" --interior "$interior" --mirrors 1 --planes "$plane" \
      --threads 1 --face-order "$order" --replay "$case_dir/input/published$n-symmetric.mesh" --output "$stem.jsonl" 2> "$stem.log"
    grep -q '^SEARCH_FINISHED candidates=1$' "$stem.log"
  done
done
echo "Speedup tests passed: independent ranks, partial-state multisets, and positive cube candidates. Logs: $out"
