#!/usr/bin/env bash
set -euo pipefail
case_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
bash "$case_dir/run.sh" build
binary="$case_dir/build/mirror_search"
out=$(mktemp -d "$case_dir/fast-tests-XXXXXX")
for mirrors in 1 2; do
  for plane in diagonal axial; do
    for fast in 0 1; do
      stem="$out/H14-m$mirrors-$plane-fast$fast"
      "$binary" --boundary-canonical 0 --input "$case_dir/input/pyramid.mesh" --cap 14 --mirrors "$mirrors" --planes "$plane" --threads 4 \
        --vertex-domains "$fast" --color-domains "$fast" --geometry-domains "$fast" \
        --state-trace "$stem.states" --output "$stem.jsonl" 2> "$stem.log"
      grep -q '^SEARCH_FINISHED candidates=0$' "$stem.log"
      LC_ALL=C sort "$stem.states" > "$stem.sorted"
    done
    cmp "$out/H14-m$mirrors-$plane-fast0.sorted" "$out/H14-m$mirrors-$plane-fast1.sorted"
    for geometry in 0 1; do
      stem="$out/cube7-m$mirrors-$plane-geometry$geometry"
      "$binary" --boundary-canonical 0 --input "$case_dir/input/cube.mesh" --cap 7 --mirrors "$mirrors" --planes "$plane" --threads 4 \
        --corner-geometry "$geometry" --location-domains "$geometry" --mirror-separation "$geometry" \
        --output "$stem.jsonl" 2> "$stem.log"
      grep -q '^SEARCH_FINISHED candidates=2$' "$stem.log"
      LC_ALL=C sort "$stem.jsonl" > "$stem.sorted"
    done
    cmp "$out/cube7-m$mirrors-$plane-geometry0.sorted" "$out/cube7-m$mirrors-$plane-geometry1.sorted"
  done
done
echo "Fast-domain and geometry candidate tests passed. Logs: $out"
