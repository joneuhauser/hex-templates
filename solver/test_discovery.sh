#!/usr/bin/env bash
set -euo pipefail
case_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
bash "$case_dir/run.sh" build
out=$(mktemp -d "$case_dir/discovery-tests-XXXXXX")
binary="$case_dir/build/mirror_search"
"$binary" --input "$case_dir/input/pyramid.mesh" --mirrors 2 --planes diagonal --cap 36 --threads 1 \
  --replay "$case_dir/input/published36-symmetric.mesh" --output "$out/known36.jsonl" 2> "$out/replay.log"
"$binary" --input "$case_dir/input/pyramid.mesh" --mirrors 2 --planes diagonal --cap 36 --threads "${THREADS:-64}" \
  --output "$out/unseeded36.jsonl" 2> "$out/search.log"
grep -q '^SEARCH_FINISHED candidates=2$' "$out/search.log"
grep -F -x -q -f "$out/known36.jsonl" "$out/unseeded36.jsonl"
echo "Unseeded discovery passed: known 36-cell connectivity/actions found. Logs: $out"
