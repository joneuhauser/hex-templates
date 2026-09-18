#!/usr/bin/env bash
set -euo pipefail
threads=${1:?Usage: $0 THREADS}
[[ $threads =~ ^[1-9][0-9]*$ ]] || exit 2
cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.."
mkdir -p runs
out=${SEARCH_OUTPUT_DIR:-$(mktemp -d runs/side-walls-Q12.XXXXXX)}
mkdir -p -- "$out"
out=$(cd -- "$out" && pwd -P)
echo "Output: $out"
pids=()
for r in {0..10}; do
  (
    "${QUAD_SOLVER:-./build/quad_disk_search}" 12 "$r" "$out/r$r-orbits.jsonl" 2> "$out/r$r-orbits.log"
    if (( r <= 8 )); then
      "${QUAD_SOLVER:-./build/quad_disk_search}" 10 "$r" "$out/r$r-individual.jsonl" --no-orbits 2> "$out/r$r-individual.log"
    fi
  ) &
  pids+=("$!")
  if (( ${#pids[@]} == threads )); then
    for pid in "${pids[@]}"; do wait "$pid"; done
    pids=()
  fi
done
for pid in "${pids[@]}"; do wait "$pid"; done
