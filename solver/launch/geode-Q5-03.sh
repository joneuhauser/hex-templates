#!/usr/bin/env bash
set -euo pipefail
threads=${1:?Usage: $0 THREADS}
cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.."
mkdir -p runs
out=${SEARCH_OUTPUT_DIR:-$(mktemp -d "runs/geode-Q5-03-H38.XXXXXX")}
mkdir -p -- "$out"
out=$(cd -- "$out" && pwd -P)
cp "input/geode-Q5-03.mesh" "$out/boundary.mesh"
echo "Output: $out"
exec "${SOLVER:-./build/mirror_search}" --input "$out/boundary.mesh" \
  --threads "$threads" --cap 38 --mirrors 2 --planes diagonal --seconds 0 \
  --boundary-canonical 0 --output "$out/candidates.jsonl" \
  --corner-geometry 0 --geometry-domains 0 --location-domains 0
