#!/usr/bin/env bash
set -euo pipefail
threads=${1:?Usage: $0 THREADS}
cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.."
mkdir -p runs
out=$(mktemp -d "runs/geode-Q7-06-flipped-H44.XXXXXX")
cp "input/geode-Q7-06-flipped.mesh" "$out/boundary.mesh"
echo "Output: $PWD/$out"
exec "${SOLVER:-./build/mirror_search}" --input "$out/boundary.mesh" \
  --threads "$threads" --cap 44 --mirrors 2 --planes diagonal --seconds 0 \
  --boundary-canonical 0 --output "$out/candidates.jsonl" \
  --corner-geometry 0 --geometry-domains 0 --location-domains 0
