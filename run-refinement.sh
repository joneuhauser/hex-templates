#!/usr/bin/env bash
set -euo pipefail
root=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
export OPENBLAS_NUM_THREADS=1
args=(search --threads "${THREADS:?Set THREADS to the number of threads}" --case "${CASE:-all}" --planes "${PLANES:-axial}")
[[ -z ${CAP:-} ]] || args+=(--cap "$CAP")
[[ -z ${OUTPUT_DIR:-} ]] || args+=(--output "$OUTPUT_DIR")
exec python3 -u "$root/tools/refinement.py" "${args[@]}"
