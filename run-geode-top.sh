#!/usr/bin/env bash
set -euo pipefail
root=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
export OPENBLAS_NUM_THREADS=1
args=(--threads "${THREADS:?Set THREADS to the number of threads}" --symmetry "${SYMMETRY:-first}")
[[ -z ${CAP:-} ]] || args+=(--cap "$CAP")
[[ -z ${OUTPUT_DIR:-} ]] || args+=(--output "$OUTPUT_DIR")
[[ -z ${HEXOPT_BINARY:-} ]] || args+=(--hexopt "$HEXOPT_BINARY")
exec python3 -u "$root/tools/geode_top.py" "${args[@]}"
