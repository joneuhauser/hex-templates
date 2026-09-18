#!/usr/bin/env bash
set -euo pipefail
root=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
mkdir -p "$root/solver/runs"
out=${OUTPUT_DIR:-$(mktemp -d "$root/solver/runs/geode-quality.XXXXXX")}
for configuration in Q2-01 Q4-01 Q5-03; do
  OUTPUT_DIR="$out/$configuration" python3 -u "$root/tools/run_template.py" "geode-$configuration"
done
