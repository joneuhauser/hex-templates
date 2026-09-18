#!/usr/bin/env bash
set -euo pipefail
root=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
for configuration in pyramid geode side-walls; do
  if [[ -n ${OUTPUT_DIR:-} ]]; then
    OUTPUT_DIR="$OUTPUT_DIR/$configuration" bash "$root/run-$configuration.sh"
  else
    bash "$root/run-$configuration.sh"
  fi
done
