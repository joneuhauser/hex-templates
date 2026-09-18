#!/usr/bin/env bash
set -euo pipefail
case_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
bash "$case_dir/run.sh" build
cmake --build "$case_dir/build" --target canonical_test -- -j "${BUILD_JOBS:-2}"
"$case_dir/build/canonical_test" "$case_dir/input"
python3 "$case_dir/test_boundary_orbits.py"
for flag in boundary-canonical leader-bound; do
  if "$case_dir/build/mirror_search" "--$flag" 2 > /dev/null 2>&1; then
    echo "Invalid $flag accepted"; exit 1
  else [[ $? == 2 ]]; fi
done
echo 'Canonical tests passed: representative retention, partial states, label invariance, bound propagation, published orbit replays, and option validation.'
