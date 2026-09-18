#!/usr/bin/env bash
set -euo pipefail
case_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
bash "$case_dir/run.sh" build
cmake --build "$case_dir/build" --target support_test -- -j "${BUILD_JOBS:-2}"
"$case_dir/build/support_test" "$case_dir/input"
"$case_dir/build/support_test" "$case_dir/input" 2
python3 "$case_dir/test_topology_symmetry.py" "$case_dir/build/mirror_search"
for flag in fast-pairs cover-closed cover-memo cover-fixed two-face-packing topology-symmetry; do
  if "$case_dir/build/mirror_search" "--$flag" 2 > /dev/null 2>&1; then
    echo "Invalid $flag accepted"; exit 1
  else [[ $? == 2 ]]; fi
done

for root in -2 16; do
  if "$case_dir/build/mirror_search" --input "$case_dir/input/pyramid.mesh" --root-face "$root" > /dev/null 2>&1; then
    echo "Invalid root face accepted"; exit 1
  else [[ $? == 2 ]]; fi
done
echo 'Support tests passed: reference comparisons for face pairs, partial-face support and cell intersections; known completions; option validation.'
