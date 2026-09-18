#!/usr/bin/env bash
set -euo pipefail
root=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd -P)
export OPENBLAS_NUM_THREADS=1
export BUILD_JOBS=${BUILD_JOBS:-2}
cmake -S "$root/solver" -B "$root/solver/build" -DCMAKE_BUILD_TYPE=Release -DNATIVE_ARCH=OFF
cmake --build "$root/solver/build" --parallel "$BUILD_JOBS"
for suite in test test_one_mirror test_published test_fast test_speedups test_cover test_global test_discovery test_zero test_mrv test_canonical test_support; do
  bash "$root/solver/$suite.sh"
done
for suite in test_one_orbits test_one_pillow test_launcher test_tree; do
  python3 "$root/solver/$suite.py"
done
python3 "$root/tests/test_published_coordinates.py"
python3 "$root/tests/test_templates.py"
python3 "$root/tests/test_site.py"
python3 "$root/tests/test_geode.py"
python3 "$root/tests/test_quad_sides.py"
python3 "$root/tests/test_geode_side_search.py"
python3 "$root/tests/test_half_turn.py"
python3 "$root/tests/test_geode_top.py"
python3 "$root/tests/test_launch_configs.py"
python3 "$root/tests/test_periodic_transition.py"

python3 "$root/tests/test_periodic_gallery.py"
