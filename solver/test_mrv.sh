#!/usr/bin/env bash
set -euo pipefail
case_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
bash "$case_dir/run.sh" build
cmake --build "$case_dir/build" --target mrv_test -- -j "${BUILD_JOBS:-2}"
"$case_dir/build/mrv_test" "$case_dir/input"
binary="$case_dir/build/mirror_search"
out=$(mktemp -d "$case_dir/mrv-tests-XXXXXX")
for n in 36 44; do
  interior=33; [[ $n != 44 ]] || interior=43
  for partial in 5 6 7; do
    for work in 0 1 1000; do
      for quick in 0 1; do
        stem="$out/replay$n-p$partial-w$work-q$quick"
        "$binary" --input "$case_dir/input/pyramid.mesh" --mirrors 0 --cap "$n" --interior "$interior" \
          --threads 1 --face-order 5 --partial-cover "$partial" --partial-work "$work" --quick-packing "$quick" \
          --boundary-canonical 0 --replay "$case_dir/input/published$n-normalized.mesh" --output "$stem.jsonl" 2> "$stem.log"
        grep -q "REPLAY_OK cells=$n cell_orbits=$n" "$stem.log"
      done
    done
  done
done
for fixture in cube skew-cube-boundary; do
  ref="$out/$fixture-reference"
  "$binary" --input "$case_dir/input/$fixture.mesh" --mirrors 0 --cap 7 --interior 8 --threads 4 \
    --face-cover 0 --face-order 3 --partial-cover 0 --quick-packing 0 --output "$ref.jsonl" 2> "$ref.log"
  grep -q '^SEARCH_FINISHED candidates=2$' "$ref.log"
  for variant in default partial-off quick-off incremental-off; do
    args=()
    case "$variant" in
      partial-off) args=(--partial-cover 0) ;;
      quick-off) args=(--quick-packing 0) ;;
      incremental-off) args=(--cover-incremental 0) ;;
    esac
    stem="$out/$fixture-$variant"
    "$binary" --input "$case_dir/input/$fixture.mesh" --mirrors 0 --cap 7 --interior 8 --threads 4 \
      "${args[@]}" --output "$stem.jsonl" 2> "$stem.log"
    grep -q '^SEARCH_FINISHED candidates=2$' "$stem.log"
    python3 "$case_dir/compare_meshes.py" 8 "$ref.jsonl" "$stem.jsonl"
  done
done
for plane in none diagonal axial; do
  stem="$out/H16-$plane"
  "$binary" --input "$case_dir/input/pyramid.mesh" --mirrors 0 --planes "$plane" --cap 16 --threads 4 \
    --state-trace "$stem.states" --output "$stem.jsonl" 2> "$stem.log"
  grep -q '^SEARCH_FINISHED candidates=0$' "$stem.log"
  LC_ALL=C sort "$stem.states" > "$stem.sorted"
  cmp "$out/H16-none.sorted" "$stem.sorted"
done
for entry in 'partial-cover 4' 'partial-cover 8' 'partial-cover -1' 'quick-packing 2' 'partial-work -1' 'partial-work 1000001' 'face-order 6'; do
  read -r flag value <<< "$entry"
  if "$binary" "--$flag" "$value" --output "$out/invalid.jsonl" 2> "$out/invalid.log"; then
    echo 'ERROR: invalid option accepted.'; exit 1
  else [[ $? == 2 ]]; fi
done
echo "MRV tests passed: partial bounds on known completions, published replays, positive candidate equivalence, orientation independence, and option validation. Logs: $out"
