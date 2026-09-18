#!/usr/bin/env bash
set -euo pipefail
case_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
action=${1:-start}
integer() { [[ $2 =~ ^[0-9]+$ && ${#2} -le 7 ]] || { echo "$1 must be a nonnegative integer." >&2; exit 2; }; }
build() {
  for command in cmake c++ make sha256sum; do
    command -v "$command" >/dev/null || { echo "Missing dependency: $command" >&2; exit 2; }
  done
  cmake -S "$case_dir" -B "$case_dir/build" -DCMAKE_BUILD_TYPE=Release
  cmake --build "$case_dir/build" -- -j "${BUILD_JOBS:-8}"
}
status_write() { printf '%s\n' "$1" > "$run_dir/status.tmp"; mv -- "$run_dir/status.tmp" "$run_dir/status"; }
worker() {
  run_dir=$2; cap=$3; threads=$4; limit=$5; progress=$6; interior=$7
  mirrors=$(cat "$run_dir/mirrors")
  solver_pid=; stopped=0; incomplete=0; failure=0
  trap 'stopped=1; touch "$run_dir/stop_requested"; if [[ -n ${solver_pid:-} ]]; then kill -TERM "$solver_pid" 2>/dev/null || true; fi' TERM INT
  printf '%s\n' "$$" > "$run_dir/launcher.pid"
  status_write running
  while read -r plane; do
    [[ ! -f $run_dir/stop_requested ]] || { stopped=1; break; }
    printf '%s\n' "$plane" > "$run_dir/current_plane"
    args=(--input "$run_dir/pyramid.mesh" --cap "$cap" --threads "$threads" --planes "$plane"
          --mirrors "$mirrors" --seconds "$limit" --progress "$progress" --output "$run_dir/$plane.candidates.jsonl")
    args+=(--cover-work "$(cat "$run_dir/cover_work")")
    while read -r option value; do args+=("$option" "$value"); done < "$run_dir/options"
    [[ $interior == full ]] || args+=(--interior "$interior")
    "$run_dir/solver" "${args[@]}" > "$run_dir/$plane.stdout" 2> "$run_dir/$plane.log" &
    solver_pid=$!; printf '%s\n' "$solver_pid" > "$run_dir/solver.pid"
    rc=0; wait "$solver_pid" || rc=$?
    if (( stopped )); then wait "$solver_pid" 2>/dev/null || true; fi
    printf '%s\n' "$rc" > "$run_dir/$plane.exit_code"
    if [[ -f $run_dir/stop_requested ]]; then stopped=1; break; fi
    if (( rc == 124 )); then incomplete=1
    elif (( rc != 0 )) || ! grep -q '^SEARCH_FINISHED ' "$run_dir/$plane.log"; then failure=1; break
    fi
  done < "$run_dir/planes"
  if (( stopped )); then status_write stopped_incomplete
  elif (( failure )); then status_write failed_incomplete
  elif (( incomplete )); then status_write time_limit_incomplete
  else status_write enumeration_finished
  fi
  date -u +%FT%TZ > "$run_dir/finished_utc"
  if (( failure )); then return 2; fi
  if (( stopped )); then return 130; fi
  if (( incomplete )); then return 124; fi
}
select_run() {
  if [[ -n ${2:-} ]]; then run_dir=$(cd -- "$2" && pwd -P)
  elif [[ -f $case_dir/runs/latest ]]; then run_dir="$case_dir/runs/$(cat "$case_dir/runs/latest")"
  else echo 'No run exists.' >&2; exit 2
  fi
}
case "$action" in
  build) build ;;
  _worker) worker "$@" ;;
  start|foreground)
    cap=${CAP:-32}; threads=${THREADS:-64}; limit=${TIME_LIMIT:-0}; progress=${PROGRESS_SECONDS:-5}
    mirrors=${MIRRORS:-1}
    [[ $mirrors == 0 || $mirrors == 1 || $mirrors == 2 ]] || { echo "MIRRORS must be 0, 1, or 2." >&2; exit 2; }
    cover_work=${COVER_WORK:-1000}; [[ $mirrors != 0 ]] || cover_work=${COVER_WORK:-100000}
    [[ $cover_work =~ ^[0-9]+$ && ${#cover_work} -le 7 ]] || { echo "COVER_WORK must be a nonnegative integer." >&2; exit 2; }; cover_work=$((10#$cover_work)); (( cover_work<=1000000 )) || { echo "COVER_WORK must be at most 1000000." >&2; exit 2; }
    partial_cover=${PARTIAL_COVER:-6}; partial_work=${PARTIAL_WORK:-1000}; [[ $mirrors != 1 ]] || partial_work=${PARTIAL_WORK:-0}; face_order=${FACE_ORDER:-auto}
    [[ $mirrors != 1 || $face_order != auto ]] || face_order=5
    integer PARTIAL_COVER "$partial_cover"; integer PARTIAL_WORK "$partial_work"
    partial_cover=$((10#$partial_cover)); partial_work=$((10#$partial_work))
    (( (partial_cover==0 || (partial_cover>=5 && partial_cover<=7)) && partial_work<=1000000 )) || { echo 'Use PARTIAL_COVER=0,5,6,7 and PARTIAL_WORK=0..1000000.' >&2; exit 2; }
    if [[ $face_order != auto ]]; then integer FACE_ORDER "$face_order"; face_order=$((10#$face_order)); (( face_order<=5 )) || { echo 'Use FACE_ORDER=auto or 0..5.' >&2; exit 2; }; fi
    root_face=${ROOT_FACE:-auto}
    if [[ $root_face != auto ]]; then
      integer ROOT_FACE "$root_face"; root_face=$((10#$root_face))
      (( root_face<=15 )) || { echo 'Use ROOT_FACE=auto or 0..15 for the supplied pyramid.' >&2; exit 2; }
    fi
    planes=${PLANES:-both}; interior=${INTERIOR_CAP:-full}
    integer CAP "$cap"; integer THREADS "$threads"; integer TIME_LIMIT "$limit"; integer PROGRESS_SECONDS "$progress"
    cap=$((10#$cap)); threads=$((10#$threads)); limit=$((10#$limit)); progress=$((10#$progress))
    (( cap>=4 && cap<=48 && threads>=1 && threads<=256 )) || { echo 'Use CAP=4..48 and THREADS=1..256.' >&2; exit 2; }
    [[ $interior == full ]] || integer INTERIOR_CAP "$interior"
    if [[ $interior != full ]]; then interior=$((10#$interior)); fi
    if [[ $interior != full ]] && (( interior > 2*cap-7 )); then echo 'INTERIOR_CAP exceeds the full bound.' >&2; exit 2; fi
    case "$planes" in none|diagonal|axial|both) ;; *) echo 'PLANES must be none, diagonal, axial, or both.' >&2; exit 2 ;; esac
    if [[ $mirrors == 0 ]]; then planes=none
    elif [[ $planes == none ]]; then echo 'PLANES=none requires MIRRORS=0.' >&2; exit 2; fi
    build
    mkdir -p -- "$case_dir/runs"
    run_dir=$(mktemp -d "$case_dir/runs/H${cap}-$(date -u +%Y%m%dT%H%M%SZ)-XXXXXX")
    printf '%s\n' "$cover_work" > "$run_dir/cover_work"
    printf -- '--partial-cover %s\n--partial-work %s\n' "$partial_cover" "$partial_work" > "$run_dir/options"
    [[ $root_face == auto ]] || printf -- '--root-face %s\n' "$root_face" >> "$run_dir/options"
    [[ $face_order == auto ]] || printf -- '--face-order %s\n' "$face_order" >> "$run_dir/options"
    printf '%s\n' "$mirrors" > "$run_dir/mirrors"
    cp -- "$case_dir/build/mirror_search" "$run_dir/solver"
    cp -- "$case_dir/input/pyramid.mesh" "$run_dir/pyramid.mesh"
    if [[ $planes == both ]]; then printf 'axial\ndiagonal\n' > "$run_dir/planes"
    else printf '%s\n' "$planes" > "$run_dir/planes"; fi
    while read -r plane; do touch "$run_dir/$plane.log"; done < "$run_dir/planes"
    {
      printf 'started_utc=%s\nhost=%s\n' "$(date -u +%FT%TZ)" "$(hostname)"
      printf 'cap=%s\nthreads=%s\nplanes=%s\ninterior_cap=%s\nfull_interior_bound=%s\ntime_limit_per_plane=%s\n' \
        "$cap" "$threads" "$planes" "$interior" "$((2*cap-7))" "$limit"
      cat "$run_dir/options"
      printf 'mirrors=%s\ncover_work=%s\n' "$mirrors" "$cover_work"
      if [[ $mirrors == 0 ]]; then printf 'enforcement=individual_cells\n'; else printf 'enforcement=cell_orbits\n'; fi
      printf 'geometry_solver=not_implemented\ncheckpointing=not_implemented\n'
      sha256sum "$run_dir/solver" "$run_dir/pyramid.mesh"
    } > "$run_dir/metadata.txt"
    status_write launching
    basename -- "$run_dir" > "$case_dir/runs/latest.tmp"; mv -- "$case_dir/runs/latest.tmp" "$case_dir/runs/latest"
    echo "Run directory: $run_dir"
    if [[ $action == foreground ]]; then worker _worker "$run_dir" "$cap" "$threads" "$limit" "$progress" "$interior"
    else
      nohup bash "$case_dir/run.sh" _worker "$run_dir" "$cap" "$threads" "$limit" "$progress" "$interior" \
        > "$run_dir/launcher.log" 2>&1 < /dev/null &
      echo "Detached launcher PID: $!"
      echo 'Use bash run.sh status / follow / stop. Zero mirrors runs once; otherwise axial and diagonal cases run sequentially.'
    fi
    ;;
  status|follow|stop)
    select_run "$@"
    if [[ $action == follow ]]; then
      logs=(); while read -r plane; do logs+=("$run_dir/$plane.log"); done < "$run_dir/planes"
      exec tail -n 3 -f "${logs[@]}"
    elif [[ $action == stop ]]; then
      [[ -f $run_dir/solver.pid ]] || { echo 'No solver PID recorded yet.' >&2; exit 2; }
      pid=$(cat "$run_dir/solver.pid")
      [[ $pid =~ ^[1-9][0-9]*$ ]] || exit 2
      [[ $(readlink "/proc/$pid/exe" 2>/dev/null || true) == "$run_dir/solver" ]] || { echo 'No matching live solver.' >&2; exit 2; }
      touch "$run_dir/stop_requested"; kill -TERM "$pid"
      echo 'Stop requested. Partial searches cannot resume.'
    else
      echo "Run directory: $run_dir"; cat "$run_dir/status"
      while read -r plane; do echo "Mirror orientation: $plane"; tail -n 2 "$run_dir/$plane.log"; done < "$run_dir/planes"
      if [[ -f $run_dir/solver.pid ]]; then
        pid=$(cat "$run_dir/solver.pid")
        if [[ $pid =~ ^[1-9][0-9]*$ && $(readlink "/proc/$pid/exe" 2>/dev/null || true) == "$run_dir/solver" ]]; then
          ps -p "$pid" -o pid,nlwp,etime,pcpu,rss,comm
        elif [[ $(cat "$run_dir/status") == running ]]; then echo 'No matching live solver; status may be stale.'; fi
      fi
    fi
    ;;
  *) echo 'Usage: bash run.sh {build|start|foreground|status|follow|stop} [run-directory]' >&2; exit 2 ;;
esac
