#!/usr/bin/env python3
"""Try diagonal-mirror or half-turn Geode fillings for catalog side patches.

Search output is topology only; any candidates still need geometric embedding
and exact validation. Vertical flips are distinct placements against fixed caps.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
from pathlib import Path
import re
import subprocess
import time

import numpy as np

from build_geode_sides import flipped, harmonic
from geode_boundary import assemble
from mesh_tools import write_mesh
from prepare_geode import face_set

ROOT = Path(__file__).resolve().parents[1]


def placements(catalog, max_quads):
    for side in catalog['patterns']:
        if len(side['quads']) > max_quads:
            continue
        yield side['id'], side, False
        if side['vertical_flip_class_size'] == 2:
            quads = flipped(side['quads'], side['boundary_vertices'])
            points = harmonic(quads, side['side_vertices']-2)
            other = dict(side, quads=quads,
                         points=[[float(x), float(y)*7/12] for x, y in points])
            yield side['id']+'-flipped', other, True


def diagonal_action(points, quads, second=False):
    reflected = points[:, [1, 0, 2]].copy()
    if second:
        reflected[:, :2] *= -1
    action = []
    for point in reflected:
        matches = np.flatnonzero(np.max(abs(points-point), axis=1) < 1e-12)
        if len(matches) != 1:
            raise ValueError('Boundary vertices do not preserve the diagonal mirror')
        action.append(int(matches[0]))
    if face_set(np.asarray(action)[quads]) != face_set(quads):
        raise ValueError('Boundary quads do not preserve the diagonal mirror')
    return action


def diagonal_actions(points, quads, mirrors):
    first = np.asarray(diagonal_action(points, quads))
    actions = [np.arange(len(points)), first]
    if mirrors == 2:
        second = np.asarray(diagonal_action(points, quads, second=True))
        if not np.array_equal(first[second], second[first]):
            raise ValueError('Diagonal reflections do not commute')
        actions.extend([second, first[second]])
    elif mirrors != 1:
        raise ValueError('Expected one or two diagonal mirrors')
    return np.asarray(actions).T.tolist()


def run_case(item, args):
    name, side, is_flipped = item
    case = args.output/name
    case.mkdir()
    points, quads = assemble(args.reference, side)
    half_turn = getattr(args, 'half_turn', False)
    actions = diagonal_actions(points, quads, 2 if half_turn else args.mirrors)
    if half_turn:
        actions = [[row[0], row[3]] for row in actions]
    write_mesh(case/'boundary.mesh', points, np.empty((0, 8), dtype=int), quads)
    # Only cap vertices are fixed. Using provisional side coordinates in
    # geometric pruning would incorrectly exclude other side embeddings.
    fixed_boundary = bool(np.all(abs(abs(points[:, 2])-7/12) < 1e-12))
    geometry = '1' if fixed_boundary else '0'
    command = [str(args.solver), '--input', str(case/'boundary.mesh'),
               *(['--half-turn', '1'] if half_turn else ['--mirrors', str(args.mirrors), '--planes', 'diagonal']), '--cap', str(args.cap),
               '--threads', str(args.threads), '--seconds', str(args.seconds),
               '--progress', '10', '--boundary-canonical', '0',
               '--corner-geometry', geometry, '--geometry-domains', geometry,
               '--location-domains', geometry, '--output', str(case/'candidates.jsonl')]
    record = dict(id=name, side_id=side['id'], flipped_placement=is_flipped,
                  cap=args.cap, side_quads=len(side['quads']), boundary_vertices=len(points),
                  boundary_quads=len(quads), full_interior_bound=2*args.cap+1-len(quads)//2,
                  mirrors=0 if half_turn else args.mirrors, half_turn=half_turn, vertex_actions=actions,
                  fixed_boundary_geometry=fixed_boundary, command=command,
                  boundary_sha256=hashlib.sha256((case/'boundary.mesh').read_bytes()).hexdigest())
    (case/'input.json').write_text(json.dumps(record, indent=2)+'\n')
    started = time.monotonic()
    with (case/'search.log').open('w') as log, (case/'stdout.log').open('w') as stdout:
        result = subprocess.run(command, stdout=stdout, stderr=log)
    log = (case/'search.log').read_text()
    progress = [line for line in log.splitlines() if line.startswith('PROGRESS ')]
    counters = dict(re.findall(r'(\w+)=([^ ]+)', progress[-1])) if progress else {}
    record.update(exit_code=result.returncode, wall_seconds=time.monotonic()-started,
                  complete=result.returncode == 0 and '\nSEARCH_FINISHED ' in log,
                  candidates=sum(bool(line.strip()) for line in (case/'candidates.jsonl').read_text().splitlines()),
                  final_counters=counters)
    (case/'result.json').write_text(json.dumps(record, indent=2)+'\n')
    print(json.dumps({k: record[k] for k in ['id', 'cap', 'complete', 'candidates', 'wall_seconds']}), flush=True)
    if result.returncode not in (0, 124):
        raise RuntimeError(f'Solver failed; see {case / "search.log"}')
    return record


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--solver', type=Path, default=ROOT/'solver/build/mirror_search')
    parser.add_argument('--reference', type=Path, default=ROOT/'docs/assets/geode/meshes/G26-Mitchell.mesh')
    parser.add_argument('--catalog', type=Path, default=ROOT/'docs/assets/geode-sides/catalog.json')
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--max-side-quads', type=int, default=6)
    parser.add_argument('--ids', nargs='+', help='Select placement IDs, e.g. Q2-01 Q5-01-flipped')
    parser.add_argument('--cap', type=int, default=26)
    parser.add_argument('--mirrors', type=int, choices=[1, 2], default=1)
    parser.add_argument('--half-turn', action='store_true', help='Constrain only the 180-degree rotation; no mirrors')
    parser.add_argument('--seconds', type=int, default=120)
    parser.add_argument('--threads', type=int, default=8)
    parser.add_argument('--jobs', type=int, default=3)
    args = parser.parse_args()
    if min(args.cap, args.threads, args.jobs) < 1 or args.seconds < 0:
        parser.error('Use positive cap/threads/jobs and nonnegative seconds')
    for key in ['solver', 'reference', 'catalog', 'output']:
        setattr(args, key, getattr(args, key).resolve())
    items = list(placements(json.loads(args.catalog.read_text()), args.max_side_quads))
    if args.ids:
        missing = set(args.ids)-{name for name, _, _ in items}
        if missing:
            parser.error(f'Unknown placement IDs: {sorted(missing)}')
        items = [item for item in items if item[0] in args.ids]
    if not items:
        parser.error('No selected side patterns')
    args.output.mkdir(parents=True, exist_ok=False)
    metadata = dict(mirrors=0 if args.half_turn else args.mirrors, half_turn=args.half_turn,
                    planes=[] if args.half_turn else ['x=y', 'x=-y'][:args.mirrors],
                    cap=args.cap, seconds_per_case=args.seconds,
                    threads_per_case=args.threads, jobs=args.jobs,
                    placements=[name for name, _, _ in items],
                    sha256={name: hashlib.sha256(path.read_bytes()).hexdigest() for name, path in
                            [('solver', args.solver), ('reference', args.reference),
                             ('catalog', args.catalog), ('driver', Path(__file__))]})
    (args.output/'metadata.json').write_text(json.dumps(metadata, indent=2)+'\n')
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        records = list(pool.map(lambda item: run_case(item, args), items))
    (args.output/'results.json').write_text(json.dumps(records, indent=2)+'\n')


if __name__ == '__main__':
    main()
