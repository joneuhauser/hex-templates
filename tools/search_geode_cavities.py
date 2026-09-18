#!/usr/bin/env python3
"""Find and embed alternative Geodes by unrestricted local cavity searches.

Enumerates fillings of three/four-cell face-adjacency cliques with caps five/four.
This finite neighborhood is not an exhaustive search over Geode templates.
"""
import argparse
from itertools import combinations
import json
from pathlib import Path
import subprocess
import sys

import numpy as np

from mesh_tools import EDGES, boundary, owners, read_mesh, topology, write_mesh, write_vtu
from geode_geometry import improve, validate
from prepare_geode import ROOT

sys.path.insert(0, str(ROOT / 'solver'))
from compare_meshes import equivalent


def splice(points, cells, quads, removed, candidate, boundary_ids):
    """Install a cavity filling, harmonically initialize new vertices, compact IDs."""
    nb = len(boundary_ids)
    count = candidate['vertices']
    local = np.asarray(candidate['cells'], dtype=int)
    p = np.vstack([points, np.zeros((count-nb, 3))])
    indices = np.asarray(boundary_ids + list(range(len(points), len(p))))
    if count > nb:
        laplacian = np.zeros((count, count))
        edges = {tuple(sorted(edge)) for cell in local for edge in cell[EDGES]}
        for i, j in edges:
            laplacian[i, i] += 1
            laplacian[j, j] += 1
            laplacian[i, j] -= 1
            laplacian[j, i] -= 1
        p[len(points):] = np.linalg.solve(laplacian[nb:, nb:], -laplacian[nb:, :nb] @ points[boundary_ids])
    h = np.vstack([np.delete(cells, removed, axis=0), indices[local]])
    used = sorted(set(h.flat))
    mapping = {v: i for i, v in enumerate(used)}
    h = np.asarray([[mapping[v] for v in cell] for cell in h])
    q = np.asarray([[mapping[v] for v in face] for face in quads])
    return p[used], h, q


def run(reference, binary, output, seconds=3, iterations=150):
    p, h, q = read_mesh(reference)
    output.mkdir(parents=True, exist_ok=True)
    adjacent = [set() for _ in h]
    for entries in owners(h).values():
        if len(entries) == 2:
            a, b = [entry[0] for entry in entries]
            adjacent[a].add(b)
            adjacent[b].add(a)
    records = []
    for size, cap in [(4, 4), (3, 5)]:
        for ids in combinations(range(len(h)), size):
            if not all(j in adjacent[i] for i, j in combinations(ids, 2)):
                continue
            old = h[list(ids)]
            faces = boundary(old)
            bv = sorted(map(int, set(faces.flat)))
            order = bv + sorted(set(map(int, old.flat)) - set(bv))
            mapping = {v: i for i, v in enumerate(order)}
            local = np.asarray([[mapping[v] for v in cell] for cell in old])
            faces = np.asarray([[mapping[v] for v in face] for face in faces])
            stem = output / '-'.join(map(str, ids))
            write_mesh(stem.with_suffix('.mesh'), p[bv], np.empty((0, 8), dtype=int), faces)
            command = [str(binary), '--input', str(stem.with_suffix('.mesh')),
                       '--mirrors', '0', '--planes', 'none', '--cap', str(cap),
                       '--threads', '1', '--seconds', str(seconds), '--boundary-canonical', '0',
                       '--corner-geometry', '0', '--geometry-domains', '0', '--location-domains', '0',
                       '--output', str(stem.with_suffix('.jsonl'))]
            with stem.with_suffix('.log').open('w') as log:
                completed = subprocess.run(command, stdout=subprocess.DEVNULL, stderr=log, check=False)
            if completed.returncode not in [0, 124]:
                raise RuntimeError(f'Solver failed; see {stem.with_suffix(".log")}')
            record = dict(cavity=list(ids), cap=cap, command=command,
                          enumeration_complete=completed.returncode == 0,
                          alternatives=0, incompatible_splices=0, accepted=[])
            original = dict(vertices=len(order), hexes=size, cells=local.tolist())
            for line in stem.with_suffix('.jsonl').read_text().splitlines():
                candidate = json.loads(line)
                if equivalent(candidate, original, len(bv)):
                    continue
                index = record['alternatives']
                record['alternatives'] += 1
                name = output / f'{stem.name}-alternative-{index}'
                name.with_suffix('.json').write_text(json.dumps(dict(candidate,
                    cavity_ids=ids, boundary_global_ids=bv), indent=2) + '\n')
                pp, hh, qq = splice(p, h, q, ids, candidate, bv)
                check = topology(hh, qq)
                if check['incompatible_pairs'] or not check['boundary_matches']:
                    record['incompatible_splices'] += 1
                    continue
                pp, history = improve(pp, hh, qq, iterations)
                write_mesh(name.with_suffix('.mesh'), pp, hh, qq)
                report = validate(name.with_suffix('.mesh'), reference)
                report['optimization'] = history
                name.with_suffix('.validation.json').write_text(json.dumps(report, indent=2) + '\n')
                if report['accepted']:
                    record['accepted'].append(name.name)
                    write_vtu(name.with_suffix('.vtu'), pp, hh)
            records.append(record)
            (output / 'summary.json').write_text(json.dumps(records, indent=2) + '\n')
            print(json.dumps({k: v for k, v in record.items() if k != 'command'}), flush=True)
    return records


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--reference', type=Path, default=ROOT / 'docs/assets/geode/meshes/G26-Mitchell.mesh')
    parser.add_argument('--binary', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--seconds', type=float, default=3)
    parser.add_argument('--iterations', type=int, default=150)
    args = parser.parse_args()
    if args.seconds <= 0 or args.iterations < 1:
        parser.error('Use positive seconds and iterations')
    run(args.reference.resolve(), args.binary.resolve(), args.output.resolve(), args.seconds, args.iterations)
