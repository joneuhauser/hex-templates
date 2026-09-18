#!/usr/bin/env python3
"""Embed Geode or pyramid candidates with HexOpt and certify their geometry."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys

import numpy as np

from geode_geometry import variables, validate
from mesh_tools import EDGES, read_mesh, write_mesh, write_vtu, sampled
from quality_metrics import dense
from prepare_geode import ROOT

sys.path.insert(0, str(ROOT/'solver'))
from compare_meshes import equivalent


def initialize(points, candidate):
    nb, n = len(points), candidate['vertices']
    cells = np.asarray(candidate['cells'], dtype=int)
    laplacian = np.zeros((n, n))
    for a, b in {tuple(sorted(e)) for c in cells for e in c[EDGES]}:
        laplacian[a, a] += 1
        laplacian[b, b] += 1
        laplacian[a, b] -= 1
        laplacian[b, a] -= 1
    interior = np.linalg.solve(laplacian[nb:, nb:], -laplacian[nb:, :nb] @ points)
    return np.vstack([points, interior]), cells


def collect(roots):
    classes = []
    for root in roots:
        for path in sorted(root.rglob('candidates.jsonl')):
            boundary_path = path.parent/'boundary.mesh'
            digest = hashlib.sha256(boundary_path.read_bytes()).hexdigest()
            points, _, quads = read_mesh(boundary_path)
            for line, row in enumerate(path.read_text().splitlines(), 1):
                candidate = json.loads(row)
                source = dict(path=str(path.relative_to(ROOT)) if path.is_relative_to(ROOT)
                              else str(path), line=line)
                found = next((c for c in classes if c['boundary_sha256'] == digest
                              and equivalent(candidate, c['candidate'], len(points))), None)
                if found is not None:
                    found['sources'].append(source)
                    continue
                classes.append(dict(side=path.parent.name, boundary_sha256=digest,
                    boundary_path=boundary_path, candidate=candidate, sources=[source]))
    classes.sort(key=lambda c: (c['side'], c['candidate']['hexes'], c['candidate']['cells']))
    counts = {}
    for c in classes:
        stem = f'{c["side"]}-H{c["candidate"]["hexes"]}'
        counts[stem] = counts.get(stem, 0)+1
        c['id'] = f'{stem}-{counts[stem]:02d}'
    return classes


def optimize(record, args):
    directory = args.output/record['id']
    directory.mkdir()
    candidate = record['candidate']
    bp, _, quads = read_mesh(record['boundary_path'])
    points, cells = initialize(bp, candidate)
    actions = np.asarray(candidate['vertex_actions'])
    pyramid = getattr(args, 'profile', 'geode') == 'pyramid'
    if pyramid:
        # Interior starting coordinates are free; project each restart onto the
        # diagonal symmetry group before using the fixed-boundary adapter.
        x0 = points[len(bp):].ravel().copy()
        matrix = np.zeros((points.size, len(x0)))
        matrix[3*len(bp):] = np.eye(len(x0))
        offset = points.ravel()-matrix@x0
        bounds = [(-1., 1.), (-1., 1.), (0., np.sqrt(2))]*(len(points)-len(bp))
    else:
        offset, matrix, x0, bounds = variables(points, quads, vertex_actions=actions)
    bounds = np.asarray(bounds)
    fixed_columns = np.any(matrix[:3*len(bp)] != 0, axis=0)
    free_sides = bool(np.any(fixed_columns))
    assert np.max(np.sum(matrix != 0, axis=1)) <= 1, 'Adapter requires disjoint affine columns'
    saved = {k: v for k, v in record.items() if k != 'boundary_path'}
    (directory/'candidate.json').write_text(json.dumps(saved, indent=2)+'\n')
    write_mesh(directory/'initial.mesh', points, cells, quads)
    attempts = []
    rng = np.random.default_rng(int(hashlib.sha256(record['id'].encode()).hexdigest()[:8], 16))
    for restart in range(args.restarts):
        for mode in (['fixed', 'free'] if free_sides else ['fixed']):
            parameters = x0.copy()
            if restart:
                noise = rng.normal(0, min(0.025*restart, 0.15), size=len(parameters))
                if mode == 'fixed':
                    noise[fixed_columns] = 0
                parameters += noise
            parameters = np.clip(parameters, bounds[:, 0], bounds[:, 1])
            start = (offset+matrix@parameters).reshape(-1, 3)
            if pyramid:
                images = []
                for g in range(4):
                    image = start[actions[:, g]].copy()
                    if g & 1:
                        image[:, :2] = image[:, [1, 0]]
                    if g & 2:
                        image[:, :2] = -image[:, [1, 0]]
                    images.append(image)
                start = np.mean(images, axis=0)
                start[:len(bp)] = bp
            stem = directory/f'{mode}-{restart}'
            with stem.with_suffix('.txt').open('w') as stream:
                stream.write(f'{len(points)} {len(cells)} {len(bp)} 0\n')
                np.savetxt(stream, start, fmt='%.17g')
                np.savetxt(stream, actions, fmt='%d')
                np.savetxt(stream, cells, fmt='%d')
            command = [str(args.binary), str(stem.with_suffix('.txt')), str(stem.with_suffix('.xyz'))]
            if mode == 'free':
                with stem.with_suffix('.constraints').open('w') as stream:
                    stream.write(f'{len(parameters)}\n')
                    for j, value in enumerate(parameters):
                        rows = np.flatnonzero(matrix[:, j])
                        stream.write(f'{value:.17g} {bounds[j, 0]:.17g} {bounds[j, 1]:.17g} {len(rows)}')
                        for row in rows:
                            stream.write(f' {row//3} {row%3} {matrix[row, j]:.17g}')
                        stream.write('\n')
                command.append(str(stem.with_suffix('.constraints')))
            try:
                with stem.with_suffix('.log').open('w') as log:
                    result = subprocess.run(command, stdout=log, stderr=log, timeout=120,
                        env=dict(os.environ, HEXOPT_STALL_LIMIT=str(args.stall),
                                 HEXOPT_MAX_ITERATIONS=str(args.iterations)))
            except subprocess.TimeoutExpired:
                attempts.append(dict(mode=mode, restart=restart, timed_out=True))
                continue
            attempt = dict(mode=mode, restart=restart, exit_code=result.returncode)
            attempts.append(attempt)
            if result.returncode:
                continue
            optimized = np.loadtxt(stem.with_suffix('.xyz'))
            with np.errstate(divide='ignore', invalid='ignore'):
                quality = sampled(optimized, cells, 9)
            attempt['sampled'] = {k: (None if isinstance(v, float) and not np.isfinite(v) else v)
                                  for k, v in quality.items()}
            if quality['min_det'] <= 0 or not np.isfinite(quality['min_scaled_jacobian']):
                continue
            trial = stem.with_suffix('.mesh')
            write_mesh(trial, optimized, cells, quads)
            if pyramid:
                from validate_mesh import validate as validate_pyramid
                validation = validate_pyramid(trial)
            else:
                validation = validate(trial, args.reference)
            mirror_error = 0.
            for g in range(4):
                image = optimized.copy()
                if g & 1:
                    image[:, :2] = image[:, [1, 0]]
                if g & 2:
                    image[:, :2] = -image[:, [1, 0]]
                mirror_error = max(mirror_error, float(np.max(abs(image-optimized[actions[:, g]]))))
            validation['diagonal_mirror_error'] = mirror_error
            validation['accepted'] &= mirror_error < 1e-12
            attempt['accepted'] = validation['accepted']
            stem.with_suffix('.validation.json').write_text(json.dumps(validation, indent=2)+'\n')
            if validation['accepted']:
                write_mesh(directory/'accepted.mesh', optimized, cells, quads)
                write_vtu(directory/'accepted.vtu', optimized, cells)
                (directory/'validation.json').write_text(json.dumps(validation, indent=2)+'\n')
                quality = dense(optimized, cells, True)
                (directory/'quality.json').write_text(json.dumps(quality, indent=2)+'\n')
                summary = dict(id=record['id'], side=record['side'], hexes=len(cells), accepted=True,
                    mode=mode, restart=restart, attempts=attempts, quality=quality)
                (directory/'result.json').write_text(json.dumps(summary, indent=2)+'\n')
                print(json.dumps(dict(id=record['id'], accepted=True,
                    min_scaled_jacobian=quality['min_scaled_jacobian'], mode=mode)), flush=True)
                return summary
    summary = dict(id=record['id'], side=record['side'], hexes=len(cells), accepted=False, attempts=attempts)
    (directory/'result.json').write_text(json.dumps(summary, indent=2)+'\n')
    print(json.dumps(dict(id=record['id'], accepted=False, attempts=len(attempts))), flush=True)
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('runs', nargs='+', type=Path)
    parser.add_argument('--binary', required=True, type=Path)
    parser.add_argument('--profile', choices=['geode', 'pyramid'], default='geode')
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--reference', type=Path, default=ROOT/'docs/assets/geode/meshes/G26-Mitchell.mesh')
    parser.add_argument('--restarts', type=int, default=6)
    parser.add_argument('--stall', type=int, default=5000)
    parser.add_argument('--iterations', type=int, default=300000)
    parser.add_argument('--only-failed', type=Path, help='Retry failed class IDs from this earlier output directory')
    parser.add_argument('--jobs', type=int, default=4)
    args = parser.parse_args()
    if min(args.restarts, args.jobs, args.stall, args.iterations) < 1:
        parser.error('Use positive restarts, jobs and iteration limits')
    args.runs = [path.resolve() for path in args.runs]
    args.binary = args.binary.resolve()
    args.output = args.output.resolve()
    args.output.mkdir(parents=True, exist_ok=False)
    classes = collect(args.runs)
    if args.only_failed:
        failed = {r['id'] for r in json.loads((args.only_failed/'results.json').read_text()) if not r['accepted']}
        classes = [c for c in classes if c['id'] in failed]
    metadata = dict(classes=len(classes), profile=args.profile, restarts=args.restarts,
        stall_limit=args.stall, max_iterations=args.iterations,
        retry_from=str(args.only_failed) if args.only_failed else None,
        optimizer='HexOpt sJGrad with repository projected-gradient adapter',
        sha256={name: hashlib.sha256(path.read_bytes()).hexdigest() for name, path in
                [('binary', args.binary), ('adapter', ROOT/'tools/hexopt_fixed.cpp'),
                 ('driver', Path(__file__)), ('reference', args.reference)]})
    (args.output/'metadata.json').write_text(json.dumps(metadata, indent=2)+'\n')
    print(f'Optimizing {len(classes)} connectivity classes', flush=True)
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        results = list(pool.map(lambda record: optimize(record, args), classes))
    (args.output/'results.json').write_text(json.dumps(results, indent=2)+'\n')
    print(json.dumps(dict(classes=len(results), accepted=sum(r['accepted'] for r in results))), flush=True)


if __name__ == '__main__':
    main()
