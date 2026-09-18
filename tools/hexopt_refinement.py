"""Optimize refinement layer height and interior coordinates with HexOpt."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
import os
from pathlib import Path
import subprocess

import numpy as np

from mesh_tools import read_mesh, write_mesh, write_vtu
from quality_metrics import dense
from refinement import INPUT, ROOT, validate

START_HEIGHTS = (.125, .25, .5, 1., 2., 4.)
HEIGHT_BOUNDS = (.05, 4.)


def save(path, data):
    path.write_text(json.dumps(data, indent=2)+'\n')


def axial_actions(points, cells):
    actions = []
    cell_sets = {tuple(sorted(c)) for c in cells}
    for g in range(4):
        reflected = points * [-1 if g & 1 else 1, -1 if g & 2 else 1, 1]
        distance = np.linalg.norm(reflected[:, None]-points[None], axis=2)
        action = distance.argmin(1)
        if distance.min(1).max() > 1e-10 or len(set(action)) != len(points):
            raise ValueError('This optimizer requires two axial mirrors')
        if {tuple(sorted(c)) for c in action[cells]} != cell_sets:
            raise ValueError('Geometry reflection does not preserve cells')
        actions.append(action)
    return np.asarray(actions).T


def height_columns(points, boundary, actions, ratio, bounds=HEIGHT_BOUNDS):
    """One shared boundary-height variable; symmetric, disjoint interior DOFs."""
    columns = [(ratio, *bounds, [(v, 2, float(z)) for v, z in enumerate(boundary[:, 2]) if z])]
    used = set()
    for v in range(len(boundary), len(points)):
        if v in used:
            continue
        used.update(actions[v])
        for axis in range(3):
            terms = {}
            for g, w in enumerate(actions[v]):
                sign = -1. if axis < 2 and g & (1 << axis) else 1.
                if w in terms and terms[w] != sign:
                    break
                terms[int(w)] = sign
            else:
                limit = bounds[1] if axis == 2 else 1.
                value = float(np.mean([points[w, axis]*sign for w, sign in terms.items()]))
                columns.append((value, -limit, limit, [(w, axis, sign) for w, sign in terms.items()]))
    return columns


def write_input(stem, points, cells, boundary, actions, ratio, bounds):
    columns = height_columns(points, boundary, actions, ratio, bounds)
    # Remove roundoff in symmetry-fixed coordinates and prescribe the boundary.
    points = np.mean([points[actions[:, g]]*[-1 if g & 1 else 1, -1 if g & 2 else 1, 1]
                      for g in range(4)], axis=0)
    points[:len(boundary)] = boundary*[1, 1, ratio]
    with stem.with_suffix('.txt').open('w') as stream:
        stream.write(f'{len(points)} {len(cells)} {len(boundary)} 1\n')
        np.savetxt(stream, points, fmt='%.17g')
        np.savetxt(stream, actions, fmt='%d')
        np.savetxt(stream, cells, fmt='%d')
    with stem.with_suffix('.constraints').open('w') as stream:
        stream.write(f'{len(columns)}\n')
        for value, lo, hi, terms in columns:
            stream.write(f'{value:.17g} {lo:.17g} {hi:.17g} {len(terms)} '+
                         ' '.join(f'{v} {axis} {coef:.17g}' for v,axis,coef in terms)+'\n')


def optimize(source, case, binary, output, heights=START_HEIGHTS, bounds=HEIGHT_BOUNDS, boundary_path=None):
    output.mkdir(parents=True, exist_ok=False)
    original, cells, quads = read_mesh(source)
    boundary_path = boundary_path or INPUT/f'{case["id"]}.mesh'
    boundary, _, _ = read_mesh(boundary_path)
    baseline_ratio = float(original[:len(boundary), 2].max())
    if not validate(source, case, baseline_ratio, boundary_path)['accepted']:
        raise ValueError('The starting mesh must be certified')
    actions = axial_actions(original, cells)
    attempts, accepted = [], []

    def consider(points, label):
        if points.shape != original.shape or not np.all(np.isfinite(points)):
            return None
        ratio = float(points[:len(boundary), 2].max())
        expected = boundary*[1, 1, ratio]
        if not bounds[0]-1e-12 <= ratio <= bounds[1]+1e-12:
            return None
        if np.max(abs(points[:len(boundary)]-expected)) > 1e-10:
            return None
        points[:len(boundary)] = expected
        path = output/f'{label}.mesh'
        write_mesh(path, points, cells, quads)
        report = validate(path, case, ratio, boundary_path)
        save(output/f'{label}.validation.json', report)
        if not report['accepted']:
            return None
        metrics = dense(points, cells, True)
        record = dict(label=label, mesh=path.name, height_ratio=ratio, metrics=metrics)
        accepted.append(record)
        return record

    baseline = consider(original.copy(), 'baseline')
    for index, ratio in enumerate(heights):
        points = original*[1, 1, ratio/baseline_ratio]
        consider(points.copy(), f'scaled-{index}')
        stem = output/f'hexopt-{index}'
        write_input(stem, points, cells, boundary, actions, ratio, bounds)
        with stem.with_suffix('.log').open('w') as log:
            try:
                run = subprocess.run([str(binary), str(stem.with_suffix('.txt')),
                    str(stem.with_suffix('.xyz')), str(stem.with_suffix('.constraints'))],
                    stdout=log, stderr=log, timeout=180,
                    env={**os.environ, 'HEXOPT_STALL_LIMIT': os.environ.get('HEXOPT_STALL_LIMIT', '50000')})
                record = consider(np.loadtxt(stem.with_suffix('.xyz')), f'optimized-{index}') if run.returncode == 0 else None
                attempts.append(dict(start_height_ratio=ratio, returncode=run.returncode,
                                     accepted=record is not None, result=record['label'] if record else None))
            except subprocess.TimeoutExpired:
                attempts.append(dict(start_height_ratio=ratio, timed_out=True))
    if not accepted:
        raise ValueError('No certified height embedding')
    by_sj = max(accepted, key=lambda r:(r['metrics']['min_scaled_jacobian'], -r['metrics']['max_condition']))
    by_condition = min(accepted, key=lambda r:(r['metrics']['max_condition'], -r['metrics']['min_scaled_jacobian']))
    for name, chosen in [('best-sj', by_sj), ('best-condition', by_condition)]:
        p, h, q = read_mesh(output/chosen['mesh'])
        write_mesh(output/f'{name}.mesh', p, h, q)
        write_vtu(output/f'{name}.vtu', p, h)
    result = dict(boundary_case=case['id'], height_bounds=list(bounds), start_heights=list(heights),
        objective='HexOpt scaled-Jacobian threshold; dense SJ and condition selections',
        max_iterations=int(os.environ.get('HEXOPT_MAX_ITERATIONS', '300000')),
        stall_limit=int(os.environ.get('HEXOPT_STALL_LIMIT', '50000')), attempt_timeout_seconds=180,
        source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
        adapter_sha256=hashlib.sha256((ROOT/'tools/hexopt_fixed.cpp').read_bytes()).hexdigest(),
        binary_sha256=hashlib.sha256(binary.read_bytes()).hexdigest(),
        baseline=baseline, best_sj=by_sj, best_condition=by_condition, attempts=attempts, candidates=accepted)
    save(output/'result.json', result)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--id', nargs='+', help='Mesh IDs; defaults to the whole refinement catalog')
    parser.add_argument('--threads', type=int, default=4, help='Parallel single-thread HexOpt jobs')
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--binary', type=Path, default=ROOT/'solver/build/hexopt_fixed')
    args = parser.parse_args()
    if args.threads < 1 or not args.binary.is_file():
        parser.error('Positive threads and a built HexOpt adapter are required')
    catalog = json.loads((ROOT/'docs/assets/refinement/catalog.json').read_text())
    templates = [t for t in catalog['templates'] if args.id is None or t['id'] in args.id]
    if args.id and set(args.id) != {t['id'] for t in templates}:
        parser.error('Unknown mesh ID')
    cases = {c['id']:c for c in json.loads((INPUT/'cases.json').read_text())}
    args.output = args.output.resolve(); args.output.mkdir(parents=True, exist_ok=False)
    args.binary = args.binary.resolve()
    def run(template):
        result = optimize(ROOT/'docs'/template['mesh'], cases[template['boundary_case']], args.binary,
                          args.output/template['id'])
        print(json.dumps(dict(id=template['id'], selections={k:dict(height=result[k]['height_ratio'],
              sj=result[k]['metrics']['min_scaled_jacobian'], condition=result[k]['metrics']['max_condition'])
              for k in ['best_sj','best_condition']})), flush=True)
        return dict(id=template['id'], **result)
    with ThreadPoolExecutor(max_workers=args.threads) as pool:
        results = list(pool.map(run, templates))
    save(args.output/'results.json', results)


if __name__ == '__main__':
    main()
