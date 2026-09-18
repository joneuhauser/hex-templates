"""Generate refinement boundaries, search their interiors and certify HexOpt meshes."""
import argparse
from collections import Counter
from fractions import Fraction as F
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import time

import numpy as np
import build_geode_sides as disks
from build_catalog import symmetries
from geode_top import embed, save
from mesh_tools import read_mesh, write_mesh, write_vtu, audit, bernstein
from mesh_homology import homology
from validate_mesh import manifold
from quality_metrics import dense

ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT/'solver/input/refinement'
SELECTED = {'9-to-1-r0-Q4-01': 13, '9-to-1-r1-Q5-01': 14,
            '16-to-4-r0-Q6-01': 28, '16-to-4-r1-Q8-03': 32,
            '16-to-4-r0-Q7-04': 34, '25-to-9-r0-Q6-01': 33}
sys.path.insert(0, str(ROOT/'solver'))
from compare_meshes import isomorphism


def rim(bottom, top, r):
    return ([(F(-1)+F(2*i, bottom), F(-1)) for i in range(bottom)]
        + [(F(1), F(-1)+F(2*i, r+1)) for i in range(r+1)]
        + [(F(1)-F(2*i, top), F(1)) for i in range(top)]
        + [(F(-1), F(1)-F(2*i, r+1)) for i in range(r+1)])


def patch_geometry(record, bottom, top, r):
    p = disks.harmonic(record['quads'], r, boundary=rim(bottom, top, r))
    scores, area = [], F(0)
    for q in record['quads']:
        for j, a in enumerate(q):
            b, c = q[(j+1)%4], q[(j-1)%4]
            u = [p[b][k]-p[a][k] for k in range(2)]
            v = [p[c][k]-p[a][k] for k in range(2)]
            det = disks.cross(u, v)
            if det <= 0:
                return None
            scores.append(2*det/sum(x*x for x in u+v))
            area += disks.cross(p[a], p[b])/2
    assert area == 4
    assert len(set(p)) == len(p)
    assert all(p[w] == (-p[v][0], p[v][1]) for v, w in enumerate(record['mirror']))
    return p, min(scores)


def boundary(bottom, top, patch, side_points):
    polygons = []
    for n, z in [(bottom, F(-1)), (top, F(1))]:
        for i in range(n):
            for j in range(n):
                x0, x1, y0, y1 = [F(-1)+F(2*k, n) for k in (i, i+1, j, j+1)]
                polygons.append([(x0,y0,z),(x1,y0,z),(x1,y1,z),(x0,y1,z)])
    for axis in (0, 1):
        for sign in (-1, 1):
            for q in patch['quads']:
                polygons.append([(F(sign),side_points[v][0],side_points[v][1]) if axis == 0
                    else (side_points[v][0],F(sign),side_points[v][1]) for v in q])
    p, lookup, faces = [], {}, []
    for poly in polygons:
        a = np.asarray(poly, float)
        if np.dot(np.cross(a[1]-a[0], a[3]-a[0]), a.mean(0)) < 0:
            poly = poly[::-1]
        face = []
        for v in poly:
            if v not in lookup:
                lookup[v] = len(p); p.append(v)
            face.append(lookup[v])
        faces.append(face)
    p, q = np.asarray(p, float), np.asarray(faces, int)
    edges = Counter((int(a),int(b)) for f in q for a,b in zip(f,np.roll(f,-1)))
    assert all(n == 1 and edges[b,a] == 1 for (a,b),n in edges.items())
    assert len(p)-len(edges)//2+len(q) == 2
    assert len(p) == len(q)+2
    for swap in (False, True):
        for sx in (-1,1):
            for sy in (-1,1):
                image = p[:,[1,0,2]] if swap else p.copy()
                image = image*[sx,sy,1]
                d = np.linalg.norm(image[:,None]-p[None],axis=2)
                permutation = d.argmin(1)
                assert d.min(1).max() < 1e-12 and len(set(permutation)) == len(p)
                assert {tuple(sorted(f)) for f in permutation[q]} == {tuple(sorted(f)) for f in q}
    for axis in (0,1):
        sides=[]
        for sign in (-1,1):
            sides.append({tuple(sorted(tuple(v[[1-axis,2]]) for v in p[f]))
                          for f in q if np.all(p[f,axis] == sign)})
        assert sides[0] == sides[1]
    return p,q


def enumerate_sides(output, grids=((1, 3), (2, 4), (3, 5)), max_quads=None):
    output.mkdir(parents=True, exist_ok=False)
    records, runs, scopes = [], [], []
    with tempfile.TemporaryDirectory() as tmp:
        for bottom, top in grids:
            cap = max_quads if max_quads is not None else (10 if (bottom, top) == (3, 5) else 8)
            scopes.append(dict(bottom_segments=bottom, top_segments=top, max_quads=cap))
            for r in (0, 1):
                rows = []
                for orbits in (True, False):
                    path = Path(tmp)/'sides.jsonl'
                    command = [str(ROOT/'solver/build/quad_disk_search'), str(cap), str(r), str(path),
                               '--bottom-segments', str(bottom), '--top-segments', str(top)]
                    if not orbits:
                        command.append('--no-orbits')
                    run = subprocess.run(command, capture_output=True, text=True, check=True)
                    rows.append([json.loads(s) for s in path.read_text().splitlines()])
                    runs.append(dict(bottom_segments=bottom, top_segments=top,
                        vertical_interior_vertices=r, orbits=orbits, log=run.stderr.strip(),
                        results_sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
                keys = lambda records: {disks.canonical(a['quads'], a['boundary_vertices']) for a in records}
                if keys(rows[0]) != keys(rows[1]):
                    raise ValueError('Side enumeration cross-check disagrees')
                valid = []
                for record in rows[0]:
                    geometry = patch_geometry(record, bottom, top, r)
                    if geometry is None:
                        raise ValueError('Side pattern needs a different embedding')
                    points, score = geometry
                    valid.append((record, points, score))
                valid.sort(key=lambda x: (len(x[0]['quads']), -x[2], x[0]['quads']))
                for i, (record, points, score) in enumerate(valid, 1):
                    name = f'{top*top}-to-{bottom*bottom}-r{r}-Q{len(record["quads"])}-{i:02d}'
                    bp, bq = boundary(bottom, top, record, points)
                    write_mesh(output/f'{name}.mesh', bp, np.empty((0, 8), int), bq)
                    records.append(dict(id=name, bottom_grid=bottom, top_grid=top,
                        side_quads=len(record['quads']), vertical_interior_vertices=r,
                        side_quality=float(score), boundary_vertices=len(bp), boundary_quads=len(bq),
                        side=record, side_points_exact=[[str(v) for v in p] for p in points]))
    save(output/'cases.json', records)
    save(output/'enumeration.json', dict(max_quads=max(s['max_quads'] for s in scopes), scopes=scopes, vertical_interior_vertices=[0, 1],
        source_sha256=hashlib.sha256((ROOT/'solver/quad_disk_search.cpp').read_bytes()).hexdigest(), runs=runs))
    print(f'Enumerated and cross-checked {len(records)} refinement side walls', flush=True)


def validate(path, case):
    p, h, q = read_mesh(path)
    bp, _, bq = read_mesh(INPUT/f'{case["id"]}.mesh')
    report = audit(path)
    report['manifold'] = manifold(h)
    report['homology'] = homology(h)
    report['fixed_boundary'] = bool(len(p) >= len(bp) and np.array_equal(p[:len(bp)], bp)
        and {tuple(sorted(f)) for f in q} == {tuple(sorted(f)) for f in bq})
    report['symmetry'] = symmetries(p, h)
    report['volume'] = sum(float(sum(bernstein(p[c]).flat))*8/27 for c in h)
    report['contained'] = bool(np.all(p >= -1-1e-12) and np.all(p <= 1+1e-12))
    report['cap_quads'] = [int(sum(np.all(p[f, 2] == z) for f in q)) for z in (-1, 1)]
    report['translation_matching'] = True
    for axis in (0, 1):
        sides = [{tuple(sorted(tuple(v[[1-axis, 2]]) for v in p[f])) for f in q
                  if np.all(p[f, axis] == sign)} for sign in (-1, 1)]
        report['translation_matching'] &= sides[0] == sides[1] and len(sides[0]) == case['side_quads']
    report['accepted'] = bool(report['all_jacobians_positive'] and report['fixed_boundary']
        and report['boundary_matches'] and report['opposite_internal_orientations']
        and set(report['face_multiplicity']) <= {1, 2} and not report['incompatible_pairs']
        and report['distinct_cell_vertices'] and report['euler'] == 1
        and report['manifold']['edge_links_valid'] and report['manifold']['vertex_links_valid']
        and report['homology']['betti_GF2'] == [1, 0, 0, 0] and report['contained']
        and abs(report['volume']-8) < 1e-9 and report['translation_matching']
        and report['cap_quads'] == [case['bottom_grid']**2, case['top_grid']**2]
        and len(report['symmetry']['mirrors']) >= 2)
    return report


def representatives(candidates, bp):
    permutations = []
    for turns in range(4):
        image = bp.copy()
        for _ in range(turns):
            image = image[:, [1, 0, 2]]*[-1, 1, 1]
        permutations.append(np.linalg.norm(image[:, None]-bp[None], axis=2).argmin(1).tolist())
    unique = []
    for candidate in sorted(candidates, key=lambda c: (c['hexes'], c['cells'])):
        if not any(isomorphism(candidate, other, len(bp), perm) is not None
                   for other in unique for perm in permutations):
            unique.append(candidate)
    return unique


def search_case(case, args, parent):
    cap = args.cap or SELECTED.get(case['id'], 28)
    directory = parent/case['id']; directory.mkdir()
    search = directory/'search'; search.mkdir()
    bp, _, q = read_mesh(INPUT/f'{case["id"]}.mesh')
    write_mesh(search/'boundary.mesh', bp, np.empty((0, 8), int), q)
    command = [str(ROOT/'solver/build/mirror_search'), '--input', str(search/'boundary.mesh'),
        '--cap', str(cap), '--threads', str(args.threads), '--mirrors', '2', '--planes', args.planes,
        '--seconds', str(args.seconds), '--boundary-canonical', '0', '--output', str(search/'candidates.jsonl')]
    status = dict(case=case['id'], cap=cap, planes=args.planes, command=command, stage='search')
    save(directory/'status.json', status)
    start = time.monotonic()
    with (search/'search.log').open('w') as log:
        result = subprocess.run(command, stdout=log, stderr=log)
    status.update(search_returncode=result.returncode, search_complete=result.returncode == 0,
                  search_seconds=time.monotonic()-start)
    if result.returncode not in (0, 124):
        save(directory/'status.json', status)
        raise RuntimeError(f'Search failed: {search}/search.log')
    candidates = [json.loads(s) for s in (search/'candidates.jsonl').read_text().splitlines()]
    unique = representatives(candidates, bp)
    status.update(stage='hexopt', candidates=len(candidates), classes=len(unique))
    save(directory/'status.json', status)
    matrix = np.eye(3)
    if args.planes == 'diagonal':
        matrix[:2, :2] = np.array([[1, -1], [1, 1]])/np.sqrt(2)
    results = []; (directory/'hexopt').mkdir()
    for i, candidate in enumerate(unique):
        output = directory/'hexopt'/f'variant-{i:03d}'
        record = embed(candidate, bp@matrix.T, q, args.hexopt, output, args.restarts)
        if record['accepted']:
            p, h, _ = read_mesh(output/'accepted.mesh')
            p = p@matrix; p[:len(bp)] = bp
            write_mesh(output/'accepted.mesh', p, h, q)
            report = validate(output/'accepted.mesh', case)
            if not report['accepted']:
                raise ValueError('Embedding did not pass physical-coordinate validation')
            write_vtu(output/'accepted.vtu', p, h)
            save(output/'validation.json', report)
            save(output/'quality.json', dense(p, h, True))
        results.append(record)
        print(json.dumps(dict(case=case['id'], variant=i, **record)), flush=True)
    save(directory/'results.json', results)
    status.update(stage='finished', accepted=sum(r['accepted'] for r in results))
    save(directory/'status.json', status)


def main():
    global INPUT
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='mode', required=True)
    sides = sub.add_parser('sides'); sides.add_argument('--output', type=Path, required=True)
    sides.add_argument('--bottom-grid', type=int)
    sides.add_argument('--top-grid', type=int)
    sides.add_argument('--max-quads', type=int)
    search = sub.add_parser('search')
    search.add_argument('--case', default='all', help='Boundary ID from cases.json; all runs the displayed cases')
    search.add_argument('--threads', type=int, required=True)
    search.add_argument('--cap', type=int)
    search.add_argument('--planes', choices=['axial', 'diagonal'], default='axial')
    search.add_argument('--seconds', type=int, default=0)
    search.add_argument('--restarts', type=int, default=int(os.environ.get('HEXOPT_RESTARTS', '6')))
    search.add_argument('--hexopt', type=Path, default=Path(os.environ.get('HEXOPT_BINARY', ROOT/'solver/build/hexopt_fixed')))
    search.add_argument('--output', type=Path)
    search.add_argument('--inputs', type=Path, default=INPUT,
                        help='Folder containing cases.json and boundary meshes')
    args = parser.parse_args()
    if args.mode == 'sides':
        if (args.bottom_grid is None) != (args.top_grid is None):
            parser.error('Specify both --bottom-grid and --top-grid')
        grids = ((args.bottom_grid, args.top_grid),) if args.bottom_grid is not None else ((1, 3), (2, 4), (3, 5))
        if (args.max_quads is not None and not 1 <= args.max_quads <= 12) or any(not 1 <= n <= 12 for pair in grids for n in pair):
            parser.error('Grid sizes and quad cap must be between 1 and 12')
        if any((bottom+top)%2 or (args.max_quads if args.max_quads is not None else (10 if (bottom,top)==(3,5) else 8)) < (bottom+top)//2+1 for bottom,top in grids):
            parser.error('Top and bottom need equal parity and a sufficient quad cap')
        enumerate_sides(args.output.resolve(), grids, args.max_quads); return
    if not 1 <= args.threads <= 256 or args.seconds < 0 or args.restarts < 1 or (args.cap is not None and not 1 <= args.cap <= 48):
        parser.error('Invalid thread count, cap, time limit or restart count')
    INPUT = args.inputs.resolve()
    cases = json.loads((INPUT/'cases.json').read_text())
    selected = [c for c in cases if (c['id'] in SELECTED if args.case == 'all' else c['id'] == args.case)]
    if not selected:
        parser.error('Unknown boundary case')
    args.hexopt = args.hexopt.resolve()
    if not args.hexopt.is_file():
        parser.error('Build the HexOpt adapter first; see tools/HEXOPT_NOTES.md')
    if args.output:
        output = args.output.resolve(); output.mkdir(parents=True, exist_ok=False)
    else:
        (ROOT/'solver/runs').mkdir(exist_ok=True)
        output = Path(tempfile.mkdtemp(prefix='refinement-', dir=ROOT/'solver/runs'))
    print(f'Output: {output}', flush=True)
    for case in selected:
        search_case(case, args, output)


if __name__ == '__main__':
    main()
