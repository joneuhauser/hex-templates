#!/usr/bin/env python3
"""Search the 6-to-8-quad upper transition with one or two axial mirrors."""
import argparse
from collections import Counter
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import time

import numpy as np

from hexopt_geode import initialize
from mesh_tools import audit, boundary, read_mesh, write_mesh, write_vtu
from mesh_homology import homology
from quality_metrics import dense
from validate_mesh import manifold

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'solver'))
from compare_meshes import isomorphism

FIXTURE = ROOT/'tests/fixtures/geode-top26.mesh'
BOUNDARY = ROOT/'solver/input/geode-top.mesh'
SITE = ROOT/'docs/assets/periodic/periodic/R52-improved.json'
COPIES = {29: (0, -1, 0), 36: (-1, 0, 0), 41: (-1, 0, 0),
          46: (-1, 0, 0), 47: (-1, 0, 0)}


def save(path, data):
    temporary = path.with_suffix('.tmp')
    temporary.write_text(json.dumps(data, indent=2)+'\n')
    temporary.replace(path)


def cut_upper(path):
    """Unwrap cells 26..51 in the documented periodic representatives."""
    data = json.loads(Path(path).read_text())
    cells = [c for b in data['blocks'] for c in b['hexes']]
    if len(cells) != 52:
        raise ValueError('Expected the 52-cell transition')
    points, lookup, hexes = [], {}, []
    for index, cell in enumerate(cells[26:], 26):
        row = []
        for point in np.asarray(cell['points'])+COPIES.get(index, (0, 0, 0)):
            key = tuple(np.round(point, 9))
            if key not in lookup:
                lookup[key] = len(points)
                points.append(point)
            row.append(lookup[key])
        hexes.append(row)
    p, h = np.asarray(points), np.asarray(hexes)
    q = boundary(h)
    order = sorted(set(q.flat))+sorted(set(range(len(p)))-set(q.flat))
    inverse = np.argsort(order)
    return p[order], inverse[h], inverse[q]


def actions(points, mirrors):
    result = []
    for g in range(1 << mirrors):
        image = points*np.array([-1 if g & 1 else 1, -1 if g & 2 else 1, 1])
        distances = np.linalg.norm(image[:, None]-points[None], axis=2)
        permutation = distances.argmin(1)
        if len(set(permutation)) != len(points) or distances.min(1).max() > 1e-12:
            raise ValueError('Boundary does not have the requested mirrors')
        result.append(permutation)
    return np.asarray(result).T


def boundary_report(points, quads, mirrors):
    """Check a spherical cut, cap counts, symmetry, and translated side pairs."""
    edges = Counter(tuple(sorted((a, b))) for f in quads for a, b in zip(f, np.roll(f, -1)))
    directed = Counter((int(a), int(b)) for f in quads for a, b in zip(f, np.roll(f, -1)))
    if set(edges.values()) != {2} or len(points)-len(edges)+len(quads) != 2:
        raise ValueError('Boundary is not a closed quadrangulated sphere')
    if any(directed[a, b] != 1 or directed[b, a] != 1 for a, b in edges):
        raise ValueError('Inconsistent boundary orientation')
    zlo, zhi = points[:, 2].min(), points[:, 2].max()
    caps = [np.flatnonzero(np.max(abs(points[quads, 2]-z), axis=1) < 1e-12)
            for z in (zlo, zhi)]
    if list(map(len, caps)) != [6, 8]:
        raise ValueError('Expected six Geode interface quads and eight top quads')
    group = actions(points, mirrors)
    faces = {tuple(sorted(f)) for f in quads}
    for g in range(group.shape[1]):
        if {tuple(sorted(f)) for f in group[:, g][quads]} != faces:
            raise ValueError('Mirror changes boundary connectivity')
    sides = sorted(set(range(len(quads)))-set(caps[0])-set(caps[1]))
    pairs = []
    for i in sides:
        matches = []
        for j in sides:
            if i == j:
                continue
            delta = points[quads[j]].mean(0)-points[quads[i]].mean(0)
            deck = np.rint([(delta[0]-delta[1])/2, (delta[0]+delta[1])/2]).astype(int)
            shift = np.array([deck.sum(), deck[1]-deck[0], 0.])
            if not deck.any() or np.max(abs(shift-delta)) > 1e-12:
                continue
            distance = np.linalg.norm((points[quads[i]]+shift)[:, None]-points[quads[j]][None], axis=2)
            if len(set(distance.argmin(1))) == 4 and distance.min(1).max() < 1e-12:
                corner_map = distance.argmin(1)
                if not np.all((np.roll(corner_map, -1)-corner_map) % 4 == 3):
                    raise ValueError('Translated side quads must have opposite orientations')
                matches.append((j, deck.tolist()))
        if len(matches) != 1:
            raise ValueError(f'Side face {i} has {len(matches)} translated partners')
        if i < matches[0][0]:
            pairs.append(dict(faces=[i, matches[0][0]], deck=matches[0][1]))
    return dict(vertices=len(points), quads=len(quads), cap_quads=[6, 8],
                mirrors=mirrors, periodic_side_pairs=pairs)


def symmetric_boundary(points, quads):
    """D2-invariant periodic cut, lofting a diamond to a 2-by-4 rectangle.

    At relative height t the domain is |Y|<=1 and
    |X|+(1-t)|Y| <= 1-t/2. Its area is two at every height.
    Four bilinear walls use 2x2 quads; two triangular walls use three quads.
    Opposite wall patches differ by integer combinations of (1,-1),(1,1).
    """
    lo, hi = points[:, 2].min(), points[:, 2].max()
    caps = [f for f in quads if min(np.max(abs(points[f, 2]-z)) for z in (lo, hi)) < 1e-12]
    polygons = [points[f].copy() for f in caps]
    A, B, C, D = np.array([[1, 0, lo], [0, 1, lo], [.5, 1, hi], [.5, 0, hi]])
    for sx in (-1, 1):
        for sy in (-1, 1):
            def at(u, v):
                return ((1-u)*(1-v)*A+u*(1-v)*B+u*v*C+(1-u)*v*D)*[sx, sy, 1]
            for i in range(2):
                for j in range(2):
                    polygons.append(np.array([at(i/2, j/2), at((i+1)/2, j/2),
                        at((i+1)/2, (j+1)/2), at(i/2, (j+1)/2)]))
    for sy in (-1, 1):
        triangle = np.array([[0, sy, lo], [.5, sy, hi], [-.5, sy, hi]])
        for i, v in enumerate(triangle):
            polygons.append(np.array([v, (v+triangle[(i+1)%3])/2,
                                      triangle.mean(0), (v+triangle[i-1])/2]))
    pp, lookup, qq = [], {}, []
    for poly in polygons:
        if np.cross(poly[1]-poly[0], poly[3]-poly[0])@poly.mean(0) < 0:
            poly = poly[::-1]
        row = []
        for v in poly:
            key = tuple(np.round(v, 10))
            if key not in lookup:
                lookup[key] = len(pp)
                pp.append(v)
            row.append(lookup[key])
        qq.append(row)
    p, q = np.asarray(pp), np.asarray(qq)
    boundary_report(p, q, 2)
    return p, q


def website_report():
    p, h, q = read_mesh(FIXTURE)
    bp, _, bq = read_mesh(BOUNDARY)
    if not np.array_equal(p[:len(bp)], bp) or not np.array_equal(q, bq):
        raise ValueError('Validation fixture and search boundary disagree')
    sp, sh, _ = cut_upper(SITE)
    target = dict(vertices=len(p), hexes=len(h), cells=h.tolist())
    witness = isomorphism(dict(vertices=len(sp), hexes=len(sh), cells=sh.tolist()), target, len(bp))
    if witness is None:
        raise ValueError('Replay fixture does not match the website upper-half connectivity')
    return dict(site=str(SITE.relative_to(ROOT)), site_sha256=hashlib.sha256(SITE.read_bytes()).hexdigest(),
                upper_cells=list(range(26, 52)), vertex_bijection=witness), target


def embed(candidate, bp, quads, binary, output, restarts=6):
    """Harmonic initialization and HexOpt with fixed boundary and axial orbits."""
    output.mkdir()
    p, h = initialize(bp, candidate)
    group = np.asarray(candidate['vertex_actions'])
    columns, used = [], set()
    for v in range(len(bp), len(p)):
        if v in used:
            continue
        used.update(group[v])
        for axis in range(3):
            terms = {}
            for g, w in enumerate(group[v]):
                sign = -1. if axis < 2 and g & (1 << axis) else 1.
                if w in terms and terms[w] != sign:
                    break  # Stabilizer fixes this coordinate to zero.
                terms[int(w)] = sign
            else:
                columns.append((p[v, axis], float(bp[:, axis].min()), float(bp[:, axis].max()),
                                [(w, axis, sign) for w, sign in terms.items()]))
    write_mesh(output/'initial.mesh', p, h, quads)
    save(output/'candidate.json', candidate)
    attempts = []
    rng = np.random.default_rng(0)
    for attempt in range(restarts):
        start = p.copy()
        parameters = [np.clip(value+(rng.normal(0, .025*attempt) if attempt else 0), lo, hi)
                      for value, lo, hi, terms in columns]
        for value, (_, _, _, terms) in zip(parameters, columns):
            for v, axis, sign in terms:
                start[v, axis] = sign*value
        stem = output/f'attempt-{attempt}'
        with stem.with_suffix('.txt').open('w') as stream:
            stream.write(f'{len(p)} {len(h)} {len(bp)} 1\n')
            np.savetxt(stream, start, fmt='%.17g')
            # The affine basis encodes symmetry; the four-action adapter table is unused.
            np.savetxt(stream, np.tile(np.arange(len(p))[:, None], (1, 4)), fmt='%d')
            np.savetxt(stream, h, fmt='%d')
        with stem.with_suffix('.constraints').open('w') as stream:
            stream.write(f'{len(columns)}\n')
            for value, (_, lo, hi, terms) in zip(parameters, columns):
                stream.write(f'{value:.17g} {lo:.17g} {hi:.17g} {len(terms)} '+
                    ' '.join(f'{v} {axis} {sign:.17g}' for v, axis, sign in terms)+'\n')
        if columns:
            with stem.with_suffix('.log').open('w') as log:
                try:
                    result = subprocess.run([str(binary), str(stem.with_suffix('.txt')),
                        str(stem.with_suffix('.xyz')), str(stem.with_suffix('.constraints'))],
                        stdout=log, stderr=log, timeout=120,
                        env=dict(os.environ, HEXOPT_STALL_LIMIT=os.environ.get('HEXOPT_STALL_LIMIT', '50000')))
                except subprocess.TimeoutExpired:
                    attempts.append(dict(attempt=attempt, timed_out=True))
                    continue
            if result.returncode:
                attempts.append(dict(attempt=attempt, exit_code=result.returncode))
                continue
            optimized = np.loadtxt(stem.with_suffix('.xyz'))
        else:
            optimized = start
        write_mesh(stem.with_suffix('.mesh'), optimized, h, quads)
        report = audit(stem.with_suffix('.mesh'))
        report['manifold'] = manifold(h)
        report['homology'] = homology(h)
        report['fixed_boundary'] = bool(np.array_equal(optimized[:len(bp)], bp))
        report['mirror_error'] = max(float(np.max(abs(optimized[group[:, g]]-
            optimized*[-1 if g & 1 else 1, -1 if g & 2 else 1, 1]))) for g in range(group.shape[1]))
        report['accepted'] = bool(report['all_jacobians_positive'] and report['boundary_matches']
            and report['opposite_internal_orientations'] and not report['incompatible_pairs']
            and report['distinct_cell_vertices'] and report['euler'] == 1
            and report['manifold']['edge_links_valid'] and report['manifold']['vertex_links_valid']
            and report['homology']['betti_GF2'] == [1, 0, 0, 0]
            and report['fixed_boundary'] and report['mirror_error'] < 1e-12)
        save(stem.with_suffix('.validation.json'), report)
        attempts.append(dict(attempt=attempt, accepted=report['accepted']))
        if report['accepted']:
            write_mesh(output/'accepted.mesh', optimized, h, quads)
            write_vtu(output/'accepted.vtu', optimized, h)
            save(output/'validation.json', report)
            quality = dense(optimized, h, True)
            save(output/'quality.json', quality)
            record = dict(hexes=len(h), accepted=True, attempts=attempts,
                          min_scaled_jacobian=quality['min_scaled_jacobian'])
            save(output/'result.json', record)
            return record
    record = dict(hexes=len(h), accepted=False, attempts=attempts)
    save(output/'result.json', record)
    return record


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--symmetry', choices=['first', 'other', 'two'], default='first')
    parser.add_argument('--cap', type=int)
    parser.add_argument('--threads', type=int, required=True)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--seconds', type=int, default=0)
    parser.add_argument('--restarts', type=int, default=int(os.environ.get('HEXOPT_RESTARTS', '6')))
    parser.add_argument('--solver', type=Path, default=ROOT/'solver/build/mirror_search')
    parser.add_argument('--hexopt', type=Path, default=ROOT/'solver/build/hexopt_fixed')
    args = parser.parse_args()
    if args.cap is None:
        args.cap = dict(first=26, other=30, two=40)[args.symmetry]
    if not 1 <= args.threads <= 256 or not 1 <= args.cap <= 48 or args.seconds < 0:
        parser.error('Use THREADS=1..256, CAP=1..48 and nonnegative seconds')
    if args.restarts < 1:
        parser.error('HexOpt restarts must be positive')
    for name in ('solver', 'hexopt'):
        path = getattr(args, name).resolve()
        if not path.is_file() or not os.access(path, os.X_OK):
            parser.error(f'Missing executable: {path}')
        setattr(args, name, path)
    bp, _, quads = read_mesh(BOUNDARY)
    if args.symmetry != 'first':
        bp, quads = symmetric_boundary(bp, quads)
    if args.symmetry == 'other':
        bp = bp[:, [1, 0, 2]]*[-1, 1, 1]  # Proper rotation maps the other mirror to X=0.
    mirrors = 2 if args.symmetry == 'two' else 1
    info = boundary_report(bp, quads, mirrors)
    if len(bp)+2*args.cap+1-len(quads)//2 >= 112:
        parser.error('Cap exceeds the solver vertex capacity for this boundary')
    if args.output:
        output = args.output.resolve()
        output.mkdir(parents=True, exist_ok=False)
    else:
        parent = ROOT/'solver/runs'
        parent.mkdir(exist_ok=True)
        output = Path(tempfile.mkdtemp(prefix=f'geode-top-{args.symmetry}-H{args.cap}-', dir=parent))
    print(f'Output: {output}', flush=True)
    search = output/'search'
    search.mkdir()
    write_mesh(search/'boundary.mesh', bp, np.empty((0, 8), int), quads)
    status = dict(symmetry=args.symmetry, cap=args.cap, threads=args.threads, boundary=info, stage='search')
    common = [str(args.solver), '--input', str(search/'boundary.mesh'), '--cap', str(args.cap),
              '--threads', str(args.threads), '--mirrors', str(mirrors), '--planes', 'axial',
              '--seconds', str(args.seconds), '--boundary-canonical', '0']
    target = None
    if args.symmetry == 'first':
        proof, target = website_report()
        save(output/'website-connectivity.json', proof)
        replay = common+['--replay', str(FIXTURE), '--output', str(output/'replay.jsonl')]
        with (output/'replay.log').open('w') as log:
            subprocess.run(replay, stdout=log, stderr=log, check=True)
        status['fixture_replay'] = True
    command = common+['--output', str(search/'candidates.jsonl')]
    status['command'] = command
    save(output/'status.json', status)
    started = time.monotonic()
    with (search/'search.log').open('w') as log:
        run = subprocess.run(command, stdout=log, stderr=log)
    status.update(search_seconds=time.monotonic()-started, search_exit_code=run.returncode,
                  search_complete=run.returncode == 0 and 'SEARCH_FINISHED ' in (search/'search.log').read_text())
    if run.returncode not in (0, 124):
        status['stage'] = 'failed'
        save(output/'status.json', status)
        raise RuntimeError(f'Search failed; see {search}')
    status['stage'] = 'hexopt'
    save(output/'status.json', status)
    (output/'hexopt').mkdir()
    results = []
    for index, line in enumerate((search/'candidates.jsonl').read_text().splitlines(), 1):
        candidate = json.loads(line)
        witness = isomorphism(candidate, target, len(bp)) if target else None
        result = embed(candidate, bp, quads, args.hexopt, output/'hexopt'/f'candidate-{index:03d}', args.restarts)
        result.update(candidate=index, website_match=witness is not None)
        if witness is not None:
            result['vertex_bijection'] = witness
        results.append(result)
        save(output/'results.json', results)
        print(json.dumps({k: v for k, v in result.items() if k not in ['attempts', 'vertex_bijection']}), flush=True)
    save(output/'results.json', results)
    status.update(stage='complete' if status['search_complete'] else 'search_incomplete',
                  candidates=len(results), accepted=sum(r['accepted'] for r in results),
                  website_matches=sum(r['website_match'] for r in results))
    save(output/'status.json', status)
    print(json.dumps(status), flush=True)
    if not status['search_complete']:
        raise SystemExit(124)


if __name__ == '__main__':
    main()
