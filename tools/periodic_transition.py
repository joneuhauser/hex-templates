"""Periodic HEX8 geometry and validation, without artificial side boundaries.

Vertices live in a quotient; each cell corner also carries a lattice offset.
All incidence tests retain those offsets, including edges that wrap a period.
"""
import argparse
from collections import Counter, defaultdict
from itertools import combinations, product
import json
import hashlib
from pathlib import Path

import numpy as np
from scipy.optimize import minimize
from scipy.sparse import csr_matrix

from mesh_tools import FACES, EDGES, derivatives, certify_hex, bernstein, topology
from geode_geometry import objective


def load(path):
    data = json.loads(Path(path).read_text())
    period = np.array([*data['period'], 0.])
    q = np.array([c['points'] for b in data['blocks'] for c in b['hexes']])
    points, lookup, cells, offsets = [], {}, [], []
    for cell in q:
        hh, oo = [], []
        for point in cell:
            off = np.r_[np.floor(point[:2]/period[:2]+1e-8).astype(int), 0]
            p = point-off*period
            key = tuple(np.round(p, 9))
            if key not in lookup:
                lookup[key] = len(points)
                points.append(p)
            hh.append(lookup[key]); oo.append(off)
        cells.append(hh); offsets.append(oo)
    p, h, o = np.array(points), np.array(cells), np.array(offsets)
    if np.max(abs(p[h]+o*period-q)) > 2e-9:
        raise ValueError('Periodic coordinate welding failed')
    return p, h, o, period


def completion(p, h, o, period, repeats=2):
    vertices, lookup, cells = [], {}, []
    q = p[h]+o*period
    decks = np.c_[-np.floor(q.mean(axis=1)[:, :2]/period[:2]).astype(int), np.zeros(len(h), int)]
    for x, y in product(range(repeats), repeat=2):
        for ci, c in enumerate(h):
            ids = []
            for v, offset in zip(c, o[ci]+decks[ci]+[x, y, 0]):
                key = (int(v), *map(int, offset))
                if key not in lookup:
                    lookup[key] = len(vertices); vertices.append(p[v]+offset*period)
                ids.append(lookup[key])
            cells.append(ids)
    return np.array(vertices), np.array(cells)


def labels(h, o):
    return np.concatenate([h[..., None], o[..., :2]], axis=-1)


def normalize(rows):
    rows = [tuple(map(int, r)) for r in rows]
    anchor = min(rows)
    return tuple((v, x-anchor[1], y-anchor[2]) for v, x, y in rows)


def key(rows):
    return tuple(sorted(normalize(rows)))


def cycle(rows):
    rows = normalize(rows)
    return min(rows[i:]+rows[:i] for i in range(len(rows)))


def incidence(h, o):
    result = defaultdict(list)
    ll = labels(h, o)
    for ci, cell in enumerate(ll):
        for fi, face in enumerate(FACES):
            result[key(cell[face])].append((ci, fi, cell[face]))
    return result


def cap_keys(p, h, o, period):
    result = [[], []]
    q = p[h]+o*period
    bounds = [q[..., 2].min(), q[..., 2].max()]
    for entries in incidence(h, o).values():
        if len(entries) != 1:
            continue
        ci, fi, _ = entries[0]
        face = q[ci, FACES[fi]].copy()
        for side, z in enumerate(bounds):
            if np.max(abs(face[:, 2]-z)) < 1e-9:
                face[:, :2] -= np.floor(face[:, :2].mean(axis=0)/period[:2]+1e-8)*period[:2]
                result[side].append(tuple(sorted(map(tuple, np.round(face, 9)))))
    return [sorted(x) for x in result]


def rank_gf2(columns):
    pivots = {}
    for column in columns:
        while column:
            pivot = column.bit_length()-1
            if pivot not in pivots:
                pivots[pivot] = column
                break
            column ^= pivots[pivot]
    return len(pivots)


def topology_report(p, h, o):
    ll = labels(h, o)
    ff = incidence(h, o)
    ee = {key(c[e]): c[e] for c in ll for e in EDGES}
    edge_index = {e: i for i, e in enumerate(ee)}
    face_index = {f: i for i, f in enumerate(ff)}
    def chain(items):
        value = 0
        for i in items:
            value ^= 1 << int(i)
        return value
    d1 = [chain(e[:, 0]) for e in ee.values()]
    d2 = [chain(edge_index[key(face[[i, (i+1)%4]])] for i in range(4))
          for entries in ff.values() for face in [entries[0][2]]]
    d3 = [chain(face_index[key(c[f])] for f in FACES) for c in ll]
    used = set(h.flat)
    ranks = [rank_gf2(d) for d in (d1, d2, d3)]
    betti = [len(used)-ranks[0], len(ee)-ranks[0]-ranks[1],
             len(ff)-ranks[1]-ranks[2], len(h)-ranks[2]]
    oriented = all(len(es) != 2 or cycle(es[0][2]) == cycle(es[1][2][::-1]) for es in ff.values())
    bad_vertices, bad_edges, incompatible = [], [], []
    cap_vertices = {int(v) for es in ff.values() if len(es) == 1 for v in es[0][2][:, 0]}
    # Build each complete lifted vertex star, even when a cell contains two
    # different periodic copies of the same quotient vertex.
    for vertex in sorted(used):
        star = []
        for ci, corner in zip(*np.where(h == vertex)):
            c = ll[ci].copy(); c[:, 1:] -= c[corner, 1:]
            star.append(c)
        points = sorted({tuple(v) for c in star for v in c})
        lookup = {v: i for i, v in enumerate(points)}
        cells = np.array([[lookup[tuple(v)] for v in c] for c in star])
        target = lookup[(int(vertex), 0, 0)]
        check = topology(cells, np.empty((0, 4), int))
        if check['incompatible_pairs'] or not check['distinct_cell_vertices']:
            incompatible.append(int(vertex))
        link_triangles = []
        for c in cells:
            link_triangles.append(tuple(int(e[0] if e[1] == target else e[1])
                for e in c[EDGES] if target in e))
        counts = Counter(tuple(sorted(e)) for t in link_triangles for e in combinations(t, 2))
        vs = set(v for t in link_triangles for v in t)
        adj = {v: set() for v in vs}
        for a, b in counts:
            adj[a].add(b); adj[b].add(a)
        seen, todo = set(), [next(iter(vs))]
        while todo:
            v = todo.pop()
            if v not in seen:
                seen.add(v); todo.extend(adj[v]-seen)
        boundary_edges = [e for e, n in counts.items() if n == 1]
        bd = Counter(v for e in boundary_edges for v in e)
        expected = 1 if vertex in cap_vertices else 2
        if (len(vs)-len(counts)+len(link_triangles) != expected or seen != vs
                or not set(counts.values()) <= {1, 2}
                or (expected == 1 and (not bd or set(bd.values()) != {2}))
                or (expected == 2 and bd)):
            bad_vertices.append(int(vertex))
        # Link of every incident edge must be one path or one cycle.
        for v in vs:
            pairs = [tuple(w for w in t if w != v) for t in link_triangles if v in t]
            degree = Counter(w for e in pairs for w in e)
            ea = {w: set() for w in degree}
            for a, b in pairs:
                ea[a].add(b); ea[b].add(a)
            seen, todo = set(), [next(iter(degree))]
            while todo:
                w = todo.pop()
                if w not in seen:
                    seen.add(w); todo.extend(ea[w]-seen)
            if seen != set(degree) or not set(degree.values()) <= {1, 2} or sum(n == 1 for n in degree.values()) not in (0, 2):
                bad_edges.append([int(vertex), int(v)])
    return dict(betti_GF2=betti, euler=len(used)-len(ee)+len(ff)-len(h),
                face_multiplicity=dict(Counter(map(len, ff.values()))),
                opposite_face_orientations=oriented, bad_vertex_links=bad_vertices,
                bad_edge_links=bad_edges, incompatible_stars=incompatible)


def quality(p, h, o, period, n=11):
    ds = derivatives(np.array(list(product(np.linspace(-1, 1, n), repeat=3))))
    jac = np.einsum('hvi,svj->hsij', p[h]+o*period, ds)
    det = np.linalg.det(jac)
    with np.errstate(divide='ignore', invalid='ignore'):
        sj = det/np.linalg.norm(jac, axis=-2).prod(-1)
    sj = np.nan_to_num(sj, nan=-1.)
    mr = 3*np.sign(det)*abs(det)**(2/3)/(jac*jac).sum(axis=(-1, -2))
    return dict(grid=n, min_scaled_jacobian=float(sj.min()), min_mean_ratio=float(mr.min()),
                max_condition=float(np.linalg.cond(jac).max()), min_det=float(det.min()))


def validate(p, h, o, period, reference=None, exact=True):
    report = topology_report(p, h, o)
    q = p[h]+o*period
    caps = cap_keys(p, h, o, period)
    bounds = [q[..., 2].min(), q[..., 2].max()]
    unmatched = []
    for es in incidence(h, o).values():
        if len(es) == 1:
            ci, fi, _ = es[0]
            if not any(np.max(abs(q[ci, FACES[fi], 2]-z)) < 1e-9 for z in bounds):
                unmatched.append([ci, fi])
    volume = sum(float(sum(bernstein(c).flat))*8/27 for c in q)
    expected_volume = np.prod(period[:2])*(bounds[1]-bounds[0])
    report.update(hexes=len(h), periodic_vertices=len(set(h.flat)), cap_faces=list(map(len, caps)),
                  unmatched_noncap_faces=unmatched, sampled=quality(p, h, o, period, 21),
                  caps_preserved=reference is None or caps == cap_keys(*reference),
                  volume_error=float(volume-expected_volume))
    if exact:
        report['jacobian_certificates'] = [certify_hex(c) for c in q]
        report['exact_positive'] = all(c['positive'] for c in report['jacobian_certificates'])
    report['accepted'] = bool(exact and report['exact_positive'] and report['betti_GF2'] == [1, 2, 1, 0]
        and report['opposite_face_orientations'] and not report['bad_vertex_links']
        and not report['bad_edge_links'] and not report['incompatible_stars']
        and set(report['face_multiplicity']) <= {1, 2} and not unmatched
        and report['caps_preserved'] and abs(report['volume_error']) < 1e-8)
    return report


def quality_penalty(jac, condition_limit=None, mean_ratio_floor=None):
    value = 0.
    dj = np.zeros_like(jac)
    if condition_limit is not None:
        u, singular, vh = np.linalg.svd(jac)
        singular = np.maximum(singular, 1e-20)
        excess = np.maximum(np.log(singular[..., 0]/singular[..., -1]/condition_limit), 0)
        # Log condition has a simple spectral derivative away from repeated
        # singular values; repeated values have zero penalty here.
        dj += (u[..., :, 0, None]*vh[..., None, 0, :]/singular[..., 0, None, None]
              -u[..., :, -1, None]*vh[..., None, -1, :]/singular[..., -1, None, None])
        dj *= (200*excess)[..., None, None]
        value += 100*np.sum(excess*excess)
    if mean_ratio_floor is not None:
        cof = np.stack([np.cross(jac[..., 1], jac[..., 2]),
                        np.cross(jac[..., 2], jac[..., 0]),
                        np.cross(jac[..., 0], jac[..., 1])], axis=-1)
        det = (jac[..., 0]*cof[..., 0]).sum(-1)
        norm = np.maximum((jac*jac).sum(axis=(-1, -2)), 1e-30)
        mr = 3*np.sign(det)*abs(det)**(2/3)/norm
        dmr = (2*cof/np.maximum(abs(det), 1e-20)[..., None, None]**(1/3)
               -2*mr[..., None, None]*jac)/norm[..., None, None]
        deficit = np.maximum(mean_ratio_floor-mr, 0)
        value += 100*np.sum(deficit*deficit)
        dj -= (200*deficit)[..., None, None]*dmr
    return value, dj


def improve(p, h, o, period, iterations=200, condition_limit=None, mean_ratio_floor=None):
    fixed = {int(v) for es in incidence(h, o).values() if len(es) == 1 for v in es[0][2][:, 0]}
    free = sorted(set(h.flat)-fixed)
    if not free:
        return p.copy(), []
    index = {v: i for i, v in enumerate(free)}
    rows, cols = [], []
    for corner, v in enumerate(h.flat):
        if v in index:
            for a in range(3):
                rows.append(3*corner+a); cols.append(3*index[v]+a)
    matrix = csr_matrix((np.ones(len(rows)), (rows, cols)), shape=(h.size*3, len(free)*3))
    x = p[free].ravel().copy()
    offset = (p[h]+o*period).ravel()-matrix@x
    local = np.arange(h.size).reshape(h.shape)
    ds = derivatives(np.array(list(product(np.linspace(-1, 1, 5), repeat=3))))
    bounds = [(None, None), (None, None), (p[:, 2].min()+1e-8, p[:, 2].max()-1e-8)]*len(free)
    best, history = p.copy(), []
    score = min(quality(p, h, o, period)[k] for k in ('min_scaled_jacobian', 'min_mean_ratio'))
    def bounded_objective(x, offset, matrix, local, ds, tau):
        value, gradient = objective(x, offset, matrix, local, ds, tau)
        if condition_limit is not None or mean_ratio_floor is not None:
            coords = (offset+matrix@x).reshape(-1, 3)[local]
            jac = np.einsum('hvi,svj->hsij', coords, ds)
        if condition_limit is not None or mean_ratio_floor is not None:
            penalty, dj = quality_penalty(jac, condition_limit, mean_ratio_floor)
            value += penalty
        if condition_limit is not None or mean_ratio_floor is not None:
            corner_gradient = np.einsum('hsij,svj->hvi', dj, ds)
            gradient += matrix.T@corner_gradient.ravel()
        return value, gradient
    for tau in (.03, .01, .003, .001, .0003):
        result = minimize(bounded_objective, x, args=(offset, matrix, local, ds, tau), jac=True,
                          bounds=bounds, method='L-BFGS-B',
                          options=dict(maxiter=iterations, ftol=1e-13, gtol=1e-8, maxls=40))
        x = result.x
        current = p.copy(); current[free] = x.reshape(-1, 3)
        metrics = quality(current, h, o, period)
        value = min(metrics[k] for k in ('min_scaled_jacobian', 'min_mean_ratio'))
        if (value > score and (condition_limit is None or metrics['max_condition'] <= condition_limit*1.001)
                and (mean_ratio_floor is None or metrics['min_mean_ratio'] >= mean_ratio_floor-1e-4)):
            best, score = current, value
        history.append(dict(tau=tau, iterations=result.nit, success=bool(result.success), **metrics))
    return best, history


def save(path, p, h, o, period, provenance):
    q = p[h]+o*period
    data = dict(description='Periodic HEX8 transition; see separate validation report.',
                canonical_cell_count=len(h), period=period[:2].tolist(), node_order='VTK_HEXAHEDRON',
                provenance=provenance, periodic_mesh=dict(points=p.tolist(), hexes=h.tolist(),
                vertex_offsets=o.tolist()), blocks=[dict(hexes=[dict(points=c.tolist()) for c in q])])
    Path(path).write_text(json.dumps(data, indent=2)+'\n')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--iterations', type=int, default=200)
    parser.add_argument('--condition-limit', type=float)
    parser.add_argument('--mean-ratio-floor', type=float)
    parser.add_argument('--audit-only', action='store_true')
    args = parser.parse_args()
    if args.iterations < 1 or (args.condition_limit is not None and args.condition_limit <= 1):
        parser.error('Use positive iterations and a condition limit greater than one')
    if args.mean_ratio_floor is not None and not 0 < args.mean_ratio_floor <= 1:
        parser.error('Mean-ratio floor must be in (0,1]')
    reference = load(args.input)
    p, h, o, period = reference
    before = quality(*reference, 21)
    if args.audit_only:
        history = []
    else:
        p, history = improve(*reference, args.iterations, args.condition_limit, args.mean_ratio_floor)
        save(args.output, p, h, o, period, dict(input=str(args.input),
             input_sha256=hashlib.sha256(args.input.read_bytes()).hexdigest(),
             iterations_per_stage=args.iterations, condition_limit=args.condition_limit,
             mean_ratio_floor=args.mean_ratio_floor, optimization=history))
    report = validate(p, h, o, period, reference)
    report['before'] = before
    args.output.with_suffix('.validation.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({k: v for k, v in report.items() if k != 'jacobian_certificates'}), flush=True)
    if not report['accepted']:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
