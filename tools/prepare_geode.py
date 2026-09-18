#!/usr/bin/env python3
"""Import Mitchell's CUBIT journal and audit a fixed-boundary replay fixture.

This prepares the reference case, not the free-side topology search.
The journal is parsed as data; CUBIT commands are never executed.
"""
import argparse
import hashlib
import json
from pathlib import Path
import sys

import numpy as np

from mesh_tools import boundary, certify_hex, sampled, topology, write_mesh
from mesh_homology import homology
from validate_mesh import manifold

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'solver'))
from compare_meshes import equivalent, isomorphism


def read_journal(path):
    points, cells = [], []
    for line in Path(path).read_text().splitlines():
        fields = line.split()
        if fields[:2] == ['create', 'node']:
            points.append([float(x) for x in fields[2:5]])
        elif fields[:3] == ['create', 'hex', 'node']:
            cells.append([int(x) - 1 for x in fields[3:11]])
    p, h = np.asarray(points), np.asarray(cells, dtype=int)
    if p.shape != (48, 3) or h.shape != (26, 8):
        raise ValueError('Expected the 48-node, 26-hex Mitchell journal')
    if not np.isfinite(p).all() or h.min() < 0 or h.max() >= len(p):
        raise ValueError('Invalid coordinates or node references')
    return p, h


def normalize(p, h):
    # Proper rotation and uniform scale: preserve orientation and shape quality.
    p = p[:, [0, 2, 1]] * np.array([1., -1., 1.]) / 12
    p[np.abs(p) < 1e-13] = 0.  # Remove journal roundoff at nominal zero only.
    q = boundary(h)
    bv = sorted(set(q.flat))
    order = bv + sorted(set(range(len(p))) - set(bv))
    inverse = np.argsort(order)
    return p[order], inverse[h], inverse[q], order, len(bv)


def cycle(face):
    face = tuple(face)
    return min(face[k:] + face[:k] for k in range(4))


def face_set(quads):
    return {min(cycle(q), cycle(tuple(q)[::-1])) for q in quads}


def translation_matching(p, quads, contract=None):
    contract = boundary_contract(p, quads) if contract is None else contract
    for first, second, shift in [('east', 'west', [-2, 0, 0]), ('north', 'south', [0, -2, 0])]:
        faces = quads[contract['patch_face_ids'][first]]
        target = quads[contract['patch_face_ids'][second]]
        ids = sorted(set(target.flat))
        mapping = {}
        for v in set(faces.flat):
            matches = [w for w in ids if np.max(abs(p[w] - p[v] - shift)) < 1e-12]
            if len(matches) != 1:
                return False
            mapping[v] = matches[0]
        if face_set([[mapping[v] for v in f] for f in faces]) != face_set(target):
            return False
    return True


def boundary_contract(p, q):
    """Record fixed caps and four identical side charts (all IDs zero-based)."""
    patches = {}
    for name, axis, value in [('bottom', 2, -7/12), ('top', 2, 7/12),
                              ('east', 0, 1), ('north', 1, 1),
                              ('west', 0, -1), ('south', 1, -1)]:
        patches[name] = np.flatnonzero(np.all(abs(p[q, axis] - value) < 1e-12, axis=1)).tolist()
    if sorted(i for ids in patches.values() for i in ids) != list(range(len(q))):
        raise ValueError('Boundary does not partition into the six box faces')
    fixed = sorted(set(q[patches['top'] + patches['bottom']].flat))
    master = sorted(set(q[patches['east']].flat))
    maps = []
    for turn, name in enumerate(['east', 'north', 'west', 'south']):
        target = sorted(set(q[patches[name]].flat))
        mapping = {}
        for v in master:
            x = p[v].copy()
            for _ in range(turn):
                x = np.array([-x[1], x[0], x[2]])
            found = [w for w in target if np.max(abs(p[w] - x)) < 1e-12]
            if len(found) != 1:
                raise ValueError('Side vertices do not match by quarter turns')
            mapping[v] = found[0]
        if face_set([[mapping[v] for v in f] for f in q[patches['east']]]) != face_set(q[patches[name]]):
            raise ValueError('Side quadrilateral connectivities differ')
        maps.append([int(mapping[v]) for v in master])
    return dict(patch_face_ids=patches, fixed_cap_vertices=[int(v) for v in fixed],
                master_side_vertices=[int(v) for v in master], side_vertex_copies=maps,
                constraint='Side copies only; no action is imposed on volume cells or interior vertices.')


def symmetry_audit(p, h, q, nb):
    """Test extensions of every cap-preserving horizontal square isometry."""
    mesh = dict(vertices=len(p), hexes=len(h), cells=h.tolist())
    result = []
    for swap in [False, True]:
        for sx in [1, -1]:
            for sy in [1, -1]:
                image = p.copy()
                image[:, :2] = p[:, [1, 0] if swap else [0, 1]] * [sx, sy]
                perm = []
                for x in image[:nb]:
                    found = np.flatnonzero(np.max(abs(p[:nb] - x), axis=1) < 1e-12)
                    if len(found) != 1:
                        break
                    perm.append(int(found[0]))
                if len(perm) != nb or face_set(np.asarray(perm)[q]) != face_set(q):
                    continue
                odd = (-1 if swap else 1) * sx * sy < 0
                cells = np.asarray(perm + list(range(nb, len(p))))[h]
                if odd:
                    cells = cells[:, [1, 0, 3, 2, 5, 4, 7, 6]]
                extends = equivalent(mesh, dict(mesh, cells=cells.tolist()), nb)
                result.append(dict(swap_xy=swap, sign_x=sx, sign_y=sy,
                                   reflection=odd, extends_to_cell_complex=extends))
    return result


def half_turn_actions(p, h, q):
    """Extend the boundary half-turn combinatorially, allowing asymmetric coordinates."""
    nb = len(set(q.flat))
    if set(q.flat) != set(range(nb)):
        raise ValueError('Expected boundary vertices first')
    image = p[:nb] * [-1, -1, 1]
    perm = []
    for point in image:
        matches = np.flatnonzero(np.max(abs(p[:nb]-point), axis=1) < 1e-12)
        if len(matches) != 1:
            raise ValueError('Boundary is not invariant under the half-turn')
        perm.append(int(matches[0]))
    mesh = dict(vertices=len(p), hexes=len(h), cells=h.tolist())
    mapping = isomorphism(mesh, mesh, nb, perm, involution=True)
    if mapping is None:
        raise ValueError('Boundary half-turn does not extend to oriented cells')
    return [[v, mapping[v]] for v in range(len(p))]


def prepare(source, output):
    raw_p, raw_h = read_journal(source)
    p, h, q, order, nb = normalize(raw_p, raw_h)
    report = topology(h, q)
    report.update(source_url='https://www.sandia.gov/files/samitch/files/geode.jou',
                  source_sha256=hashlib.sha256(Path(source).read_bytes()).hexdigest(),
                  transform='(x,y,z) -> (x,-z,y)/12; abs(coordinate)<1e-13 snapped to zero',
                  source_node_ids=[int(v) + 1 for v in order], boundary_vertices=nb,
                  interior_vertices=len(p)-nb, contract=boundary_contract(p, q),
                  sampled=sampled(p, h, 9), manifold=manifold(h), homology=homology(h),
                  jacobian_certificates=[certify_hex(p[c], exact=True) for c in h],
                  boundary_isometry_extensions=symmetry_audit(p, h, q, nb))
    if not all(c['positive'] for c in report['jacobian_certificates']):
        raise ValueError('Journal fails whole-cell positive-Jacobian certification')
    output.mkdir(parents=True, exist_ok=True)
    write_mesh(output / 'geode-boundary.mesh', p[:nb], np.empty((0, 8), dtype=int), q)
    write_mesh(output / 'geode-seed.mesh', p, h, q)
    np.savetxt(output / 'geode-half-turn.actions', half_turn_actions(p, h, q), fmt='%d')
    (output / 'geode-reference.json').write_text(json.dumps(report, indent=2) + '\n')
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, default=ROOT / 'docs/assets/geode/geode.jou')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    report = prepare(args.source, args.output)
    print(f"GEODE_REFERENCE hexes={report['hexes']} boundary_vertices={report['boundary_vertices']} "
          f"boundary_quads={report['boundary_quads']} exact_positive=26/26")
