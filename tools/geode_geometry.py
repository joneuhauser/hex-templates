"""Geode embedding with fixed caps, tied side coordinates and free interiors."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.optimize import minimize
from scipy.special import logsumexp

from mesh_tools import (audit, bernstein, boundary, derivatives,
                        read_mesh, sampled, write_mesh)
from mesh_homology import homology
from prepare_geode import boundary_contract, face_set, translation_matching
from validate_mesh import manifold


def variables(points, quads, side_reversal=True, vertex_actions=None):
    """p = offset + A*x, tying side copies and optional diagonal vertex orbits."""
    contract = boundary_contract(points, quads)
    fixed = set(contract['fixed_cap_vertices'])
    boundary_vertices = set(quads.flat)
    columns, values, bounds = [], [], []
    master = contract['master_side_vertices']
    seen = set()
    for index, v in enumerate(master):
        if v in fixed:
            continue
        on_edge = abs(abs(points[v, 1]) - 1) < 1e-12
        if v in seen:
            continue
        partners = [(index, 1)]
        if side_reversal:
            reflected = points[v].copy()
            reflected[1] *= -1
            matches = [j for j, w in enumerate(master) if np.max(abs(points[w]-reflected)) < 1e-12]
            if len(matches) != 1:
                raise ValueError('Translation matching requires a reversal-symmetric initial side patch')
            partner = matches[0]
            if partner != index:
                partners.append((partner, -1))
        seen.update(master[j] for j, _ in partners)
        for axis in ([2] if on_edge else [1, 2]):
            if side_reversal and len(partners) == 1 and axis == 1:
                continue  # Center-line side vertex stays on u=0.
            column = np.zeros_like(points)
            for member, sign in partners:
                direction = np.eye(3)[axis] * (sign if axis == 1 else 1)
                for turn in range(4):
                    w = contract['side_vertex_copies'][turn][member]
                    column[w] = direction
                    direction = np.array([-direction[1], direction[0], direction[2]])
            columns.append(column.ravel())
            values.append(points[v, axis])
            bound = 1. if axis == 1 else 7/12
            bounds.append((-bound + 1e-8, bound - 1e-8))
    rotations = np.asarray([np.eye(3), [[0, 1, 0], [1, 0, 0], [0, 0, 1]],
                            [[0, -1, 0], [-1, 0, 0], [0, 0, 1]], np.diag([-1, -1, 1])])
    if vertex_actions is None:
        actions = np.arange(len(points))[:, None]
    else:
        actions = np.asarray(vertex_actions, dtype=int)
        if actions.shape not in [(len(points), 2), (len(points), 4)]:
            raise ValueError('Expected one- or two-diagonal-mirror vertex actions')
        if actions.min() < 0 or actions.max() >= len(points):
            raise ValueError('Invalid symmetry vertex index')
        if not np.array_equal(actions[:, 0], np.arange(len(points))):
            raise ValueError('Invalid symmetry identity action')
        for g in range(actions.shape[1]):
            for h in range(actions.shape[1]):
                if not np.array_equal(actions[actions[:, g], h], actions[:, g ^ h]):
                    raise ValueError('Invalid symmetry group action')
            if not np.allclose(points[actions[:, g]], points @ rotations[g].T, atol=1e-12, rtol=0):
                raise ValueError('Initial coordinates do not realize the diagonal action')
    used = set()
    for v in sorted(set(range(len(points))) - boundary_vertices):
        if v in used:
            continue
        orbit = set(actions[v])
        if orbit & boundary_vertices:
            raise ValueError('Symmetry mixes interior and boundary vertices')
        used.update(orbit)
        stabilizer = [g for g, w in enumerate(actions[v]) if w == v]
        projection = rotations[stabilizer].mean(axis=0)
        directions = []
        for direction in projection.T:
            nonzero = np.flatnonzero(abs(direction) > 1e-12)
            if not len(nonzero):
                continue
            axis = int(nonzero[0])
            direction = direction/direction[axis]
            if np.linalg.matrix_rank(np.asarray(directions+[direction])) <= len(directions):
                continue
            directions.append(direction)
            column = np.zeros_like(points)
            for g, w in enumerate(actions[v]):
                column[w] = rotations[g] @ direction
            columns.append(column.ravel())
            values.append(points[v, axis])
            bound = 1. if axis < 2 else 7/12
            bounds.append((-bound + 1e-8, bound - 1e-8))
    matrix = np.stack(columns, axis=1) if columns else np.empty((points.size, 0))
    x = np.asarray(values)
    offset = points.ravel() - matrix @ x
    return offset, matrix, x, bounds


def objective(x, offset, matrix, cells, ds, tau):
    p = (offset + matrix @ x).reshape(-1, 3)
    jac = np.einsum('hvi,svj->hsij', p[cells], ds)
    cof = np.stack([np.cross(jac[..., 1], jac[..., 2]),
                    np.cross(jac[..., 2], jac[..., 0]),
                    np.cross(jac[..., 0], jac[..., 1])], axis=-1)
    det = np.sum(jac[..., 0] * cof[..., 0], axis=-1)
    norms = np.sum(jac * jac, axis=2) + 1e-30
    lengths = np.sqrt(np.prod(norms, axis=-1))
    sj = det / lengths
    dsj = cof / lengths[..., None, None] - sj[..., None, None] * jac / norms[:, :, None, :]
    # Smooth positive determinant barrier remains finite while untangling.
    eps = 1e-8
    root = np.sqrt(det * det + eps * eps)
    chi = np.where(det >= 0, (det + root) / 2, eps * eps / (2 * (root - np.minimum(det, 0))))
    frobenius = norms.sum(axis=-1)
    mean_ratio = 3 * chi**(2/3) / frobenius
    dmean = mean_ratio[..., None, None] * ((2/3)*cof/root[..., None, None] - 2*jac/frobenius[..., None, None])
    z = -np.stack([sj, mean_ratio]) / tau
    logz = logsumexp(z)
    weights = np.exp(z - logz)
    barrier = -0.005 * np.mean(np.log(chi))
    dj = -weights[0, ..., None, None] * dsj - weights[1, ..., None, None] * dmean - 0.005 * cof / root[..., None, None] / det.size
    local = np.einsum('hsij,svj->hvi', dj, ds)
    gradient = np.zeros_like(p)
    np.add.at(gradient, cells.ravel(), local.reshape(-1, 3))
    return tau * logz + barrier, matrix.T @ gradient.ravel()


def improve(points, cells, quads, iterations=250, side_reversal=True, vertex_actions=None):
    from itertools import product
    offset, matrix, x, bounds = variables(points, quads, side_reversal, vertex_actions)
    ds = derivatives(np.asarray(list(product(np.linspace(-1, 1, 5), repeat=3))))
    def score(p):
        jac = np.einsum('hvi,svj->hsij', p[cells], ds)
        det = np.linalg.det(jac)
        sj = det/(np.prod(np.linalg.norm(jac, axis=2), axis=2)+1e-30)
        mr = 3*np.sign(det)*abs(det)**(2/3)/(np.sum(jac*jac, axis=(2, 3))+1e-30)
        return float(min(sj.min(), mr.min()))
    history = []
    best = points.copy()
    best_quality = score(points)
    for tau in [0.08, 0.03, 0.01, 0.003, 0.001]:
        result = minimize(objective, x, args=(offset, matrix, cells, ds, tau),
                          jac=True, bounds=bounds, method='L-BFGS-B',
                          options=dict(maxiter=iterations, ftol=1e-13, gtol=1e-8, maxls=40))
        x = result.x
        current = (offset + matrix @ x).reshape(-1, 3)
        quality = score(current)
        if quality > best_quality:
            best, best_quality = current.copy(), quality
        history.append(dict(tau=tau, iterations=result.nit, success=bool(result.success),
                            min_sampled_shape_score=quality, message=result.message))
    return best, history


def validate(path, reference, require_translation=True):
    """Validate these candidates against the reference cap geometry/connectivity.

    Boundary and reference IDs need not agree. Side equality checks vertices AND
    quad cycles; Jacobians are certified over entire trilinear cells.
    """
    p, h, q = read_mesh(path)
    rp, _, rq = read_mesh(reference)
    report = audit(path)
    report['manifold'] = manifold(h)
    report['homology'] = homology(h)
    report['sampled'] = sampled(p, h, 15)
    try:
        contract = boundary_contract(p, boundary(h))
        reference_contract = boundary_contract(rp, rq)
        caps_match = True
        b = boundary(h)
        for cap in ['top', 'bottom']:
            faces = b[contract['patch_face_ids'][cap]]
            target = rq[reference_contract['patch_face_ids'][cap]]
            mapping = {}
            for v in set(faces.flat):
                matches = np.flatnonzero(np.max(abs(rp - p[v]), axis=1) < 1e-12)
                if len(matches) != 1:
                    caps_match = False
                    break
                mapping[v] = int(matches[0])
            if caps_match:
                caps_match &= face_set([[mapping[v] for v in f] for f in faces]) == face_set(target)
        report['caps_match'] = bool(caps_match)
        report['equal_sides'] = True
        report['translation_matching'] = translation_matching(p, b, contract)
    except ValueError as error:
        report.update(caps_match=False, equal_sides=False, translation_matching=False, boundary_error=str(error))
    volume = sum(float(sum(bernstein(p[c]).flat)) * 8/27 for c in h)
    report['volume_error'] = volume - 14/3
    report['inside_box'] = bool(np.all(abs(p[:, :2]) <= 1 + 1e-12) and np.all(abs(p[:, 2]) <= 7/12 + 1e-12))
    # Positive planar quad corner determinants give a nonfolded bilinear side.
    planar_positive = True
    for f in boundary(h):
        xy = p[f]
        normal = np.cross(xy[1] - xy[0], xy[3] - xy[0])
        planar_positive &= all(np.dot(np.cross(xy[(i+1)%4] - xy[i], xy[(i-1)%4] - xy[i]), normal) > 0 for i in range(4))
    report['planar_boundary_quads_positive'] = bool(planar_positive)
    report['accepted'] = bool(report['all_jacobians_positive'] and report['boundary_matches']
        and report['opposite_internal_orientations'] and not report['incompatible_pairs']
        and report['distinct_cell_vertices'] and report['euler'] == 1
        and set(report['face_multiplicity']) <= {1, 2}
        and report['manifold']['edge_links_valid'] and report['manifold']['vertex_links_valid']
        and report['homology']['betti_GF2'] == [1, 0, 0, 0]
        and report['caps_match'] and report['equal_sides'] and report['inside_box']
        and (report['translation_matching'] or not require_translation)
        and report['planar_boundary_quads_positive'] and abs(report['volume_error']) < 1e-10)
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mesh', type=Path)
    parser.add_argument('--reference', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--iterations', type=int, default=250)
    parser.add_argument('--rotated-sides-only', action='store_true')
    args = parser.parse_args()
    p, h, q = read_mesh(args.mesh)
    p, history = improve(p, h, q, args.iterations, not args.rotated_sides_only)
    write_mesh(args.output, p, h, q)
    report = validate(args.output, args.reference, not args.rotated_sides_only)
    report['optimization'] = history
    args.output.with_suffix('.validation.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({k: report[k] for k in ['accepted', 'hexes', 'sampled', 'caps_match', 'equal_sides']}))
    raise SystemExit(0 if report['accepted'] else 1)
