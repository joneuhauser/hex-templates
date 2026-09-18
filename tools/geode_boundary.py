#!/usr/bin/env python3
"""Assemble fixed Geode caps with four copies of a supplied planar side mesh.

Side JSON: {"points": [[u,z], ...], "quads": [[i,j,k,l], ...]} with zero-based
counterclockwise quads, u in [-1,1], z in [-7/12,7/12]. This is a boundary
builder, not an exhaustive planar-quadrangulation enumerator.
"""
import argparse
from collections import Counter, defaultdict
import json
from pathlib import Path

import numpy as np

from mesh_tools import read_mesh, write_mesh
from prepare_geode import ROOT, boundary_contract, translation_matching, cycle


TWO_QUADS = dict(points=[[-1, -7/12], [0, -7/12], [1, -7/12],
                        [1, 7/12], [0, 7/12], [-1, 7/12]],
                 quads=[[0, 1, 4, 5], [1, 2, 3, 4]])


def assemble(reference, side):
    p, _, q = read_mesh(reference)
    contract = boundary_contract(p, q)
    coords = np.asarray(side['points'], dtype=float)
    quads = np.asarray(side['quads'], dtype=int)
    if coords.ndim != 2 or coords.shape[1] != 2 or not np.isfinite(coords).all():
        raise ValueError('Expected finite (u,z) coordinates')
    if quads.ndim != 2 or quads.shape[1] != 4 or not len(quads) or quads.min() < 0 or quads.max() >= len(coords):
        raise ValueError('Expected indexed quadrilaterals')
    if np.any(abs(coords) > np.array([1, 7/12]) + 1e-13):
        raise ValueError('Side vertices leave the side rectangle')
    for f in quads:
        if len(set(f)) != 4:
            raise ValueError('Repeated quad vertex')
        xy = coords[f]
        for i in range(4):
            a, b = xy[(i+1)%4] - xy[i], xy[(i-1)%4] - xy[i]
            if a[0]*b[1] - a[1]*b[0] <= 0:
                raise ValueError('Side quad must be strictly convex and counterclockwise')
    points, faces = [], []
    def vertex(x):
        found = [i for i, v in enumerate(points) if np.max(abs(np.asarray(v) - x)) < 1e-13]
        if found:
            return found[0]
        points.append(list(x))
        return len(points)-1
    for f in q[contract['patch_face_ids']['bottom'] + contract['patch_face_ids']['top']]:
        faces.append([vertex(p[v]) for v in f])
    for turn in range(4):
        indices = []
        for u, z in coords:
            x = np.array([1., u, z])
            for _ in range(turn):
                x = np.array([-x[1], x[0], x[2]])
            indices.append(vertex(x))
        faces.extend([[indices[v] for v in f] for f in quads])
    edges = Counter()
    oriented = Counter()
    adjacency = defaultdict(set)
    for f in faces:
        for i in range(4):
            a, b = f[i], f[(i+1)%4]
            edges[tuple(sorted([a, b]))] += 1
            oriented[a, b] += 1
            adjacency[a].add(b)
            adjacency[b].add(a)
    todo, seen = [0], set()
    while todo:
        v = todo.pop()
        if v not in seen:
            seen.add(v)
            todo.extend(adjacency[v] - seen)
    if (set(edges.values()) != {2} or any(oriented[a, b] != 1 or oriented[b, a] != 1 for a, b in edges)
            or len(points)-len(edges)+len(faces) != 2 or len(seen) != len(points)):
        raise ValueError('Side mesh does not glue to the fixed caps as an oriented quadrangulated sphere')
    points, faces = np.asarray(points), np.asarray(faces)
    boundary_contract(points, faces)
    if not translation_matching(points, faces):
        raise ValueError('Opposite side meshes must match by translation')
    return points, faces


def boundary_mapping(p, q, target_p, target_q):
    """Oriented surface isomorphism fixing cap coordinates, with free side nodes."""
    n = len(p)
    if n != len(target_p) or len(q) != len(target_q):
        return None
    def relations(faces):
        r = np.zeros((n, n), dtype=int)
        for f in faces:
            for i in range(4):
                for j in range(4):
                    if i != j:
                        r[f[i], f[j]] |= 1 if (i-j)%2 else 2
        return r
    a, b = relations(q), relations(target_q)
    mapping = {}
    for v in range(n):
        if abs(abs(p[v, 2])-7/12) < 1e-12:
            matches = np.flatnonzero(np.max(abs(target_p-p[v]), axis=1) < 1e-12)
            if len(matches) != 1:
                return None
            mapping[v] = int(matches[0])
    target = {cycle(tuple(f)) for f in target_q}
    def solve():
        if len(mapping) == n:
            return {cycle(tuple(mapping[v] for v in f)) for f in q} == target
        available = {v: [w for w in range(n) if w not in mapping.values()
                        and sorted(a[v]) == sorted(b[w])
                        and all(a[v, x] == b[w, y] for x, y in mapping.items())]
                     for v in range(n) if v not in mapping}
        v = min(available, key=lambda x: (len(available[x]), x))
        for w in available[v]:
            mapping[v] = w
            if solve():
                return True
            del mapping[v]
        return False
    return [mapping[v] for v in range(n)] if solve() else None

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--reference', type=Path, default=ROOT / 'docs/assets/geode/meshes/G26-Mitchell.mesh')
    parser.add_argument('--side', type=Path, help='Side JSON; default is the two-quad side')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    p, q = assemble(args.reference, json.loads(args.side.read_text()) if args.side else TWO_QUADS)
    write_mesh(args.output, p, np.empty((0, 8), dtype=int), q)
    print(f'GEODE_BOUNDARY vertices={len(p)} quads={len(q)} side_quads={(len(q)-10)//4}')
