#!/usr/bin/env python3
"""Compare tiny candidate sets modulo interior labels, fixing boundary vertices.

This checks all cube faces and their orientations, not just cell vertex sets.
Uses only the standard library; intended for the small positive fixtures.
"""
import json
import sys

FACES = [(0, 3, 2, 1), (4, 5, 6, 7), (0, 1, 5, 4),
         (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7)]
BITS = [0, 1, 3, 2, 4, 5, 7, 6]


def relations(mesh):
    n = mesh['vertices']
    result = [[0] * n for _ in range(n)]
    for cell in mesh['cells']:
        for i, a in enumerate(cell):
            for j, b in enumerate(cell):
                if i != j:
                    result[a][b] |= 1 << (bin(BITS[i] ^ BITS[j]).count("1") - 1)
    return result


def oriented_cells(mesh, mapping):
    def cycle(values):
        return min(values[i:] + values[:i] for i in range(4))
    return sorted(tuple(sorted(cycle(tuple(mapping[cell[j]] for j in face))
                               for face in FACES)) for cell in mesh['cells'])


def isomorphism(a, b, boundary, boundary_mapping=None, involution=False):
    """Return an oriented mesh isomorphism extending the prescribed boundary."""
    if (a['vertices'], a['hexes']) != (b['vertices'], b['hexes']):
        return None
    n = a['vertices']
    ar, br = relations(a), relations(b)
    mapping = dict(enumerate(range(boundary) if boundary_mapping is None else boundary_mapping))
    if set(mapping) != set(range(boundary)) or set(mapping.values()) != set(range(boundary)):
        raise ValueError('Boundary mapping must be a permutation')
    if any(ar[v][w] != br[mapping[v]][mapping[w]] for v in mapping for w in mapping):
        return None
    if involution and any(mapping[mapping[v]] != v for v in mapping):
        return None

    def signature(r, v, order):
        return tuple(r[v][w] for w in order), tuple(sorted(r[v]))

    choices = {v: [w for w in range(boundary, n)
                   if signature(ar, v, range(boundary)) == signature(br, w, mapping.values())]
               for v in range(boundary, n)}
    target = oriented_cells(b, dict(enumerate(range(n))))

    def solve():
        if len(mapping) == n:
            return oriented_cells(a, mapping) == target
        used = set(mapping.values())
        available = {v: [w for w in choices[v] if w not in used and
                         all(ar[v][x] == br[w][y] for x, y in mapping.items())]
                     for v in choices if v not in mapping}
        v = min(available, key=lambda x: (len(available[x]), x))
        for w in available[v]:
            if involution and ((w in mapping and mapping[w] != v) or
                               (w not in mapping and v not in available[w])):
                continue
            mapping[v] = w
            paired = involution and w != v and w not in mapping
            if paired:
                mapping[w] = v
            if solve():
                return True
            if paired:
                del mapping[w]
            del mapping[v]
        return False

    return dict(mapping) if solve() else None


def equivalent(a, b, boundary):
    return isomorphism(a, b, boundary) is not None


def main():
    boundary = int(sys.argv[1])
    left, right = [[json.loads(line) for line in open(path) if line.strip()]
                   for path in sys.argv[2:]]
    assert len(left) == len(right), 'Candidate counts differ'
    for mesh in left:
        index = next((i for i, other in enumerate(right)
                      if equivalent(mesh, other, boundary)), None)
        assert index is not None, 'Candidate connectivity or orientation differs'
        right.pop(index)
    print(f'CANDIDATES_EQUIVALENT count={len(left)} boundary_fixed={boundary}')


if __name__ == '__main__':
    main()
