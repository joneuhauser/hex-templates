"""Regression controls for periodic validation, caps and quality gradients."""
import json
from pathlib import Path
import sys
import tempfile
import unittest
from itertools import product
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'tools'))
from mesh_tools import SIGNS
from periodic_transition import load, validate, topology_report, cap_keys, key, quality_penalty
from draw_periodic_boundary import caps


def grid():
    q = [(SIGNS+1)/4+np.array([i/2, j/2, k/2]) for i, j, k in product(range(2), repeat=3)]
    with tempfile.TemporaryDirectory() as tmp:
        path = Path(tmp)/'grid.json'
        path.write_text(json.dumps(dict(period=[1, 1], blocks=[dict(hexes=[dict(points=c.tolist()) for c in q])])))
        return load(path)


class PeriodicTests(unittest.TestCase):
    def test_quality_penalty_gradient(self):
        jac = np.array([[[[.2, .03, -.1], [.02, 1.1, .05], [.01, .1, .7]]]])
        value, gradient = quality_penalty(jac, 4., .9)
        self.assertGreater(value, 0)
        for i in range(3):
            for j in range(3):
                delta = np.zeros_like(jac); delta[0, 0, i, j] = 1e-6
                difference = (quality_penalty(jac+delta, 4., .9)[0]-quality_penalty(jac-delta, 4., .9)[0])/2e-6
                self.assertAlmostEqual(difference, gradient[0, 0, i, j], places=5)

    def test_both_rotated_cap_resolutions(self):
        for count in (2, 8):
            data = caps(count)
            p, faces, o = [np.array(data[k]) for k in ('points', 'quads', 'vertex_offsets')]
            q = p[faces]+o
            self.assertEqual(data['patches'].count('bottom'), 4)
            self.assertEqual(data['patches'].count('top'), count)
            for patch in ('bottom', 'top'):
                area = 0.
                edges = {}
                for i, face in enumerate(faces):
                    if data['patches'][i] != patch:
                        continue
                    area += abs(np.cross(q[i, 1]-q[i, 0], q[i, 3]-q[i, 0])[2])
                    for j in range(4):
                        edge = np.c_[face[[j, (j+1)%4]], o[i, [j, (j+1)%4], :2]]
                        k = key(edge); edges[k] = edges.get(k, 0)+1
                self.assertAlmostEqual(area, 1)
                self.assertEqual(set(edges.values()), {2})

    def test_torus_times_interval(self):
        mesh = grid()
        report = validate(*mesh, mesh)
        self.assertTrue(report['accepted'])
        self.assertEqual(report['betti_GF2'], [1, 2, 1, 0])
        self.assertEqual(report['cap_faces'], [4, 4])
        self.assertAlmostEqual(report['sampled']['min_scaled_jacobian'], 1)

    def test_independent_vertex_and_cell_gauge_changes(self):
        p, h, o, period = grid()
        rng = np.random.default_rng(4)
        shift = rng.integers(-3, 4, (len(p), 3)); shift[:, 2] = 0
        deck = rng.integers(-3, 4, (len(h), 1, 3)); deck[:, :, 2] = 0
        moved = (p+shift*period, h, o-shift[h]+deck, period)
        self.assertEqual(topology_report(p, h, o), topology_report(*moved[:3]))
        self.assertEqual(cap_keys(p, h, o, period), cap_keys(*moved))
        self.assertTrue(validate(*moved, (p, h, o, period))['accepted'])

    def test_broken_seam_and_missing_cell_rejected(self):
        p, h, o, period = grid()
        broken = o.copy(); broken[0, 0, 0] += 1
        self.assertFalse(validate(p, h, broken, period)['accepted'])
        self.assertFalse(validate(p, h[1:], o[1:], period)['accepted'])

    def test_negative_orientation_rejected(self):
        p, h, o, period = grid()
        p[:, 2] *= -1
        self.assertFalse(validate(p, h, o, period)['accepted'])

if __name__ == '__main__':
    unittest.main()
