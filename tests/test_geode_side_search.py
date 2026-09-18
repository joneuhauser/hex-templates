"""Boundary placements and diagonal symmetry for the Geode search batch."""
import json
from pathlib import Path
import sys
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'tools'))
from geode_boundary import assemble
from search_geode_sides import diagonal_action, diagonal_actions, placements
from geode_geometry import variables
from prepare_geode import boundary_contract, translation_matching


class GeodeSideSearchTests(unittest.TestCase):
    def test_movable_edge_heights_and_interior_orbits(self):
        catalog = json.loads((ROOT/'docs/assets/geode-sides/catalog.json').read_text())
        side = next(s for s in catalog['patterns'] if s['id'] == 'Q6-01')
        bp, q = assemble(ROOT/'docs/assets/geode/meshes/G26-Mitchell.mesh', side)
        # Generic, first-mirror, second-mirror and axis vertex orbits.
        p = np.vstack([bp, [[.25, .35, .1], [.35, .25, .1], [-.35, -.25, .1], [-.25, -.35, .1],
                           [.2, .2, 0], [-.2, -.2, 0], [.2, -.2, -.1], [-.2, .2, -.1], [0, 0, .2]]])
        actions = np.asarray(diagonal_actions(p, q, 2))
        offset, matrix, x, bounds = variables(p, q, vertex_actions=actions)
        np.testing.assert_allclose(offset+matrix@x, p.ravel(), atol=1e-15)
        self.assertLessEqual(int(np.max(np.sum(matrix != 0, axis=1))), 1)
        moved = (offset+matrix@(x+0.002*np.sin(np.arange(len(x))+1))).reshape(-1, 3)
        fixed = boundary_contract(bp, q)['fixed_cap_vertices']
        np.testing.assert_array_equal(moved[fixed], p[fixed])
        self.assertTrue(translation_matching(moved, q))
        self.assertEqual(diagonal_actions(moved, q, 2), actions.tolist())
        edges = [v for v in range(len(bp)) if v not in fixed and abs(bp[v, 0]) == abs(bp[v, 1]) == 1]
        self.assertTrue(np.any(moved[edges, 2] != p[edges, 2]))
        np.testing.assert_array_equal(moved[edges, :2], p[edges, :2])
        self.assertTrue(np.all(abs(moved[:, :2]) <= 1))
        self.assertTrue(np.all(abs(moved[:, 2]) <= 7/12))

    def test_all_small_placements_preserve_diagonal_and_translation(self):
        catalog = json.loads((ROOT/'docs/assets/geode-sides/catalog.json').read_text())
        cases = list(placements(catalog, 6))
        self.assertEqual([name for name, _, _ in cases],
                         ['Q2-01', 'Q4-01', 'Q5-01', 'Q5-01-flipped', 'Q5-02', 'Q5-03', 'Q6-01'])
        self.assertEqual(sum(flipped for _, _, flipped in cases), 1)
        geometries = {}
        for name, side, _ in cases:
            # assemble independently checks cap matching and translated side quads.
            p, q = assemble(ROOT/'docs/assets/geode/meshes/G26-Mitchell.mesh', side)
            action = np.asarray(diagonal_action(p, q))
            np.testing.assert_array_equal(action[action], np.arange(len(p)))
            np.testing.assert_allclose(p[action], p[:, [1, 0, 2]], atol=1e-14)
            group = np.asarray(diagonal_actions(p, q, 2))
            for g in range(4):
                for h in range(4):
                    np.testing.assert_array_equal(group[group[:, g], h], group[:, g ^ h])
            np.testing.assert_array_equal(group[:, 1], action)
            np.testing.assert_allclose(p[group[:, 2]], p[:, [1, 0, 2]]*[-1, -1, 1], atol=1e-14)
            np.testing.assert_allclose(p[group[:, 3]], p*[-1, -1, 1], atol=1e-14)
            geometries[name] = {tuple(x) for x in p}
            all_fixed = np.all(abs(abs(p[:, 2])-7/12) < 1e-12)
            self.assertEqual(bool(all_fixed), name == 'Q2-01')
            broken = p.copy()
            broken[0, 0] += 0.01
            with self.assertRaises(ValueError):
                diagonal_action(broken, q)
        self.assertNotEqual(geometries['Q5-01'], geometries['Q5-01-flipped'])


if __name__ == '__main__':
    unittest.main()
