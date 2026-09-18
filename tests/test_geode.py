"""Reference provenance, boundary constraints, derivatives and solver replay."""
import copy
from itertools import product
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'tools'))
from prepare_geode import prepare
from geode_boundary import assemble, TWO_QUADS
from geode_geometry import objective, variables, validate
from mesh_tools import derivatives, read_mesh, write_mesh


class GeodeTests(unittest.TestCase):
    reference = ROOT / 'docs/assets/geode/meshes/G26-Mitchell.mesh'

    def test_reference_and_provenance(self):
        with tempfile.TemporaryDirectory() as tmp:
            report = prepare(ROOT / 'docs/assets/geode/geode.jou', Path(tmp))
            self.assertEqual((report['hexes'], report['boundary_vertices'], report['boundary_quads']), (26, 32, 30))
            self.assertTrue(all(x['positive'] for x in report['jacobian_certificates']))
            self.assertEqual([len(x) for x in report['contract']['patch_face_ids'].values()], [4, 6, 5, 5, 5, 5])
            self.assertEqual(Path(tmp, 'geode-seed.mesh').read_bytes(), self.reference.read_bytes())
        report = validate(self.reference, self.reference)
        self.assertTrue(report['accepted'])
        self.assertTrue(report['translation_matching'])

    def test_fixed_caps_and_side_constraints(self):
        p, h, q = read_mesh(self.reference)
        offset, matrix, x, _ = variables(p, q)
        self.assertEqual(len(x), 51)  # 48 independent interior coordinates + 3 side parameters.
        moved = (offset + matrix @ (x + np.linspace(-0.001, 0.001, len(x)))).reshape(-1, 3)
        fixed = np.flatnonzero(abs(abs(p[:, 2]) - 7/12) < 1e-12)
        np.testing.assert_array_equal(moved[fixed], p[fixed])
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp, 'candidate.mesh')
            write_mesh(path, moved, h, q)
            result = validate(path, self.reference)
            self.assertTrue(result['accepted'])
            self.assertTrue(result['translation_matching'])
            moved[fixed[0], 0] += 1e-3
            write_mesh(path, moved, h, q)
            self.assertFalse(validate(path, self.reference)['accepted'])

    def test_optimizer_gradient(self):
        p, h, q = read_mesh(self.reference)
        offset, matrix, x, _ = variables(p, q)
        ds = derivatives(np.asarray(list(product([-1., 0., 1.], repeat=3))))
        _, gradient = objective(x, offset, matrix, h, ds, 0.03)
        for i in range(len(x)):
            delta = np.zeros_like(x)
            delta[i] = 1e-6
            numeric = (objective(x+delta, offset, matrix, h, ds, 0.03)[0]
                       - objective(x-delta, offset, matrix, h, ds, 0.03)[0]) / 2e-6
            self.assertAlmostEqual(gradient[i], numeric, delta=1e-6)

    def test_rotated_congruence_is_not_translation_matching(self):
        p, h, q = read_mesh(self.reference)
        offset, matrix, x, _ = variables(p, q, side_reversal=False)
        x[0] += 0.001
        moved = (offset + matrix @ x).reshape(-1, 3)
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp, 'rotated-only.mesh')
            write_mesh(path, moved, h, q)
            report = validate(path, self.reference)
            self.assertTrue(report['equal_sides'])
            self.assertFalse(report['translation_matching'])
            self.assertFalse(report['accepted'])

    def test_side_builder(self):
        p, q = assemble(self.reference, TWO_QUADS)
        self.assertEqual((len(p), len(q)), (20, 18))
        invalid = copy.deepcopy(TWO_QUADS)
        invalid['points'][1][0] += 0.1  # Mismatched cap midpoint must not be accepted.
        with self.assertRaises(ValueError):
            assemble(self.reference, invalid)
        invalid = copy.deepcopy(TWO_QUADS)
        invalid['quads'][0].reverse()
        with self.assertRaises(ValueError):
            assemble(self.reference, invalid)

    def test_candidates(self):
        meshes = sorted((ROOT/'docs/assets/geode/meshes').glob('G28-*.mesh'))
        self.assertEqual(len(meshes), 2)
        for mesh in meshes:
            report = validate(mesh, self.reference)
            self.assertTrue(report['accepted'], mesh.name)
            self.assertTrue(report['translation_matching'], mesh.name)
            self.assertEqual(report['hexes'], 28)

    def test_solver_replays(self):
        binary = Path(os.environ.get('GEODE_SOLVER', ROOT / 'solver/build/mirror_search'))
        if not binary.is_file():
            self.skipTest('Build solver or set GEODE_SOLVER to run seeded replay checks')
        paths = [self.reference] + list((ROOT/'docs/assets/geode/meshes').glob('G28-*.mesh'))
        # Input cavity shells contain no cells; only replay complete candidate meshes.
        for path in paths:
            p, h, q = read_mesh(path)
            if not len(h):
                continue
            nb = len(set(q.flat))
            self.assertEqual(set(q.flat), set(range(nb)))
            with tempfile.TemporaryDirectory() as tmp:
                shell = Path(tmp, 'boundary.mesh')
                write_mesh(shell, p[:nb], np.empty((0, 8), dtype=int), q)
                result = subprocess.run([str(binary), '--input', str(shell), '--replay', str(path),
                    '--mirrors', '0', '--planes', 'none', '--cap', str(len(h)), '--threads', '1',
                    '--boundary-canonical', '0', '--output', str(Path(tmp, 'replay.jsonl'))],
                    capture_output=True, text=True, timeout=30)
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertIn(f'REPLAY_OK cells={len(h)} cell_orbits={len(h)}', result.stderr)


if __name__ == '__main__':
    unittest.main()
