"""Independent discovery, replay and orientation regressions for C2 search."""
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT/'tools'), str(ROOT/'solver')]
from compare_meshes import equivalent, oriented_cells
from mesh_tools import read_mesh, write_mesh
from prepare_geode import half_turn_actions
from geode_boundary import boundary_mapping
from search_geode_sides import placements
from geode_boundary import assemble

SOLVER = Path(os.environ.get('GEODE_SOLVER', ROOT/'solver/build/mirror_search'))


class HalfTurnTests(unittest.TestCase):
    def run_search(self, directory, name, boundary, cap, flags=(), code=0):
        output = directory/(name+'.jsonl')
        result = subprocess.run([str(SOLVER), '--input', str(boundary), '--cap', str(cap),
            '--threads', '1', '--seconds', '30', '--boundary-canonical', '0',
            '--output', str(output), *map(str, flags)], capture_output=True, text=True, timeout=40)
        self.assertEqual(result.returncode, code, result.stderr)
        return [json.loads(line) for line in output.read_text().splitlines()], result.stderr

    def assert_rotation(self, candidate):
        action = dict(enumerate(row[1] for row in candidate['vertex_actions']))
        self.assertTrue(all(action[action[v]] == v for v in action))
        self.assertEqual(oriented_cells(candidate, action), oriented_cells(candidate, dict(enumerate(range(candidate['vertices'])))))

    def test_unseeded_cube_with_and_without_pruning(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            flags = ['--half-turn', '1']
            fast, _ = self.run_search(tmp, 'fast', ROOT/'solver/input/cube.mesh', 7, flags)
            slow, _ = self.run_search(tmp, 'slow', ROOT/'solver/input/cube.mesh', 7,
                flags+['--early-overlap', '0', '--overlap-domains', '0', '--face-cover', '0',
                       '--component-bound', '0', '--early-budget', '0', '--leader-bound', '0'])
            zero, _ = self.run_search(tmp, 'zero', ROOT/'solver/input/cube.mesh', 7, ['--mirrors', '0'])
            self.assertEqual(sorted(c['hexes'] for c in fast), [1, 7])
            self.assertEqual(len(fast), len(slow))
            self.assertEqual(len(fast), len(zero))
            for candidate in fast:
                self.assert_rotation(candidate)
                self.assertTrue(any(equivalent(candidate, other, 8) for other in slow))
                self.assertTrue(any(equivalent(candidate, other, 8) for other in zero))

    def test_half_turn_without_diagonal_mirrors(self):
        p, _, q = read_mesh(ROOT/'solver/input/cube.mesh')
        p[:, :2] = p[:, :2] @ np.array([[1, .2], [.3, 1]])
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            write_mesh(tmp/'skew.mesh', p, np.empty((0, 8), dtype=int), q)
            found, _ = self.run_search(tmp, 'rotation', tmp/'skew.mesh', 1, ['--half-turn', '1'])
            self.assertEqual(len(found), 1)
            self.assert_rotation(found[0])
            # A diagonal reflection is not even a symmetry of this boundary.
            result = subprocess.run([str(SOLVER), '--input', str(tmp/'skew.mesh'),
                '--mirrors', '1', '--cap', '1', '--output', str(tmp/'mirror.jsonl')],
                capture_output=True, text=True, timeout=10)
            self.assertEqual(result.returncode, 2)

    def test_published_geode_replay_and_cap(self):
        ref = ROOT/'docs/assets/geode/meshes/G26-Mitchell.mesh'
        p, h, q = read_mesh(ref)
        action = half_turn_actions(p, h, q)
        self.assertEqual(sum(v == w for v, w in action), 4)
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            np.savetxt(tmp/'rotation.actions', action, fmt='%d')
            flags = ['--half-turn', '1', '--replay', ref, '--replay-actions', tmp/'rotation.actions']
            boundary = ROOT/'solver/input/geode-boundary.mesh'
            for label, extra in [('default', []), ('free-side', ['--corner-geometry', '0', '--geometry-domains', '0', '--location-domains', '0']),
                                 ('no-cover', ['--face-cover', '0', '--early-overlap', '0', '--overlap-domains', '0'])]:
                found, log = self.run_search(tmp, label, boundary, 26, flags+extra)
                self.assertIn('REPLAY_OK cells=26 cell_orbits=13', log)
                self.assert_rotation(found[0])
                self.assertTrue(equivalent(found[0], dict(vertices=48, hexes=26, cells=h.tolist()), 32))
            found, log = self.run_search(tmp, 'too-small', boundary, 25, flags, code=3)
            self.assertEqual(found, [])
            self.assertIn('REPLAY_REJECTED', log)
            for face in range(len(q)):
                found, _ = self.run_search(tmp, f'root-{face}', boundary, 26,
                                           flags+['--root-face', face])
                self.assertEqual(len(found), 1)

    def test_catalog_placement(self):
        ref = ROOT/'docs/assets/geode/meshes/G26-Mitchell.mesh'
        p, h, q = read_mesh(ref)
        catalog = json.loads((ROOT/'docs/assets/geode-sides/catalog.json').read_text())
        for name, side, flipped in placements(catalog, 5):
            if name.startswith('Q5-01'):
                bp, bq = assemble(ref, side)
                mapping = boundary_mapping(p[:32], q, bp, bq)
                self.assertEqual(mapping is None, flipped)

    def test_existing_mirror_replays(self):
        with tempfile.TemporaryDirectory() as tmp:
            for mirrors in (1, 2):
                found, _ = self.run_search(Path(tmp), f'm{mirrors}', ROOT/'solver/input/pyramid.mesh', 36,
                    ['--mirrors', mirrors, '--replay', ROOT/'solver/input/published36-symmetric.mesh'])
                self.assertEqual(len(found), 1)


if __name__ == '__main__':
    unittest.main()
