"""Upper transition: website identity, mirror scope, periodic caps and replay."""
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'tools'))
import geode_top as g


class GeodeTopTests(unittest.TestCase):
    def test_website_connectivity_and_mirror_orbits(self):
        proof, target = g.website_report()
        self.assertEqual(len(proof['vertex_bijection']), 60)
        self.assertEqual(target['hexes'], 26)
        p, h, _ = g.read_mesh(g.FIXTURE)
        action = g.actions(p, 1)[:, 1]
        cells = {tuple(sorted(c)) for c in h}
        self.assertEqual(cells, {tuple(sorted(c)) for c in action[h]})
        fixed = sum(tuple(sorted(c)) == tuple(sorted(action[c])) for c in h)
        self.assertEqual(fixed, 4)
        self.assertEqual(fixed+(len(h)-fixed)//2, 15)

    def test_periodic_boundary_and_alternative_cut(self):
        p, h, q = g.read_mesh(g.BOUNDARY)
        self.assertEqual(len(h), 0)  # Search input contains no cells.
        report = g.boundary_report(p, q, 1)
        self.assertEqual(len(report['periodic_side_pairs']), 18)
        with self.assertRaises(ValueError):
            g.boundary_report(p, q, 2)
        pp, qq = g.symmetric_boundary(p, q)
        report = g.boundary_report(pp, qq, 2)
        self.assertEqual((len(pp), len(qq)), (38, 36))
        self.assertEqual(len(report['periodic_side_pairs']), 11)

        def caps(points, faces):
            return {tuple(sorted(tuple(np.round(v, 12)) for v in points[f]))
                    for f in faces if any(np.max(abs(points[f, 2]-z)) < 1e-12
                        for z in (points[:, 2].min(), points[:, 2].max()))}

        self.assertEqual(len(caps(p, q)), 14)
        self.assertEqual(caps(p, q), caps(pp, qq))
        other = pp[:, [1, 0, 2]]*[-1, 1, 1]
        g.boundary_report(other, qq, 1)
        broken = pp.copy()
        broken[0, 0] += .01
        with self.assertRaises(ValueError):
            g.boundary_report(broken, qq, 2)

    def test_native_replay_and_cap(self):
        solver = os.environ.get('GEODE_SOLVER', str(ROOT/'solver/build/mirror_search'))
        _, target = g.website_report()
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary)/'candidate.jsonl'
            for cap, code in [(26, 0), (25, 3)]:
                result = subprocess.run([solver, '--input', str(g.BOUNDARY),
                    '--replay', str(g.FIXTURE), '--cap', str(cap), '--threads', '1',
                    '--mirrors', '1', '--planes', 'axial', '--boundary-canonical', '0',
                    '--output', str(output)], capture_output=True, text=True, timeout=30)
                self.assertEqual(result.returncode, code, result.stderr)
                if cap == 26:
                    self.assertIn('REPLAY_OK cells=26 cell_orbits=15', result.stderr)
                    candidate = json.loads(output.read_text())
                    self.assertIsNotNone(g.isomorphism(candidate, target, 52))
                else:
                    self.assertIn('REPLAY_REJECTED', result.stderr)


if __name__ == '__main__':
    unittest.main()
