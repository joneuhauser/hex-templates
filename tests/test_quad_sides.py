"""Complete small-disk searches, flip equivalence and independent exact checks."""
from collections import Counter
from fractions import Fraction as F
import hashlib
import json
import os
from pathlib import Path
import random
import subprocess
import sys
import tempfile
import unittest
import xml.etree.ElementTree as ET
import zipfile

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'tools'))
from build_geode_sides import canonical, flipped, harmonic, certify, published_side_key
from geode_boundary import assemble


class QuadSidesTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.assets = ROOT/'docs/assets/geode-sides'
        cls.catalog = json.loads((cls.assets/'catalog.json').read_text())

    def test_counts_and_exact_embeddings(self):
        patterns = self.catalog['patterns']
        self.assertEqual(self.catalog['max_quads'], 12)
        self.assertEqual(len(patterns), 432)
        self.assertEqual(Counter(len(p['quads']) for p in patterns),
                         {2: 1, 4: 1, 5: 3, 6: 1, 7: 6, 8: 11, 9: 22, 10: 55, 11: 91, 12: 241})
        self.assertEqual(sum(p['vertical_flip_class_size'] for p in patterns), 811)
        keys = set()
        for p in patterns:
            r = p['side_vertices']-2
            points = [tuple(F(x) for x in v) for v in p['square_points_exact']]
            self.assertEqual(certify(points, p['quads'], r), p['certificate'])
            key = canonical(p['quads'], p['boundary_vertices'])
            self.assertLessEqual(key, flipped(key, p['boundary_vertices']))
            keys.add((r, key))
            # Every downloadable patch must assemble into a usable Geode shell.
            vertices, faces = assemble(ROOT/'docs/assets/geode/meshes/G26-Mitchell.mesh', p)
            self.assertEqual(len(faces), 10+4*len(p['quads']))
            self.assertEqual(len(vertices), len(faces)+2)
        self.assertEqual(len(keys), 432)

    def test_relabeling_and_vertical_flips(self):
        rng = random.Random(9721)
        for p in self.catalog['patterns']:
            q, nb = p['quads'], p['boundary_vertices']
            key = canonical(q, nb)
            self.assertEqual(flipped(flipped(q, nb), nb), key)
            n = len(p['points'])
            interior = list(range(nb, n))
            rng.shuffle(interior)
            mapping = list(range(nb))+interior
            trial = [[mapping[v] for v in face] for face in q]
            rng.shuffle(trial)
            trial = [face[j:]+face[:j] for face, j in zip(trial, [rng.randrange(4) for _ in trial])]
            self.assertEqual(canonical(trial, nb), key)

    def test_published_side_present(self):
        found = [p for p in self.catalog['patterns'] if p['published_side_topology']]
        self.assertEqual(len(found), 1)
        self.assertEqual(len(found[0]['quads']), 5)
        self.assertEqual(found[0]['side_vertices'], 2)
        self.assertEqual(canonical(found[0]['quads'], 6), published_side_key())

    def test_both_enumerations(self):
        binary = Path(os.environ.get('QUAD_SOLVER', ROOT/'solver/build/quad_disk_search'))
        if not binary.is_file():
            self.skipTest('Build quad_disk_search or set QUAD_SOLVER')
        evidence = ROOT/'tests/fixtures/side-quadrangulations'
        with tempfile.TemporaryDirectory() as tmp:
            for r in range(11):
                outputs = []
                # Full orbit search; keep the live individual-cell regression small.
                for option in ([[], ['--no-orbits']] if r <= 6 else [[]]):
                    path = Path(tmp)/'quads.jsonl'
                    result = subprocess.run([str(binary), '8' if option else '12', str(r), str(path), *option],
                                            capture_output=True, text=True, check=True, timeout=30)
                    self.assertIn('QUAD_SEARCH_FINISHED', result.stderr)
                    outputs.append({canonical(a['quads'], a['boundary_vertices']) for a in map(json.loads, path.read_text().splitlines())})
                if len(outputs) == 2:
                    self.assertEqual({k for k in outputs[0] if len(k) <= 8}, outputs[1])
                recorded = [json.loads(s) for s in (evidence/f'r{r}-orbits.jsonl').read_text().splitlines()]
                self.assertEqual(outputs[0], {canonical(a['quads'], a['boundary_vertices']) for a in recorded})

    def test_recorded_crosscheck(self):
        evidence = ROOT/'tests/fixtures/side-quadrangulations'
        record = json.loads((self.assets/'enumeration.json').read_text())
        self.assertEqual(record['verification_max_quads'], 10)
        self.assertEqual(record['source_sha256'], hashlib.sha256((ROOT/'solver/quad_disk_search.cpp').read_bytes()).hexdigest())
        outputs = {}
        for run in record['runs']:
            r, orbits = run['vertical_interior_vertices'], run['symmetry_orbits']
            path = evidence/f'r{r}-{"orbits" if orbits else "individual"}.jsonl'
            self.assertEqual(run['results_sha256'], hashlib.sha256(path.read_bytes()).hexdigest())
            self.assertTrue(run['completed'])
            self.assertEqual(run['max_quads'], 12 if orbits else 10)
            self.assertTrue(run['log'].startswith(f'QUAD_SEARCH_FINISHED max_quads={run["max_quads"]} '))
            keys = {canonical(a['quads'], a['boundary_vertices']) for a in map(json.loads, path.read_text().splitlines())}
            self.assertEqual(len(keys), run['candidates'])
            outputs[r, orbits] = keys
        self.assertEqual(set(outputs), {(r, True) for r in range(11)} | {(r, False) for r in range(9)})
        for r in range(9):
            self.assertEqual({k for k in outputs[r, True] if len(k) <= 10}, outputs[r, False])

    def test_downloads_and_figures(self):
        with zipfile.ZipFile(self.assets/'patterns.zip') as archive:
            self.assertIsNone(archive.testzip())
            for p in self.catalog['patterns']:
                for extension in ['svg', 'json']:
                    name = f'{p["id"]}.{extension}'
                    self.assertEqual(archive.read(name), (self.assets/name).read_bytes())
                ET.parse(self.assets/f'{p["id"]}.svg')
        page = (ROOT/'docs/geode-sides.html').read_text()
        self.assertEqual(page.count('class="side-card"'), 432)
        self.assertEqual(page.count('loading="lazy"'), 432)
        positions = [page.index(f'id="quads-{q}"') for q in [2, 4, 5, 6, 7, 8, 9, 10, 11, 12]]
        self.assertEqual(positions, sorted(positions))


if __name__ == '__main__':
    unittest.main()
