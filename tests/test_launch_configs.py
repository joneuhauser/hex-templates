"""Fixed-cap launchers admit the smallest published meshes and reproduce 2D results."""
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'tools'))
from mesh_tools import read_mesh, write_mesh
from prepare_geode import half_turn_actions
from geode_boundary import boundary_mapping
from search_geode_sides import diagonal_actions


class ConfigurationLaunchers(unittest.TestCase):
    def test_hex_caps_and_known_completions(self):
        catalog = json.loads((ROOT/'docs/assets/geode/catalog.json').read_text())
        selected = {}
        for entry in catalog['templates']:
            name = 'geode-'+(entry.get('side') or 'Q5-01')
            if name not in selected or entry['hexes'] < selected[name]['hexes']:
                selected[name] = entry
        selected['pyramid'] = dict(hexes=36, mesh='../solver/input/published36-symmetric.mesh')
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            shutil.copytree(ROOT/'solver/launch', tmp/'launch')
            (tmp/'input').symlink_to(ROOT/'solver/input', target_is_directory=True)
            capture = tmp/'args.json'
            stub = tmp/'capture'
            stub.write_text('#!/usr/bin/env python3\nimport json,os,sys\n'
                            'open(os.environ["CAPTURE"],"w").write(json.dumps(sys.argv[1:]))\n')
            stub.chmod(0o755)
            for name, entry in selected.items():
                with self.subTest(configuration=name):
                    subprocess.run(['bash', str(tmp/'launch'/f'{name}.sh'), '3'], check=True,
                        env=dict(os.environ, SOLVER=str(stub), CAPTURE=str(capture)), capture_output=True)
                    args = json.loads(capture.read_text())
                    options = dict(zip(args[::2], args[1::2]))
                    self.assertEqual(int(options['--threads']), 3)
                    self.assertEqual(int(options['--cap']), entry['hexes'])
                    self.assertEqual(options['--seconds'], '0')
                    bp, empty, bq = read_mesh(tmp/options['--input'])
                    self.assertEqual(len(empty), 0)
                    seed = ROOT/'docs'/entry['mesh']
                    extra = []
                    if name != 'pyramid':
                        p, h, q = read_mesh(seed)
                        nb = len(bp)
                        mapping = boundary_mapping(p[:nb], q, bp, bq)
                        self.assertIsNotNone(mapping)
                        perm = np.asarray(mapping+list(range(nb, len(p))))
                        actions = np.asarray(half_turn_actions(p, h, q) if name == 'geode-Q5-01'
                                             else diagonal_actions(p, q, 2))
                        remapped = np.empty_like(actions)
                        remapped[perm] = perm[actions]
                        seed = tmp/'seed.mesh'
                        write_mesh(seed, np.vstack([bp, p[nb:]]), perm[h], bq)
                        np.savetxt(tmp/'seed.actions', remapped, fmt='%d')
                        extra = ['--replay-actions', str(tmp/'seed.actions')]
                    replay = subprocess.run([str(ROOT/'solver/build/mirror_search'), *args,
                        '--replay', str(seed), *extra], cwd=tmp, capture_output=True, text=True, timeout=30)
                    self.assertEqual(replay.returncode, 0, replay.stderr)
                    self.assertIn(f'REPLAY_OK cells={entry["hexes"]}', replay.stderr)

    def test_parallel_side_wall_enumeration(self):
        record = json.loads((ROOT/'docs/assets/geode-sides/enumeration.json').read_text())
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            (tmp/'launch').mkdir()
            shutil.copyfile(ROOT/'solver/launch/side-walls.sh', tmp/'launch/side-walls.sh')
            run = subprocess.run(['bash', str(tmp/'launch/side-walls.sh'), '2'], check=True,
                env=dict(os.environ, QUAD_SOLVER=str(ROOT/'solver/build/quad_disk_search')),
                capture_output=True, text=True, timeout=300)
            output = Path(run.stdout.strip().removeprefix('Output: '))
            for case in record['runs']:
                mode = 'orbits' if case['symmetry_orbits'] else 'individual'
                path = output/f'r{case["vertical_interior_vertices"]}-{mode}.jsonl'
                self.assertEqual(hashlib.sha256(path.read_bytes()).hexdigest(), case['results_sha256'])


if __name__ == '__main__':
    unittest.main()
