"""Recheck the published assets and replay them through each search symmetry mode."""
from pathlib import Path
import argparse, hashlib, json, subprocess, sys, tempfile
import numpy as np
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'tools'))
from build_catalog import symmetries, TRANSFORMS
from mesh_tools import read_mesh, write_mesh
from validate_mesh import validate
from quality_metrics import dense

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--skip-replay', action='store_true')
    args = parser.parse_args()
    catalog = json.loads((ROOT / 'docs/assets/catalog.json').read_text())
    rows = catalog['templates']
    assert len(rows) == catalog['mesh_count']
    assert sum(r['hexes'] <= 48 for r in rows) == catalog['study_mesh_count']
    assert len({r['id'] for r in rows}) == len(rows)
    replays = 0
    with tempfile.TemporaryDirectory(prefix='hex-atlas-test-') as tmp:
        for r in rows:
            mesh = ROOT / 'docs' / r['mesh']
            assert hashlib.sha256(mesh.read_bytes()).hexdigest() == r['sha256'], r['id']
            p, h, q = read_mesh(mesh,infer_boundary=True)
            assert (len(p), len(h), len(q)) == (r['vertices'], r['hexes'], 16)
            result = validate(mesh,r.get('validation_profile','prescribed'))
            assert result['accepted'] and result['all_jacobians_positive'], r['id']
            stored = json.loads((ROOT / 'docs' / r['certificate']).read_text())
            assert stored['accepted'] and stored['all_jacobians_positive']
            assert result['exact_bernstein_bounds'] == stored['exact_bernstein_bounds'], r['id']
            measured = symmetries(p @ np.array(r.get('display_matrix',np.eye(3))).T, h)
            assert measured == r['symmetry'], r['id']
            if len(h) <= 48 and not r.get('validation_profile'): assert measured['group'] == 'C2v'
            if r.get('validation_profile') == 'published44':
                assert r['sha256'] == r['provenance']['source_sha256'] == '8c61a75d882c0935e1e1412cb7deb582cf7ed3e06d774c3e02a7baf0240eccb1'
            assert r['metrics']['min_scaled_jacobian'] > 0
            metrics = dense(p, h, True)
            for key, value in metrics.items():
                assert np.allclose(value, r['metrics'][key], rtol=1e-10, atol=1e-13), (r['id'], key)
            if not args.skip_replay and len(h) <= 48 and not r.get('validation_profile'):
                plane = 'axial' if 'mx' in [a['id'] for a in measured['actions']] else 'diagonal'
                for mirrors in (0, 1, 2):
                    command = [str(ROOT / 'solver/build/mirror_search'), '--input', str(ROOT / 'solver/input/pyramid.mesh'), '--cap', str(len(h)), '--threads', '1', '--mirrors', str(mirrors), '--planes', plane, '--replay', str(mesh), '--output', str(Path(tmp) / 'candidate.jsonl')]
                    if mirrors == 1:
                        command += ['--face-order', '5', '--partial-cover', '6', '--partial-work', '0']
                    # A seed may be a noncanonical boundary image. At least one
                    # admissible image must survive with canonical pruning enabled.
                    reflection = np.diag([-1,1,1]) if plane == 'axial' else np.array([[0,1,0],[1,0,0],[0,0,1]])
                    for _, _, raw_matrix in TRANSFORMS:
                        matrix = np.array(raw_matrix)
                        if mirrors and not np.array_equal(matrix @ reflection, reflection @ matrix): continue
                        transformed = p @ matrix.T
                        permutation = np.arange(len(p))
                        distance = np.linalg.norm(transformed[:18,None,:]-p[None,:18,:],axis=2)
                        permutation[:18] = distance.argmin(1)
                        assert len(set(permutation[:18])) == 18
                        points = np.empty_like(p); points[permutation] = transformed
                        cells = permutation[h].copy()
                        if np.linalg.det(matrix) < 0:
                            cells[:,[1,3]] = cells[:,[3,1]]; cells[:,[5,7]] = cells[:,[7,5]]
                        seed = Path(tmp) / 'image.mesh'; write_mesh(seed, points, cells, q)
                        command[command.index('--replay')+1] = str(seed)
                        run = subprocess.run(command, text=True, capture_output=True, timeout=120)
                        if run.returncode == 0 and 'SEARCH_FINISHED candidates=1' in run.stderr: break
                        assert run.returncode == 3, (r['id'], mirrors, run.stderr)
                    else:
                        raise AssertionError((r['id'], mirrors, 'All admissible images rejected', run.stderr))
                    replays += 1
            print('PASS', r['id'], flush=True)
    print(f'TEMPLATES_OK {len(rows)} exact validations; {replays} seeded replays')

if __name__ == '__main__':
    main()
