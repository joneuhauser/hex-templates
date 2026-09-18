"""Check compact tracing against the existing full text trace and a reference run."""
from collections import Counter
import json
from pathlib import Path
import re
import subprocess
import tempfile

import numpy as np

root = Path(__file__).resolve().parent
(root/'experiments').mkdir(exist_ok=True)
directory = Path(tempfile.mkdtemp(prefix='tree-validation-', dir=root/'experiments'))
for fixture, cap in [('pyramid', 16), ('pyramid', 20), ('cube', 7)]:
    path = directory/f'{fixture}-{cap}'
    path.mkdir()
    common = [str(root/'build/mirror_search'), '--input', str(root/'input'/f'{fixture}.mesh'),
              '--mirrors', '0', '--cap', str(cap), '--threads', '4', '--cover-work', '10000']
    if fixture == 'cube':
        common += ['--interior', '8', '--face-order', '3']
    command = common + ['--tree-trace', str(path), '--state-trace', str(path/'states.txt'),
                        '--output', str(path/'candidates.jsonl')]
    run = subprocess.run(command, capture_output=True, text=True)
    (path/'search.log').write_text(run.stderr)
    assert run.returncode == 0, run.stderr
    assert json.loads((path/'format.json').read_text()) == json.loads((root/'tree_trace_format.json').read_text())
    records = {}
    for file in path.glob('states-*.bin'):
        for row in np.fromfile(file, dtype='<u8').reshape(-1, 6):
            records[int(row[0])] = list(map(int, row))
    actual = Counter()
    for row in records.values():
        cells, cursor = [], row
        while cursor[1]:
            cells.append(tuple((cursor[3] >> (8*j) & 255)-1 for j in range(8)))
            cursor = records[cursor[1]]
        actual[(row[2] >> 8 & 255, tuple(cells[::-1]))] += 1
    expected = Counter()
    for line in (path/'states.txt').read_text().splitlines():
        vertices, hexes = line.split(' h')
        values = list(map(int, hexes.split()))
        expected[(len(vertices.split()[1:])//4,
                  tuple(tuple(values[i:i+8]) for i in range(0, len(values), 8)))] += 1
    assert actual == expected, 'Compact parent links do not reconstruct the text trace'
    if fixture == 'pyramid' and cap == 20:
        assert any((row[2] >> 40 & 255) == 5 for row in records.values()), 'Two-face packing trace reason was not exercised'
    analyzed = subprocess.run(['python3', str(root/'analyze_tree.py'), str(path), '--csv'],
                              capture_output=True, text=True)
    assert analyzed.returncode == 0, analyzed.stderr
    reference = subprocess.run(common + ['--output', str(path/'reference.jsonl')],
                               capture_output=True, text=True)
    assert reference.returncode == 0, reference.stderr
    assert sorted((path/'reference.jsonl').read_text().splitlines()) == sorted((path/'candidates.jsonl').read_text().splitlines())
    def counters(log):
        line = next(x for x in reversed(log.splitlines()) if x.startswith('PROGRESS'))
        return dict(re.findall(r'(\w+)=([^ ]+)', line))
    observed, baseline = counters(run.stderr), counters(reference.stderr)
    for key in ['states', 'branches', 'orbit_trials', 'candidates', 'cover_work',
                'partial_cover_rejections', 'reject_placed_cells', 'canonical_rejections']:
        assert observed[key] == baseline[key], (fixture, key, observed[key], baseline[key])
    summary = json.loads((path/'summary.json').read_text())
    queried = subprocess.run(['python3', str(root/'analyze_tree.py'), str(path), '--state',
                             str(summary['root']['id'])], capture_output=True, text=True)
    assert queried.returncode == 0, queried.stderr
    assert json.loads(queried.stdout)['local_expansion_counters']['children'] == summary['root']['children']
print(f'TREE_TEST_OK all prefixes reconstructed; counters, candidates, subtree totals and queries agree. Logs: {directory}')
