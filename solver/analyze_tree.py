#!/usr/bin/env python3
"""Index a completed --tree-trace run; retain every state and summarize subtrees.

Usage: analyze_tree.py RUN_DIRECTORY [--state ID] [--csv]
State IDs belong to this run. Parent links reconstruct the complete ordered mesh
prefix; no isomorphism or hash quotient is applied to the trace.
"""
import argparse
import gzip
import json
from pathlib import Path
import re

import numpy as np


def build_index(path):
    files = sorted(path.glob('states-*.bin'), key=lambda f: int(f.stem.split('-')[1]))
    if not files:
        raise ValueError('No state trace files')
    counts = np.zeros(256, dtype=np.int64)
    for f in files:
        if f.stat().st_size % 48:
            raise ValueError(f'Truncated record: {f}')
        counts[int(f.stem.split('-')[1])] = f.stat().st_size // 48
    offsets = np.r_[0, np.cumsum(counts)[:-1]]
    n = int(counts.sum())
    nodes = np.lib.format.open_memmap(path/'nodes.npy', mode='w+', dtype='<u8', shape=(n, 6))
    for f in files:
        worker = int(f.stem.split('-')[1])
        rows = np.fromfile(f, dtype='<u8').reshape(-1, 6)
        serial = (rows[:, 0] >> 8).astype(np.int64)
        if not np.all((rows[:, 0] & 255) == worker):
            raise ValueError('Worker ID mismatch')
        if not np.array_equal(np.sort(serial), np.arange(1, len(rows)+1)):
            raise ValueError('Missing or duplicate state ID: run may be incomplete')
        nodes[offsets[worker]+serial-1] = rows
    nodes.flush()
    parents = np.full(n, -1, dtype=np.int64)
    nonroot = nodes[:, 1] != 0
    pid = nodes[nonroot, 1]
    workers = (pid & 255).astype(np.int64)
    serial = (pid >> 8).astype(np.int64)
    if np.any(serial < 1) or np.any(serial > counts[workers]):
        raise ValueError('Missing parent')
    parents[nonroot] = offsets[workers]+serial-1
    if np.count_nonzero(~nonroot) != 1:
        raise ValueError('Trace must have exactly one root')
    depth = (nodes[:, 2] & 255).astype(np.uint8)
    if not np.all(depth[nonroot] == depth[parents[nonroot]]+1):
        raise ValueError('Parent does not precede child by one cell')
    children = np.bincount(parents[nonroot], minlength=n)
    outcome = ((nodes[:, 2] >> 32) & 255).astype(np.uint8)
    reason = ((nodes[:, 2] >> 40) & 255).astype(np.uint8)
    if np.any(outcome == 6):
        raise ValueError('Interrupted expansion: not a complete enumeration')
    terminal_names = ['quick_packing', 'greedy_packing', 'cover_proof',
                      'cell_or_face_budget', 'local_expansion_exhausted',
                      'candidate', 'topology_rejected', 'two_face_packing']
    terminal = np.zeros((n, len(terminal_names)), dtype=np.uint32)
    for k in range(3):
        terminal[:, k] = (outcome == 1) & (reason == k+1)
    terminal[:, 3] = outcome == 2
    terminal[:, 4] = (outcome == 3) & (children == 0)
    terminal[:, 5] = outcome == 4
    terminal[:, 6] = outcome == 5
    terminal[:, 7] = (outcome == 1) & (reason == 5)
    if np.any((outcome == 1) & ~np.isin(reason, [1, 2, 3, 5])):
        raise ValueError('Unclassified cover rejection')
    subtree = np.ones(n, dtype=np.uint64)
    work = np.array(nodes[:, 4])
    branches = np.array(nodes[:, 5])
    for d in range(int(depth.max()), 0, -1):
        ix = np.flatnonzero(depth == d)
        target = parents[ix]
        np.add.at(subtree, target, subtree[ix])
        np.add.at(work, target, work[ix])
        np.add.at(branches, target, branches[ix])
        np.add.at(terminal, target, terminal[ix])
    root = int(np.flatnonzero(~nonroot)[0])
    if subtree[root] != n:
        raise ValueError('Disconnected trace')
    for name, value in [('parents', parents), ('children', children),
                        ('subtree_states', subtree), ('subtree_cover_work', work),
                        ('subtree_choose_calls', branches), ('terminal_counts', terminal)]:
        np.save(path/(name+'.npy'), value)
    (path/'index.json').write_text(json.dumps(dict(states=n, offsets=offsets.tolist(),
        worker_counts=counts.tolist(), root_id=int(nodes[root, 0]),
        terminal_names=terminal_names), indent=2)+'\n')
    return nodes, parents, children, subtree, work, branches, terminal


def load_index(path):
    return tuple(np.load(path/(name+'.npy'), mmap_mode='r') for name in
        ['nodes', 'parents', 'children', 'subtree_states', 'subtree_cover_work',
         'subtree_choose_calls', 'terminal_counts'])


def describe(path, data, ix, cells=True, local_details=False):
    nodes, parents, children, subtree, work, branches, terminal = data
    fmt = json.loads((path/'format.json').read_text())
    index = json.loads((path/'index.json').read_text())
    row = list(map(int, nodes[ix]))
    result = dict(id=row[0], parent_id=row[1], placed_cells=row[2] & 255,
        vertices=row[2] >> 8 & 255, front_faces=row[2] >> 16 & 65535,
        outcome=fmt['outcomes'][str(row[2] >> 32 & 255)],
        cover_reason=fmt['cover_reasons'][str(row[2] >> 40 & 255)],
        children=int(children[ix]), subtree_states=int(subtree[ix]),
        fraction_of_run=float(subtree[ix]/len(nodes)), local_cover_work=row[4],
        subtree_cover_work=int(work[ix]), local_choose_calls=row[5],
        subtree_choose_calls=int(branches[ix]),
        terminal_reasons=dict(zip(index['terminal_names'], map(int, terminal[ix]))))
    if local_details and result['outcome'] in ['exhausted_expansion', 'interrupted']:
        fields = fmt['expansion_words']
        records = np.fromfile(path/f'expansions-{row[0] & 255}.bin', dtype='<u8').reshape(-1, len(fields))
        found = records[records[:, 0] == row[0]]
        if len(found) != 1:
            raise ValueError('Missing expansion detail')
        result['local_expansion_counters'] = {name:int(value) for name, value in zip(fields[1:], found[0, 1:]) if value}
    if cells:
        mesh = []
        cursor = ix
        while parents[cursor] >= 0:
            packed = int(nodes[cursor, 3])
            mesh.append([(packed >> (8*j) & 255)-1 for j in range(8)])
            cursor = int(parents[cursor])
        result['hexes'] = mesh[::-1]
    return result


def summarize(path, data):
    nodes, parents, children, subtree, work, branches, terminal = data
    index = json.loads((path/'index.json').read_text())
    fmt = json.loads((path/'format.json').read_text())
    root = int(np.flatnonzero(parents < 0)[0])
    depth = (nodes[:, 2] & 255).astype(np.uint8)
    by_depth = []
    for d in range(int(depth.max())+1):
        ix = np.flatnonzero(depth == d)
        if not len(ix):
            continue
        top = ix[np.argsort(subtree[ix])[-5:][::-1]]
        by_depth.append(dict(placed_cells=d, states=len(ix),
            expanded=int(np.count_nonzero(children[ix])),
            local_choose_calls=int(nodes[ix, 5].sum()),
            local_cover_work=int(nodes[ix, 4].sum()),
            top_subtrees=[describe(path, data, int(i)) for i in top]))
    fields = fmt['expansion_words']
    totals = np.zeros(len(fields)-1, dtype=np.uint64)
    expansion_states = 0
    for file in path.glob('expansions-*.bin'):
        records = np.fromfile(file, dtype='<u8').reshape(-1, len(fields))
        expansion_states += len(records)
        totals += records[:, 1:].sum(axis=0)
    local = dict(zip(fields[1:], map(int, totals)))
    progress = [x for x in (path/'search.log').read_text().splitlines() if x.startswith('PROGRESS')][-1]
    counters = dict(re.findall(r'(\w+)=([^ ]+)', progress))
    if int(counters['states']) != len(nodes) or int(counters['branches']) != int(branches[root]):
        raise ValueError('Trace totals disagree with solver progress')
    if int(counters['cover_work']) != int(work[root]):
        raise ValueError('Cover work not attributed exactly once')
    comparisons = {'trials':'orbit_trials', 'partial_cover_calls':'partial_cover_calls',
        'partial_cover_rejections':'partial_cover_rejections', 'early_domain':'early_domain',
        'early_geometry':'early_geometry', 'early_budget':'early_budget',
        'early_parity':'early_parity', 'reject_location':'reject_location',
        'canonical_rejections':'canonical_rejections', 'insert_placed_cells':'reject_placed_cells'}
    for field, counter in comparisons.items():
        if local[field] != int(counters[counter]):
            raise ValueError(f'Incomplete attribution: {field}: {local[field]} != {counters[counter]}')
    result = dict(states=len(nodes), expanded_records=expansion_states,
        complete='SEARCH_FINISHED' in (path/'search.log').read_text(),
        root=describe(path, data, root), by_depth=by_depth,
        expansion_local_counters=local, final_progress=counters,
        interpretation='States are recursive visits, not distinct isomorphism classes. '
        'Cover work counts auxiliary search nodes. Option-removal counters and '
        'branch-rejection counters have different units and must not be added as leaf counts.')
    (path/'summary.json').write_text(json.dumps(result, indent=2)+'\n')
    return result


def export_csv(path, data):
    nodes, parents, children, subtree, work, branches, terminal = data
    with gzip.open(path/'all-states.csv.gz', 'wt', compresslevel=1) as out:
        out.write('id,parent_id,placed_cells,vertices,front_faces,outcome,cover_reason,'
                  'children,subtree_states,local_choose_calls,subtree_choose_calls,'
                  'local_cover_work,subtree_cover_work\n')
        for start in range(0, len(nodes), 100000):
            stop = min(start+100000, len(nodes))
            a = nodes[start:stop]
            rows = np.column_stack([a[:, 0], a[:, 1], a[:, 2] & 255,
                a[:, 2] >> 8 & 255, a[:, 2] >> 16 & 65535,
                a[:, 2] >> 32 & 255, a[:, 2] >> 40 & 255, children[start:stop].astype(np.uint64),
                subtree[start:stop], a[:, 5], branches[start:stop], a[:, 4], work[start:stop]])
            np.savetxt(out, rows, fmt='%u', delimiter=',')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('directory', type=Path)
    parser.add_argument('--state', type=int)
    parser.add_argument('--csv', action='store_true')
    args = parser.parse_args()
    path = args.directory
    if 'SEARCH_FINISHED' not in (path/'search.log').read_text():
        raise ValueError('Wait for a completed enumeration before indexing')
    if json.loads((path/'format.json').read_text())['version'] != 2:
        raise ValueError('Unsupported trace format version')
    data = load_index(path) if (path/'index.json').exists() else build_index(path)
    if args.state is not None:
        info = json.loads((path/'index.json').read_text())
        worker, serial = args.state & 255, args.state >> 8
        if not 0 < serial <= info['worker_counts'][worker]:
            raise ValueError('Unknown state ID')
        print(json.dumps(describe(path, data, info['offsets'][worker]+serial-1, local_details=True), indent=2))
        return
    summary = summarize(path, data)
    if args.csv:
        export_csv(path, data)
    print(json.dumps({k:summary[k] for k in ['states', 'complete', 'root']}, indent=2))


if __name__ == '__main__':
    main()
