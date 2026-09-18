#!/usr/bin/env python3
"""Run a published configuration, then embed its distinct hex connectivities."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import json
import os
from pathlib import Path
import subprocess
import tempfile
import time
import traceback
from types import SimpleNamespace

ROOT = Path(__file__).resolve().parents[1]
CAPS = {'pyramid': 36, 'geode-Q2-01': 34, 'geode-Q4-01': 38,
        'geode-Q5-03': 38, 'side-walls': 12}


def save(path, data):
    temporary = path.with_suffix('.tmp')
    temporary.write_text(json.dumps(data, indent=2)+'\n')
    temporary.replace(path)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('case', choices=CAPS)
    args = parser.parse_args()

    def positive(name, default=None):
        raw = os.environ.get(name, default)
        if raw is None or not str(raw).isdigit() or int(raw) < 1:
            parser.error(f'Set {name} to a positive integer')
        return int(raw)

    threads = positive('THREADS')
    side_walls = args.case == 'side-walls'
    key = 'QUAD_SOLVER' if side_walls else 'SOLVER'
    binary = 'quad_disk_search' if side_walls else 'mirror_search'
    solver = Path(os.environ.get(key, ROOT/'solver/build'/binary)).resolve()
    hexopt = Path(os.environ.get('HEXOPT_BINARY', ROOT/'solver/build/hexopt_fixed')).resolve()
    for path in ([solver] if side_walls else [solver, hexopt]):
        if not path.is_file() or not os.access(path, os.X_OK):
            parser.error(f'Missing executable: {path}. See README.md and tools/HEXOPT_NOTES.md')
    restarts = positive('HEXOPT_RESTARTS', '6')
    stall = positive('HEXOPT_STALL_LIMIT', '50000')
    iterations = positive('HEXOPT_MAX_ITERATIONS', '300000')
    if max(stall, iterations) > 10000000:
        parser.error('HexOpt iteration and stall limits must be at most 10000000')

    # Avoid BLAS thread pools multiplying the per-variant worker count.
    os.environ['OPENBLAS_NUM_THREADS'] = '1'
    if not side_walls:
        from hexopt_geode import collect, optimize

    if 'OUTPUT_DIR' in os.environ:
        output = Path(os.environ['OUTPUT_DIR']).resolve()
        output.mkdir(parents=True, exist_ok=False)
    else:
        parent = ROOT/'solver/runs'
        parent.mkdir(exist_ok=True)
        output = Path(tempfile.mkdtemp(prefix=f'{args.case}-quality-', dir=parent))
    search = output/'search'
    search.mkdir()
    status = dict(case=args.case, cap=CAPS[args.case], threads=threads,
                  stage='search', started_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()),
                  solver=str(solver), hexopt=None if side_walls else str(hexopt),
                  restarts=restarts, stall=stall, iterations=iterations)
    save(output/'status.json', status)
    print(f'Output: {output}', flush=True)
    try:
        with (search/'stdout.log').open('w') as stdout, (search/'search.log').open('w') as stderr:
            result = subprocess.run(['bash', str(ROOT/'solver/launch'/f'{args.case}.sh'), str(threads)],
                stdout=stdout, stderr=stderr, env=dict(os.environ, **{key: str(solver)},
                                                     SEARCH_OUTPUT_DIR=str(search)))
        status['search_exit_code'] = result.returncode
        if result.returncode:
            raise RuntimeError(f'Search exited with {result.returncode}; see {search}')
        status['search_complete'] = True
        if not side_walls:
            status['stage'] = 'hexopt'
            classes = collect([search])
            # Names describe the case rather than the generic search folder.
            for record in classes:
                record['side'] = args.case
                record['id'] = args.case+record['id'].removeprefix('search')
            status['classes'] = len(classes)
            save(output/'status.json', status)
            (output/'hexopt').mkdir()
            options = SimpleNamespace(output=output/'hexopt', binary=hexopt,
                reference=ROOT/'docs/assets/geode/meshes/G26-Mitchell.mesh',
                profile='pyramid' if args.case == 'pyramid' else 'geode',
                restarts=restarts, stall=stall, iterations=iterations)

            def embed(record):
                try:
                    return optimize(record, options)
                except Exception as error:
                    directory = options.output/record['id']
                    directory.mkdir(exist_ok=True)
                    failure = dict(id=record['id'], accepted=False,
                                   error=f'{type(error).__name__}: {error}')
                    save(directory/'result.json', failure)
                    (directory/'error.log').write_text(traceback.format_exc())
                    return failure

            results = []
            with ThreadPoolExecutor(max_workers=threads) as pool:
                for result in pool.map(embed, classes):
                    results.append(result)
                    save(output/'results.json', results)
                    status['processed'] = len(results)
                    save(output/'status.json', status)
            save(output/'results.json', results)
            status['accepted'] = sum(r['accepted'] for r in results)
            status['errors'] = sum('error' in r for r in results)
            if status['errors']:
                raise RuntimeError('Some variants failed to process; see results.json and per-variant error.log')
        status['stage'] = 'complete'
        save(output/'status.json', status)
        print(json.dumps(status), flush=True)
    except BaseException as error:
        status.update(stage='failed', error=f'{type(error).__name__}: {error}')
        save(output/'status.json', status)
        raise


if __name__ == '__main__':
    main()
