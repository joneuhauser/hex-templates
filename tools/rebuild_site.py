#!/usr/bin/env python3
"""Rebuild the final site from its published meshes and compact source records."""
import os
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]


def main():
    env = dict(os.environ, OPENBLAS_NUM_THREADS='1')
    commands = [
        ['build_catalog.py'],
        ['build_geode_catalog.py'],
        ['build_periodic_catalog.py'],
        ['build_geode_sides.py', '--search-results',
         str(ROOT/'tests/fixtures/side-quadrangulations')],
        ['make_figures.py'],
        ['build_home.py'],
    ]
    for script, *args in commands:
        subprocess.run([sys.executable, str(ROOT/'tools'/script), *args],
                       cwd=ROOT, env=env, check=True)


if __name__ == '__main__':
    main()
