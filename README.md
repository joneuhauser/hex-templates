# Hexahedral meshing: Template library

> [!CAUTION]
> This is an experimental library that I created by discussion with GPT-6 Astra for two days. The generated meshes check out and I'm using them in production for fluid dynamics simulations, but I'm not claiming that they are optimal. Use at your own risk. In particular I do not assert that there are no meshes with smaller element count for the same boundaries.

> [!TIP]
> Check out the ![documentation and mesh library](http://joneuhauser.github.io/hex-templates).  

A static mesh library with a connectivity search solver, validation tests,
and geometry optimization tools:

- [Schneiders’ pyramid](docs/pyramid.html): 56 certified meshes.
- [Mitchell’s Geode](docs/geode.html): 29 certified meshes, including the published
  26-cell reference and 28 alternatives.
- [Geode side walls](docs/geode-sides.html): 432 symmetric quadrangulations through
  12 quads, identifying vertical flips.
- [Periodic grid rotation](docs/periodic.html): seven certified periodic
  connectivity classes, with 38–110 cells per tile.

Every displayed mesh has an exact whole-cell positive-Jacobian certificate.
Sampled quality measurements are separate from that certification. Search
completion and symmetry restrictions are explained on the site; no unrestricted
optimality claim is made.

## Preview

```sh
python3 -m http.server 8085 --directory docs
```

Open <http://localhost:8085>. No frontend build is required. The interactive
viewers load Three.js from jsDelivr and require JavaScript and WebGL; static
measurements and downloads remain available independently.

## Rebuild the current site

```sh
python3 -m pip install -r tools/requirements.txt
OPENBLAS_NUM_THREADS=1 MPLCONFIGDIR=/tmp/hex-atlas-matplotlib python3 tools/rebuild_site.py
```

The checked-in catalogs and mesh files are the authoritative site inputs.
The rebuild validates Geode and periodic geometry, regenerates exports and
archives, refreshes symmetry information, rebuilds the side-wall catalog from
its enumeration fixtures, and regenerates figures and navigation. All source
meshes and enumeration fixtures are included.

`tools/data/` contains the completed pyramid enumeration summary.
`tests/fixtures/` contains the published connectivity comparison and side-wall
enumeration cross-checks.

## Solver and tests

Requirements: a C++17 compiler, CMake, OpenMP, Boost.Multiprecision headers and
Bash, plus the Python dependencies above. Debian/Ubuntu packages are
`build-essential cmake libboost-dev`.

```sh
cmake -S solver -B solver/build -DCMAKE_BUILD_TYPE=Release -DNATIVE_ARCH=OFF
cmake --build solver/build --parallel 2
THREADS=4 BUILD_JOBS=2 bash tests/run.sh
```

The tests cover the search constraints, published-connectivity replays,
independent topology and exact geometry checks, side-wall enumeration, periodic
identification, relative website links, and download consistency.

```sh
# Example pyramid search; MIRRORS may be 0, 1 or 2.
MIRRORS=2 CAP=36 THREADS=4 PLANES=diagonal bash solver/run.sh foreground
```

`mirror_search` supports a separate `--half-turn 1` constraint for Geodes.
`tools/search_geode_sides.py` assembles matching side walls from the
catalog. The root launchers below run search, HexOpt embedding, and validation.
`tools/search_geode_cavities.py` implements the local replacements used for the
28-cell alternatives. Their command-line help describes the available options.
New run output belongs in `solver/runs/` or outside the repository.

Pyramid, Geode and side-wall searches have root launchers. The hex launchers complete an
unseeded search, deduplicate its connectivities, then run HexOpt and exact geometry
validation on every variant. Build the solver above and the
[HexOpt adapter](tools/HEXOPT_NOTES.md) before running them:

```sh
THREADS=8 bash run-pyramid.sh
THREADS=8 bash run-geode.sh
THREADS=8 bash run-side-walls.sh
```

`THREADS` is required: it controls solver threads, then concurrent HexOpt variants
(one optimizer thread per variant). The 2D side-wall finder uses it for concurrent
boundary cases and needs no HexOpt. `THREADS=8 bash run.sh` runs these five cases
sequentially. `run-geode.sh` runs Q2-01 (cap 34), Q4-01 (cap 38), and Q5-03
(cap 38), with HexOpt after each search.

Each invocation prints a fresh folder under `solver/runs/`. To choose a new output
folder explicitly:

```sh
THREADS=8 OUTPUT_DIR=/tmp/geode-results bash run-geode.sh
```

Existing case folders are never overwritten. `run-geode.sh` creates `Q2-01/`,
`Q4-01/`, and `Q5-03/` inside its output folder. With `run.sh`, `OUTPUT_DIR` is a
parent containing `pyramid/`, `geode/`, and `side-walls/`, with the three Geode cases
inside `geode/`. Hex case folders contain:

- `search/`: boundary, candidate connectivities, solver stdout and search log.
- `hexopt/<variant>/`: initialization, optimizer inputs, coordinates and logs,
  plus `accepted.mesh`, `accepted.vtu`, validation and quality reports when certified.
- `status.json` and `results.json`: progress and per-variant outcomes, including
  failed embedding attempts. Failure to embed is not an impossibility proof.

The side-wall folder contains enumeration JSONL files and logs under `search/`,
which can be passed to `tools/build_geode_sides.py --search-results DIRECTORY`.
Optional environment variables are `HEXOPT_BINARY` (default
`solver/build/hexopt_fixed`), `HEXOPT_RESTARTS` (6), `HEXOPT_STALL_LIMIT` (50000),
and `HEXOPT_MAX_ITERATIONS` (300000). Search has no time limit; individual HexOpt
attempts have a 120-second limit and timed-out attempts are recorded before retrying.

The lower-level, search-only launchers take the thread count as their argument:

```sh
bash solver/launch/pyramid.sh 8
bash solver/launch/geode-Q4-01.sh 8
bash solver/launch/side-walls.sh 8
```

| Script in `solver/launch/` | Fixed cap | Search constraint | Included in default |
| --- | ---: | --- | --- |
| `pyramid.sh` | 36 hexes | Two diagonal mirrors | Yes |
| `geode-Q2-01.sh` | 34 hexes | Two diagonal mirrors | Yes |
| `geode-Q4-01.sh` | 38 hexes | Two diagonal mirrors | Yes |
| `geode-Q5-03.sh` | 38 hexes | Two diagonal mirrors | Yes |
| `side-walls.sh` | 12 quads | Left/right symmetry; cross-check through 10 | Yes |
| `geode-Q5-01.sh` | 26 hexes | 180° rotation | No |
| `geode-Q6-01.sh` | 42 hexes | Two diagonal mirrors | No |
| `geode-Q7-06.sh` | 42 hexes | Two diagonal mirrors | No |
| `geode-Q7-06-flipped.sh` | 44 hexes | Two diagonal mirrors | No |

Hex caps equal the smallest displayed mesh for each boundary configuration.
These scripts launch unseeded searches with no time limit and print a fresh
output directory. They require the binaries built above. For the serial 2D
enumerator, the thread argument limits concurrent boundary cases. Its output
can be passed to `tools/build_geode_sides.py --search-results DIRECTORY`.
The upper-transition search fills the gap between the six-quad Geode interface
and the eight-quad rotated top grid:

```sh
THREADS=8 bash run-geode-top.sh
THREADS=8 SYMMETRY=other CAP=30 bash run-geode-top.sh
THREADS=8 SYMMETRY=two CAP=40 bash run-geode-top.sh
```

The default uses an unseeded cap-26 search with only
`solver/input/geode-top.mesh`: six interface quads, eight top quads, and 36
prescribed side quads in 18 translated pairs. One mirror constrains the interior.
The search finds the 26-cell upper connectivity displayed in the site's
52-cell template, with four fixed cells and eleven exchanged pairs (15 cell
orbits). Each discovered connectivity gets a fresh
harmonic initialization, HexOpt optimization with fixed boundary and mirror
constraints, and exact whole-cell Jacobian validation. An oriented connectivity
comparison identifies matches to cells 26–51 of `R52-improved.json`. The displayed
coordinates and the search's fixed-boundary embedding can differ.

This cut admits only this mirror. `SYMMETRY=other` and `SYMMETRY=two`
preserve the same six-quad interface and eight-quad top, replacing the sides
with a doubly symmetric periodic cut containing 22 quads in eleven translated
pairs. They impose the other diagonal mirror or both mirrors, respectively.
These are searches on prescribed finite boundaries, not an unrestricted
enumeration of periodic transitions.

`run-geode-top.sh` is separate from `run.sh`. It accepts `OUTPUT_DIR`, `CAP`,
`HEXOPT_BINARY`, and the HexOpt settings above. `THREADS` controls its search;
embedding runs sequentially after search completion. Output follows the same
`search/`, `hexopt/`, `status.json`, and `results.json` layout, with an additional
fixture validation log and website connectivity bijection in the default case.
All required inputs are in the repository. The fixture provenance is in
`tests/fixtures/geode-top26.json`.

The search frame is `X=x+y-1`, `Y=y-x`, `Z=z-(zmin+zmax)/2`, with mirror `X=0`
(physical plane `x+y=1`); the other diagonal is `Y=0`. This scales horizontal
lengths by √2 while preserving vertical lengths, so quality values in this
frame are not directly comparable to the website's geometry.

The HexOpt gradient kernel is an external dependency needed only to repeat
coordinate optimization, not to rebuild the site or run the standard tests.
See [the adapter instructions](tools/HEXOPT_NOTES.md).

## Layout and license

- `docs/`: the final website, published geometry, provenance, and downloads.
- `solver/`: search engines, launchers, input fixtures, and solver regressions.
- `tests/`: independent validation and website tests with compact fixtures.
- `tools/`: publication generators, geometry utilities, optimizers, and their inputs.

GitHub Actions tests the sources and publishes only `docs/` to Pages.
The solver, tools, and website code are [MIT licensed](LICENSE). Published mesh
data retain their attribution and upstream notices in
[the data provenance](docs/assets/DATA_SOURCES.md). The independent solver follows
the earlier Hextreme investigations; it is not an upstream Hextreme release.
