# HexOpt adapter

The adapter uses the `sJGrad` kernel from [CMU-CBML/HexOpt](https://github.com/CMU-CBML/HexOpt) at commit `5f46bf4f1c8dbff6cbfbec7ab2c2295ffd4d16b9`. Obtain that upstream checkout separately; its source is not vendored here. This adapter uses fixed-step projected gradients; the full augmented-Lagrangian/L-BFGS procedure is described in [Tong and Zhang's paper](https://arxiv.org/abs/2410.11656).

From the repository root, with `HEXOPT_SOURCE` pointing to that checkout:

```bash
mkdir -p solver/build
c++ -O3 -std=c++17 -fopenmp -include cmath -include cfloat -I "$HEXOPT_SOURCE" tools/hexopt_fixed.cpp \
  "$HEXOPT_SOURCE/meshQuality.cpp" "$HEXOPT_SOURCE/geometry.cpp" \
  -o solver/build/hexopt_fixed
```

The forced standard includes supply headers omitted by this upstream checkout.

After building the search solver too, `THREADS=8 bash run-pyramid.sh` and
`THREADS=8 bash run-geode.sh` run search, HexOpt, and exact validation in sequence.
The Geode launcher processes Q2-01, Q4-01, and Q5-03 in separate folders.
See the root README for output locations and environment settings. The shared
Python driver also supports `--profile pyramid` for existing pyramid candidates;
it fixes every boundary vertex, projects each restart onto the two diagonal
mirrors, and applies the pyramid-specific geometry validator.

The executable takes an input text file and an output XYZ file. The input contains:

1. `vertex_count cell_count fixed_vertex_count axial_flag` (1 for axial mirrors, 0 for diagonal).
2. One XYZ row per vertex. Fixed boundary vertices must come first.
3. One row of four zero-based vertex indices per vertex: identity, first reflection, second reflection, product.
4. One row of eight zero-based vertex indices per hex, in the `.mesh` ordering used by `mesh_tools.py`.

Run `solver/build/hexopt_fixed input.txt output.xyz`. It preserves connectivity; replace the source vertex coordinates with the returned XYZ rows, then run `validate_mesh.py` and recompute `quality_metrics.dense(points, cells, True)`. Never publish optimizer output based on the reported threshold alone. Retain the initial embedding if it is certified and has the better dense sampled minimum SJ.

For each iteration, cell gradients are accumulated at shared vertices, averaged under the two-reflection group, and zeroed by retaining fixed boundary coordinates. The step scale is `0.001 / max(1, gradient_norm)`. Threshold increments are 0.01, with at most 300,000 iterations and a 5,000-iteration stall limit. The implementation returns the last threshold-feasible coordinates, or the input if no threshold was attained. Other starts can converge to different local optima.

The 88-cell literature asset is the normalized published geometry, without this optimization. Its exported geometry has only one measured mirror; the two-reflection adapter must not be applied to it without a separate justified constraint choice.

For Geodes, an optional third argument supplies bounded affine coordinates:
`hexopt_fixed input.txt output.xyz constraints.txt`. The file starts with the
number of degrees of freedom. Each subsequent row is `initial lower upper
term_count`, followed by that many `vertex coordinate coefficient` triples.
Indices start at zero. Columns must have disjoint coordinate support. The
adapter projects gradients onto these columns and coordinates onto their bounded
affine span. This mode supersedes the fixed-count and four-action projection;
the supplied basis must encode all required constraints.

`tools/hexopt_geode.py` builds this basis with fixed caps, shared side coordinates,
shared vertical-edge heights, and the two diagonal interior actions. It retries
deterministic perturbations, certifies whole-cell Jacobians and all Geode boundary
conditions, and exports only accepted meshes. `HEXOPT_MAX_ITERATIONS` and
`HEXOPT_STALL_LIMIT` override the defaults, accepting integers from 1 to 10,000,000.
The root launchers use a 50,000-iteration stall limit. Failed attempts do not
establish geometric impossibility.

H44-class16 also retains the authors’ original coordinates byte-for-byte. Its gallery asset is not an output of this adapter. Use `validate_mesh.py --profile published44` to account explicitly for the rounded published boundary; Jacobian positivity is still certified on the original coordinates.

For periodic rotated-grid slabs, `tools/hexopt_periodic.py` uses the same affine
adapter. It ties every periodic copy of each interior quotient vertex, fixes
the bottom grid, and gives the entire top grid two shared x/y translation
variables. Cap heights and the top grid’s relative vertex positions are fixed.
The adapter uses broad finite translation/displacement bounds of ±2 periods.
Output is checked against a reference with precisely that rigid top shift,
then receives independent periodic topology and exact Jacobian checks. Dense
21³ quality, rather than the optimizer’s corner threshold, determines whether
to adopt the geometry.
