# Search inputs and validation fixtures

## Published pyramid connectivities

The normalized fixtures preserve the published cell connectivity, with a vertex
renumbering that places the solver's 18 boundary vertices first. No cell is
removed, added, split, or reconnected.

- `published36-*`: Xiang–Liu's 36-cell connectivity, from the mesh accompanying
  https://kilian.ac/posts/meshing-game/ at
  https://kilian.ac/posts/meshing-game/meshes/pyramid.mesh . The original local
  input is `published36.mesh`. Coordinates are transformed by
  `(x,y,z) -> (x/2,-z/2,(y+1.5)*sqrt(2)/3)`.
- `published44-*`: Verhetsel–Pellerin–Remacle's 44-cell mesh, from
  `results/schneiders-44.mesh` in the topological-hex 0.2.0 distribution,
  https://www.hextreme.eu/Download/topological-hex-0.2.0.tar.gz . Coordinates are
  transformed by `(x,y,z) -> (x,-z,y*sqrt(2)/1.41421)`.
  These 44-cell data derive from the GPL-distributed upstream files; their
  upstream license is included in `COPYING-44.txt`. The standalone solver's
  source license does not replace that data provenance.

For both, rounded boundary coordinates are restored to the exact prescribed
pyramid coordinates. Interior coordinates in `*-normalized.mesh` are otherwise
untouched after the affine transformation. Because these exported positions are
not exactly reflection symmetric, `*.actions` supplies the exact combinatorial
reflection permutations, verified on the boundary, edges, faces, and cells.
The four columns are identity, first reflection, second reflection, and product;
indices are zero-based. The solver additionally checks the group law and
oriented cell preservation before using these actions for replay.

The `*-symmetric.mesh` files average the same interior vertex orbits under the
verified reflection action. They preserve cell connectivity and boundary
coordinates and allow automatic coordinate-based symmetry detection. Maximum
coordinate changes from the normalized fixtures are 0.0001210315 for 36 cells
and 0.00135305327599 for 44 cells. Both normalized and symmetrized fixtures pass
independent topology checks and exact-rational whole-cell Jacobian positivity
certification locally; results are recorded in the accompanying `*.validation.json` files.

The 36-cell topology admits the diagonal pair x=y, x=-y. The 44-cell topology
admits the axial pair x=0, y=0 in this common coordinate convention. They have
12 and 14 cell orbits respectively. These are seeded branch checks, not
unseeded discoveries or exhaustive enumerations at caps 36 or 44.

## One-plane regression derivatives

`one-axial-*` and `one-diagonal-*` translate the supplied grid2 fixture by
(0,3,0) and (3,3,0), respectively. They preserve the first reflection and remove
the second from the coordinate set. `published36-diagonal-one.actions` and
`published44-axial-one.actions` contain the first two columns of the existing
verified action files: identity and the first reflection.

## Skew cube fixture

`skew-cube-boundary.mesh` and `skew-cube-seed.mesh` are generated from the
MIT-licensed cube fixture using the linear map
`(x,y,z) -> (x+y/4+z/8, 3*y/2+3*z/8, 2*z)`. Its determinant is 3, so the
single affine hex has strictly positive Jacobian everywhere. Its three edge
directions have distinct lengths and nonzero pairwise inner products,
precluding a plane reflection symmetry. It exercises zero-mirror input and
positive unseeded one-/seven-cell searches.

## Geode upper transition

`geode-top.mesh` contains only the 52 boundary vertices and 50 quads of a
26-cell upper transition. Its replay mesh and provenance are
`tests/fixtures/geode-top26.mesh` and `tests/fixtures/geode-top26.json`.
Cells 26–51 of the regularized 52-cell transition use the periodic
representatives recorded in the provenance. The regularized coordinates are
stored with the Geode52 ParaView viewer. The launcher checks their oriented connectivity
against the repository's `docs/assets/periodic/periodic/R52-improved.json`.

Coordinates are projected onto the physical mirror `x+y=1` (maximum source
reflection discrepancy about 5.47e-8), then transformed to
`(X,Y,Z)=(x+y-1,y-x,z-(zmin+zmax)/2)`. Horizontal distances scale by √2.
Boundary faces comprise six Geode interface quads, eight top quads, and 36
artificial side quads paired by translations in the lattice `(1,-1,0)`, `(1,1,0)`.
The 26 cells form 15 mirror orbits, including four fixed cells.
`tools/geode_top.py` also constructs a two-mirror periodic cut with identical
caps and 22 side quads, for the alternative-symmetry searches.

## Refinement boundaries

`refinement/` contains boundary-only cubes for 9-to-1, 16-to-4, 25-to-9 and 25-to-1 grid
transitions. `cases.json` records exact rational side-wall coordinates and
connectivity; the top and bottom use regular grids on `[-1,1]²` at z = ±1.
All four walls are identical, left/right symmetric, and translation matching.

Generate these inputs into a fresh folder with
`python3 tools/refinement.py sides --output DIRECTORY`. This enumerates side
disks through eight quads for 9-to-1 and 16-to-4, and ten for 25-to-9 and 25-to-1,
with zero or one extra vertex per vertical edge,
and cross-checks orbit insertion against individual-quad insertion.
`enumeration.json` records the source hash and output hashes for both modes.
No volume filling is supplied as search input.
