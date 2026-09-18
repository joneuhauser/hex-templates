# Refinement templates

6 certified hexahedral meshes for 9-to-1, 16-to-4 and 25-to-9 regular quad-grid
transitions, with 28 enumerated side-wall patterns. Coordinates use
[-1,1]³; the fine grid is on top. Four identical walls match by translation.
All displayed interiors have two axial mirrors. JSON connectivity indices
start at zero; Medit mesh indices start at one.

The catalog and mesh files are the authoritative geometry inputs. Certificates
bind exact Jacobian and topology checks to mesh hashes. Quality measurements
use 21³ samples per cell and are not certified extrema. The 14- and 34-cell
meshes have relatively poor condition numbers despite their positive Jacobians.

Side JSON files contain the quads and reflection under `side`, and rational
coordinates under `side_points_exact`. Boundary-only Medit files are supplied
for all 28 cases. Enumeration and 3D search records state the finite
bounds and distinguish complete searches from time-limited searches.

From the repository root, `python3 tools/build_refinement_catalog.py` rebuilds
the collection. `THREADS=8 bash run-refinement.sh` finds volume connectivities
from the displayed boundaries, then runs HexOpt and exact validation.
The meshes are independent search results. Schneiders, Refining Quadrilateral
and Hexahedral Element Meshes (1996), figures 10 and 16a/b, gives closely related
3-refinement and 2-refinement constructions. F16-4-H28 has four seven-cell
quarters; F9-1-H13 is comparable to the transition below figure 10's fine-grid
layer. These are figure-based comparisons, not graph-isomorphism checks against
published mesh data. See the collection page for references. No global minimum
or novelty is claimed.
