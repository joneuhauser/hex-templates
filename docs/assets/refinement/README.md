# Refinement templates

Five certified hexahedral meshes for 9-to-1 and 16-to-4 regular quad-grid
transitions, with sixteen enumerated side-wall patterns. Coordinates use
[-1,1]³; the fine grid is on top. Four identical walls match by translation.
All displayed interiors have two axial mirrors. JSON connectivity indices
start at zero; Medit mesh indices start at one.

The catalog and mesh files are the authoritative geometry inputs. Certificates
bind exact Jacobian and topology checks to mesh hashes. Quality measurements
use 21³ samples per cell and are not certified extrema. The 14- and 34-cell
meshes have relatively poor condition numbers despite their positive Jacobians.

Side JSON files contain the quads and reflection under `side`, and rational
coordinates under `side_points_exact`. Boundary-only Medit files are supplied
for all sixteen cases. Enumeration and 3D search records state the finite
bounds and distinguish complete searches from time-limited searches.

From the repository root, `python3 tools/build_refinement_catalog.py` rebuilds
the collection. `THREADS=8 bash run-refinement.sh` finds volume connectivities
from the five displayed boundaries, then runs HexOpt and exact validation.
The meshes are independent search results; equivalence to published refinement
template sets has not been established. No global minimum or novelty is claimed.
