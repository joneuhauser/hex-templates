# Published periodic grid-rotation meshes

This collection contains the seven connectivity classes displayed on the site.
Four rectangular bottom quads connect to eight rotated top quads in an x/y
periodic slab. Each entry retains its stored height and coordinate geometry.

`current-catalog.json` and the seven `periodic/*.json` meshes are the authoritative
inputs. The JSON meshes preserve quotient vertex indices and per-corner integer
lattice offsets; a corner is its vertex coordinate plus offset times period.
Source labels and hashes identify the source of each mesh.

Medit/VTU downloads show a representative tile or a 2×2 completion. Only top and
bottom are physical boundaries; display cuts are artificial periodic seams.
Certificates check quotient incidence, links, homology, cap matching, volume,
and exact whole-cell positive Jacobians. Quality values use 21³ samples per
cell and are not certified extrema. Connectivity classes are compared by exact
colored cubical-flag isomorphism, preserving top/bottom roles and attachments.

Rebuild: `python3 tools/build_periodic_catalog.py` from the repository root.
The code is MIT licensed; source attribution accompanies the mesh data. No complete enumeration or minimum cell count is claimed.
