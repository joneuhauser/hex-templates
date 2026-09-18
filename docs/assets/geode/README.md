# Geode mesh collection

The archive contains Mitchell's published 26-cell geometry, two alternative
28-cell connectivities, and twenty-six certified embeddings from two-mirror searches.
All have fixed published caps and identical sides matching directly by translation.
JSON connectivity indices are zero-based; Medit .mesh indices are one-based.
Mesh hashes bind certificates and connectivity exports to exact coordinate files.

The 34–44-cell meshes use the repository's projected-gradient HexOpt adapter
(upstream sJGrad kernel, commit 5f46bf4f1c8dbff6cbfbec7ab2c2295ffd4d16b9).
This is a fixed-step projected-gradient adapter. Failed embedding does not prove
geometric impossibility. The four Q7 meshes have highly stretched cells:
worst sampled condition numbers range from about 13,000 to 201 million.

Every exported cell has an exact positive-Jacobian certificate throughout
its reference cube. Reports also check incidence, orientation, edge/vertex
links, homology, cap preservation, side matching, containment, and total volume.
Quality values are measurements on 21³ samples per cell, not certified minima.
The symmetry display measures exported coordinates and oriented cells; the
published 26-cell connectivity has a half-turn even though its coordinates do
not realize it exactly. No complete enumeration of Geode fillings is claimed.

Source journal: https://www.sandia.gov/files/samitch/files/geode.jou
Rebuild from the repository root: python3 tools/build_geode_catalog.py
