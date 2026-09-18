#!/usr/bin/env python3
"""Publish certified Geode connectivities, freshly auditing every exported mesh."""
from collections import Counter
import csv
import hashlib
import json
from pathlib import Path
import shutil
import zipfile

from site_layout import header, search_summary

from build_catalog import symmetries
from geode_geometry import validate
from mesh_tools import read_mesh, write_vtu
from quality_metrics import dense

ROOT = Path(__file__).resolve().parents[1]
DOCS = ROOT/'docs'
OUT = DOCS/'assets/geode'
REFERENCE = ROOT/'docs/assets/geode/meshes/G26-Mitchell.mesh'


def inputs():
    # The published meshes are the final source inventory, as for the pyramid.
    catalog = json.loads((OUT/'catalog.json').read_text())
    for entry in catalog['templates']:
        yield entry['id'], DOCS/entry['mesh'], entry['origin'], entry.get('side')


def build():
    for folder in ('meshes', 'certificates', 'connectivities', 'vtu'):
        (OUT/folder).mkdir(parents=True, exist_ok=True)
    templates = []
    for index, (name, path, origin, side) in enumerate(inputs()):
        p, h, q = read_mesh(path)
        report = validate(path, REFERENCE)
        assert report['accepted'], name
        files = dict(mesh=f'assets/geode/meshes/{name}.mesh',
                     certificate=f'assets/geode/certificates/{name}.json',
                     connectivity=f'assets/geode/connectivities/{name}.json',
                     vtu=f'assets/geode/vtu/{name}.vtu')
        if path.resolve() != (DOCS/files['mesh']).resolve():
            shutil.copyfile(path, DOCS/files['mesh'])
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        report.update(id=name, mesh_sha256=digest)
        (DOCS/files['certificate']).write_text(json.dumps(report, indent=2)+'\n')
        (DOCS/files['connectivity']).write_text(json.dumps(dict(id=name, index_base=0,
            points=p.tolist(), cells=h.tolist(), boundary_quads=q.tolist(),
            mesh_sha256=digest, source=str(path.relative_to(ROOT))), indent=2)+'\n')
        write_vtu(DOCS/files['vtu'], p, h)
        metrics = dense(p, h, True)
        entry = dict(id=name, candidate_id=index, hexes=len(h), vertices=len(p),
            interior_vertices=len(p)-len(set(q.flat)), origin=origin, side=side,
            sha256=digest, metrics=metrics, symmetry=symmetries(p, h),
            camera_target=[0, 0, 0], plane_range=[-.62, .62], plane_radius=1.5, **files)
        if name == 'G26-Mitchell':
            entry['reference_url'] = 'https://www.sandia.gov/files/samitch/files/geode.jou'
        templates.append(entry)
        print(name, len(h), 'certified; min SJ', metrics['min_scaled_jacobian'], flush=True)
    templates.sort(key=lambda t: (t['hexes'], t['id']))
    catalog = dict(mesh_count=len(templates), partial_enumeration_counts=[],
        count_notes={'26': 'published reference'}, templates=templates)
    (OUT/'catalog.json').write_text(json.dumps(catalog, indent=2)+'\n')
    with (OUT/'quality.csv').open('w') as stream:
        writer = csv.writer(stream)
        writer.writerow(['id', 'hexes', 'side', 'min_scaled_jacobian', 'max_condition', 'mesh'])
        for t in templates:
            writer.writerow([t['id'], t['hexes'], t['side'] or '', t['metrics']['min_scaled_jacobian'], t['metrics']['max_condition'], t['mesh']])
    (OUT/'README.md').write_text('''# Geode mesh collection

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
''')
    (DOCS/'geode.html').write_text(page(templates))
    with zipfile.ZipFile(OUT/'geodes.zip', 'w', zipfile.ZIP_DEFLATED) as archive:
        for t in templates:
            for key in ('mesh', 'certificate', 'connectivity', 'vtu'):
                archive.write(DOCS/t[key], t[key])
        for name in ('catalog.json', 'quality.csv', 'README.md', 'geode.jou'):
            archive.write(OUT/name, name)


def page(templates):
    counts = Counter(t['hexes'] for t in templates)
    tabs = ''.join(f'<button data-count="{n}" aria-pressed="false">{n} <span>{count}</span></button>' for n, count in sorted(counts.items()))
    rows = ''.join(f'<tr><td>{t["id"]}</td><td>{t["hexes"]}</td><td><a href="geode-sides.html#{(t["side"] or "Q5-01").removesuffix("-flipped")}">{t["side"] or "Q5-01"}</a></td><td>{t["metrics"]["min_scaled_jacobian"]:.6f}</td><td><a href="{t["mesh"]}" download>Mesh</a> · <a href="{t["connectivity"]}" download>Connectivity</a></td></tr>' for t in templates)
    return f'''<!doctype html>
<html lang="en"><head>
<meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>Mitchell’s Geode — Hexahedral meshes</title>
<meta name="description" content="Certified alternative hexahedral connectivities for Mitchell’s Geode, with interactive views and downloads.">
<link rel="stylesheet" href="assets/site.css"><link rel="icon" href="assets/figures/mark.svg" type="image/svg+xml">
<script type="importmap">{{"imports":{{"three":"https://cdn.jsdelivr.net/npm/three@0.180.0/build/three.module.js","three/addons/":"https://cdn.jsdelivr.net/npm/three@0.180.0/examples/jsm/"}}}}</script>
<script defer src="assets/catalog.js"></script></head>
<body data-catalog="assets/geode/catalog.json">
<a class="skip-link" href="#templates">Skip to meshes</a>
{header('geode.html')}
<main><section class="introduction"><div>
<h1>Hexahedral meshes of Mitchell’s Geode</h1>
<p>A hexahedral transition between a quadrilateral grid and a quadrangulated tetrahedral-mesh boundary.</p>
<p>{len(templates)-1} alternative connectivities alongside Mitchell’s published 26-cell reference. The top and bottom caps remain fixed. All four side walls are equal and opposite sides match directly by translation.</p>
<p>Every mesh below has strictly positive Jacobian determinants throughout all its cells, certified with exact arithmetic. The alternatives have 28–44 cells; none improves on the published cell count.</p>
<p><a href="assets/geode/geodes.zip" download>All meshes and connectivities (.zip)</a> · <a href="assets/geode/quality.csv" download>Quality table (.csv)</a> · <a href="#search">Search and validation</a></p>
</div></section>
<aside class="collection-highlight" aria-labelledby="highlight-title">
<h2 id="highlight-title">Simpler side walls with usable element quality</h2>
<p><a href="geode.html?cells=34#G-Q2-01-H34-01">G-Q2-01-H34-01</a> uses 34 hexes and just <a href="geode-sides.html#Q2-01">two quads per side wall</a>, compared with five in Mitchell’s 26-cell Geode. These simpler interfaces make the template easier to integrate into surrounding meshes while retaining a minimum sampled scaled Jacobian of 0.250.</p>
</aside>
{search_summary([
    ('What is searched', 'Hexahedral fillings of the whole Geode, with Mitchell’s top and bottom caps prescribed. Each run fixes one <a href="geode-sides.html">side-wall connectivity</a>; all four walls are identical and match directly by translation.'),
    ('Starting point', 'The whole-volume solver starts from the prescribed boundary and an empty interior. With <a href="geode-sides.html#Q5-01">Q5-01 side walls</a>, a 26-cell cap, and 180° rotational symmetry, it independently finds Mitchell’s published 26-cell connectivity without a supplied interior mesh.'),
    ('Symmetry', 'The 34–44-cell alternatives impose both vertical diagonal mirror planes, x = y and x = −y. A separate rotation-only mode imposes a 180° turn about the vertical axis, allowing interiors without mirrors; Mitchell’s connectivity admits this constraint.'),
    ('Search limits', 'The default <code>run-geode.sh</code> searches Q2-01 through 34 cells, Q4-01 through 38, and Q5-03 through 38, with no time limit. Set <code>THREADS</code> to the worker count. The gallery includes additional side-wall cases through 44 cells; it is not a complete enumeration of Geode fillings.'),
    ('Embedding', 'HexOpt fixes the cap coordinates and may move side and interior vertices while preserving the imposed symmetry and translation matching. Each displayed mesh passes exact positive-Jacobian certification.'),
    ('Literature', '<a href="https://www.sandia.gov/files/samitch/files/geode.jou">Mitchell’s published 26-cell Geode</a> is the reference and has side pattern <a href="geode-sides.html#Q5-01">Q5-01</a>. The 28–44-cell alternatives are results of this project; none has fewer cells than the reference.'),
])}
<section class="collection" id="templates"><h2>Geode meshes</h2>
<p class="collection-description">Grouped by cell count; sorting applies within each group. Drag to rotate; scroll to zoom. Shrink separates the cells. Symmetry is measured on the exported coordinates. SJ means scaled Jacobian: larger values indicate better cell shape. Smaller condition numbers indicate less distortion. Both are sampled measurements; see the <a href="methodology.html#quality">quality definitions</a>.</p>
<div class="collection-toolbar"><div class="count-tabs" role="group" aria-label="Filter by element count"><button class="active" data-count="all" aria-pressed="true">All <span>{len(templates)}</span></button>{tabs}</div>
<label class="sort-control">Sort by <select id="sort"><option value="quality">Minimum SJ (descending)</option><option value="condition">Condition number (ascending)</option><option value="id">Mesh ID</option></select></label></div>
<p id="collection-status" class="collection-status" role="status">Loading the Geode catalog…</p><div id="gallery"></div>
<noscript><p>The interactive views require JavaScript. All downloads are available in the table below.</p></noscript></section>
<section class="collection" id="downloads"><h2>Connectivity downloads</h2>
<p>JSON files include coordinates, eight vertex indices per cell, and boundary quads. Indices start at zero. The archive also includes VTU files and validation reports.</p>
<p>Side-wall IDs link to the pattern catalog: for example, Q5-03 is a five-quad pattern. A <code>-flipped</code> placement reverses top and bottom. These placements share one side-pattern ID, but can give different Geode fillings because the two caps differ.</p>
<div style="overflow-x:auto"><table class="metrics"><thead><tr><th>Mesh</th><th>Cells</th><th>Side wall</th><th>Minimum sampled SJ</th><th>Downloads</th></tr></thead><tbody>{rows}</tbody></table></div></section>
<section class="collection" id="search"><h2>Search and validation</h2>
<p>The whole-volume search independently finds Mitchell’s 26-cell connectivity from its boundary alone. The twenty-six 34–44-cell alternatives likewise start from empty interiors, using both diagonal mirrors and HexOpt embedding. The two 28-cell alternatives come from a separate local cavity-replacement tool: it takes the published Geode and replaces three cells with five. That local operation is separate from the whole-volume solver’s starting condition. The two 28-cell meshes are included specifically for their improved minimum sampled scaled Jacobian: 0.306 and 0.310 versus 0.256 for Mitchell’s mesh, at the cost of two additional cells.</p>
<p>The Q7-06 meshes and their flipped placements contain highly stretched cells, with worst sampled condition numbers from about 13,000 to 201 million. Positive Jacobians certify valid geometry; condition numbers describe the remaining distortion.</p>
<p>The rotation-only constraint is <code>(x,y,z) → (−x,−y,z)</code>. Mitchell’s 26-cell connectivity admits this symmetry as thirteen cell pairs. Its exported interior coordinates do not realize that symmetry exactly, so the viewer’s measured geometry symmetry differs from the connectivity symmetry.</p>
<p>The collection is not an unrestricted enumeration of Geode fillings. Failed embedding attempts do not prove geometric impossibility.</p>
<p>The top and bottom cap coordinates stay fixed during embedding; side vertices may move. Embedding ties side copies and edge heights, preserving direct translation matching. HexOpt uses the repository’s projected-gradient adapter and the upstream scaled-Jacobian gradient. Each final mesh passes whole-cell Jacobian certification and the boundary and topology checks described in its report. Sampled quality values do not replace those certificates.</p>
<p><a href="geode-sides.html">All 432 side-wall patterns through twelve quads</a> · <a href="assets/geode/README.md">Data provenance</a> · <a href="assets/geode/geode.jou" download>Published source journal</a> · <a href="methodology.html#quality">Quality definitions</a></p></section></main>
<footer><span>Geode mesh collection</span><a href="assets/geode/geodes.zip" download>Mesh archive</a></footer></body></html>
'''


if __name__ == '__main__':
    build()
