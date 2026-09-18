"""Revalidate and rebuild the seven displayed periodic connectivity classes.

The current catalog and its periodic JSON meshes are the source inventory.
No development runs, discarded embeddings or research directories are needed.
"""
import csv
import hashlib
import json
from pathlib import Path
import zipfile
from collections import Counter

import numpy as np
from site_layout import header, search_summary
from mesh_tools import FACES, write_mesh, write_vtu
from periodic_transition import completion, load, validate
from quality_metrics import dense
from periodic_topology_classes import classify
from draw_periodic_boundary import sketch

ROOT = Path(__file__).resolve().parents[1]
DOCS = ROOT/'docs'
OUT = DOCS/'assets/periodic'


def dump(path, obj):
    path.write_text(json.dumps(obj, indent=2, allow_nan=False)+'\n')


def physical_caps(p, h):
    zmin, zmax = p[:, 2].min(), p[:, 2].max()
    return np.array([c[f] for c in h for f in FACES
                     if np.all(abs(p[c[f], 2]-zmin)<1e-9) or np.all(abs(p[c[f], 2]-zmax)<1e-9)])


def build():
    catalog = json.loads((OUT/'current-catalog.json').read_text())
    templates = catalog['templates']
    for entry in templates:
        source = DOCS/entry['connectivity']
        p, h, offsets, period = load(source)
        report = validate(p, h, offsets, period)
        assert report['accepted'], entry['id']
        digest = hashlib.sha256(source.read_bytes()).hexdigest()
        report.update(id=entry['id'], periodic_json_sha256=digest)
        dump(DOCS/entry['certificate'], report)
        pp, hh = completion(p, h, offsets, period, 1)
        write_mesh(DOCS/entry['mesh'], pp, hh, physical_caps(pp, hh))
        write_vtu(DOCS/entry['vtu'], pp, hh)
        rp, rh = completion(p, h, offsets, period, 2)
        write_mesh(DOCS/entry['repeated_mesh'], rp, rh, physical_caps(rp, rh))
        write_vtu(DOCS/entry['periodic_vtu'], rp, rh)
        metrics = dense(pp, hh, True)
        metrics['min_mean_ratio'] = report['sampled']['min_mean_ratio']
        entry.update(metrics=metrics, sha256=digest)
        print(entry['id'], len(h), 'CERTIFIED', metrics['min_scaled_jacobian'], flush=True)
    groups, witnesses = classify([(t['id'], load(DOCS/t['connectivity'])) for t in templates])
    assert len(groups) == len(templates), 'Published representatives must have distinct connectivities'
    dump(OUT/'topology-classes.json', dict(groups=groups, witnesses=witnesses,
        definition='Exact isomorphism of colored periodic cubical flags, preserving top/bottom roles and face attachments; coordinates and vertex labels ignored.'))
    current_downloads(catalog)
    (DOCS/'periodic.html').write_text(page(catalog))
    print('GALLERY', len(templates), 'certified periodic connectivity classes')


def current_downloads(catalog):
    templates=[{key:value for key,value in t.items() if key!='variants'} for t in catalog['templates']]
    current=dict(mesh_count=len(templates),geometry_count=len(templates),partial_enumeration_counts=[],count_notes={},templates=templates)
    dump(OUT/'current-catalog.json',current)
    with (OUT/'current-quality.csv').open('w') as stream:
        w=csv.writer(stream);w.writerow(['class','geometry','cells_per_tile','height','min_sampled_sj','max_sampled_condition','periodic_json'])
        for t in templates:w.writerow([t['topology_class'],t['id'],t['hexes'],t['height'],t['metrics']['min_scaled_jacobian'],t['metrics']['max_condition'],t['connectivity']])
    with zipfile.ZipFile(OUT/'current-meshes.zip','w',zipfile.ZIP_DEFLATED) as archive:
        for t in templates:
            for key in ['mesh','repeated_mesh','connectivity','certificate','vtu','periodic_vtu']:archive.write(DOCS/t[key],t[key])
        archive.write(OUT/'current-catalog.json','catalog.json')
        archive.write(OUT/'current-quality.csv','quality.csv')
    sketch(OUT/'boundary.svg')


def page(catalog):
    templates=catalog['templates'];counts=Counter(t['hexes'] for t in templates)
    tabs=''.join(f'<button data-count="{n}" aria-pressed="false">{n} <span>{count}</span></button>' for n,count in sorted(counts.items()))
    rows=''.join(f'<tr><th scope="row">{t["topology_class"]}</th><td>{t["hexes"]}</td><td>{t["metrics"]["min_scaled_jacobian"]:.6f}</td><td>{t["metrics"]["max_condition"]:.3f}</td><td><a href="{t["connectivity"]}" download>Periodic JSON</a> · <a href="{t["periodic_vtu"]}" download>2×2 VTU</a> · <a href="{t["certificate"]}">Certificate</a></td></tr>' for t in templates)
    return f'''<!doctype html>
<html lang="en"><head>
<meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>Periodic grid rotation — Hexahedral meshes</title>
<meta name="description" content="Known hexahedral fillings between rectangular and rotated periodic quad grids, with interactive periodic completions, quality measurements and downloads.">
<link rel="stylesheet" href="assets/site.css"><link rel="icon" href="assets/figures/mark.svg" type="image/svg+xml">
<script type="importmap">{{"imports":{{"three":"https://cdn.jsdelivr.net/npm/three@0.180.0/build/three.module.js","three/addons/":"https://cdn.jsdelivr.net/npm/three@0.180.0/examples/jsm/"}}}}</script>
<script defer src="assets/catalog.js"></script></head>
<body data-catalog="assets/periodic/current-catalog.json">
<a class="skip-link" href="#templates">Skip to meshes</a>
{header('periodic.html')}
<main><section class="introduction"><div>
<h1>From a rectangular grid to a rotated grid</h1>
<p>Hexahedral fillings of a slab periodic in x and y. Four rectangular quads on the bottom connect to eight quads in the rotated top grid. Only the two end grids are physical boundaries: cells may cross the sides of the drawn unit tile.</p>
<p>This collection contains <strong>{len(templates)} distinct periodic connectivity classes</strong> with {min(counts)}–{max(counts)} cells per tile. Each viewer shows the selected validated geometry for one class.</p>
<p><a href="assets/periodic/current-meshes.zip" download>All meshes (.zip)</a> · <a href="assets/periodic/current-quality.csv" download>Quality table (.csv)</a> · <a href="#downloads">Download table</a></p>
</div><figure><img src="assets/periodic/boundary.svg" alt="A 1×1×1 wireframe cube with the rectangular bottom grid and rotated top grid; the side walls are empty."><figcaption>1×1×1 periodic tile. Only the end grids are prescribed.</figcaption></figure></section>
<aside class="collection-highlight" aria-labelledby="highlight-title">
<h2 id="highlight-title">Featured periodic transition</h2>
<p><a href="periodic.html?cells=52#P52-01">P52-01</a> combines a 26-cell Geode with a 26-cell upper transition. Its 52 cells per tile achieve a minimum sampled scaled Jacobian of 0.360, the highest in this collection. The upper-half connectivity can be found with the one-mirror search described below.</p>
</aside>
{search_summary([
    ('Problem and scope', 'The gallery contains full periodic transitions from four bottom quads to eight rotated top quads per tile. The runnable search targets the upper half of the 52-cell template P52-01: a six-quad Geode interface to an eight-quad top. The gallery is not a complete periodic enumeration.'),
    ('Prescribed lower half', 'The lower 26-cell Geode is prescribed, not searched. Only its upper interface enters the upper-half search boundary. The default finite cut also fixes 36 side quads in 18 translated pairs; all boundary coordinates stay fixed during HexOpt embedding.'),
    ('Starting point', 'The upper interior starts empty, with no prescribed upper cells. The default search finds a 26-cell upper filling, giving 52 cells when joined to the lower Geode.'),
    ('Symmetry', 'The default upper search imposes the diagonal mirror x + y = 1. The lower Geode connectivity has a 180° rotation; a common mirror of the whole 52-cell assembly is not required. The other diagonal mirror, y = x, and both mirrors together are also supported for the upper half.'),
    ('Search limits', '<code>run-geode-top.sh</code> defaults to a cap of 26 upper cells. <code>SYMMETRY=other</code> defaults to 30 and <code>SYMMETRY=two</code> to 40; these use a different symmetric side cut with 22 quads in 11 translated pairs, preserving the six- and eight-quad caps. All caps count upper cells only. Set <code>THREADS</code> for workers and optionally <code>CAP</code>; there is no default time limit.'),
    ('Literature', '<a href="https://www.sandia.gov/files/samitch/files/geode.jou">Mitchell’s 26-cell Geode</a> supplies the lower building block of P52-01. The full periodic transitions and upper fillings shown here are results of this project; no external literature solution is used as a reference for the full transition.'),
])}
<section class="collection" id="templates"><h2>Distinct periodic connectivities</h2>
<p class="collection-description">Classes are compared after periodic identification, ignoring coordinates and vertex labels. Each class is shown with one certified embedding; no globally optimal quality is claimed. SJ means scaled Jacobian: larger values indicate better cell shape. Smaller condition numbers indicate less distortion. Both are sampled measurements; see the <a href="methodology.html#quality">quality definitions</a>. Sorting applies within each cell-count group. Drag to rotate; scroll to zoom. Switch between one tile and a 2×2 completion. Shrink exposes cell interfaces. Counts and quality values always refer to one periodic tile; cap outlines show the physical boundaries.</p>
<div class="collection-toolbar"><div class="count-tabs" role="group" aria-label="Filter by element count"><button class="active" data-count="all" aria-pressed="true">All <span>{len(templates)}</span></button>{tabs}</div>
<label class="sort-control">Sort by <select id="sort"><option value="quality">Minimum SJ (descending)</option><option value="condition">Condition number (ascending)</option><option value="id">Class ID</option></select></label></div>
<p id="collection-status" class="collection-status" role="status">Loading the periodic mesh catalog…</p><div id="gallery"></div>
<noscript><p>The interactive views require JavaScript. Every mesh and certificate is available in the download table below.</p></noscript></section>
<section class="collection" id="downloads"><h2>Mesh downloads</h2>
<p>The periodic JSON preserves quotient vertex indices and lattice offsets. VTU files show unshrunk 2×2 completions. Heights and quality are measured in the stored flat geometry, before any cylindrical or flow-domain mapping.</p>
<div style="overflow-x:auto"><table class="metrics"><thead><tr><th>Class</th><th>Cells / tile</th><th>Min. sampled SJ ↑</th><th>Max. condition ↓</th><th>Downloads</th></tr></thead><tbody>{rows}</tbody></table></div></section>
<section class="collection" id="validation"><h2>Validation</h2>
<p>Every displayed mesh passes periodic topology checks and an exact positive-Jacobian certificate throughout every trilinear cell. Quality measurements use 21³ samples per cell and are not certified extrema.</p>
<p>HexOpt keeps the bottom grid fixed and permits the top grid to translate rigidly in x and y, preserving its shape and height. Periodic vertex copies move together.</p>
<h3>Searching the upper transition</h3>
<p>The 52-cell template consists of a 26-cell Geode below a 26-cell upper transition. The shared interface has six quads; the upper transition connects it to the eight-quad top grid.</p>
<p>The repository launcher <code>THREADS=8 bash run-geode-top.sh</code> searches from a boundary-only input at cap 26 and finds the upper connectivity of the displayed 52-cell template. A diagonal mirror constrains the search; the 26-cell solution has four fixed cells and eleven exchanged pairs, giving 15 cell orbits. HexOpt optimizes each candidate, followed by exact cell validation.</p>
<p>The boundary consists of the six-quad Geode interface, eight-quad top, and a periodic cut with 36 side quads in translated pairs. The displayed coordinates and the search's fixed-boundary embedding can differ. Alternative one- and two-mirror searches preserve the caps and use a different periodic side cut. These finite-boundary searches do not enumerate all periodic transitions.</p>
<p><a href="methodology.html#quality">Quality definitions</a></p></section></main>
<footer><span>Periodic rotation collection</span><a href="assets/periodic/current-meshes.zip" download>Mesh archive</a></footer></body></html>
'''


if __name__=='__main__':
    build()
