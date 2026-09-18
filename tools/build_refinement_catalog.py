"""Validate and publish refinement meshes and their enumerated side-wall patterns."""
from collections import Counter
import csv
import hashlib
import json
from pathlib import Path
import shutil
import zipfile

from build_geode_sides import diagram
from mesh_tools import read_mesh, write_vtu
from quality_metrics import dense
from refinement import INPUT, validate
from site_layout import header, search_summary

ROOT = Path(__file__).resolve().parents[1]
DOCS = ROOT/'docs'
OUT = DOCS/'assets/refinement'


def dump(path, data):
    path.write_text(json.dumps(data, indent=2)+'\n')


def page(catalog, cases, searches):
    templates = catalog['templates']
    counts = Counter(t['hexes'] for t in templates)
    tabs = ''.join(f'<button data-count="{n}" aria-pressed="false">{n} <span>{k}</span></button>' for n,k in sorted(counts.items()))
    rows = ''.join(f'<tr><th><a href="?cells={t["hexes"]}#{t["id"]}">{t["id"]}</a></th>'
        f'<td>{t["fine_quads"]} → {t["coarse_quads"]}</td><td>{t["hexes"]}</td>'
        f'<td><a href="#side-{t["boundary_case"]}">{t["side_quads"]}</a></td>'
        f'<td>{t["metrics"]["min_scaled_jacobian"]:.3f}</td><td>{t["metrics"]["max_condition"]:.2f}</td>'
        f'<td><a href="{t["mesh"]}" download>Mesh</a> · <a href="{t["connectivity"]}" download>JSON</a></td></tr>' for t in templates)
    walls = []
    for fine,coarse in sorted({(c['top_grid']**2,c['bottom_grid']**2) for c in cases}):
        selected = [c for c in cases if c['top_grid']**2 == fine]
        walls.append(f'<h3>{fine} → {coarse}: {len(selected)} side-wall patterns</h3>')
        for count in sorted({c['side_quads'] for c in selected}):
            walls.append(f'<section class="side-group"><h4>{count} quads per wall</h4><div class="side-grid">')
            for c in sorted((c for c in selected if c['side_quads'] == count), key=lambda c:(-c['side_quality'],c['id'])):
                name = c['id']; found = [t for t in templates if t['boundary_case']==name]
                links = ' · '.join(f'<a href="?cells={t["hexes"]}#{t["id"]}">{t["hexes"]}-cell mesh</a>' for t in found)
                status = links or 'No certified volume mesh in this collection.'
                limits = []
                for planes in ['axial', 'diagonal']:
                    runs = [r for r in searches['runs'] if r['case']==name and r['planes']==planes]
                    complete = max((r['cap'] for r in runs if r['returncode']==0), default=0)
                    partial = max((r['cap'] for r in runs if r['returncode']==124), default=0)
                    label = f'{planes.capitalize()}: complete through {complete}' if complete else f'{planes.capitalize()}: no completed search'
                    if partial > complete:
                        label += f'; cap {partial} incomplete'
                    if runs and not any(r['candidates'] for r in runs):
                        label += ', no candidates'
                    limits.append(label)
                walls.append(f'<article class="side-card" id="side-{name}"><h3>{name}</h3>'
                    f'<img src="assets/refinement/sides/{name}.svg" width="280" height="278" loading="lazy" alt="{fine} to {coarse} side wall with {count} quads">'
                    f'<p>{c["vertical_interior_vertices"]+2} vertices per vertical edge · minimum corner mean ratio {c["side_quality"]:.3f}</p>'
                    f'<p>{status}</p><p>{"<br>".join(limits)}</p><div class="side-downloads"><a href="assets/refinement/sides/{name}.json" download>Side JSON</a>'
                    f'<a href="assets/refinement/boundaries/{name}.mesh" download>3D boundary</a></div></article>')
            walls.append('</div></section>')
    return f'''<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Refinement templates — Hexahedral meshes</title>
<meta name="description" content="Certified 9-to-1, 16-to-4 and 25-to-9 hexahedral refinement templates and symmetric side-wall quadrangulations.">
<link rel="stylesheet" href="assets/site.css"><link rel="stylesheet" href="assets/geode-sides.css">
<link rel="icon" href="assets/figures/mark.svg" type="image/svg+xml">
<script type="importmap">{{"imports":{{"three":"https://cdn.jsdelivr.net/npm/three@0.180.0/build/three.module.js","three/addons/":"https://cdn.jsdelivr.net/npm/three@0.180.0/examples/jsm/"}}}}</script>
<script defer src="assets/catalog.js"></script></head>
<body data-catalog="assets/refinement/catalog.json"><a class="skip-link" href="#templates">Skip to meshes</a>
{header('refinement.html')}
<main><section class="introduction"><div><h1>Hexahedral refinement templates</h1>
<p>Transitions from a 3×3 quad grid to one quad (9 → 1), from a 4×4 grid to a 2×2 grid (16 → 4), and from a 5×5 grid to a 3×3 grid (25 → 9). The fine grid is on top and the coarse grid below; the same connectivity can be used in either direction.</p>
<p>{len(templates)} certified meshes and all {len(cases)} tested side-wall patterns. Every mesh has strictly positive Jacobian determinants throughout its cells.</p>
<p><a href="assets/refinement/templates.zip" download>Meshes and side walls (.zip)</a> · <a href="#side-walls">Side-wall patterns</a> · <a href="assets/refinement/quality.csv" download>Quality table (.csv)</a></p>
</div><figure><img src="assets/previews/refinement.svg" alt="The 32-cell 16-to-4 refinement mesh, with cells separated for inspection."><figcaption>32-cell 16 → 4 transition.</figcaption></figure></section>
<aside class="collection-highlight" aria-labelledby="highlight-title"><h2 id="highlight-title">Compact transitions and a quality alternative</h2>
<p><a href="?cells=13#F9-1-H13">F9-1-H13</a>, <a href="?cells=28#F16-4-H28">F16-4-H28</a>, and <a href="?cells=33#F25-9-H33">F25-9-H33</a> are the smallest fillings found for the three transitions, respectively. The 33-cell 25 → 9 mesh uses six quads per side wall and has minimum sampled scaled Jacobian 0.260.
<a href="?cells=32#F16-4-H32">F16-4-H32</a> improves both minimum sampled scaled Jacobian (0.511 versus 0.200) and maximum sampled condition number (4.79 versus 7.12), using four extra cells.</p></aside>
{search_summary([
('What is searched', 'Conforming hex fillings of a cube with prescribed regular top and bottom grids. Each run fixes one of the <a href="#side-walls">side-wall patterns below</a>. All four walls are identical and opposite walls match by direct translation.'),
('Starting point', 'The solver starts with the boundary and an empty interior. No existing volume connectivity or lower submesh is prescribed.'),
('Symmetry', 'Both axial mirrors (x = 0, y = 0) or both diagonal mirrors (x = y, x = −y) constrain each search. All displayed fillings have the axial pair. The completed diagonal searches yield no candidates within the stated scope.'),
('Search limits', 'The six 9 → 1 boundaries are completely enumerated through 28 hexes in both orientations. The ten 16 → 4 boundaries are completely enumerated through 28 hexes in both orientations and through 36 with diagonal mirrors; some axial cap-36 searches are incomplete. The twelve 25 → 9 boundaries were tested at cap 36, with higher-cap trials on selected walls; some runs are incomplete. The six-quad 25 → 9 wall has a complete axial cap-36 search, yielding the 33-cell connectivity. Per-wall limits are listed below. These are fixed-boundary results, not unrestricted minimum cell counts. <a href="assets/refinement/searches.json">Search records</a> identify completed and time-limited runs.'),
('Side-wall search', 'Up to eight quads for 9 → 1 and 16 → 4, and ten for 25 → 9, with zero or one additional vertex on each vertical edge and left–right reflection. Top/bottom edges have 3/1, 4/2, or 5/3 segments. Reflection-orbit enumeration agrees with single-quad enumeration. Top–bottom flips are distinct because the cap subdivisions differ.'),
('Embedding', 'HexOpt fixes the entire boundary and preserves the axial mirrors. All coordinates use the cube [−1,1]³. Exact cell positivity, boundary matching, topology, volume, and symmetry are checked independently. Sampled quality describes the stored embedding.'),
('Literature', '<a href="https://www.academia.edu/75511990/Refining_Quadrilateral_and_Hexahedral_Element_Meshes">Schneiders, <em>Refining Quadrilateral and Hexahedral Element Meshes</em> (1996)</a>, gives directly relevant transitions in figures 10 and 16a/b; see the <a href="#literature">comparison below</a>. Other template sets include <a href="https://doi.org/10.1016/j.proeng.2015.10.120">Owen and Shih’s two-refinement templates</a> and <a href="https://arxiv.org/abs/2512.14862">Tong and Zhang’s three-refinement templates</a>. The solver finds the meshes here from their boundaries; that does not establish new connectivity or literature-wide optimality.'),
])}
<section class="collection" id="comparison"><h2>Cell count and element quality</h2>
<p>Higher minimum sampled SJ is better; lower maximum condition number is better. The 14-cell 9 → 1 mesh improves SJ over the 13-cell mesh but has a worse condition number. The 34-cell 16 → 4 mesh is included for its seven-quad side wall; the 32-cell mesh has better quality and fewer cells, with eight quads per wall.</p>
<div style="overflow-x:auto"><table class="metrics"><thead><tr><th>Mesh</th><th>Transition</th><th>Hexes</th><th>Quads / wall</th><th>Min. SJ ↑</th><th>Max. condition ↓</th><th>Downloads</th></tr></thead><tbody>{rows}</tbody></table></div></section>
<section class="collection" id="literature"><h2>Relation to Schneiders’ refinement templates</h2>
<p><a href="?cells=28#F16-4-H28">F16-4-H28</a> divides along its axial mirrors into four seven-cell quarters, consistent with the four-block face transition in Schneiders’ <a href="https://figures.academia-assets.com/83531426/figure_015.jpg">figures 16a/b</a>. This is the closest displayed match to that construction.</p>
<p><a href="?cells=13#F9-1-H13">F9-1-H13</a> is comparable to the transition beneath the regular fine-grid layer in <a href="https://figures.academia-assets.com/83531426/figure_010.jpg">figure 10</a>. The paper calls this 3-refinement by its edge ratio; the collection uses 9 → 1 to count quads on the two caps. These are figure-based comparisons; exact connectivity equivalence to a published mesh file has not been checked. The displayed coordinates are optimized here.</p></section>
<section class="collection" id="templates"><h2>Refinement meshes</h2>
<p class="collection-description">Quality is sampled on 21³ points per cell; exact positivity is certified separately. Drag to rotate, scroll to zoom, and use shrink to inspect cell interfaces. <a href="methodology.html#quality">Quality definitions</a>.</p>
<div class="collection-toolbar"><div class="count-tabs" role="group" aria-label="Filter by cell count"><button class="active" data-count="all" aria-pressed="true">All <span>{len(templates)}</span></button>{tabs}</div>
<label class="sort-control">Sort by <select id="sort"><option value="quality">Minimum SJ (descending)</option><option value="condition">Condition number (ascending)</option><option value="id">Mesh ID</option></select></label></div>
<p id="collection-status" class="collection-status" role="status">Loading refinement meshes…</p><div id="gallery"></div>
<noscript><p>Interactive views require JavaScript; downloads are available in the comparison table above.</p></noscript></section>
<section class="collection" id="side-walls"><h2>Side-wall quadrangulations</h2>
<p>All {len(cases)} patterns from the 2D search, grouped by transition and quad count. Within each group, the best minimum corner mean ratio comes first. The score is <code>2 det(u,v) / (|u|² + |v|²)</code> for outgoing corner edges u and v, minimized over all corners; it penalizes skew and stretching. These are the actual square side walls of [−1,1]³. A valid wall alone does not establish a valid volume filling.</p>
{''.join(walls)}
<p><a href="assets/refinement/side-catalog.json" download>Side-wall catalog</a> · <a href="assets/refinement/enumeration.json" download>2D enumeration record</a></p></section>
<section class="collection"><h2>Running the searches</h2>
<p><code>THREADS=8 bash run-refinement.sh</code> searches the six displayed boundaries at caps 13, 14, 28, 32, 33, and 34, then runs HexOpt and exact validation. It starts from the boundary alone. Optimized coordinates can differ from the displayed embeddings. See the repository README for output folders and optional settings.</p></section>
</main><footer><span>Refinement templates</span><a href="assets/refinement/templates.zip" download>Mesh archive</a></footer></body></html>'''


def build():
    cases = json.loads((INPUT/'cases.json').read_text())
    lookup = {c['id']: c for c in cases}
    templates = json.loads((OUT/'catalog.json').read_text())['templates']
    for folder in ['certificates', 'connectivities', 'vtu', 'sides', 'boundaries']:
        (OUT/folder).mkdir(parents=True, exist_ok=True)
    for i,t in enumerate(templates):
        path = DOCS/t['mesh']; case = lookup[t['boundary_case']]
        p,h,q = read_mesh(path); report = validate(path, case)
        if not report['accepted']:
            raise ValueError(f'Invalid refinement mesh: {t["id"]}')
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        files = {key:f'assets/refinement/{folder}/{t["id"]}.{ext}' for key,folder,ext in
                 [('certificate','certificates','json'),('connectivity','connectivities','json'),('vtu','vtu','vtu')]}
        report.update(mesh_sha256=digest, id=t['id'], boundary_case=case['id'])
        dump(DOCS/files['certificate'], report)
        dump(DOCS/files['connectivity'], dict(points=p.tolist(), cells=h.tolist(), boundary_quads=q.tolist(), index_base=0, mesh_sha256=digest))
        write_vtu(DOCS/files['vtu'], p, h)
        t.update(candidate_id=i, hexes=len(h), vertices=len(p), interior_vertices=len(p)-len(set(q.flat)),
            coarse_quads=case['bottom_grid']**2, fine_quads=case['top_grid']**2, side_quads=case['side_quads'],
            origin=f'{case["top_grid"]**2} → {case["bottom_grid"]**2} refinement · {case["side_quads"]} quads per side wall',
            sha256=digest, symmetry=report['symmetry'], metrics=dense(p,h,True),
            camera_target=[0,0,0], camera_position=[4,-5.5,4], plane_range=[-1,1], plane_radius=1.5, **files)
        print(t['id'], 'certified', t['metrics']['min_scaled_jacobian'], flush=True)
    templates.sort(key=lambda t: (t['hexes'], t['id']))
    catalog = dict(mesh_count=len(templates), partial_enumeration_counts=[], count_notes={}, templates=templates)
    enumeration = json.loads((INPUT/'enumeration.json').read_text())
    dump(OUT/'catalog.json', catalog); dump(OUT/'side-catalog.json', dict(max_quads=enumeration['max_quads'], scopes=enumeration['scopes'], patterns=cases))
    shutil.copyfile(INPUT/'enumeration.json', OUT/'enumeration.json')
    for c in cases:
        record = dict(id=c['id'], quads=c['side']['quads'], square_points_exact=c['side_points_exact'],
                      side_vertices=c['vertical_interior_vertices']+2)
        svg = diagram(record).replace('>2 edges</text>', f'>{c["top_grid"]} fine edges</text>', 1)
        svg = svg.replace('>2 edges</text>', f'>{c["bottom_grid"]} coarse edges</text>', 1)
        (OUT/'sides'/f'{c["id"]}.svg').write_text(svg)
        dump(OUT/'sides'/f'{c["id"]}.json', c)
        shutil.copyfile(INPUT/f'{c["id"]}.mesh', OUT/'boundaries'/f'{c["id"]}.mesh')
    with (OUT/'quality.csv').open('w') as stream:
        writer=csv.writer(stream);writer.writerow(['id','fine_quads','coarse_quads','hexes','side_quads','min_sampled_sj','max_sampled_condition'])
        for t in templates:writer.writerow([t[k] for k in ['id','fine_quads','coarse_quads','hexes','side_quads']]+[t['metrics']['min_scaled_jacobian'],t['metrics']['max_condition']])
    (OUT/'README.md').write_text(f'''# Refinement templates

{len(templates)} certified hexahedral meshes for 9-to-1, 16-to-4 and 25-to-9 regular quad-grid
transitions, with {len(cases)} enumerated side-wall patterns. Coordinates use
[-1,1]³; the fine grid is on top. Four identical walls match by translation.
All displayed interiors have two axial mirrors. JSON connectivity indices
start at zero; Medit mesh indices start at one.

The catalog and mesh files are the authoritative geometry inputs. Certificates
bind exact Jacobian and topology checks to mesh hashes. Quality measurements
use 21³ samples per cell and are not certified extrema. The 14- and 34-cell
meshes have relatively poor condition numbers despite their positive Jacobians.

Side JSON files contain the quads and reflection under `side`, and rational
coordinates under `side_points_exact`. Boundary-only Medit files are supplied
for all {len(cases)} cases. Enumeration and 3D search records state the finite
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
''')
    shutil.copyfile(ROOT/'LICENSE', OUT/'LICENSE.txt')
    (DOCS/'refinement.html').write_text(page(catalog,cases,json.loads((OUT/'searches.json').read_text())))
    with zipfile.ZipFile(OUT/'templates.zip','w',zipfile.ZIP_DEFLATED) as archive:
        for path in sorted(OUT.rglob('*')):
            if path.is_file() and path.suffix!='.zip':archive.write(path,path.relative_to(OUT))


if __name__=='__main__':build()
