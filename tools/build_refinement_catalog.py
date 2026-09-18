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
from refinement import INPUT, SELECTED, validate
from site_layout import header, search_summary

ROOT = Path(__file__).resolve().parents[1]
DOCS = ROOT/'docs'
OUT = DOCS/'assets/refinement'


def dump(path, data):
    path.write_text(json.dumps(data, indent=2)+'\n')


def height_catalogs(catalog, cases):
    records = json.loads((OUT/'height/optimization.json').read_text())
    lookup = {c['id']: c for c in cases}
    catalogs = {}
    for mode in ['sj', 'condition']:
        templates = []
        for baseline in catalog['templates']:
            record = records['templates'][baseline['id']]
            if record['source_sha256'] != baseline['sha256']:
                raise ValueError('Height comparison baseline changed')
            selection = record[f'best_{mode}']
            ratio = selection['height_ratio']
            stem = OUT/'height'/f'{baseline["id"]}-{mode}'
            p,h,q = read_mesh(stem.with_suffix('.mesh'))
            report = validate(stem.with_suffix('.mesh'), lookup[baseline['boundary_case']], ratio)
            if not report['accepted']:
                raise ValueError(f'Invalid height embedding: {stem.name}')
            digest = hashlib.sha256(stem.with_suffix('.mesh').read_bytes()).hexdigest()
            report.update(mesh_sha256=digest, id=baseline['id'], boundary_case=baseline['boundary_case'])
            dump(stem.with_suffix('.certificate.json'), report)
            dump(stem.with_suffix('.json'), dict(points=p.tolist(), cells=h.tolist(), boundary_quads=q.tolist(), index_base=0, mesh_sha256=digest))
            write_vtu(stem.with_suffix('.vtu'),p,h)
            metrics = dense(p,h,True)
            scale = max(1.,ratio)
            template = dict(baseline, height_ratio=ratio, geometry_mode=mode, sha256=digest,
                metrics=metrics, symmetry=report['symmetry'], plane_range=[-ratio,ratio],
                camera_position=[4*scale,-5.5*scale,4*scale])
            for key, suffix in [('mesh','.mesh'),('certificate','.certificate.json'),('connectivity','.json'),('vtu','.vtu')]:
                template[key] = str(stem.with_suffix(suffix).relative_to(DOCS))
            templates.append(template)
        catalogs[mode] = dict(mesh_count=len(templates), partial_enumeration_counts=[], count_notes={}, templates=templates)
        dump(OUT/'height'/f'{mode}-catalog.json',catalogs[mode])
    with (OUT/'height/quality.csv').open('w') as stream:
        writer=csv.writer(stream)
        writer.writerow(['id','geometry','hexes','height_width_ratio','min_sampled_sj','max_sampled_condition'])
        for mode,data in [('fixed',catalog),*catalogs.items()]:
            for t in data['templates']:
                writer.writerow([t['id'],mode,t['hexes'],t['height_ratio'],t['metrics']['min_scaled_jacobian'],t['metrics']['max_condition']])
    return catalogs


def page(catalog, cases, searches, height):
    templates = catalog['templates']
    counts = Counter(t['hexes'] for t in templates)
    tabs = ''.join(f'<button data-count="{n}" aria-pressed="false">{n} <span>{k}</span></button>' for n,k in sorted(counts.items()))
    rows = ''.join(f'<tr><th><a href="?geometry=fixed&amp;cells={t["hexes"]}#{t["id"]}">{t["id"]}</a></th>'
        f'<td>{t["fine_quads"]} → {t["coarse_quads"]}</td><td>{t["hexes"]}</td>'
        f'<td><a href="#side-{t["boundary_case"]}">{t["side_quads"]}</a></td>'
        f'<td>{t["metrics"]["min_scaled_jacobian"]:.3f}</td><td>{t["metrics"]["max_condition"]:.2f}</td>'
        f'<td><a href="{t["mesh"]}" download>Mesh</a> · <a href="{t["connectivity"]}" download>JSON</a></td></tr>' for t in templates)
    height_rows = []
    for t in templates:
        row = f'<tr><th>{t["id"]}</th><td>{t["hexes"]}</td>'
        for mode in ['sj','condition']:
            v = next(v for v in height[mode]['templates'] if v['id']==t['id'])
            row += (f'<td>{v["height_ratio"]:.3f}</td><td>{v["metrics"]["min_scaled_jacobian"]:.3f}</td>'
                    f'<td>{v["metrics"]["max_condition"]:.2f}</td><td><a href="?geometry={mode}&amp;cells={t["hexes"]}#{t["id"]}">View</a></td>')
        height_rows.append(row+'</tr>')
    walls = []
    for fine,coarse in sorted({(c['top_grid']**2,c['bottom_grid']**2) for c in cases}):
        selected = [c for c in cases if (c['top_grid']**2,c['bottom_grid']**2) == (fine,coarse)]
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
<meta name="description" content="Certified 9-to-1, 16-to-4, 25-to-9 and 25-to-1 hexahedral refinement templates and symmetric side-wall quadrangulations.">
<link rel="stylesheet" href="assets/site.css"><link rel="stylesheet" href="assets/geode-sides.css">
<link rel="icon" href="assets/figures/mark.svg" type="image/svg+xml">
<script type="importmap">{{"imports":{{"three":"https://cdn.jsdelivr.net/npm/three@0.180.0/build/three.module.js","three/addons/":"https://cdn.jsdelivr.net/npm/three@0.180.0/examples/jsm/"}}}}</script>
<script defer src="assets/catalog.js"></script></head>
<body data-catalog="assets/refinement/catalog.json" data-catalog-sj="assets/refinement/height/sj-catalog.json" data-catalog-condition="assets/refinement/height/condition-catalog.json"><a class="skip-link" href="#templates">Skip to meshes</a>
{header('refinement.html')}
<main><section class="introduction"><div><h1>Hexahedral refinement templates</h1>
<p>Transitions from a 3×3 quad grid to one quad (9 → 1), from a 4×4 grid to a 2×2 grid (16 → 4), and from a 5×5 grid to either a 3×3 grid (25 → 9) or one quad (25 → 1). The fine grid is on top and the coarse grid below; the same connectivity can be used in either direction.</p>
<p>{len(templates)} connectivities with selectable height-optimized and fixed-height embeddings, and {len(cases)} enumerated side-wall patterns. Every displayed embedding has strictly positive Jacobian determinants throughout its cells.</p>
<p><a href="assets/refinement/templates.zip" download>Meshes and side walls (.zip)</a> · <a href="#side-walls">Side-wall patterns</a> · <a href="assets/refinement/height/quality.csv" download>Quality comparison (.csv)</a></p>
</div><figure><img src="assets/previews/refinement.svg" alt="The 32-cell 16-to-4 refinement mesh, with cells separated for inspection."><figcaption>32-cell 16 → 4 transition.</figcaption></figure></section>
<aside class="collection-highlight" aria-labelledby="highlight-title"><h2 id="highlight-title">Compact transitions and a quality alternative</h2>
<p><a href="?cells=13#F9-1-H13">F9-1-H13</a>, <a href="?cells=28#F16-4-H28">F16-4-H28</a>, and <a href="?cells=33#F25-9-H33">F25-9-H33</a> are the smallest displayed fillings for those three transitions. The 33-cell 25 → 9 mesh uses six quads per side wall.</p>
<p>Removing the four bottom cells from <a href="?cells=32#F16-4-H32">F16-4-H32</a> gives the H28 connectivity, verified by oriented mesh isomorphism. The <a href="#comparison">height comparison</a> uses the same horizontal footprint and height freedom for every template.</p></aside>
<p><a href="?cells=37#F25-1-H37-02">F25-1-H37-02</a> gives a direct 25 → 1 transition with 37 cells and seven quads per wall. Three distinct 37-cell fillings are displayed, together with three 38-cell alternatives on an eight-quad wall.</p>
{search_summary([
('What is searched', 'Conforming hex fillings of a cube with prescribed regular top and bottom grids. Each run fixes one of the <a href="#side-walls">side-wall patterns below</a>. All four walls are identical and opposite walls match by direct translation.'),
('Starting point', 'The solver starts with the boundary and an empty interior. No existing volume connectivity or lower submesh is prescribed.'),
('Symmetry', 'Both axial mirrors (x = 0, y = 0) or both diagonal mirrors (x = y, x = −y) constrain each search. All displayed fillings have the axial pair. The completed diagonal searches yield no candidates within the stated scope.'),
('Search limits', 'The six 9 → 1 boundaries are completely enumerated through 28 hexes in both orientations. The ten 16 → 4 boundaries are completely enumerated through 28 hexes in both orientations and through 36 with diagonal mirrors; some axial cap-36 searches are incomplete. The 25 → 9 search scope includes twelve boundaries at cap 36 and selected walls through cap 44; some searches are incomplete. The six-quad 25 → 9 wall has a complete axial cap-36 search, yielding the 33-cell connectivity. The 25 → 1 search scope includes seven boundaries at cap 40 with time limits; complete axial searches at caps 37 and 38 give three classes each for the displayed seven- and eight-quad walls. Per-wall limits are listed below. These are fixed-boundary results, not unrestricted minimum cell counts. <a href="assets/refinement/searches.json">Search records</a> identify completed and time-limited runs.'),
('Side-wall search', 'Up to eight quads for 9 → 1 and 16 → 4, and ten for 25 → 9 and 25 → 1, with zero or one additional vertex on each vertical edge and left–right reflection. Top/bottom edges have 3/1, 4/2, 5/3, or 5/1 segments. Reflection-orbit enumeration agrees with single-quad enumeration. Top–bottom flips are distinct because the cap subdivisions differ.'),
('Embedding', 'Searches and fixed-height references use the cube [−1,1]³. Height optimization keeps the horizontal footprint [−1,1]² and the normalized side-wall coordinates fixed, while one shared variable scales all boundary heights. HexOpt also moves interior vertices and preserves both axial mirrors. Exact cell positivity, the scaled boundary, topology, volume, and symmetry are checked independently.'),
('Literature', '<a href="https://www.academia.edu/75511990/Refining_Quadrilateral_and_Hexahedral_Element_Meshes">Schneiders, <em>Refining Quadrilateral and Hexahedral Element Meshes</em> (1996)</a>, gives directly relevant transitions in figures 10 and 16a/b; see the <a href="#literature">comparison below</a>. Other template sets include <a href="https://doi.org/10.1016/j.proeng.2015.10.120">Owen and Shih’s two-refinement templates</a> and <a href="https://arxiv.org/abs/2512.14862">Tong and Zhang’s three-refinement templates</a>. The solver finds the meshes here from their boundaries; that does not establish new connectivity or literature-wide optimality.'),
])}
<section class="collection" id="comparison"><h2>Cell count and element quality</h2>
<p>Layer height is optimized independently for each connectivity. Height/width is the full layer height divided by the footprint width; 1 is the cube reference. Every mesh uses starts 0.125, 0.25, 0.5, 1, 2, and 4, with the height variable bounded to [0.05, 4]. HexOpt optimizes its scaled-Jacobian threshold. The table selects the best dense sampled SJ and, separately, the lowest worst condition number among the same certified trials, including simply scaled reference meshes. The second selection is not a direct condition-number optimization. These are local optimization results, not proven optima.</p>
<p>Higher SJ is better; lower condition number is better. SJ measures angular distortion and can miss stretching, so both are reported. Each column group may select a different height and interior embedding. <a href="assets/refinement/height/quality.csv" download>Comparison CSV</a> · <a href="assets/refinement/height/optimization.json">Optimization records</a>.</p>
<div style="overflow-x:auto"><table class="metrics"><thead><tr><th rowspan="2">Mesh</th><th rowspan="2">Hexes</th><th colspan="4">Best sampled SJ</th><th colspan="4">Lowest condition among trials</th></tr><tr><th>Height/width</th><th>Min. SJ ↑</th><th>Max. condition ↓</th><th>Mesh</th><th>Height/width</th><th>Min. SJ ↑</th><th>Max. condition ↓</th><th>Mesh</th></tr></thead><tbody>{''.join(height_rows)}</tbody></table></div>
<details><summary>Fixed-height reference meshes (height/width = 1)</summary>
<p>Higher minimum sampled SJ is better; lower maximum condition number is better. The 14-cell 9 → 1 mesh improves SJ over the 13-cell mesh but has a worse condition number. The 34-cell 16 → 4 mesh is included for its seven-quad side wall; the 32-cell mesh has better quality and fewer cells, with eight quads per wall.</p>
<p>For 25 → 1, the 38-cell variants improve minimum SJ over the 37-cell variants but have worse maximum condition numbers. Strict positivity alone does not imply good element quality.</p>
<div style="overflow-x:auto"><table class="metrics"><thead><tr><th>Mesh</th><th>Transition</th><th>Hexes</th><th>Quads / wall</th><th>Min. SJ ↑</th><th>Max. condition ↓</th><th>Downloads</th></tr></thead><tbody>{rows}</tbody></table></div></details></section>
<section class="collection" id="literature"><h2>Relation to Schneiders’ refinement templates</h2>
<p><a href="?cells=28#F16-4-H28">F16-4-H28</a> divides along its axial mirrors into four seven-cell quarters, consistent with the four-block face transition in Schneiders’ <a href="https://figures.academia-assets.com/83531426/figure_015.jpg">figures 16a/b</a>. This is the closest displayed match to that construction.</p>
<p><a href="?cells=13#F9-1-H13">F9-1-H13</a> is comparable to the transition beneath the regular fine-grid layer in <a href="https://figures.academia-assets.com/83531426/figure_010.jpg">figure 10</a>. The paper calls this 3-refinement by its edge ratio; the collection uses 9 → 1 to count quads on the two caps. These are figure-based comparisons; exact connectivity equivalence to a published mesh file has not been checked. The displayed coordinates are optimized here.</p></section>
<section class="collection" id="templates"><h2>Refinement meshes</h2>
<label>Geometry <select id="geometry"><option value="sj">Variable height: best sampled SJ</option><option value="condition">Variable height: lowest condition among trials</option><option value="fixed">Fixed height/width = 1</option></select></label>
<p class="collection-description">Quality is sampled on 21³ points per cell; exact positivity is certified separately. Drag to rotate, scroll to zoom, and use shrink to inspect cell interfaces. <a href="methodology.html#quality">Quality definitions</a>.</p>
<div class="collection-toolbar"><div class="count-tabs" role="group" aria-label="Filter by cell count"><button class="active" data-count="all" aria-pressed="true">All <span>{len(templates)}</span></button>{tabs}</div>
<label class="sort-control">Sort by <select id="sort"><option value="quality">Minimum SJ (descending)</option><option value="condition">Condition number (ascending)</option><option value="id">Mesh ID</option></select></label></div>
<p id="collection-status" class="collection-status" role="status">Loading refinement meshes…</p><div id="gallery"></div>
<noscript><p>Interactive views require JavaScript; downloads are available in the comparison table above.</p></noscript></section>
<section class="collection" id="side-walls"><h2>Side-wall quadrangulations</h2>
<p>All {len(cases)} patterns from the 2D search, grouped by transition and quad count. Within each group, the best minimum corner mean ratio comes first. The score is <code>2 det(u,v) / (|u|² + |v|²)</code> for outgoing corner edges u and v, minimized over all corners; it penalizes skew and stretching. These patterns use the square side walls of the height/width = 1 reference; other heights scale their vertical coordinates. A valid wall alone does not establish a valid volume filling.</p>
{''.join(walls)}
<p><a href="assets/refinement/side-catalog.json" download>Side-wall catalog</a> · <a href="assets/refinement/enumeration.json" download>2D enumeration record</a></p></section>
<section class="collection"><h2>Running the searches</h2>
<p><code>THREADS=8 bash run-refinement.sh</code> searches the {len(SELECTED)} displayed boundaries at caps {', '.join(map(str,sorted(set(SELECTED.values()))))}, then runs HexOpt and exact validation for the cube and variable-height embeddings. It starts from the boundary alone. <code>python3 tools/hexopt_refinement.py --threads 4 --output DIRECTORY</code> runs the height comparison from the fixed-height assets. Optimized coordinates can differ between runs. See the repository README for output folders and optional settings.</p></section>
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
            height_ratio=1., geometry_mode='fixed', camera_target=[0,0,0], camera_position=[4,-5.5,4], plane_range=[-1,1], plane_radius=1.5, **files)
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

{len(templates)} certified hexahedral meshes for 9-to-1, 16-to-4, 25-to-9 and 25-to-1 regular quad-grid
transitions, with {len(cases)} enumerated side-wall patterns. Fixed-height references
use [-1,1]³; variable-height meshes use [-1,1]² × [-h,h], where h is height/width.
The fine grid is on top. Four identical walls match by translation.
All displayed interiors have two axial mirrors. JSON connectivity indices
start at zero; Medit mesh indices start at one.

The catalog and mesh files are the authoritative geometry inputs. Certificates
bind exact Jacobian and topology checks to mesh hashes. Quality measurements
use 21³ samples per cell and are not certified extrema.

`height/sj-catalog.json` selects the highest sampled minimum scaled Jacobian;
`height/condition-catalog.json` selects the lowest worst sampled condition number
among the same trials. All meshes use six starts (0.125, 0.25, 0.5, 1, 2, 4),
with height/width bounded to [0.05, 4]. HexOpt varies height and interior coordinates
using its scaled-Jacobian threshold objective. Trials also include simply scaled
references. Condition selection is not direct condition-number optimization;
neither selection proves an optimum. See `height/optimization.json` for records
and `height/quality.csv` for comparisons. Connectivity is unchanged.

Side JSON files contain the quads and reflection under `side`, and rational
coordinates under `side_points_exact`. Boundary-only Medit files are supplied
for all {len(cases)} cases. Enumeration and 3D search records state the finite
bounds and distinguish complete searches from time-limited searches.

From the repository root, `python3 tools/build_refinement_catalog.py` rebuilds
the collection. `THREADS=8 bash run-refinement.sh` finds volume connectivities
from the displayed boundaries, then runs HexOpt and exact validation. Axial searches
also optimize height. `python3 tools/hexopt_refinement.py --threads 4 --output DIRECTORY`
runs height optimization on all fixed-height assets into a fresh directory.
The meshes are independent search results. Schneiders, Refining Quadrilateral
and Hexahedral Element Meshes (1996), figures 10 and 16a/b, gives closely related
3-refinement and 2-refinement constructions. F16-4-H28 has four seven-cell
quarters; F9-1-H13 is comparable to the transition below figure 10's fine-grid
layer. These are figure-based comparisons, not graph-isomorphism checks against
published mesh data. See the collection page for references. No global minimum
or novelty is claimed.
''')
    shutil.copyfile(ROOT/'LICENSE', OUT/'LICENSE.txt')
    height = height_catalogs(catalog,cases)
    (DOCS/'refinement.html').write_text(page(catalog,cases,json.loads((OUT/'searches.json').read_text()),height))
    with zipfile.ZipFile(OUT/'templates.zip','w',zipfile.ZIP_DEFLATED) as archive:
        for path in sorted(OUT.rglob('*')):
            if path.is_file() and path.suffix!='.zip':archive.write(path,path.relative_to(OUT))


if __name__=='__main__':build()
