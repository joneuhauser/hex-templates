#!/usr/bin/env python3
"""Enumerate, identify vertical flips, certify and publish Geode side disks."""
import argparse
from collections import Counter, defaultdict, deque
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import subprocess
import time
import zipfile

from site_layout import header, search_summary

ROOT = Path(__file__).resolve().parents[1]


def canonical(quads, boundary):
    """Oriented disk traversal, fixing the boundary labels pointwise."""
    directed = {}
    for i, face in enumerate(quads):
        q = tuple(face)
        for j in range(4):
            rotated = q[j:] + q[:j]
            directed[rotated[:2]] = (i, rotated)
    labels = {v: v for v in range(boundary)}
    seen, result, queue = set(), [], deque([(0, 1)])
    while queue:
        found = directed.get(queue.popleft())
        if found is None or found[0] in seen:
            continue
        i, q = found
        seen.add(i)
        for v in q:
            if v not in labels:
                labels[v] = len(labels)
        result.append(tuple(labels[v] for v in q))
        queue.extend((q[(j+1)%4], q[j]) for j in range(1, 4))
    if len(seen) != len(quads):
        raise ValueError('Disconnected quad adjacency')
    return tuple(result)


def flipped(quads, boundary):
    r = (boundary-6)//2
    def image(v):
        return (r+5-v) % boundary if v < boundary else v
    return canonical([[image(q[i]) for i in [0, 3, 2, 1]] for q in quads], boundary)


def boundary_points(r):
    return [(F(-1), F(-1)), (F(0), F(-1)), (F(1), F(-1))] + [
        (F(1), -1+F(2*j, r+1)) for j in range(1, r+1)] + [
        (F(1), F(1)), (F(0), F(1)), (F(-1), F(1))] + [
        (F(-1), 1-F(2*j, r+1)) for j in range(1, r+1)]


def harmonic(quads, r, boundary=None):
    """Solve the unweighted graph Laplacian in exact rational arithmetic."""
    points = boundary_points(r) if boundary is None else list(boundary)
    nb, n = len(points), 1+max(v for q in quads for v in q)
    adjacent = defaultdict(set)
    for q in quads:
        for j in range(4):
            a, b = q[j], q[(j+1)%4]
            adjacent[a].add(b)
            adjacent[b].add(a)
    matrix = []
    for v in range(nb, n):
        row = [F(len(adjacent[v]) if v == w else -int(w in adjacent[v])) for w in range(nb, n)]
        row += [sum((points[w][axis] for w in adjacent[v] if w < nb), F(0)) for axis in range(2)]
        matrix.append(row)
    for j in range(n-nb):
        pivot = next(i for i in range(j, n-nb) if matrix[i][j])
        matrix[j], matrix[pivot] = matrix[pivot], matrix[j]
        scale = matrix[j][j]
        matrix[j] = [x/scale for x in matrix[j]]
        for i in range(n-nb):
            if i != j:
                scale = matrix[i][j]
                matrix[i] = [x-scale*y for x, y in zip(matrix[i], matrix[j])]
    return points + [tuple(row[-2:]) for row in matrix]


def cross(a, b):
    return a[0]*b[1]-a[1]*b[0]


def certify(points, quads, r):
    """Independent exact straight-line, incidence and symmetry checks."""
    nb = 6+2*r
    if len(set(points)) != len(points) or len(points) != len(quads)+r+4:
        raise ValueError('Coincident vertices or incorrect Euler vertex count')
    if points[:nb] != boundary_points(r):
        raise ValueError('Boundary coordinates changed')
    edges, oriented, determinants = Counter(), Counter(), []
    area = F(0)
    for q in quads:
        if len(set(q)) != 4:
            raise ValueError('Repeated quad vertex')
        for j in range(4):
            a, b, c = points[q[j]], points[q[(j+1)%4]], points[q[(j-1)%4]]
            det = cross((b[0]-a[0], b[1]-a[1]), (c[0]-a[0], c[1]-a[1]))
            if det <= 0:
                raise ValueError('Nonpositive quad corner determinant')
            determinants.append(det)
            u, v = q[j], q[(j+1)%4]
            edges[tuple(sorted((u, v)))] += 1
            oriented[u, v] += 1
            area += cross(a, b)/2
    rim = {tuple(sorted((v, (v+1)%nb))) for v in range(nb)}
    if {e for e, count in edges.items() if count == 1} != rim or any(count not in [1, 2] for count in edges.values()):
        raise ValueError('Incorrect boundary or face incidence')
    if any(oriented[v, (v+1)%nb] != 1 for v in range(nb)):
        raise ValueError('Boundary orientation disagrees with the square')
    for a, b in edges:
        if (a, b) not in rim and (oriented[a, b], oriented[b, a]) != (1, 1):
            raise ValueError('Internal edge orientations disagree')
    if area != 4 or len(points)-len(edges)+len(quads) != 1:
        raise ValueError('Incorrect area or Euler characteristic')
    # Distinct straight edges cannot cross or touch except at a shared endpoint.
    es = list(edges)
    for i, (a, b) in enumerate(es):
        for c, d in es[i+1:]:
            if len({a, b, c, d}) != 4:
                continue
            def orient(a, b, c):
                x, y, z = points[a], points[b], points[c]
                return cross((y[0]-x[0], y[1]-x[1]), (z[0]-x[0], z[1]-x[1]))
            if orient(a, b, c)*orient(a, b, d) <= 0 and orient(c, d, a)*orient(c, d, b) <= 0:
                if (max(min(points[a][0], points[b][0]), min(points[c][0], points[d][0])) <= min(max(points[a][0], points[b][0]), max(points[c][0], points[d][0]))
                    and max(min(points[a][1], points[b][1]), min(points[c][1], points[d][1])) <= min(max(points[a][1], points[b][1]), max(points[c][1], points[d][1]))):
                    raise ValueError('Disjoint edges intersect')
    lookup = {p: i for i, p in enumerate(points)}
    mirror = [lookup[(-x, y)] for x, y in points]
    key = lambda f: min(tuple(f[j:]+f[:j]) for j in range(4))
    if {key(list(q)) for q in quads} != {key([mirror[q[i]] for i in [0, 3, 2, 1]]) for q in quads}:
        raise ValueError('Reflection does not preserve oriented quads')
    return dict(arithmetic='exact rational', minimum_corner_determinant=str(min(determinants)),
                area=str(area), convex=True, nonintersecting_edges=True, reflection=mirror)


def published_side_key():
    from mesh_tools import read_mesh
    from prepare_geode import boundary_contract
    p, _, q = read_mesh(ROOT/'docs/assets/geode/meshes/G26-Mitchell.mesh')
    contract = boundary_contract(p, q)
    side = q[contract['patch_face_ids']['east']]
    ids = sorted(set(side.flat))
    labels = {}
    for v, (x, y) in enumerate(boundary_points(0)):
        matches = [w for w in ids if abs(p[w, 1]-float(x)) < 1e-12 and abs(p[w, 2]-float(y*F(7, 12))) < 1e-12]
        if len(matches) != 1:
            raise ValueError('Published side rim changed')
        labels[matches[0]] = v
    for w in ids:
        if w not in labels:
            labels[w] = len(labels)
    key = canonical([[labels[v] for v in face] for face in side], 6)
    return min(key, flipped(key, 6))


def diagram(record):
    points = [(float(F(x)), float(F(y))) for x, y in record['square_points_exact']]
    xy = lambda v: (140+105*points[v][0], 137-105*points[v][1])
    svg = ['<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 280 278" role="img" '
           f'aria-label="{record["id"]}: {len(record["quads"])} quads, {record["side_vertices"]} vertices on each vertical edge">']
    for q in record['quads']:
        svg.append('<polygon points="'+' '.join(f'{xy(v)[0]:.6f},{xy(v)[1]:.6f}' for v in q)+'" fill="#e7eff4" stroke="none"/>')
    edges = sorted({tuple(sorted((q[j], q[(j+1)%4]))) for q in record['quads'] for j in range(4)})
    for a, b in edges:
        svg.append(f'<path d="M {xy(a)[0]:.6f},{xy(a)[1]:.6f} L {xy(b)[0]:.6f},{xy(b)[1]:.6f}" stroke="#245a82" stroke-width="1.8" fill="none"/>')
    svg.append('<path class="reflection-axis" d="M140,20 V254" stroke="#bc7648" stroke-width="1" stroke-dasharray="4 4"/>')
    for v in range(len(points)):
        svg.append(f'<circle cx="{xy(v)[0]:.6f}" cy="{xy(v)[1]:.6f}" r="2.9" fill="#202428"/>')
    svg.append('<text x="140" y="14" text-anchor="middle" fill="#5d646b" font-family="system-ui,sans-serif" font-size="10">2 edges</text>')
    svg.append('<text x="140" y="268" text-anchor="middle" fill="#5d646b" font-family="system-ui,sans-serif" font-size="10">2 edges</text></svg>')
    return ''.join(svg)


def corner_quality(record):
    """Minimum corner mean ratio on the Geode rectangle, computed exactly."""
    points = [(F(x), F(y)*F(7, 12)) for x, y in record['square_points_exact']]
    values = []
    for q in record['quads']:
        for j in range(4):
            a, b, c = (points[q[k % 4]] for k in (j, j+1, j-1))
            u = (b[0]-a[0], b[1]-a[1])
            v = (c[0]-a[0], c[1]-a[1])
            values.append(2*cross(u, v)/sum(t*t for t in (*u, *v)))
    return min(values)


def page(catalog):
    records = catalog['patterns']
    quality = {r['id']: corner_quality(r) for r in records}
    max_quads = catalog['max_quads']
    side_counts = sorted({r['side_vertices'] for r in records})
    counts = Counter(len(r['quads']) for r in records)
    groups = []
    for count in sorted(counts):
        cards = []
        for r in sorted(records, key=lambda r: (-quality[r['id']], r['id'])):
            if len(r['quads']) != count:
                continue
            name = r['id']
            badge = '<span class="published-side">Mitchell side topology</span>' if r['published_side_topology'] else ''
            cards.append(f'<article class="side-card" data-side-vertices="{r["side_vertices"]}" id="{name}">'
                f'<h3>{name}{badge}</h3><img src="assets/geode-sides/{name}.svg" '
                f'loading="lazy" width="280" height="278" alt="{name}: symmetric square mesh with {count} quads">'
                f'<p>{r["side_vertices"]} vertices per vertical edge · '
                f'{r["interior_vertices"]} interior vertices</p>'
                f'<p>Minimum corner mean ratio: <strong>{float(quality[name]):.3f}</strong></p>'
                f'<div class="side-downloads"><a href="assets/geode-sides/{name}.json" download>Side mesh JSON</a>'
                f'<a href="assets/geode-sides/{name}.svg" download>SVG</a></div></article>')
        groups.append(f'<section class="side-group" id="quads-{count}" data-quads="{count}">'
                      f'<h2>{count} quads <span>{counts[count]} '+('pattern' if counts[count] == 1 else 'patterns')+'</span></h2>'
                      '<div class="side-grid">'+''.join(cards)+'</div></section>')
    rows = ''.join('<tr><th>'+str(q)+'</th>'+''.join(f'<td>{sum(len(r["quads"]) == q and r["side_vertices"] == v for r in records)}</td>' for v in side_counts)+f'<td>{counts[q]}</td></tr>' for q in range(2, max_quads+1))
    jumps = ''.join(f'<a href="#quads-{q}">{q} quads <span>{counts[q]}</span></a>' for q in sorted(counts))
    return f'''<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Geode side-wall quadrangulations</title><meta name="description" content="All symmetric square side-wall quad meshes through {max_quads} quads, identifying vertical flips.">
<link rel="stylesheet" href="assets/site.css"><link rel="stylesheet" href="assets/geode-sides.css">
<link rel="icon" href="assets/figures/mark.svg" type="image/svg+xml"><script defer src="assets/geode-sides.js"></script></head>
<body><a class="skip-link" href="#patterns">Skip to patterns</a>{header('geode-sides.html')}
<main><section class="side-intro"><p class="side-eyebrow">GEODE BOUNDARY TEMPLATES</p><h1>Square side-wall quadrangulations</h1>
<p><strong>{len(records)} distinct patterns through {max_quads} quads.</strong> Each has left–right reflection symmetry, exactly three vertices on its top and bottom edges, and matching subdivisions on its vertical edges. Vertically flipped patterns count as the same pattern.</p>
<p class="side-note">These are connectivity types, each shown with one convex symmetric placement. Interior vertices may have any valid valence. Conforming quads share at most one complete edge or one vertex. The orange dashed line marks the required reflection.</p>
<p class="side-note">Within each quad-count group, patterns are ordered by minimum corner mean ratio, best first. This score penalizes skew and stretching: 1 is ideal and values near 0 indicate poor shape. It is measured on the Geode side rectangle, before the figures are rescaled to a square.</p>
<p><a href="assets/geode-sides/patterns.zip" download>Download all patterns (.zip)</a> · <a href="assets/geode-sides/catalog.json" download>Catalog and exact coordinates</a> · <a href="#enumeration">Enumeration details</a></p></section>
{search_summary([
    ('What is searched', 'Conforming quadrangulations of a square disk for use as Geode side walls. The top and bottom edges each have exactly two corners and one midpoint. The two vertical edges have equal numbers of vertices and matching subdivisions.'),
    ('Starting point', 'For each prescribed boundary cycle, the interior starts empty. The enumerator inserts quads and introduces interior vertices; it does not start from an existing side mesh.'),
    ('Symmetry', 'Left–right reflection is required for direct translation matching of assembled side walls. A pattern and its top–bottom reflection count as one connectivity; top–bottom symmetry itself is not required.'),
    ('Search limits', f'Complete through {max_quads} quads, enumerating 0–{max_quads-2} additional vertices on each vertical edge. Reflection-orbit insertion is independently cross-checked by single-quad insertion through {catalog["verification_max_quads"]} quads. <code>run-side-walls.sh</code> uses <code>THREADS</code> to limit concurrent boundary cases.'),
    ('Geometry', 'Every pattern has an exactly checked convex harmonic embedding. A valid side pattern alone does not establish that the assembled Geode boundary has a hexahedral filling.'),
    ('Literature', '<a href="#Q5-01">Q5-01</a> is the side connectivity of <a href="https://www.sandia.gov/files/samitch/files/geode.jou">Mitchell’s published 26-cell Geode</a>. It serves as the literature reference within this enumerated pattern collection.'),
])}
<section id="patterns"><div class="side-toolbar"><nav aria-label="Quad counts">{jumps}</nav><label>Vertices per vertical edge <select id="side-vertices"><option value="all">All</option>{''.join(f'<option value="{v}">{v}</option>' for v in side_counts)}</select></label></div>
<p id="side-status" class="collection-status" role="status">Showing all {len(records)} patterns, ordered by quad count, then best quality first. There are no one- or three-quad patterns.</p>{''.join(groups)}</section>
<section id="enumeration" class="side-method"><h2>Enumeration and checks</h2><p>If each vertical edge has <em>r</em> additional vertices, the boundary has <em>B = 6 + 2r</em> vertices. A disk with <em>Q</em> quads has <em>I = Q − r − 2</em> interior vertices. For <em>Q ≤ {max_quads}</em>, only <em>r = 0,…,{max_quads-2}</em> need to be enumerated. No patterns occur with more than {max(side_counts)} vertices per vertical edge under this quad limit.</p>
<div class="side-table"><table><caption>Distinct patterns after identifying vertical flips; columns give vertices per vertical edge</caption><thead><tr><th>Quads</th>{''.join(f'<th>{v}</th>' for v in side_counts)}<th>Total</th></tr></thead><tbody>{rows}</tbody></table></div>
<p>The completed search through {max_quads} quads inserts reflection orbits of quads and produces {catalog['oriented_count']} patterns before identifying vertical flips. A second search inserts individual quads and tests reflection only on completed disks; it agrees exactly through {catalog['verification_max_quads']} quads. Both modes share the local validity checks. Interior vertex labels are removed by an oriented traversal from a fixed boundary edge. A second canonicalization identifies each disk with its vertical reflection.</p>
<p>Every displayed mesh has an exact rational harmonic embedding. Exact checks verify positive quad corner determinants, edge nonintersection, area, boundary incidence and reflection symmetry. The coordinates in each JSON file use the Geode side rectangle <code>[-1,1] × [-7/12,7/12]</code>; the figures rescale it to a square. These are side meshes, not claims that each assembled boundary has a hexahedral filling.</p>
<p>For outgoing corner edges u and v in counterclockwise order, the mean ratio is <code>2 det(u,v) / (|u|² + |v|²)</code>. The card score is the minimum over all corners of all quads in the stored embedding. Sorting uses exact rational scores, with pattern ID breaking ties; it does not change pattern IDs or claim optimal embeddings.</p>
<p><a href="assets/geode-sides/enumeration.json" download>Enumeration record</a> · <a href="assets/geode-sides/quad_disk_search.cpp" download>Enumerator source</a> · <a href="assets/geode-sides/LICENSE.txt">Code license</a></p></section></main>
<footer><span>Geode side-wall patterns</span><a href="pyramid.html">Pyramid templates</a><a href="assets/geode-sides/patterns.zip" download>Side meshes</a></footer></body></html>
'''


def build(binary, max_quads=12, search_results=None, verification_max_quads=10):
    if not 2 <= max_quads <= 12:
        raise ValueError('max_quads must be between 2 and 12')
    verification_max_quads = min(max_quads, verification_max_quads)
    if verification_max_quads < 2:
        raise ValueError('verification_max_quads must be at least 2')
    evidence = ROOT / 'solver/runs/side-quadrangulations'
    assets = ROOT / 'docs/assets/geode-sides'
    evidence.mkdir(parents=True, exist_ok=True)
    assets.mkdir(parents=True, exist_ok=True)
    patterns, runs, oriented_count = [], [], 0
    published = published_side_key()
    for r in range(max_quads-1):
        outputs = []
        for orbits in ([True, False] if r <= verification_max_quads-2 else [True]):
            search_cap = max_quads if orbits else verification_max_quads
            stem = evidence / f'r{r}-{"orbits" if orbits else "individual"}'
            command = [str(binary), str(search_cap), str(r), str(stem.with_suffix('.jsonl'))]
            if not orbits:
                command.append('--no-orbits')
            started = time.monotonic()
            if search_results is None:
                run = subprocess.run(command, capture_output=True, text=True, check=True)
                log = run.stderr
            else:
                source = search_results/stem.name
                log = source.with_suffix('.log').read_text()
                stem.with_suffix('.jsonl').write_bytes(source.with_suffix('.jsonl').read_bytes())
            marker = (f'QUAD_SEARCH_FINISHED max_quads={search_cap} '
                      f'vertical_interior_vertices={r} symmetry_orbits={int(orbits)} ')
            if not log.startswith(marker):
                raise RuntimeError('Enumeration did not finish')
            stem.with_suffix('.log').write_text(log)
            data = [json.loads(line) for line in stem.with_suffix('.jsonl').read_text().splitlines()]
            keys = {canonical(row['quads'], row['boundary_vertices']) for row in data}
            outputs.append(keys)
            runs.append(dict(max_quads=search_cap, vertical_interior_vertices=r, symmetry_orbits=orbits, completed=True,
                             candidates=len(keys), log=log.strip(),
                             results_sha256=hashlib.sha256(stem.with_suffix('.jsonl').read_bytes()).hexdigest()))
            print(f'r={r} orbits={int(orbits)}: {len(keys)} candidates '
                  f'({time.monotonic()-started:.2f}s)', flush=True)
        if len(outputs) == 2 and {k for k in outputs[0] if len(k) <= verification_max_quads} != outputs[1]:
            raise RuntimeError('Orbit and individual-cell enumerations disagree')
        oriented_count += len(outputs[0])
        groups = defaultdict(list)
        nb = 6+2*r
        for key in outputs[0]:
            flip = flipped(key, nb)
            if flip not in outputs[0]:
                raise RuntimeError('Enumeration is missing a vertical reflection')
            groups[min(key, flip)].append(key)
        for key, members in sorted(groups.items()):
            points = harmonic(key, r)
            certificate = certify(points, key, r)
            patterns.append(dict(quads=key, points=[[float(x), float(y*F(7, 12))] for x, y in points],
                square_points_exact=[[str(x), str(y)] for x, y in points], boundary_vertices=nb,
                side_vertices=r+2, interior_vertices=len(points)-nb, mirror=certificate['reflection'],
                vertical_flip_class_size=len(members), published_side_topology=r == 0 and key == published,
                certificate=certificate))
    patterns.sort(key=lambda a: (len(a['quads']), a['side_vertices'], a['quads']))
    counters = Counter()
    for record in patterns:
        count = len(record['quads'])
        counters[count] += 1
        record['id'] = f'Q{count}-{counters[count]:02d}'
        (assets / f'{record["id"]}.json').write_text(json.dumps(record, indent=2)+'\n')
        (assets / f'{record["id"]}.svg').write_text(diagram(record)+'\n')
    catalog = dict(max_quads=max_quads, verification_max_quads=verification_max_quads,
                   vertical_flips_identified=True, oriented_count=oriented_count,
                   counts=dict(sorted(counters.items())), patterns=patterns)
    record = dict(max_quads=max_quads, verification_max_quads=verification_max_quads,
        runs=runs, oriented_count=oriented_count,
        classes=len(patterns), counts=catalog['counts'],
        source_sha256=hashlib.sha256((ROOT/'solver/quad_disk_search.cpp').read_bytes()).hexdigest())
    for target in [assets, evidence]:
        (target/'catalog.json').write_text(json.dumps(catalog, indent=2)+'\n')
        (target/'enumeration.json').write_text(json.dumps(record, indent=2)+'\n')
    (assets/'quad_disk_search.cpp').write_bytes((ROOT/'solver/quad_disk_search.cpp').read_bytes())
    (assets/'LICENSE.txt').write_bytes((ROOT/'LICENSE').read_bytes())
    with zipfile.ZipFile(assets/'patterns.zip', 'w', compression=zipfile.ZIP_DEFLATED) as archive:
        for record in patterns:
            for suffix in ['json', 'svg']:
                path = assets/f'{record["id"]}.{suffix}'
                archive.write(path, path.name)
        for name in ['catalog.json', 'enumeration.json', 'quad_disk_search.cpp', 'LICENSE.txt']:
            archive.write(assets/name, name)
    (ROOT/'docs/geode-sides.html').write_text(page(catalog))
    print(json.dumps(dict(classes=len(patterns), counts=catalog['counts'], oriented_count=oriented_count)))
    return catalog


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--solver', type=Path, default=ROOT/'solver/build/quad_disk_search')
    parser.add_argument('--max-quads', type=int, default=12, choices=range(2, 13))
    parser.add_argument('--verification-max-quads', type=int, default=10, choices=range(2, 13),
                        help='Cap for the slower individual-quad cross-check (default: 10)')
    parser.add_argument('--search-results', type=Path,
                        help='Reuse completed rN-orbits/individual JSONL and log files from this directory')
    args = parser.parse_args()
    build(args.solver.resolve(), args.max_quads, args.search_results, args.verification_max_quads)
