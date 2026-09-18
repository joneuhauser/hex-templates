"""Build a static four-collection landing page and data-derived SVG previews."""
from html import escape
import json
from pathlib import Path
import re

import numpy as np

from mesh_tools import read_mesh,FACES
from site_layout import header

ROOT=Path(__file__).resolve().parents[1]
DOCS=ROOT/'docs'
OUT=DOCS/'assets/previews'


def write_svg(name,body,title):
    (OUT/name).write_text(f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 620 310" role="img" aria-label="{escape(title)}"><title>{escape(title)}</title>'+body+'</svg>\n')


def mesh_preview(name,path,title):
    p,h,_=read_mesh(path)
    q=p[h];center=q.mean(axis=1,keepdims=True);q=center+.86*(q-center)
    direction=np.array([3.,-5.,3.]);direction/=np.linalg.norm(direction)
    right=np.cross([0,0,1],direction);right/=np.linalg.norm(right);up=np.cross(direction,right)
    xy=np.stack([q@right,-q@up],axis=-1)
    lo=xy.min(axis=(0,1));hi=xy.max(axis=(0,1));mid=(lo+hi)/2
    scale=min(530/(hi[0]-lo[0]),264/(hi[1]-lo[1]))
    light=np.array([-2.,-3.,6.]);light/=np.linalg.norm(light)
    patches=[]
    for cell in q:
        for face in cell[FACES]:
            face_patches=[]
            def at(u,v):return (1-u)*(1-v)*face[0]+u*(1-v)*face[1]+u*v*face[2]+(1-u)*v*face[3]
            for i in range(3):
                for j in range(3):
                    corners=np.array([at(i/3,j/3),at((i+1)/3,j/3),at((i+1)/3,(j+1)/3),at(i/3,(j+1)/3)])
                    normal=np.cross(corners[1]-corners[0],corners[3]-corners[0]);norm=np.linalg.norm(normal)
                    if norm<1e-15:continue
                    # Opaque bilinear patches, depth sorted for a compact preview.
                    brightness=.82+.18*abs(normal@light/norm)
                    rgb=np.rint(np.array([116,161,177])*brightness+25*max(0,normal@light/norm)).astype(int)
                    color='#'+''.join(f'{min(255,n):02x}' for n in rgb)
                    projected=(np.stack([corners@right,-corners@up],axis=-1)-mid)*scale+[310,155]
                    points=' '.join(f'{x:.2f},{y:.2f}' for x,y in projected)
                    face_patches.append(f'<polygon points="{points}" fill="{color}" stroke="{color}" stroke-width=".45"/>')
            outline=(np.stack([face@right,-face@up],axis=-1)-mid)*scale+[310,155]
            points=' '.join(f'{x:.2f},{y:.2f}' for x,y in outline)
            face_patches.append(f'<polygon points="{points}" fill="none" stroke="#36576a" stroke-width=".75"/>')
            patches.append((float(face.mean(axis=0)@direction),''.join(face_patches)))
    write_svg(name,''.join(body for _,body in sorted(patches,key=lambda x:x[0])),title)


def side_preview(patterns):
    chosen=[next(p for p in patterns if p['id']==name) for name in ['Q2-01','Q5-01','Q7-06']]
    body=''
    for i,pattern in enumerate(chosen):
        p=np.asarray(pattern['points']);lo=p.min(axis=0);hi=p.max(axis=0)
        xy=(p-lo)/(hi-lo)*[156,-176]+[32+i*204,235]
        for face in pattern['quads']:
            points=' '.join(f'{x:.2f},{y:.2f}' for x,y in xy[face])
            body+=f'<polygon points="{points}" fill="#e4eff2" stroke="#46748a" stroke-width="1.5"/>'
        for x,y in xy:body+=f'<circle cx="{x:.2f}" cy="{y:.2f}" r="2.6" fill="#36576a"/>'
        body+=f'<text x="{110+i*204}" y="270" text-anchor="middle" font-family="system-ui,sans-serif" font-size="14" fill="#5d646b">{len(pattern["quads"])} quads</text>'
    write_svg('side-walls.svg',body,'Examples of Geode side-wall patterns with two, five and seven quadrilaterals')


def build():
    OUT.mkdir(parents=True,exist_ok=True)
    pyramid=json.loads((DOCS/'assets/catalog.json').read_text())
    geode=json.loads((DOCS/'assets/geode/catalog.json').read_text())
    sides=json.loads((DOCS/'assets/geode-sides/catalog.json').read_text())
    periodic=json.loads((DOCS/'assets/periodic/current-catalog.json').read_text())
    mesh_preview('pyramid.svg',DOCS/pyramid['templates'][0]['mesh'],'A 36-cell hexahedral filling of Schneiders’ pyramid, with cells separated for inspection')
    mesh_preview('geode.svg',DOCS/geode['templates'][0]['mesh'],'Mitchell’s 26-cell Geode mesh, with cells separated for inspection')
    side_preview(sides['patterns'])
    # Reuse the exact cap sketch; the nested viewport keeps all previews equal.
    sketch=(DOCS/'assets/periodic/boundary.svg').read_text()
    sketch=re.sub(r'<svg\b[^>]*>', '<svg x="163" y="5" width="294" height="300" viewBox="0 0 450 478">',sketch,count=1)
    write_svg('periodic.svg',sketch,'A periodic unit cube with rectangular bottom grid, rotated top grid, and unmeshed side walls')
    cards=[
        ('pyramid.html','pyramid.svg','Schneiders’ pyramid','Hexahedral fillings of the prescribed pyramid boundary, with interactive geometry and quality measurements.',f'{len(pyramid["templates"])} meshes','36–48 and 88 cells','A hexahedral pyramid mesh with separated cells'),
        ('geode.html','geode.svg','Mitchell’s Geode','A hexahedral transition between a quadrilateral grid and a quadrangulated tetrahedral-mesh boundary.',f'{len(geode["templates"])} meshes','26–44 cells','The 26-cell Geode mesh with separated cells'),
        ('geode-sides.html','side-walls.svg','Geode side walls','Quadrilateral boundary patterns for assembling Geode cells with matching sides.',f'{len(sides["patterns"])} patterns',f'Up to {sides["max_quads"]} quads','Three quadrilateral side-wall patterns'),
        ('periodic.html','periodic.svg','Periodic grid rotation','Transitions between rectangular and rotated quad grids, grouped by periodic connectivity.',f'{len(periodic["templates"])} connectivities','38–110 cells per tile','A unit cube with meshed top and bottom and empty side walls')]
    tiles=''.join(f'''<a class="collection-card" href="{url}">
<div class="collection-preview"><img src="assets/previews/{preview}" alt="{alt}" width="620" height="310"></div>
<div class="collection-card-body"><div class="collection-card-title"><h2>{title}</h2><span class="collection-arrow" aria-hidden="true">↗</span></div>
<p>{description}</p><div class="collection-meta"><span>{count}</span><span>{size}</span></div></div></a>''' for url,preview,title,description,count,size,alt in cards)
    home=f'''<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>Hex meshing templates through symmetry-constrained exhaustive search</title><meta name="description" content="Hex meshing templates from symmetry-constrained exhaustive connectivity searches, with optimized geometry, exact validation, and mesh downloads.">
<link rel="stylesheet" href="assets/site.css"><link rel="icon" href="assets/figures/mark.svg" type="image/svg+xml"></head>
<body class="atlas-home"><a class="skip-link" href="#collections">Skip to collections</a>{header('index.html')}
<main><section class="atlas-intro"><p class="eyebrow">Mesh collections</p><h1>Hex meshing templates through symmetry-constrained exhaustive search</h1>
<p>Symmetry constraints make exhaustive connectivity searches tractable within prescribed boundaries and cell-count limits. Explore the resulting templates, optimized geometry, and certified mesh downloads.</p>
</section>
<aside class="site-caution" aria-label="Experimental library caution"><h2>Caution</h2><p>This is an experimental library that I created by discussion with GPT-6 Astra for two days. The generated meshes check out and I'm using them in production for fluid dynamics simulations, but I'm not claiming that they are optimal. Use at your own risk. In particular I do not assert that there are no meshes with smaller element count for the same boundaries.</p></aside>
<section id="collections" class="collection-links" aria-label="Mesh collections">{tiles}</section>
<section class="atlas-method"><div><h2>How the meshes are checked</h2><p>Read about connectivity search, quality measurements, and exact geometric validation.</p></div><a href="methodology.html">Explore the methodology <span aria-hidden="true">→</span></a></section>
</main><footer><span>Hex meshing templates</span><a href="methodology.html">Methodology</a><a href="assets/DATA_SOURCES.md">Data sources</a></footer></body></html>
'''
    (DOCS/'index.html').write_text(home)
    for page in DOCS.glob('*.html'):
        if page.name=='index.html':continue
        text=page.read_text();text=re.sub(r'<header class="site-header">.*?</header>',header(page.name),text,count=1,flags=re.S)
        text=text.replace('<a href="index.html">Pyramid templates</a>','<a href="pyramid.html">Pyramid templates</a>').replace('<a href="index.html">Templates</a>','<a href="pyramid.html">Pyramid meshes</a>')
        page.write_text(text)
    print('HOME_OK: four collection previews; shared navigation on all pages')


if __name__=='__main__':build()
