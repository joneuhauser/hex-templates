"""Vector sketch of the prescribed caps on a unit periodic cube."""
from itertools import product
from pathlib import Path
import numpy as np


def caps(top_quads=2, height=1.):
    if top_quads not in (2, 8) or height <= 0:
        raise ValueError('Use 2 or 8 top quads and positive height')
    points, lookup, faces, offsets, patches = [], {}, [], [], []
    def add(coords, patch):
        indices, shifts = [], []
        for p in coords:
            shift = np.r_[np.floor(p[:2]+1e-9).astype(int), 0]
            base = p-shift
            k = tuple(base)
            if k not in lookup:
                lookup[k] = len(points); points.append(base)
            indices.append(lookup[k]); shifts.append(shift)
        faces.append(indices); offsets.append(shifts); patches.append(patch)
    for i in range(2):
        for j in range(2):
            add(np.array([[i/2, j/2, 0], [i/2, (j+1)/2, 0],
                          [(i+1)/2, (j+1)/2, 0], [(i+1)/2, j/2, 0]]), 'bottom')
    n = 1 if top_quads == 2 else 2
    r = .5/n
    diamond = np.array([[0, -r, 0], [r, 0, 0], [0, r, 0], [-r, 0, 0]])
    for i in range(n):
        for j in range(n):
            for phase in (0, .5):
                add(diamond+[(i+phase)/n, (j+phase)/n, height], 'top')
    return dict(period=[1, 1], height=height, points=np.array(points).tolist(),
                quads=faces, vertex_offsets=np.array(offsets).tolist(), patches=patches,
                half_turn=dict(center=[0, 0], linear_action=[[-1, 0, 0], [0, -1, 0], [0, 0, 1]]),
                domain='T^2 x [0,height]; no physical side boundary',
                note='Boundary specification only, not a hexahedral filling or sphere-solver input.')


def clip_edge(a,b):
    lo,hi=0.,1.;delta=b-a
    for axis in (0,1):
        if abs(delta[axis])<1e-12:
            if a[axis]<-1e-12 or a[axis]>1+1e-12:return None
        else:
            limits=sorted((-a[axis]/delta[axis],(1-a[axis])/delta[axis]))
            lo=max(lo,limits[0]);hi=min(hi,limits[1])
    if hi-lo<1e-10:return None
    return a+lo*delta,a+hi*delta


def sketch(path):
    data=caps(8);p=np.asarray(data['points']);q=np.asarray(data['quads']);o=np.asarray(data['vertex_offsets'])
    def project(v):
        x,y,z=v
        return 72+232*x+88*y,408-110*y-232*z
    def line(a,b,color,width=1.5,dash=''):
        x,y=project(a);xx,yy=project(b)
        return f'<line x1="{x:.2f}" y1="{y:.2f}" x2="{xx:.2f}" y2="{yy:.2f}" stroke="{color}" stroke-width="{width}"'+(f' stroke-dasharray="{dash}"' if dash else '')+'/>'
    body=['<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 450 478" role="img" aria-labelledby="title description">',
          '<title id="title">Rectangular and rotated end grids on a unit cube</title>',
          '<desc id="description">A one by one by one wireframe cube. Four rectangular bottom quads and eight rotated top quads per periodic tile. All four side walls are unmeshed.</desc>',
          '<g fill="none" stroke-linecap="round" stroke-linejoin="round">']
    # Only cube edges connect the two caps; no side-face subdivisions or fill.
    for x,y in product((0,1),repeat=2):
        body.append(line((x,y,0),(x,y,1),'#9ba5ab',1.25,'5 5' if (x,y)==(0,1) else ''))
    for z,color in [(0,'#587885'),(1,'#245a82')]:
        seen=set();body.append(f'<g id="{"bottom" if z==0 else "top"}-grid">')
        for face,offset in zip(q,o):
            corners=p[face]+offset
            if abs(corners[0,2]-z)>1e-9:continue
            for tx,ty in product(range(-2,3),repeat=2):
                points=corners+[tx,ty,0]
                for j in range(4):
                    segment=clip_edge(points[j],points[(j+1)%4])
                    if segment is None:continue
                    key=tuple(sorted(tuple(np.round(v,9)) for v in segment))
                    if key in seen:continue
                    seen.add(key);body.append(line(*segment,color,1.65))
        body.append('</g>')
        corners=[(0,0,z),(1,0,z),(1,1,z),(0,1,z)]
        for j in range(4):body.append(line(corners[j],corners[(j+1)%4],'#526570',1.5))
    body+=['</g>',
        '<g font-family="system-ui, sans-serif" fill="#5d646b" font-size="14" text-anchor="middle">',
        '<text x="232" y="31" fill="#245a82">Rotated top grid</text>',
        '<text x="232" y="51" font-size="12">8 quads per period</text>',
        '<text x="213" y="462" fill="#587885">Rectangular bottom grid · 4 quads</text>',
        '<text x="43" y="298">1</text><text x="188" y="434">1</text><text x="371" y="365">1</text>',
        '</g></svg>']
    Path(path).write_text('\n'.join(body)+'\n')


if __name__=='__main__':
    sketch(Path(__file__).resolve().parents[1]/'docs/assets/periodic/boundary.svg')
