"""Schneiders pyramid I/O, topology checks, and trilinear Jacobian bounds."""
from collections import Counter, defaultdict
from fractions import Fraction as Q
from itertools import product, combinations
from pathlib import Path
import json
import numpy as np
FACES=np.array([[0,3,2,1],[4,5,6,7],[0,1,5,4],[1,2,6,5],[2,3,7,6],[3,0,4,7]])
SIGNS=np.array([[-1,-1,-1],[1,-1,-1],[1,1,-1],[-1,1,-1],[-1,-1,1],[1,-1,1],[1,1,1],[-1,1,1]])
EDGES=np.array([[0,1],[1,2],[2,3],[3,0],[4,5],[5,6],[6,7],[7,4],[0,4],[1,5],[2,6],[3,7]])
def read_mesh(path,infer_boundary=False):
    lines=Path(path).read_text().splitlines(); out={}
    for name,width,typ in [('Vertices',3,float),('Hexahedra',8,int),('Quadrilaterals',4,int)]:
        i=next((i for i,s in enumerate(lines) if s.strip()==name),None)
        if i is None: out[name]=np.empty((0,width),dtype=typ); continue
        n=int(lines[i+1]); a=np.array([[typ(x) for x in l.split()[:width]] for l in lines[i+2:i+2+n]])
        out[name]=a if typ==float else a-1
    if infer_boundary and not len(out['Quadrilaterals']):out['Quadrilaterals']=boundary(out['Hexahedra'])
    return out['Vertices'],out['Hexahedra'],out['Quadrilaterals']
def write_mesh(path,p,h,q=None):
    if q is None: q=boundary(h)
    with open(path,'w') as f:
        f.write('MeshVersionFormatted 2\nDimension 3\nVertices\n%d\n'%len(p))
        for v in p:f.write(' '.join(format(x,'.17g') for x in v)+' 0\n')
        for name,a in [('Hexahedra',h),('Quadrilaterals',q)]:
            f.write(name+'\n'+str(len(a))+'\n')
            for v in a:f.write(' '.join(str(int(x)+1) for x in v)+' 0\n')
        f.write('End\n')
def write_vtu(path,p,h):
    from xml.etree.ElementTree import Element,SubElement,ElementTree
    root=Element('VTKFile',type='UnstructuredGrid',version='0.1',byte_order='LittleEndian')
    piece=SubElement(SubElement(root,'UnstructuredGrid'),'Piece',NumberOfPoints=str(len(p)),NumberOfCells=str(len(h)))
    a=SubElement(SubElement(piece,'Points'),'DataArray',type='Float64',NumberOfComponents='3',format='ascii');a.text=' '.join(format(x,'.17g') for x in p.flat)
    c=SubElement(piece,'Cells')
    for name,a in [('connectivity',h.ravel()),('offsets',8*np.arange(1,len(h)+1)),('types',np.full(len(h),12))]:
        SubElement(c,'DataArray',type='Int32',Name=name,format='ascii').text=' '.join(map(str,a))
    ElementTree(root).write(path,encoding='utf-8',xml_declaration=True)
def owners(h):
    d=defaultdict(list)
    for i,c in enumerate(h):
        for f in c[FACES]:d[tuple(sorted(f))].append((i,tuple(f)))
    return d
def boundary(h):return np.array([a[0][1] for a in owners(h).values() if len(a)==1])
def topology(h,q):
    d=owners(h); b=boundary(h)
    bad=[]; hs=[set(c) for c in h]
    for i,j in combinations(range(len(h)),2):
        s=hs[i]&hs[j]
        if len(s)<=1:continue
        if len(s)==2 and all(any(set(h[k][e])==s for e in EDGES) for k in (i,j)):continue
        if len(s)==4 and all(any(set(h[k][f])==s for f in FACES) for k in (i,j)):continue
        bad.append([i,j])
    opposite=all(len(a)!=2 or any(a[0][1]==tuple(np.roll(a[1][1][::-1],r)) for r in range(4)) for a in d.values())
    edges={tuple(sorted(e)) for c in h for e in c[EDGES]}; verts=set(h.flat)
    return dict(hexes=len(h),vertices=len(verts),boundary_quads=len(b),boundary_matches={tuple(sorted(f)) for f in b}=={tuple(sorted(f)) for f in q},face_multiplicity=dict(Counter(map(len,d.values()))),opposite_internal_orientations=opposite,incompatible_pairs=bad,distinct_cell_vertices=all(len(set(c))==8 for c in h),euler=len(verts)-len(edges)+len(d)-len(h))
def derivatives(r):
    a=1+r[:,None,:]*SIGNS[None,:,:]
    return np.stack([SIGNS[None,:,k]*np.prod(a[:,:,np.arange(3)!=k],axis=2)/8 for k in range(3)],axis=2)
def sampled(p,h,n=7):
    r=np.array(list(product(np.linspace(-1,1,n),repeat=3)));J=np.einsum('hvi,svj->hsij',p[h],derivatives(r));det=np.linalg.det(J)
    sj=det/np.prod(np.linalg.norm(J,axis=2),axis=2)
    return dict(grid=n,min_det=float(det.min()),min_scaled_jacobian=float(sj.min()),nonpositive_hexes=int(np.sum(det.min(axis=1)<=0)))
def det3(a):
    return a[0,0]*(a[1,1]*a[2,2]-a[1,2]*a[2,1])-a[0,1]*(a[1,0]*a[2,2]-a[1,2]*a[2,0])+a[0,2]*(a[1,0]*a[2,1]-a[1,1]*a[2,0])
def bernstein(p,exact=False):
    """Degree (2,2,2) determinant coefficients, reference cube [-1,1]^3."""
    r=np.array(list(product([-1.,0.,1.],repeat=3)))
    if exact:
        # Treat each stored float as its exact binary rational value.
        p=np.array([[Q(float(x)) for x in row] for row in p],dtype=object)
        ds=np.array([[[Q(float(x)) for x in row] for row in d] for d in derivatives(r)],dtype=object)
        vals=np.array([det3(p.T@d) for d in ds],dtype=object).reshape(3,3,3)
    else:vals=np.linalg.det(np.einsum('vi,svj->sij',p,derivatives(r))).reshape(3,3,3)
    for ax in range(3):
        v=np.moveaxis(vals,ax,0);v[1]=2*v[1]-(v[0]+v[2])/2
    return vals
def split(b,ax):
    v=np.moveaxis(b,ax,0);m=(v[0]+2*v[1]+v[2])/4
    l=np.stack([v[0],(v[0]+v[1])/2,m]);r=np.stack([m,(v[1]+v[2])/2,v[2]])
    return np.moveaxis(l,0,ax),np.moveaxis(r,0,ax)
def certify_hex(p,exact=True,max_depth=24):
    todo=[(bernstein(p,exact),0)];lower=None;leaves=0;depth=0
    while todo:
        b,d=todo.pop();lo=min(b.flat);hi=max(b.flat);depth=max(depth,d)
        if lo>0:
            lower=lo if lower is None else min(lower,lo);leaves+=1;continue
        if hi<=0:return dict(positive=False,reason='nonpositive subbox',depth=depth)
        corners=b[::2,::2,::2]
        if min(corners.flat)<=0:return dict(positive=False,reason='nonpositive point',depth=depth)
        if d>=max_depth:return dict(positive=False,reason='unresolved',depth=depth)
        todo.extend((c,d+1) for c in split(b,d%3))
    return dict(positive=True,lower_bound=float(lower),depth=depth,leaves=leaves)
def audit(path,exact=True,infer_boundary=False):
    p,h,q=read_mesh(path,infer_boundary=infer_boundary);t=topology(h,q);t['sampled']=sampled(p,h)
    t['jacobian_certificates']=[certify_hex(p[c],exact) for c in h]
    t['all_jacobians_positive']=all(c['positive'] for c in t['jacobian_certificates'])
    t['arithmetic']='exact rational from binary64 coordinates' if exact else 'float64'
    return t
if __name__=='__main__':
    import sys
    a=audit(sys.argv[1]);Path(sys.argv[1]+'.audit.json').write_text(json.dumps(a,indent=2));print(json.dumps({k:v for k,v in a.items() if k!='jacobian_certificates'},indent=2))
