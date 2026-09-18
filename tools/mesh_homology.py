"""Cellular homology over GF(2), a necessary topology check for a solid pyramid."""
from mesh_tools import *
def rank2(columns):
    pivots={}
    for col in columns:
        while col:
            p=col.bit_length()-1
            if p in pivots:col^=pivots[p]
            else:pivots[p]=col;break
    return len(pivots)
def homology(h):
    verts=sorted(set(h.flat));vi={v:i for i,v in enumerate(verts)}
    ee=sorted({tuple(sorted(e)) for c in h for e in c[EDGES]});ei={e:i for i,e in enumerate(ee)}
    ff=owners(h);fi={f:i for i,f in enumerate(ff)}
    d1=[(1<<vi[a])^(1<<vi[b]) for a,b in ee]
    d2=[]
    for key,own in ff.items():
        f=own[0][1];col=0
        for k in range(4):col^=1<<ei[tuple(sorted((f[k],f[(k+1)%4]))) ]
        d2.append(col)
    d3=[]
    for c in h:
        col=0
        for f in c[FACES]:col^=1<<fi[tuple(sorted(f))]
        d3.append(col)
    r1,r2,r3=rank2(d1),rank2(d2),rank2(d3)
    return dict(betti_GF2=[len(verts)-r1,len(ee)-r1-r2,len(ff)-r2-r3,len(h)-r3],cell_counts=[len(verts),len(ee),len(ff),len(h)],boundary_ranks=[r1,r2,r3])
if __name__=='__main__':
    import sys
    for name in sys.argv[1:]:
        p,h,q=read_mesh(name);print(name,homology(h))
