"""Oriented connectivity equivalence under boundary D4 and interior relabeling."""
from pathlib import Path
from collections import Counter
import sys,itertools,numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'solver'))
from compare_meshes import relations,oriented_cells
from mesh_tools import read_mesh
def mapping(a,b):
 if (a['vertices'],a['hexes'])!=(b['vertices'],b['hexes']):return None
 n=a['vertices'];ar,br=relations(a),relations(b)
 ls=[list(range(18))+[18]*(n-18) for _ in range(2)]
 for _ in range(n):
  signatures=[[(labels[v],tuple(sorted((rel[v][w],labels[w]) for w in range(n) if rel[v][w]))) for v in range(n)] for labels,rel in zip(ls,[ar,br])]
  ids={s:i for i,s in enumerate(sorted(set(signatures[0]+signatures[1])))}
  new=[[ids[s] for s in ss] for ss in signatures]
  if Counter(new[0])!=Counter(new[1]):return None
  stable=len(set(new[0]))==len(set(ls[0]));ls=new
  if stable:break
 choices={v:[w for w in range(n) if ls[0][v]==ls[1][w]] for v in range(n)}
 target=oriented_cells(b,dict(enumerate(range(n))))
 m={v:ws[0] for v,ws in choices.items() if len(ws)==1}
 def solve():
  if len(m)==n:return dict(m) if oriented_cells(a,m)==target else None
  used=set(m.values());avail={v:[w for w in ws if w not in used and all(ar[v][x]==br[w][y] for x,y in m.items())] for v,ws in choices.items() if v not in m}
  v=min(avail,key=lambda x:(len(avail[x]),x))
  for w in avail[v]:
   m[v]=w;r=solve()
   if r is not None:return r
   del m[v]
  return None
 return solve()

bp,_,_=read_mesh(ROOT/'solver/input/pyramid.mesh')
transforms=[]
for swap,sx,sy in itertools.product([0,1],[1,-1],[1,-1]):
 matrix=np.array([[0,sx,0],[sy,0,0],[0,0,1]]) if swap else np.diag([sx,sy,1])
 perm=np.argmin(np.linalg.norm((bp@matrix.T)[:,None,:]-bp[None,:,:],axis=2),axis=1).tolist()
 assert np.max(abs(bp[perm]-bp@matrix.T))<1e-12
 transforms.append((perm,int(round(np.linalg.det(matrix))),matrix))

def transform(a,perm,det):
 perm=perm+list(range(18,a['vertices']));b=dict(a);b['cells']=[[perm[v] for v in h] for h in a['cells']]
 if det<0:
  for h in b['cells']:h[1],h[3]=h[3],h[1];h[5],h[7]=h[7],h[5]
 return b

def match(a,b):
 for perm,det,matrix in transforms:
  m=mapping(transform(a,perm,det),b)
  if m is not None:
   full=perm+list(range(18,a['vertices']))
   return {v:m[full[v]] for v in range(a['vertices'])},matrix
 return None
