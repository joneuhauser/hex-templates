"""Build the atlas from the certified .mesh assets; never infer symmetry from run labels."""
from pathlib import Path
import argparse, hashlib, json, zipfile
import numpy as np
from mesh_tools import read_mesh, FACES, EDGES

ROOT=Path(__file__).resolve().parents[1]
TRANSFORMS=[('identity','Identity',[[1,0,0],[0,1,0],[0,0,1]]),('r90','90° about z',[[0,-1,0],[1,0,0],[0,0,1]]),('r180','180° about z',[[-1,0,0],[0,-1,0],[0,0,1]]),('r270','270° about z',[[0,1,0],[-1,0,0],[0,0,1]]),('mx','x = 0',[[-1,0,0],[0,1,0],[0,0,1]]),('my','y = 0',[[1,0,0],[0,-1,0],[0,0,1]]),('md','x = y',[[0,1,0],[1,0,0],[0,0,1]]),('ma','x = −y',[[0,-1,0],[-1,0,0],[0,0,1]])]

def oriented_cell(h):
 return tuple(sorted(min(tuple(q[i:]+q[:i]) for i in range(4)) for q in ([int(h[j]) for j in f] for f in FACES)))

def symmetries(points,cells):
 tolerance=1e-10*max(1.,float(np.ptp(points,axis=0).max()))
 original={oriented_cell(h) for h in cells};cell_ids={tuple(sorted(map(int,h))):i for i,h in enumerate(cells)}
 actions=[]
 for key,label,matrix in TRANSFORMS:
  matrix=np.array(matrix);image=points@matrix.T
  distance=np.linalg.norm(image[:,None,:]-points[None,:,:],axis=2);perm=distance.argmin(axis=1)
  residual=float(distance[np.arange(len(points)),perm].max())
  if residual>tolerance or len(set(perm))!=len(points):continue
  transformed=perm[cells].copy()
  if np.linalg.det(matrix)<0:transformed[:,[1,3]]=transformed[:,[3,1]];transformed[:,[5,7]]=transformed[:,[7,5]]
  if {oriented_cell(h) for h in transformed}!=original:continue
  mapping=[cell_ids[tuple(sorted(map(int,h)))] for h in transformed]
  actions.append(dict(id=key,label=label,kind='mirror' if key.startswith('m') else 'rotation' if key!='identity' else 'identity',matrix=matrix.tolist(),max_residual=residual,cell_permutation=mapping,fixed_cells=sum(i==j for i,j in enumerate(mapping))))
 mirrors=[a for a in actions if a['kind']=='mirror'];rotations=[a for a in actions if a['kind']=='rotation']
 parent=list(range(len(cells)))
 def root(i):
  while parent[i]!=i:i=parent[i]
  return i
 for a in actions:
  for i,j in enumerate(a['cell_permutation']):parent[root(j)]=root(i)
 representatives=sorted({root(i) for i in range(len(cells))});orbit=[representatives.index(root(i)) for i in range(len(cells))]
 group='C4v' if len(actions)==8 else 'C2v' if len(actions)==4 and len(mirrors)==2 else 'Cs' if len(actions)==2 and mirrors else 'C2' if len(actions)==2 else 'C1'
 return dict(group=group,order=len(actions),mirrors=[a['label'] for a in mirrors],rotations=[a['label'] for a in rotations],actions=actions,cell_orbits=orbit,orbit_count=len(representatives),tolerance=tolerance,method='Bijection of coordinates within tolerance, followed by exact oriented-cell connectivity comparison under all eight D4 boundary isometries.')

def refresh():
 path=ROOT/'docs/assets/catalog.json'
 data=json.loads(path.read_text())
 for r in data['templates']:
  mesh=ROOT/'docs'/r['mesh'];p,h,q=read_mesh(mesh,infer_boundary=True)
  assert hashlib.sha256(mesh.read_bytes()).hexdigest()==r['sha256'],r['id']
  cert=json.loads((ROOT/'docs'/r['certificate']).read_text())
  assert cert['accepted'] and cert['all_jacobians_positive'],r['id']
  r['symmetry']=symmetries(p@np.array(r.get('display_matrix',np.eye(3))).T,h)
 data['mesh_count']=len(data['templates'])
 path.write_text(json.dumps(data,indent=2)+'\n')
 with zipfile.ZipFile(ROOT/'docs/assets/templates.zip','w',zipfile.ZIP_DEFLATED) as archive:
  for r in data['templates']:
   for key in ['mesh','certificate']:archive.write(ROOT/'docs'/r[key],r[key])
  for name in ['catalog.json','quality.csv','DATA_SOURCES.md','COPYING-HEXTREME.txt']:archive.write(ROOT/'docs/assets'/name,name)
 print(f"Refreshed {data['mesh_count']} templates and archive")

if __name__=='__main__':
 parser=argparse.ArgumentParser(description=__doc__);parser.parse_args();refresh()
