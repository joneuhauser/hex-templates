"""Published coordinates retain their bytes and exact determinant requirement."""
from pathlib import Path
import hashlib,json,sys,tempfile
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'tools'))
from mesh_tools import read_mesh,write_mesh
from validate_mesh import validate
from connectivity import match
mesh=ROOT/'docs/assets/meshes/H44-class16.mesh'
assert hashlib.sha256(mesh.read_bytes()).hexdigest()=='8c61a75d882c0935e1e1412cb7deb582cf7ed3e06d774c3e02a7baf0240eccb1'
p,h,q=read_mesh(mesh,infer_boundary=True)
assert len(q)==16
report=validate(mesh,'published44');assert report['accepted'] and report['all_jacobians_positive']
# Identify boundary labels after the viewer's rigid rotation. Connectivity must
# remain the class-16 complex; rounding is used only to identify boundary labels.
bp,_,_=read_mesh(ROOT/'solver/input/pyramid.mesh');rotated=p[:,[0,2,1]]*[1,-1,1]
indices=np.linalg.norm(bp[:,None,:]-rotated[None,:,:],axis=2).argmin(1)
assert len(set(indices))==18
perm=np.empty(len(p),dtype=int);perm[indices]=np.arange(18)
interior=sorted(set(range(len(p)))-set(indices));perm[interior]=np.arange(18,len(p))
a=dict(vertices=len(p),hexes=len(h),cells=perm[h].tolist())
reference=json.loads((ROOT/'tests/fixtures/published44-connectivity.json').read_text())
assert match(a,reference) is not None
with tempfile.TemporaryDirectory(prefix='published44-test-') as tmp:
 path=Path(tmp)/'mesh.mesh'
 # Even after the rigid rotation and explicit boundary extraction, the strict
 # prescribed-boundary profile must reject the publication's rounded values.
 write_mesh(path,rotated,h,q);assert not validate(path)['accepted']
 changed=p.copy();changed[indices[0],0]+=1e-3
 write_mesh(path,changed,h,q);assert not validate(path,'published44')['accepted']
 inverted=h.copy();inverted[:,[1,3]]=inverted[:,[3,1]];inverted[:,[5,7]]=inverted[:,[7,5]]
 write_mesh(path,p,inverted)
 result=validate(path,'published44');assert not result['accepted'] and not result['all_jacobians_positive']
print('PUBLISHED_COORDINATES_OK unchanged source, class-16 isomorphism, strict boundary and inversion rejection')
