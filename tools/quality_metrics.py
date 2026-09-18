"""Dense quality measurements shared by study and literature templates."""
from mesh_tools import *
def dense(points,h,condition=False):
 r=np.array(list(product(np.linspace(-1,1,21),repeat=3)))
 J=np.einsum('hvi,svj->hsij',points[h],derivatives(r));det=np.linalg.det(J)
 sj=det/np.prod(np.linalg.norm(J,axis=2),axis=2);cellmin=sj.min(1)
 rr=np.concatenate([SIGNS,np.zeros((1,3))]);JJ=np.einsum('hvi,svj->hsij',points[h],derivatives(rr))
 nine=np.linalg.det(JJ)/np.prod(np.linalg.norm(JJ,axis=2),axis=2)
 result=dict(grid=21,min_scaled_jacobian=float(sj.min()),min_det=float(det.min()),nonpositive_cells=int(np.count_nonzero(det.min(1)<=0)),nine_point_min_sj=float(nine.min()),cell_min_sj=cellmin.tolist(),cell_min_sj_p05=float(np.quantile(cellmin,.05)),cell_min_sj_median=float(np.median(cellmin)))
 if condition:
  singular=np.linalg.svd(J,compute_uv=False);kappa=singular[:,:,0]/singular[:,:,2]
  gauss=np.array(list(product([-1/np.sqrt(3),1/np.sqrt(3)],repeat=3)))
  volume=np.linalg.det(np.einsum('hvi,svj->hsij',points[h],derivatives(gauss))).sum(1)
  edge=points[h[:,EDGES[:,1]]]-points[h[:,EDGES[:,0]]];length=np.linalg.norm(edge,axis=2)
  result.update(max_condition=float(kappa.max()),cell_max_condition=kappa.max(1).tolist(),volume_min=float(volume.min()),volume_max=float(volume.max()),volume_ratio=float(volume.max()/volume.min()),volume_cv=float(volume.std()/volume.mean()),cell_volumes=volume.tolist(),worst_cell_edge_ratio=float((length.max(1)/length.min(1)).max()))
 return result
