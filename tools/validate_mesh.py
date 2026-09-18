"""Validate the exported pyramid, including exact whole-cell determinant bounds."""
from mesh_tools import *
from mesh_homology import homology

def connected(nodes,adj):
    if not nodes:return False
    seen=set();todo=[next(iter(nodes))]
    while todo:
        v=todo.pop()
        if v in seen:continue
        seen.add(v);todo.extend(adj[v]-seen)
    return seen==set(nodes)

def manifold(h):
    fs=owners(h);b=boundary(h);bv=set(b.flat)
    es=defaultdict(list)
    for i,c in enumerate(h):
        for e in c[EDGES]:es[tuple(sorted(e))].append(i)
    be={tuple(sorted((f[k],f[(k+1)%4]))) for f in b for k in range(4)}
    bad_edges=[]
    for e,cs in es.items():
        adjacent=defaultdict(set);degree=Counter();nodes=set()
        for i in cs:
            ff=[f for f in fs if set(e)<=set(f) and i in [a[0] for a in fs[f]]]
            assert len(ff)==2
            a,bf=ff;nodes.update(ff);adjacent[a].add(bf);adjacent[bf].add(a);degree.update(ff)
        expected=Counter({2:len(nodes)-2,1:2}) if e in be else Counter({2:len(nodes)})
        expected=+expected
        if Counter(degree.values())!=expected or not connected(nodes,adjacent):bad_edges.append(e)
    bad_vertices=[]
    for v in set(h.flat):
        cells={i for i,c in enumerate(h) if v in c};edges=[e for e in es if v in e];faces=[f for f in fs if v in f]
        euler=len(edges)-len(faces)+len(cells);adj=defaultdict(set)
        for f in faces:
            for a,bf in combinations([x[0] for x in fs[f]],2):adj[a].add(bf);adj[bf].add(a)
        if euler!=(1 if v in bv else 2) or not connected(cells,adj):bad_vertices.append(int(v))
    return dict(edge_links_valid=not bad_edges,vertex_links_valid=not bad_vertices,bad_edges=bad_edges,bad_vertices=bad_vertices)

def validate(path,profile='prescribed'):
    if profile not in ['prescribed','published44']:raise ValueError('Unknown boundary profile')
    published=profile=='published44'
    p,h,q=read_mesh(path,infer_boundary=published);a=audit(path,infer_boundary=published);a['manifold']=manifold(h);a['homology']=homology(h)
    boundary_tolerance=5e-6 if published else 1e-13
    volume_tolerance=1e-5 if published else 1e-12
    original_points=p
    if published:
        if len(h)!=44 or len(p)!=61:raise ValueError('Published-44 profile requires 44 cells and 61 vertices')
        # A rigid coordinate permutation for boundary comparison only.
        # Determinant certificates below use the original file coordinates.
        p=p[:,[0,2,1]]*np.array([1.,-1.,1.])
        a['boundary_profile']=dict(name=profile,coordinate_transform=[[1,0,0],[0,0,-1],[0,1,0]],boundary_tolerance=boundary_tolerance,volume_tolerance=volume_tolerance,inferred_boundary=True)
    height=np.sqrt(2);side_counts=[];b=boundary(h);points=p[b]
    base=np.all(np.abs(points[:,:,2])<boundary_tolerance,axis=1)
    for dim,sgn in [(0,1),(0,-1),(1,1),(1,-1)]:
        side_counts.append(int(np.sum(np.all(np.abs(sgn*points[:,:,dim]+points[:,:,2]/height-1)<boundary_tolerance,axis=1))))
    bpoints=p[sorted(set(q.flat))]
    expected=np.array([[x,y,0.] for x in [-1.,0.,1.] for y in [-1.,0.,1.]]+[[0.,0.,height]]+[[x/2,y/2,height/2] for x in [-1.,1.] for y in [-1.,1.]]+[[x,y,height/3] for x,y in [(2/3,0.),(-2/3,0.),(0.,2/3),(0.,-2/3)]])
    err=max(min(np.linalg.norm(v-bpoints,axis=1)) for v in expected)
    volume=sum(float(sum(bernstein(p[c]).flat))*8/27 for c in h)
    a['pyramid']=dict(base_quads=int(sum(base)),side_quads=side_counts,boundary_vertices=len(bpoints),max_prescribed_vertex_error=float(err),max_outside_halfspace=float(max(np.max(np.abs(p[:,:2])+p[:,2,None]/height-1),np.max(-p[:,2]))),integrated_volume=volume,expected_volume=4*height/3,volume_error=volume-4*height/3)
    coeff=[]
    for i,c in enumerate(h):
        bb=bernstein(original_points[c],True);lo=min(bb.flat)
        coeff.append(dict(cell=i,all_27_bernstein_coefficients_positive=bool(lo>0),min_coefficient_exact=str(lo),min_coefficient_float=float(lo)))
    a['exact_bernstein_bounds']=coeff
    a['accepted']=bool(a['homology']['betti_GF2']==[1,0,0,0] and a['all_jacobians_positive'] and a['boundary_matches'] and a['opposite_internal_orientations'] and not a['incompatible_pairs'] and a['distinct_cell_vertices'] and a['euler']==1 and a['manifold']['edge_links_valid'] and a['manifold']['vertex_links_valid'] and a['pyramid']['base_quads']==4 and side_counts==[3]*4 and err<boundary_tolerance and a['pyramid']['max_outside_halfspace']<boundary_tolerance and abs(a['pyramid']['volume_error'])<volume_tolerance)
    return a
if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser();parser.add_argument('mesh');parser.add_argument('output',nargs='?');parser.add_argument('--profile',choices=['prescribed','published44'],default='prescribed')
    args=parser.parse_args();a=validate(args.mesh,args.profile);dest=args.output or args.mesh+'.validation.json'
    Path(dest).write_text(json.dumps(a,indent=2)+'\n');print(json.dumps({k:v for k,v in a.items() if k not in ['exact_bernstein_bounds','jacobian_certificates']},indent=2));raise SystemExit(0 if a['accepted'] else 1)
