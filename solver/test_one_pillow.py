from pathlib import Path
import subprocess,json
p=Path(__file__).resolve().parent
# Import only the mesh reader; the published-orbit test also has a main loop.
ns={};exec((p/'test_boundary_orbits.py').read_text().split('records=[]')[0],dict(__file__=str(p/'test_boundary_orbits.py')),ns)
read=ns['read'];mesh=read(p/'input/one-pillow20-seed.mesh')
out=p/'experiments/pillow-orbits';out.mkdir(parents=True,exist_ok=True)
points=mesh['Vertices'];bv=len(set(v for q in mesh['Quadrilaterals'] for v in q));accepted=0;records=[]
for sx in [1,-1]:
 for sy in [1,-1]:
  transform=lambda q:[sx*q[0],sy*q[1],q[2]]
  perm=[points[:bv].index(transform(q)) for q in points[:bv]]+list(range(bv,len(points)))
  mapped=[None]*len(points)
  for v,q in enumerate(points):mapped[perm[v]]=transform(q)
  cells=[[perm[v] for v in h] for h in mesh['Hexahedra']]
  if sx*sy<0:
   for h in cells:h[1],h[3]=h[3],h[1];h[5],h[7]=h[7],h[5]
  stem=out/f'x{sx}-y{sy}';lines=['MeshVersionFormatted 2','Dimension 3','Vertices',str(len(points))]
  lines+=[' '.join(map(str,q))+' 0' for q in mapped];lines+=['Hexahedra',str(len(cells))];lines+=[' '.join(str(v+1) for v in h)+' 0' for h in cells]
  lines+=['Quadrilaterals',str(len(mesh['Quadrilaterals']))];lines+=[' '.join(str(v+1) for v in q)+' 0' for q in mesh['Quadrilaterals']];lines+=['End'];Path(str(stem)+'.mesh').write_text('\n'.join(lines)+'\n')
  cmd=[str(p/'build/mirror_search'),'--input',str(p/'input/grid2-boundary.mesh'),'--mirrors','1','--planes','axial','--cap','20','--threads','1','--replay',str(stem)+'.mesh','--leader-bound','0','--partial-cover','0','--output',str(stem)+'.jsonl']
  with open(str(stem)+'.log','w') as log:r=subprocess.run(cmd,stderr=log)
  log=Path(str(stem)+'.log').read_text();assert 'boundary_group=4 ' in log,log
  assert r.returncode in [0,3],log
  records.append(dict(sx=sx,sy=sy,exit=r.returncode))
  if r.returncode==0:
   accepted+=1
   for partial in [5,6,7]:
    cmd[cmd.index('--leader-bound')+1]='1';cmd[cmd.index('--partial-cover')+1]=str(partial)
    with open(str(stem)+f'.partial{partial}.log','w') as log:r2=subprocess.run(cmd,stderr=log)
    assert r2.returncode==0,Path(str(stem)+f'.partial{partial}.log').read_text()
assert accepted==2,records
(out/'results.json').write_text(json.dumps(records,indent=2)+'\n')
print('ONE_PILLOW_OK 20 cells; 2 of 4 orientations retained; all partial stages preserve both representatives')
