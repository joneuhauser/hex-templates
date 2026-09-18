from pathlib import Path
import subprocess,json
p=Path(__file__).resolve().parent;out=p/'experiments/one-orbit-seeds';out.mkdir(parents=True,exist_ok=True)
def read(path):
 t=iter(path.read_text().split());data={}
 for k in t:
  if k in ['MeshVersionFormatted','Dimension']:next(t)
  elif k in ['Vertices','Hexahedra','Quadrilaterals']:
   n=int(next(t));width={'Vertices':3,'Hexahedra':8,'Quadrilaterals':4}[k];rows=[]
   for _ in range(n):
    rows.append([float(next(t)) for _ in range(width)] if k=='Vertices' else [int(next(t))-1 for _ in range(width)]);next(t)
   data[k]=rows
  elif k=='End':break
 return data
records=[]
for n in [36,44]:
 mesh=read(p/f'input/published{n}-symmetric.mesh');points=mesh['Vertices'];boundary=points[:18];accepted=0
 for swap in [0,1]:
  for sx in [1,-1]:
   for sy in [1,-1]:
    if (n==44 and swap) or (n==36 and sx!=sy):continue
    transform=lambda q:[sx*q[swap],sy*q[1-swap],q[2]]
    perm=[boundary.index(transform(q)) for q in boundary]+list(range(18,len(points)))
    mapped=[None]*len(points)
    for v,q in enumerate(points):mapped[perm[v]]=transform(q)
    cells=[[perm[v] for v in h] for h in mesh['Hexahedra']]
    if (-1 if swap else 1)*sx*sy<0:
     for h in cells:h[1],h[3]=h[3],h[1];h[5],h[7]=h[7],h[5]
    stem=out/f'H{n}-s{swap}-x{sx}-y{sy}';lines=['MeshVersionFormatted 2','Dimension 3','Vertices',str(len(mapped))]
    lines+=[' '.join(map(repr,q))+' 0' for q in mapped];lines+=['Hexahedra',str(len(cells))];lines+=[' '.join(str(v+1) for v in h)+' 0' for h in cells]
    lines+=['Quadrilaterals',str(len(mesh['Quadrilaterals']))];lines+=[' '.join(str(v+1) for v in q)+' 0' for q in mesh['Quadrilaterals']];lines+=['End'];Path(str(stem)+'.mesh').write_text('\n'.join(lines)+'\n')
    cmd=[str(p/'build/mirror_search'),'--input',str(p/'input/pyramid.mesh'),'--mirrors','1','--planes','axial' if n==44 else 'diagonal','--cap',str(n),'--interior',str(33 if n==36 else 43),'--threads','1','--leader-bound','0','--partial-cover','0','--replay',str(stem)+'.mesh','--output',str(stem)+'.jsonl']
    with open(str(stem)+'.log','w') as log:r=subprocess.run(cmd,stderr=log)
    log=Path(str(stem)+'.log').read_text();record=dict(cells=n,swap=swap,sx=sx,sy=sy,exit=r.returncode,log=log);records.append(record)
    if r.returncode==0:
     accepted+=1
     cmd[cmd.index('--leader-bound')+1]='1';cmd[cmd.index('--partial-cover')+1]='6'
     with open(str(stem)+'.all.log','w') as full:r2=subprocess.run(cmd,stderr=full)
     assert r2.returncode==0,Path(str(stem)+'.all.log').read_text()
    elif r.returncode!=3:raise RuntimeError(log)
    else:
     progress=[x for x in log.splitlines() if x.startswith('PROGRESS')][-1];d=dict(x.split('=') for x in progress.split()[1:]);assert int(d['canonical_rejections'])>0,log
 print(n,'accepted orientations',accepted,flush=True);assert accepted
(p/'experiments/one-orbit-replays.json').write_text(json.dumps(records,indent=2))
