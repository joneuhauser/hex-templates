"""Optimize each periodic quotient class with HexOpt and a rigid movable top.

Uses the existing bounded affine adapter for the upstream sJGrad kernel.
Bottom coordinates and cap heights are fixed. Two shared DOFs translate the
whole top in x/y; all periodic copies of an interior vertex share 3 DOFs.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess

import numpy as np

from periodic_transition import load, save, validate, incidence

ROOT=Path(__file__).resolve().parents[1]


def boundary_sets(p,h,o):
    caps=[set(),set()]
    for entries in incidence(h,o).values():
        if len(entries)==1:
            vertices=entries[0][2][:,0]
            side=0 if np.max(abs(p[vertices,2]-p[:,2].min()))<1e-9 else 1
            caps[side].update(map(int,vertices))
    return [sorted(v) for v in caps]


def affine_problem(p,h,o,period):
    bottom,top=boundary_sets(p,h,o);fixed=set(bottom+top)
    points=[];indices={};cells=[];copies={v:[] for v in range(len(p))}
    for cell,offsets in zip(h,o):
        ids=[]
        for v,offset in zip(cell,offsets):
            k=(int(v),*map(int,offset))
            if k not in indices:
                indices[k]=len(points);copies[int(v)].append(len(points));points.append(p[v]+offset*period)
            ids.append(indices[k])
        cells.append(ids)
    columns=[]
    # Finite broad bounds keep affine inputs well-defined; enough room for
    # two full phase periods in either direction, without a shape constraint.
    for axis in (0,1):
        columns.append((0.,-2*period[axis],2*period[axis],[(i,axis,1.) for v in top for i in copies[v]]))
    height=p[:,2].max()-p[:,2].min()
    for v in range(len(p)):
        if v in fixed:continue
        for axis in range(3):
            lower,upper=(-2*period[axis],2*period[axis]) if axis<2 else (p[:,2].min()+1e-9*height-p[v,2],p[:,2].max()-1e-9*height-p[v,2])
            columns.append((0.,lower,upper,[(i,axis,1.) for i in copies[v]]))
    return np.asarray(points),np.asarray(cells),copies,columns


def write_problem(directory,p,h,o,period):
    points,cells,copies,columns=affine_problem(p,h,o,period)
    with (directory/'input.txt').open('w') as f:
        f.write(f'{len(points)} {len(cells)} 0 0\n');np.savetxt(f,points,fmt='%.17g')
        np.savetxt(f,np.repeat(np.arange(len(points))[:,None],4,axis=1),fmt='%d')
        np.savetxt(f,cells,fmt='%d')
    with (directory/'constraints.txt').open('w') as f:
        f.write(f'{len(columns)}\n')
        for value,lo,hi,terms in columns:
            f.write(f'{value:.17g} {lo:.17g} {hi:.17g} {len(terms)}')
            for v,k,a in terms:f.write(f' {v} {k} {a:.17g}')
            f.write('\n')
    return points,copies


def optimize(entry,args):
    name=entry['topology_class'];directory=args.output/name;directory.mkdir()
    source=ROOT/'docs'/entry['connectivity'];shutil.copy2(source,directory/'initial.json')
    p,h,o,period=load(source);original,copies=write_problem(directory,p,h,o,period)
    bottom,top=boundary_sets(p,h,o)
    command=[str(args.output/'hexopt'),str(directory/'input.txt'),str(directory/'optimized.xyz'),str(directory/'constraints.txt')]
    status=dict(topology_class=name,initial_geometry=entry['id'],status='running',command=command,
                optimizer='HexOpt sJGrad / bounded affine projected-gradient adapter',
                bottom_fixed=True,top_rigid_translation=True,heights_fixed=True,periodic_copies_tied=True,
                initial_quality=entry['metrics'],iterations=args.iterations,stall=args.stall)
    (directory/'status.json').write_text(json.dumps(status,indent=2)+'\n')
    with (directory/'optimizer.log').open('w') as log:
        try:
            result=subprocess.run(command,stdout=log,stderr=log,timeout=args.seconds,
                env=dict(os.environ,HEXOPT_MAX_ITERATIONS=str(args.iterations),HEXOPT_STALL_LIMIT=str(args.stall)))
            code=result.returncode
        except subprocess.TimeoutExpired:code=124
    status['exit_code']=code;status['status']='finished';status['adopted']=False
    if code==0:
        optimized=np.loadtxt(directory/'optimized.xyz');displacement=optimized-original
        new=p.copy();tie_error=0.
        for v,ids in copies.items():
            new[v]+=displacement[ids].mean(axis=0)
            tie_error=max(tie_error,float(np.max(abs(displacement[ids]-displacement[ids].mean(axis=0)))))
        delta=new[top]-p[top];shift=delta.mean(axis=0)
        fixed_error=float(np.max(abs(new[bottom]-p[bottom])))
        rigid_error=float(np.max(abs(delta-shift)))
        reference=p.copy();reference[top]+=shift
        report=validate(new,h,o,period,reference=(reference,h,o,period))
        report.update(bottom_error=fixed_error,top_rigidity_error=rigid_error,periodic_copy_error=tie_error,top_translation=shift.tolist())
        report['accepted'] &= max(fixed_error,rigid_error,tie_error,abs(float(shift[2])))<1e-10
        status['validated']=report['accepted'];status['top_translation']=shift.tolist();status['quality']=report['sampled']
        provenance=dict(source=entry['connectivity'],optimizer=status['optimizer'],bottom_fixed=True,
                        top_translation=shift.tolist(),top_rigid=True,heights_fixed=True)
        save(directory/'optimized.json',new,h,o,period,provenance)
        report['periodic_json_sha256']=hashlib.sha256((directory/'optimized.json').read_bytes()).hexdigest()
        (directory/'validation.json').write_text(json.dumps(report,indent=2)+'\n')
        status['adopted']=bool(report['accepted'] and report['sampled']['min_scaled_jacobian']>entry['metrics']['min_scaled_jacobian']+1e-9)
    (directory/'status.json').write_text(json.dumps(status,indent=2)+'\n')
    print(json.dumps({k:status[k] for k in ['topology_class','exit_code','adopted']}|{k:status[k] for k in ['top_translation','quality'] if k in status}),flush=True)
    return status


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--catalog',type=Path,default=ROOT/'docs/assets/periodic/current-catalog.json')
    parser.add_argument('--binary',type=Path,default=ROOT/'solver/build/hexopt_periodic')
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--iterations',type=int,default=300000)
    parser.add_argument('--stall',type=int,default=15000)
    parser.add_argument('--seconds',type=float,default=900)
    parser.add_argument('--jobs',type=int,default=4)
    args=parser.parse_args()
    if min(args.iterations,args.stall,args.seconds,args.jobs)<1:parser.error('Limits must be positive')
    args.output=args.output.resolve();args.output.mkdir(parents=True,exist_ok=False)
    shutil.copy2(args.binary,args.output/'hexopt')
    catalog=json.loads(args.catalog.read_text());entries=catalog['templates']
    if not all('topology_class' in t for t in entries):raise ValueError('Expected a catalog deduplicated by periodic topology')
    shutil.copy2(args.catalog,args.output/'input-catalog.json')
    metadata=dict(classes=len(entries),constraints='fixed bottom; rigid xy top translation; fixed cap heights; periodic copies tied',
        binary_sha256=hashlib.sha256((args.output/'hexopt').read_bytes()).hexdigest(),
        source_sha256={name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in ['tools/hexopt_fixed.cpp','tools/hexopt_periodic.py']},
        upstream_commit='5f46bf4f1c8dbff6cbfbec7ab2c2295ffd4d16b9')
    (args.output/'metadata.json').write_text(json.dumps(metadata,indent=2)+'\n')
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        results=list(pool.map(lambda entry:optimize(entry,args),entries))
    (args.output/'results.json').write_text(json.dumps(results,indent=2)+'\n')


if __name__=='__main__':main()
