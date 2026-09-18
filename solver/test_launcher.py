from pathlib import Path
import os,subprocess,json
p=Path(__file__).resolve().parent
base=dict(os.environ,CAP='12',THREADS='4',TIME_LIMIT='5',MIRRORS='1')
for key in ['PARTIAL_COVER','PARTIAL_WORK','FACE_ORDER','ROOT_FACE','COVER_WORK','PLANES','INTERIOR_CAP']:base.pop(key,None)
records=[]
for key,value in [('PARTIAL_COVER','4'),('PARTIAL_COVER','8'),('PARTIAL_WORK','-1'),('PARTIAL_WORK','1000001'),('FACE_ORDER','6'),('ROOT_FACE','-1'),('ROOT_FACE','16')]:
 r=subprocess.run(['bash',str(p/'run.sh'),'foreground'],env=dict(base,**{key:value}),capture_output=True,text=True)
 assert r.returncode==2,(key,value,r.stdout,r.stderr);records.append(dict(option=key,value=value,rejected=True))
for mirrors in [0,1,2]:
 r=subprocess.run(['bash',str(p/'run.sh'),'foreground'],env=dict(base,MIRRORS=str(mirrors),ROOT_FACE='2' if mirrors==0 else 'auto'),capture_output=True,text=True)
 assert r.returncode==0,r.stdout+r.stderr
 line=next(x for x in r.stdout.splitlines() if x.startswith('Run directory: '));run=Path(line.split(': ',1)[1])
 assert (run/'status').read_text().strip()=='enumeration_finished'
 opts=(run/'options').read_text()
 assert f'--partial-work {0 if mirrors==1 else 1000}\n' in opts
 assert int((run/'cover_work').read_text())==(100000 if mirrors==0 else 1000)
 if mirrors==0:assert '--root-face 2\n' in opts
 if mirrors==1:assert '--face-order 5\n' in opts
 planes=(run/'planes').read_text().splitlines();assert planes==(['none'] if mirrors==0 else ['axial','diagonal'])
 for plane in planes:
  assert (run/(plane+'.exit_code')).read_text().strip()=='0'
  log=(run/(plane+'.log')).read_text();assert 'SEARCH_FINISHED candidates=0' in log
 records.append(dict(mirrors=mirrors,run=str(run.relative_to(p)),complete=True,options=opts))
(p/'experiments/validation').mkdir(parents=True,exist_ok=True)
(p/'experiments/validation/launcher.json').write_text(json.dumps(records,indent=2)+'\n')
print('LAUNCHER_OK invalid settings rejected; all three modes completed with correct frozen options')
