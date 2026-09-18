#!/usr/bin/env python3
"""Check the scope and elementary results of topological boundary quotients."""
from pathlib import Path
import subprocess
import sys
import tempfile

case=Path(__file__).resolve().parent
binary=Path(sys.argv[1]).resolve() if len(sys.argv)>1 else case/'build/mirror_search'
common=['--mirrors','0','--threads','1','--corner-geometry','0','--location-domains','0',
        '--geometry-domains','0','--boundary-canonical','1','--topology-symmetry','1','--progress','0']
with tempfile.TemporaryDirectory(prefix='hex-topology-symmetry-') as work:
    output=Path(work)/'candidates.jsonl'
    for fixture,cap,group,candidates in [('cube',1,48,1),('pyramid',10,16,0)]:
        cmd=[str(binary),'--input',str(case/f'input/{fixture}.mesh'),'--cap',str(cap),*common,'--output',str(output)]
        r=subprocess.run(cmd,capture_output=True,text=True)
        assert r.returncode==0,(cmd,r.stderr)
        assert f'boundary_group={group}' in r.stderr
        assert r.stderr.rstrip().endswith(f'SEARCH_FINISHED candidates={candidates}')
    for override in [['--corner-geometry','1'],['--location-domains','1'],['--geometry-domains','1'],
                     ['--boundary-canonical','0'],['--mirrors','1']]:
        cmd=[str(binary),'--input',str(case/'input/cube.mesh'),'--cap','1',*common,*override,'--output',str(output)]
        r=subprocess.run(cmd,capture_output=True,text=True)
        assert r.returncode==2,(cmd,r.stderr)
        assert 'Topological boundary automorphisms require' in r.stderr
print('TOPOLOGY_SYMMETRY_OK cube and pyramid quotients; five scope guards')
