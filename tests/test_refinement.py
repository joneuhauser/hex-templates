"""Refinement boundary enumeration, independent discovery and certified exports."""
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
import zipfile

import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/'tools'),str(ROOT/'solver')]
from refinement import INPUT, boundary, enumerate_sides, patch_geometry, validate, representatives
from mesh_tools import read_mesh, write_mesh
from compare_meshes import isomorphism


class RefinementTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.cases=json.loads((INPUT/'cases.json').read_text())
        cls.catalog=json.loads((ROOT/'docs/assets/refinement/catalog.json').read_text())

    def test_side_enumeration_and_fixed_boundaries(self):
        with tempfile.TemporaryDirectory() as tmp:
            output=Path(tmp)/'sides';enumerate_sides(output, grids=((1,3),(2,4)))
            self.assertEqual(json.loads((output/'cases.json').read_text()),[c for c in self.cases if c['top_grid']<5])
            large=Path(tmp)/'large';enumerate_sides(large, grids=((3,5),), max_quads=8)
            self.assertEqual(json.loads((large/'cases.json').read_text()),[c for c in self.cases if c['top_grid']==5 and c['side_quads']<=8])
            self.assertEqual(len(self.cases),28)
            self.assertEqual([sum(c['top_grid']==n for c in self.cases) for n in (3,4,5)],[6,10,12])
            for c in self.cases:
                generated=(large if c['top_grid']==5 else output)/f'{c["id"]}.mesh'
                if generated.exists():
                    self.assertEqual(generated.read_bytes(),(INPUT/f'{c["id"]}.mesh').read_bytes())
                p,score=patch_geometry(c['side'],c['bottom_grid'],c['top_grid'],c['vertical_interior_vertices'])
                bp,q=boundary(c['bottom_grid'],c['top_grid'],c['side'],p)
                self.assertEqual(len(q),c['top_grid']**2+c['bottom_grid']**2+4*c['side_quads'])
                self.assertAlmostEqual(float(score),c['side_quality'])
                self.assertEqual(len(bp),len(q)+2)

    def test_odd_boundary_parity_and_bad_options(self):
        binary=ROOT/'solver/build/quad_disk_search'
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'quads.jsonl'
            command=[str(binary),'4','0',str(path),'--bottom-segments','1','--top-segments','3']
            result=subprocess.run(command,capture_output=True,text=True,check=True)
            rows=[json.loads(s) for s in path.read_text().splitlines()]
            self.assertEqual(len(rows),1)
            self.assertEqual(len(rows[0]['quads']),4)
            self.assertTrue(all(v!=w for v,w in enumerate(rows[0]['mirror'])))
            for options in [['--top-segments','2'],['--bottom-segments','0'],['--top-segments']]:
                result=subprocess.run(command+options,capture_output=True,text=True)
                self.assertEqual(result.returncode,2)

    def test_unseeded_smallest_connectivities(self):
        with tempfile.TemporaryDirectory() as tmp:
            for name in ['F9-1-H13','F16-4-H28','F25-9-H33']:
                t=next(t for t in self.catalog['templates'] if t['id']==name)
                bp,_,_=read_mesh(INPUT/f'{t["boundary_case"]}.mesh')
                p,h,_=read_mesh(ROOT/'docs'/t['mesh'])
                path=Path(tmp)/f'{name}.jsonl'
                result=subprocess.run([str(ROOT/'solver/build/mirror_search'), '--input',str(INPUT/f'{t["boundary_case"]}.mesh'),
                    '--cap',str(t['hexes']),'--mirrors','2','--planes','axial','--threads','2','--seconds','20',
                    '--boundary-canonical','0','--output',str(path)],capture_output=True,text=True,check=True,timeout=30)
                self.assertIn('SEARCH_FINISHED',result.stderr)
                candidates=[json.loads(s) for s in path.read_text().splitlines()]
                unique=representatives(candidates,bp)
                self.assertEqual(len(unique),1)
                target=dict(vertices=len(p),hexes=len(h),cells=h.tolist())
                self.assertTrue(any(isomorphism(c,target,len(bp)) is not None for c in candidates))

    def test_large_boundary_vertex_allowance(self):
        t=next(t for t in self.catalog['templates'] if t['id']=='F25-9-H33')
        with tempfile.TemporaryDirectory() as tmp:
            result=subprocess.run([str(ROOT/'solver/build/mirror_search'),
                '--input',str(INPUT/f'{t["boundary_case"]}.mesh'),'--cap','46',
                '--threads','1','--mirrors','2','--planes','axial',
                '--replay',str(ROOT/'docs'/t['mesh']),'--output',str(Path(tmp)/'replay.jsonl')],
                capture_output=True,text=True,check=True,timeout=30)
            self.assertIn('REPLAY_OK cells=33',result.stderr)

    def test_exact_geometry_and_archive(self):
        with zipfile.ZipFile(ROOT/'docs/assets/refinement/templates.zip') as archive:
            self.assertIsNone(archive.testzip())
            for t in self.catalog['templates']:
                case=next(c for c in self.cases if c['id']==t['boundary_case'])
                report=validate(ROOT/'docs'/t['mesh'],case)
                self.assertTrue(report['accepted'],t['id'])
                self.assertEqual(report['cap_quads'],[t['coarse_quads'],t['fine_quads']])
                digest=hashlib.sha256((ROOT/'docs'/t['mesh']).read_bytes()).hexdigest()
                self.assertEqual(digest,t['sha256'])
                self.assertEqual(json.loads((ROOT/'docs'/t['certificate']).read_text())['mesh_sha256'],digest)
                for key in ['mesh','certificate','connectivity','vtu']:
                    self.assertEqual(archive.read(str(Path(t[key]).relative_to('assets/refinement'))),(ROOT/'docs'/t[key]).read_bytes())
            for c in self.cases:
                for path in [f'sides/{c["id"]}.json',f'sides/{c["id"]}.svg',f'boundaries/{c["id"]}.mesh']:
                    self.assertEqual(archive.read(path),(ROOT/'docs/assets/refinement'/path).read_bytes())

    def test_inversion_and_boundary_motion_rejected(self):
        t=self.catalog['templates'][0]
        c=next(c for c in self.cases if c['id']==t['boundary_case'])
        p,h,q=read_mesh(ROOT/'docs'/t['mesh'])
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'bad.mesh'
            moved=p.copy();moved[0,0]+=.001
            write_mesh(path,moved,h,q)
            self.assertFalse(validate(path,c)['accepted'])
            inverted=h.copy();inverted[:,[1,3]]=inverted[:,[3,1]];inverted[:,[5,7]]=inverted[:,[7,5]]
            write_mesh(path,p,inverted,q)
            self.assertFalse(validate(path,c)['accepted'])


if __name__=='__main__':unittest.main()
