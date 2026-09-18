"""Periodic gallery: quotient classes, complete exports, and rigid-cap DOFs."""
import hashlib
import json
from pathlib import Path
import sys
import unittest
import zipfile

import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'tools'))
from periodic_transition import load, validate
from periodic_topology_classes import quotient_graph,isomorphism
from hexopt_periodic import affine_problem,boundary_sets
from mesh_tools import read_mesh


class PeriodicGalleryTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.catalog=json.loads((ROOT/'docs/assets/periodic/current-catalog.json').read_text())
        cls.meshes={t['id']:load(ROOT/'docs'/t['connectivity']) for t in cls.catalog['templates']}

    def test_periodic_exports_and_certificates(self):
        with zipfile.ZipFile(ROOT/'docs/assets/periodic/current-meshes.zip') as archive:
            self.assertIsNone(archive.testzip())
            for t in self.catalog['templates']:
                report=json.loads((ROOT/'docs'/t['certificate']).read_text())
                self.assertTrue(report['accepted']);self.assertTrue(report['exact_positive'])
                self.assertEqual(report['betti_GF2'],[1,2,1,0])
                self.assertEqual(report['cap_faces'],[4,8])
                self.assertEqual(report['periodic_json_sha256'],hashlib.sha256((ROOT/'docs'/t['connectivity']).read_bytes()).hexdigest())
                for key in ['mesh','repeated_mesh','connectivity','certificate','vtu','periodic_vtu']:
                    self.assertEqual(archive.read(t[key]),(ROOT/'docs'/t[key]).read_bytes())
                for key,repeats in [('mesh',1),('repeated_mesh',4)]:
                    points,cells,caps=read_mesh(ROOT/'docs'/t[key])
                    self.assertEqual(len(cells),repeats*t['hexes'])
                    self.assertEqual(len(caps),repeats*12)
                    for face in points[caps]:
                        edges=np.roll(face,-1,axis=0)-face
                        bottom=abs(face[:,2].mean()-points[:,2].min())<1e-8
                        expected=0.5 if bottom else np.sqrt(2)/4
                        np.testing.assert_allclose(np.linalg.norm(edges,axis=1),expected,atol=1e-8)
                        np.testing.assert_allclose((edges*np.roll(edges,1,axis=0)).sum(axis=1),0,atol=1e-8)
                    self.assertTrue(np.all((abs(points[caps,2]-points[:,2].min()).max(axis=1)<1e-8)|(abs(points[caps,2]-points[:,2].max()).max(axis=1)<1e-8)))

    def test_every_group_has_verified_flag_isomorphisms(self):
        evidence=json.loads((ROOT/'docs/assets/periodic/topology-classes.json').read_text())
        graphs={name:quotient_graph(*mesh[:3]) for name,mesh in self.meshes.items()}
        for name,witness in evidence['witnesses'].items():
            ak,aa=graphs[name];bk,ba=graphs[witness['reference']];mapping=witness['node_mapping']
            self.assertEqual(sorted(mapping),list(range(len(ak))))
            self.assertEqual(len(bk),len(ak))
            for v,w in enumerate(mapping):
                self.assertEqual(ak[v],bk[w])
                self.assertEqual({mapping[x]:label for x,label in aa[v].items()},ba[w])
        # There must be one viewer per actual class, not per coordinate fit.
        self.assertEqual(len(self.catalog['templates']),len(evidence['groups']))
        self.assertEqual(len({t['topology_class'] for t in self.catalog['templates']}),len(evidence['groups']))
        reps=self.catalog['templates']
        for i,a in enumerate(reps):
            for b in reps[:i]:
                if a['hexes']==b['hexes']:
                    self.assertIsNone(isomorphism(graphs[a['id']],graphs[b['id']]))

    def test_vertex_labels_and_deck_gauge_do_not_change_class(self):
        mesh=self.meshes['R38-improved'];p,h,o,period=mesh
        rng=np.random.default_rng(123);perm=rng.permutation(len(p));shift=np.c_[rng.integers(-2,3,(len(p),2)),np.zeros(len(p),int)]
        pp=np.empty_like(p);pp[perm]=p+shift*period
        hh=perm[h];oo=o-shift[h]
        self.assertIsNotNone(isomorphism(quotient_graph(p,h,o),quotient_graph(pp,hh,oo)))

    def test_current_meshes_pass_independent_exact_validation(self):
        for name, mesh in self.meshes.items():
            report = validate(*mesh)
            self.assertTrue(report['accepted'], name)
            self.assertTrue(report['exact_positive'], name)
            self.assertEqual(report['betti_GF2'], [1, 2, 1, 0])

    def test_affine_space_preserves_caps_and_periodicity(self):
        for t in self.catalog['templates']:
            p,h,o,period=self.meshes[t['id']];points,cells,copies,columns=affine_problem(p,h,o,period)
            change=np.zeros_like(points)
            rng=np.random.default_rng(42)
            for _,lo,hi,terms in columns:
                delta=rng.uniform(max(lo,-.03),min(hi,.03))
                for v,k,a in terms:change[v,k]+=a*delta
            bottom,top=boundary_sets(p,h,o)
            quotient=np.array([change[copies[v][0]] for v in range(len(p))])
            self.assertTrue(np.all(quotient[bottom]==0))
            np.testing.assert_allclose(quotient[top],np.tile(quotient[top[0]],(len(top),1)))
            self.assertEqual(quotient[top[0],2],0)
            self.assertTrue(any(abs(quotient[top[0],:2])>1e-6))
            for v,ids in copies.items():np.testing.assert_allclose(change[ids],np.tile(quotient[v],(len(ids),1)))
            self.assertEqual(cells.shape,(len(h),8))


if __name__=='__main__':unittest.main()
