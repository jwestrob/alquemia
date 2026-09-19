import sys
from pathlib import Path
import unittest
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import second_shell_hydrogen as h
from affordable_common import read_json,xyz,verify,InvalidArtifact

class TestPhysicalH(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.root=Path(__file__).resolve().parents[1]
        cls.manifest=cls.root/'workspaces/second_shell_20260919/pqq_h_v1/manifest.json'
        cls.m=read_json(cls.manifest)

    def test_pinned_four_real_contexts(self):
        self.assertEqual(h.validate(self.manifest)['tasks'],4)
        for t in self.m['tasks']:
            prep=read_json(verify(t['source_preparation']))
            caps={a['qm_index'] for a in prep['mapping']['source_to_qm'] if a['kind']=='sigma_link_H'}
            rows=xyz(verify(t['xyz']))
            self.assertFalse(caps&set(t['mobile']))
            self.assertTrue(all(rows[i][0]=='H' for i in t['mobile']))

    def test_fixed_geometry_and_corrupted_heavy_copy(self):
        t=self.m['tasks'][0];rows=xyz(verify(t['xyz']));coords=np.array([a[1:] for a in rows])
        self.assertEqual(float(h.chemical_checks(rows,coords,t['mobile'],t['hydrogen_parents']).max()),0.)
        corrupted=coords.copy();corrupted[0,0]+=.001
        with self.assertRaises(InvalidArtifact):h.chemical_checks(rows,corrupted,t['mobile'],t['hydrogen_parents'])

    def test_declared_domain_on_real_H(self):
        t=self.m['tasks'][0];rows=xyz(verify(t['xyz']));coords=np.array([a[1:] for a in rows])
        i=t['mobile'][0];parent=t['hydrogen_parents'][0];axis=coords[i]-coords[parent];axis/=np.linalg.norm(axis)
        perturbed=coords.copy();perturbed[i]+=.01*axis
        self.assertAlmostEqual(float(h.chemical_checks(rows,perturbed,t['mobile'],t['hydrogen_parents']).max()),.01)
        perturbed[i]=coords[i]+.36*axis
        with self.assertRaises(InvalidArtifact):h.chemical_checks(rows,perturbed,t['mobile'],t['hydrogen_parents'])

    def test_actual_prepared_core_and_DFT_algebra(self):
        path=self.root/'workspaces/second_shell_20260919/pqq_h_dft_v1/collection_1202101.json'
        if not path.exists():self.skipTest('actual DFT outputs unavailable; no substitute')
        result=read_json(path)
        for row in result['rows']:
            task=next(t for t in self.m['tasks'] if t['case']==row['case'] and t['metal']==row['metal'])
            initial=xyz(verify(task['core']['xyz']));final=xyz(verify(row['proposal']['core_xyz']))
            mobile=set(row['proposal']['transferred_core_H_indices'])
            self.assertTrue(all(initial[i][0]=='H' for i in mobile))
            self.assertTrue(all(initial[i]==final[i] for i in range(len(initial)) if i not in mobile))
        before={r['case']:r['before_R_hartree'] for r in result['cases']}
        after={r['case']:r['after_R_hartree'] for r in result['cases']}
        self.assertAlmostEqual((after['4MAE']-after['1H4I'])*h.HA_TO_KCAL,result['PQQ_gap_after_kcal_mol'])
        self.assertAlmostEqual(((after['4MAE']-before['4MAE'])-(after['1H4I']-before['1H4I']))*h.HA_TO_KCAL,
                               result['PQQ_gap_after_kcal_mol']-result['PQQ_gap_before_kcal_mol'])

if __name__=='__main__':unittest.main()
