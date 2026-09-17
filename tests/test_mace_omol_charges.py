"""Charge projection and guard tests on actual normalized quantum artifacts."""
import copy
from pathlib import Path
import sys
import unittest
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from affordable_common import InvalidArtifact,read_json,verify,xyz
from mace_omol_charges import projection,validate,quality
W=ROOT/'workspaces/mace_omol_20260917';P=W/'matched_H_prepared_v1/preparation.json';M=W/'normalized_charge_v1/manifest.json'

@unittest.skipUnless(P.exists(),'real normalized preparation missing')
class ChargeProjection(unittest.TestCase):
    def test_projection_reproduces_actual_source_geometry_and_overlap(self):
        for pin in read_json(P)['cases'].values():
            case=read_json(verify(pin));atoms=xyz(verify(case['endpoints']['La']['xyz']));p=projection(case,atoms)
            matrix=np.array(p['weights']);coords=np.array(p['coordinates_A'])
            self.assertEqual(len(p['physical_ids']),len(set(p['physical_ids'])))
            np.testing.assert_allclose(matrix.sum(0),1,rtol=0,atol=1e-12)
            np.testing.assert_allclose(coords.T@matrix,np.array([a[1:] for a in atoms]).T,atol=1e-9,rtol=0)
            self.assertTrue(any(np.count_nonzero(row)>1 for row in matrix))
            self.assertTrue(all(0<c['lambda']<1 for c in p['caps']))

    def test_corrupted_actual_cap_or_source_is_rejected(self):
        c=read_json(verify(read_json(P)['cases']['GGR_extended']));atoms=xyz(verify(c['endpoints']['La']['xyz']))
        altered=copy.deepcopy(c);cap=next(a for a in altered['mapping'] if a['kind']=='sigma_link_H');cap['length_A']=10.
        with self.assertRaisesRegex(InvalidArtifact,'cap anchors/geometry'):projection(altered,atoms)
        corrupted=list(atoms);corrupted[1]=(corrupted[1][0],corrupted[1][1]+.01,*corrupted[1][2:])
        with self.assertRaisesRegex(InvalidArtifact,'mapped source coordinates'):projection(c,corrupted)

    def test_real_manifest_and_paired_probes(self):
        if not M.exists():self.skipTest('charge preparation not complete')
        self.assertEqual(validate(M)['utility_calls'],16)

    def test_actual_charge_potential_outputs(self):
        p=W/'normalized_charge_report_v1/result.json'
        if not p.exists():self.skipTest('scientific utility outputs not yet executed/reported')
        r=read_json(p);a=read_json(verify(r['potential_arrays']))
        self.assertEqual(r['status'],'complete');self.assertEqual(len(r['rows']),8)
        for key,v in r['rows'].items():
            self.assertEqual(v['status'],'computed');self.assertTrue(v['projection_conservation_pass'])
            self.assertLessEqual(abs(v['charge_sum_e']-(-1 if key.endswith('_Ca') else 0)),5e-5)
            self.assertEqual(v['fit_quality'],quality(np.array(a[key]['fitted']),np.array(a[key]['exact'])))
        self.assertEqual(len(r['paired_Ca_minus_La']),4)
        self.assertIsNone(r['solvent_correction']);self.assertIsNone(r['affinity_score']);self.assertFalse(r['baseline_changed'])

if __name__=='__main__':unittest.main()
