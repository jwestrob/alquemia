"""Real pinned sources; corrupted real inputs are failure tests, not scientific data."""
import copy
from pathlib import Path
import sys
import tempfile
import unittest
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from affordable_common import InvalidArtifact, read_json, verify, write_new, xyz
from mace_explicit_field_short import validate,probes,key,collect,TOL
W=ROOT/'workspaces/mace_omol_20260917';M=W/'explicit_field_short_v1/manifest.json'

@unittest.skipUnless(M.exists(),'real preparation unavailable')
class ExplicitField(unittest.TestCase):
    def test_exact_eight_pairs_and_native_sources(self):
        self.assertEqual(validate(M)['tasks'],8)
        m=read_json(M)
        for t in m['tasks']:
            s=read_json(verify(t['state']));p,w=probes(s,xyz(verify(t['xyz'])))
            self.assertGreater(w['nearest_QM_nucleus_or_cap_A'],1e-6)
            self.assertEqual(len(p),np.count_nonzero(s['environment_charges_e']))
            self.assertEqual(t['charge'],-1 if t['metal']=='Ca' else 0)

    def test_corrupted_real_probe_overlap_rejected(self):
        t=read_json(M)['tasks'][0];s=read_json(verify(t['state']));atoms=xyz(verify(t['xyz']))
        i=next(i for i,q in enumerate(s['environment_charges_e']) if q)
        s['physical_atoms'][i]['xyz_A']=list(atoms[0][1:])
        with self.assertRaisesRegex(InvalidArtifact,'coincides'):probes(s,atoms)

    def test_scientific_cache_and_state_corruption(self):
        m=read_json(M);t=m['tasks'][0]
        for field in ('GB','charges','reuse_audit','short_reference'):
            altered=copy.deepcopy(m);altered[field]['sha256']='corrupted_real_fixture'
            self.assertNotEqual(key(t,m),key(t,altered))
        altered=copy.deepcopy(m);altered['tasks'][0]['charge']+=1
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'manifest.json';write_new(p,altered)
            with self.assertRaisesRegex(InvalidArtifact,'source changed'):validate(p)

    def test_partial_results_do_not_become_zero(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'manifest.json';write_new(p,read_json(M));r=collect(p)
            self.assertEqual(r['status'],'incomplete');self.assertIsNone(r['partition_shift_kcal_scale'])
            self.assertFalse(r['ordering_gate_pass']);self.assertFalse(r['charge_representation_gate_pass'])
            self.assertTrue(all(v['short'] is None and v['potential_receipt'] is None for v in r['rows'].values()))

    @unittest.skipUnless((W/'explicit_field_short_result_v1.json').exists(),'actual native execution unavailable')
    def test_actual_endpoint_algebra_and_unavailable_reference(self):
        r=read_json(W/'explicit_field_short_result_v1.json');self.assertEqual(r['status'],'complete')
        self.assertTrue(r['numerical_gate_pass']);self.assertEqual(len(r['coupling_quality']),5)
        for c in r['cases'].values():
            self.assertAlmostEqual(sum(c['components_R'].values()),c['hybrid_R_kcal_scale'],places=8)
            a,b=(c['endpoints'][m]['hybrid_kcal_scale'] for m in ('Ca','La'))
            self.assertLess(abs(a-b-c['hybrid_R_kcal_scale']),TOL['algebra_kcal_scale'])
        for k in ('reference','calibrated_class','combined_gradient','relaxation_correction'):self.assertIsNone(r[k])

if __name__=='__main__':unittest.main()
