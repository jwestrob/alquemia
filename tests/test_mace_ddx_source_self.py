"""Pinned real protein inputs for the resolved-boundary component diagnostic."""
from pathlib import Path
import sys
import tempfile
import unittest
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from affordable_common import InvalidArtifact,read_json,verify,write_new
from mace_ddx_source_self import validate
W=ROOT/'workspaces/mace_omol_20260917/ddx_source_self_v1'

class DDxSource(unittest.TestCase):
    def test_actual_full_cavity_and_source_only_inventory(self):
        m=validate(W/'manifest.json');self.assertEqual(len(m['groups']),8)
        self.assertEqual(sum(len(g['states']) for g in m['groups']),20)
        for g in m['groups']:
            d=read_json(verify(g['source']));n=len(d['physical_ids']);self.assertGreater(n,4000)
            mask=np.ones(n,dtype=bool);mask[d['source_indices']]=False
            for metal,e in d['endpoints'].items():
                q=np.array(e['charge_e']);self.assertTrue(np.all(q[mask]==0))
                self.assertLessEqual(abs(sum(q)-(-1 if metal=='Ca' else 0)),5e-5)
            self.assertEqual(d['radii_A'][d['physical_ids'].index('metal')],1.82485)

    def test_corrupted_real_model_rejected(self):
        m=read_json(W/'manifest.json');m['model']['solvent_epsilon']=4.
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'corrupted_real_pcm_manifest.json';write_new(p,m)
            with self.assertRaisesRegex(InvalidArtifact,'scope differs'):validate(p)

    def test_corrupted_real_source_support_rejected(self):
        m=read_json(W/'manifest.json');m['groups'][0]['states'].append('La_repeat')
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'corrupted_real_pcm_inventory.json';write_new(p,m)
            with self.assertRaisesRegex(InvalidArtifact,'scope differs'):validate(p)

    def test_direct_source_variant_preserves_all_physical_inputs(self):
        old=validate(W/'manifest.json');new=validate(W.parent/'ddx_source_self_v2/manifest.json')
        self.assertEqual(new['source_potential_policy'],'direct_Coulomb_fsum_v1')
        self.assertEqual(new['model'],old['model']);self.assertEqual(new['tolerances'],old['tolerances'])
        for a,b in zip(old['groups'],new['groups']):
            self.assertEqual(a['source']['sha256'],b['source']['sha256'])
            self.assertEqual(a['states'],b['states']);self.assertNotEqual(a['cache_key'],b['cache_key'])

    def test_actual_scientific_output_contract(self):
        paths=list(W.glob('collection_job_*.json'))
        if not paths:self.skipTest('native ddX scientific integration not yet executed')
        r=read_json(paths[-1]);self.assertEqual(r['new_DFT_calls'],0);self.assertEqual(r['new_MACE_calls'],0);self.assertIsNone(r['new_score'])
        for g in r['groups'].values():
            verify(g['arrays'])
            for row in g['rows'].values():
                if row['status']=='computed':
                    self.assertTrue(row['forward_solve_started']);self.assertTrue(np.isfinite(row['energy_hartree']))
                    self.assertEqual(row['contraction_pass'],row['contraction_error_kcal']<=1e-8)
                else:self.assertIn('failure_reason',row)

if __name__=='__main__':unittest.main()
