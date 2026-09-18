"""Actual-source conductor model checks; native results are never fabricated."""
from pathlib import Path
import sys
import tempfile
import unittest
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from affordable_common import HA_TO_KCAL,InvalidArtifact,cache_key,read_json,verify,write_new
from mace_ddx_cpcm import preflight
from mace_ddx_recovery import retained_state
W=ROOT/'workspaces/mace_omol_20260917'


class DDXCPCM(unittest.TestCase):
    def test_refined_manifest_reuses_only_matching_actual_grid(self):
        m=preflight(W/'ddx_cpcm_source_v2/manifest.json')
        self.assertEqual(m['requested_new_forward_solves'],16);self.assertEqual(m['reused_state_roles'],4)
        for g in m['groups']:
            for label,pin in g['reuse'].items():
                self.assertEqual((g['variant'],g['grid']),('primary','coarse'))
                row,a=retained_state(pin)
                self.assertEqual(row['actual_parameters']['lmax'],m['grids']['coarse']['lmax'])
                self.assertEqual(row['actual_parameters']['n_lebedev'],m['grids']['coarse']['n_lebedev'])
                self.assertEqual(row['energy_prefactor'],m['energy_prefactor'])
                self.assertFalse(row['new_forward_solve_started'])

    def test_real_sources_preserved_and_no_pcm_cache_reuse(self):
        m=preflight(W/'ddx_cpcm_source_v1/manifest.json');p=read_json(verify(m['parent']))
        self.assertEqual(m['requested_new_forward_solves'],20);self.assertEqual(m['reused_state_roles'],0)
        self.assertEqual(m['energy_prefactor'],(78.3-1)/78.3)
        for g,old in zip(m['groups'],p['groups']):
            self.assertEqual(g['source'],old['source']);self.assertEqual(g['states'],old['states'])
            self.assertFalse(g['reuse']);self.assertNotEqual(g['cache_key'],old['cache_key'])

    def test_corrupted_real_energy_scale_rejected(self):
        m=read_json(W/'ddx_cpcm_source_v1/manifest.json');m['energy_prefactor']=1.
        m['cache_key']=cache_key({k:v for k,v in m.items() if k!='cache_key'})
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'corrupted_real_conductor_scale.json';write_new(p,m)
            with self.assertRaises(InvalidArtifact):preflight(p)

    def test_actual_native_energy_scaling_and_failure_visibility(self):
        paths=list((W/'ddx_cpcm_source_v1').glob('collection_job_*.json'))
        if not paths:self.skipTest('actual conductor integration has not completed')
        r=read_json(paths[-1]);self.assertEqual(r['protocol'],'fixed_source_full_protein_ddCPCM_component_v1')
        self.assertIsNone(r['new_score']);self.assertFalse(r['baseline_changed'])
        for group in r['groups'].values():
            with np.load(verify(group['arrays']),allow_pickle=False) as a:
                for label,row in group['rows'].items():
                    if row['status']!='computed':
                        self.assertIsNone(row['energy_kcal_mol']);self.assertIn('failure_reason',row);continue
                    raw=.5*np.sum(a[label+'_psi']*a[label+'_x'])
                    self.assertLess(abs(raw-row['raw_native_energy_hartree'])*HA_TO_KCAL,1e-8)
                    self.assertEqual(row['energy_hartree'],row['raw_native_energy_hartree']*r['energy_prefactor'])
                    self.assertEqual(row['energy_kcal_mol'],row['energy_hartree']*HA_TO_KCAL)


if __name__=='__main__':unittest.main()
