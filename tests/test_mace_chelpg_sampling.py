"""Fixed-orbital replay contracts pinned to real GGR wavefunctions."""
import copy
import json
from pathlib import Path
import sys
import tempfile
import unittest
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from affordable_common import InvalidArtifact,read_json,verify
from mace_chelpg_sampling import validate,native_replay,validate_identity
W=ROOT/'workspaces/mace_omol_20260917/chelpg_sampling_qualification_v1'

class Sampling(unittest.TestCase):
    def test_real_paired_state_and_byte_exact_orbital_import(self):
        r=validate(W/'manifest.json');self.assertEqual(r['tasks'],2)
        m=read_json(W/'manifest.json')
        for t in m['tasks']:
            self.assertEqual(verify(t['seed_gbw']).read_bytes(),verify(t['source_gbw']).read_bytes())
            inp=verify(t['input']).read_text()
            self.assertIn('CalcGuessEnergy NoIter CHELPG',inp);self.assertNotIn('EnGrad',inp)
            self.assertEqual(t['sampling'],{'grid_A':.3,'extent_A':2.8})

    def test_undeclared_sampling_change_rejected(self):
        m=read_json(W/'manifest.json');m['tasks'][0]['sampling']['grid_A']=.2
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'corrupted_real_sampling_manifest.json';p.write_text(json.dumps(m))
            with self.assertRaisesRegex(InvalidArtifact,'changed method or sampling'):validate(p)

    def test_actual_NoIter_receipts_and_exact_default_charge_identity(self):
        m=read_json(W/'manifest.json')
        if any(not Path(t['output_path']).exists() for t in m['tasks']):self.skipTest('native NoIter replay unavailable')
        for t in m['tasks']:
            r=native_replay(t,W/'manifest.json')
            self.assertTrue(r['preliminary_identity_pass'])
            self.assertFalse(r['original_generic_SCF_converged'])
            self.assertFalse(r['SCF_optimization_performed'])
            self.assertEqual(r['default_charge_max_error_e'],0.)
            self.assertIsNone(r['qualification_pass'])

    def test_real_density_identity_query_scope(self):
        root=W.parent/'chelpg_sampling_density_identity_v1';m=validate_identity(root/'manifest.json')
        self.assertEqual(m['requested_queries'],2);self.assertEqual(m['requested_DFT_calls'],0)
        for t in m['tasks']:
            expected=read_json(verify(t['expected']))
            self.assertEqual(int(verify(t['points']).read_text().splitlines()[0]),len(expected['potential_au']))

    def test_actual_density_identity_results(self):
        paths=list((W.parent/'chelpg_sampling_density_identity_v1').glob('collection_job_*.json'))
        if not paths:self.skipTest('new native density queries not yet executed')
        r=read_json(paths[-1]);self.assertTrue(r['complete']);self.assertTrue(r['qualification_pass'])
        for row in r['rows'].values():
            self.assertLessEqual(row['maximum_identity_error_au'],1e-8)
            verify(row['receipt']['log']);verify(row['potential'])

if __name__=='__main__':unittest.main()
