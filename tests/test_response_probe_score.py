"""Prepared-core scoring checks using completed real continuation artifacts."""
import copy
import json
from pathlib import Path
import sys
import tempfile
import unittest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
import response_probe_score as scorer


class PreparedResponseScoringTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.base=ROOT/'workspaces/response_probe_20260919/continuation_v2'
        if not (cls.base/'GGR_2FW0_score_v1.json').exists():
            raise unittest.SkipTest('actual completed response continuation unavailable')
        cls.prep=ROOT/'workspaces/ggr_mechanism_20260915/stage_b_prepared_v1/2fw0/alpha_caps/preparation_manifest.json'
        cls.collection=cls.base/'mace/collection_job_1202089.json'
        cls.evaluation=ROOT/'workspaces/response_probe_20260919/evaluation_v2.json'
        cls.release=ROOT/'diagnostics/pqq_pmdh_fixed_core_calibration_20260914/result.json'

    def test_real_prepared_pair_and_completed_receipts(self):
        p,states,*_=scorer.prepared_states(self.prep)
        pair=scorer.dft_pair(self.base/'dft/manifest.json','GGR_2FW0',states,p['protocol_id'])
        saved=scorer.read_json(self.base/'GGR_2FW0_score_v1.json')
        self.assertEqual(saved['DFT_endpoints'],pair)
        self.assertEqual(saved['response_heads']['status'],'absolute_cross_target_classification_unavailable')
        self.assertIsNone(saved['response_heads']['class'])
        self.assertIsNone(saved['baseline']['S_kcal_mol'])

    def test_corrupt_real_MACE_artifact_cannot_fall_back_to_baseline_success(self):
        data=scorer.read_json(self.collection)
        data['rows']['GGR_2FW0_Ca']['forces']['sha256']='corrupted_real_fixture'
        with tempfile.TemporaryDirectory() as tmp:
            collection=Path(tmp)/'collection.json';collection.write_text(json.dumps(data))
            r=scorer.score(self.prep,collection,self.base/'dft/manifest.json','GGR_2FW0',
                self.evaluation,self.release,Path(tmp)/'result.json')
        self.assertEqual(r['status'],'unavailable')
        self.assertIsNotNone(r['baseline']['R_kcal_mol'])
        self.assertIsNone(r['response_heads'])

    def test_actual_failed_ORCA_launch_remains_unavailable(self):
        prep=ROOT/'workspaces/benchmark_augmentation_20260918/prepared_8gy2_v1/01_8GY2/pmdh_fc_holdout_8gy2_carve_manifest.json'
        with tempfile.TemporaryDirectory() as tmp:
            r=scorer.score(prep,self.collection,self.base/'dft/manifest.json','8GY2',
                self.evaluation,self.release,Path(tmp)/'failed.json')
        self.assertEqual(r['status'],'unavailable')
        self.assertIsNone(r['baseline'])
        self.assertIsNone(r['response_heads'])

    def test_runtime_input_final_newline_on_actual_repaired_8GY2(self):
        for metal in ('Ca','La'):
            path=self.base/f'dft_8gy2_retry_v1/8GY2_{metal}/endpoint.inp'
            self.assertTrue(path.read_bytes().endswith(b'\n'))
            original=self.base/f'dft/8GY2_{metal}/endpoint.inp'
            self.assertEqual(path.read_bytes(),original.read_bytes()+b'\n')


if __name__=='__main__':unittest.main()
