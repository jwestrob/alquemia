"""Shared-feature boundaries on actual source states and native references."""
from pathlib import Path
import copy
import sys
import tempfile
import unittest
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from affordable_common import InvalidArtifact,read_json,write_new
from mace_omol_neutral_run import validate,collect
from mace_omol_neutral import INDEX,RAW_SHA

class SharedNeutral(unittest.TestCase):
    def setUp(self):
        self.work=ROOT/'workspaces/mace_omol_20260917'
        self.manifest=self.work/'shared_neutral_core_v2/manifest.json'
        if not self.manifest.exists():self.skipTest('requires exact real core manifest and native references')

    def test_actual_charges_stay_distinct_and_native_references_are_neutral(self):
        self.assertEqual(validate(self.manifest)['actual_native_references'],4)
        m=read_json(self.manifest);self.assertEqual(len(m['tasks']),8)
        self.assertEqual({(t['metal'],t['charge']) for t in m['tasks']},{('Ca',-1),('La',0)})
        audit=read_json(m['embedding_audit']['path']);zero=next(r for r in audit['rows'] if r['charge']==0)
        self.assertEqual(audit['specs']['total_charge']['offset'],INDEX)
        self.assertEqual(zero['raw_sha256'],RAW_SHA)
        self.assertGreater(zero['raw_norm'],0)
        for ref in m['native_references'].values():
            self.assertEqual(ref['result']['input_state_check']['charge'],0)
            self.assertIsNotNone(ref['result']['forces'])

    def test_corrupted_real_manifest_cannot_neutralize_physical_charge(self):
        m=copy.deepcopy(read_json(self.manifest));m['tasks'][0]['charge']=0
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'corrupted_real_charge.json';write_new(p,m)
            with self.assertRaisesRegex(InvalidArtifact,'neutral task differs'):validate(p)

    def test_actual_numeric_qualification_and_partition_decision_are_separate(self):
        if not (self.work/'shared_neutral_core_report_v1/result.json').exists():
            self.skipTest('requires completed actual neutral-feature inference')
        r=collect(self.manifest);self.assertEqual(r['status'],'complete')
        self.assertTrue(r['numerical_gate_pass']);self.assertEqual(len(r['checks']),16)
        self.assertEqual(r['whole_inference_eligible'],r['partition']['pass'])
        self.assertIsNone(r['relaxation_correction_kcal_mol']);self.assertIsNone(r['calibrated_class'])
        self.assertFalse(r['physical_charge_modified'])
        for row in r['rows'].values():
            feature=row['charge_feature_adapter']
            self.assertEqual(feature['conditioned_table_index'],100)
            self.assertEqual(feature['original_index_min'],row['input_state_check']['charge']+100)
            self.assertTrue(feature['all_rows_match_learned_feature'])

if __name__=='__main__':unittest.main()
