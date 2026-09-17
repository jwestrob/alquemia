"""Ablation inventory and failure handling against real, pinned prior calculations."""
from pathlib import Path
import copy
import sys
import tempfile
import unittest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from affordable_common import InvalidArtifact, cache_key, read_json, verify, write_new, xyz
from mace_file_checks import cached_file_checks
from mace_hybrid import check_atoms
from mace_omol import accepted_state
from mace_omol_ablation import ADAPTER,COMPONENT
from mace_omol_ablation_run import validate,collect,execution_gate

MANIFEST=ROOT/'workspaces/mace_omol_20260917/charge_ablation_development_v2/manifest.json'


@unittest.skipUnless(MANIFEST.exists(),'requires real 42-task ablation preparation')
class AblationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.m=read_json(MANIFEST)

    @cached_file_checks
    def test_all_real_coordinates_charges_and_electron_states_preserved(self):
        self.assertEqual(validate(MANIFEST)['tasks'],42)
        for t in self.m['tasks']:
            parent=read_json(verify(t['original_endpoint']['manifest']))
            old=next(r for r in parent['tasks'] if r['task_id']==t['original_endpoint']['task_id'])
            for field in ('xyz','source_xyz','charge','spin_multiplicity','state','metal_index','position','variant'):
                self.assertEqual(t[field],old[field])
            self.assertEqual(check_atoms(xyz(verify(t['xyz'])),t['charge']),t['state'])
            self.assertEqual(t['charge_feature_adapter'],ADAPTER)
            self.assertEqual(t['energy_component'],COMPONENT)
            self.assertNotEqual(t['cache_key'],old['cache_key'])

    @cached_file_checks
    def test_native_results_cannot_satisfy_modified_descriptor(self):
        source=read_json(verify(self.m['sources']['core_collection']))
        for t in self.m['tasks'][:8]:
            old=source['rows'][t['original_endpoint']['task_id']]
            self.assertFalse(accepted_state(old,t))

    @cached_file_checks
    def test_actual_masked_forward_requires_complete_adapter_receipt(self):
        t=self.m['tasks'][0]
        path=MANIFEST.parent/'execution'/t['task_id']/'attempt_001/result.json'
        if not path.exists():self.skipTest('masked scientific forward has not run')
        result=read_json(path)
        self.assertTrue(accepted_state(result,t))
        self.assertEqual(result['charge_feature_adapter']['rows'],result['input_state_check']['atoms'])
        wrong=copy.deepcopy(result);wrong['charge_feature_adapter']['rows']-=1
        self.assertFalse(accepted_state(wrong,t))
        missing=copy.deepcopy(result);missing.pop('charge_feature_adapter')
        self.assertFalse(accepted_state(missing,t))
        wrong_component=copy.deepcopy(result);wrong_component['energy_component']='MACE_OMOL_total_vacuum_energy'
        self.assertFalse(accepted_state(wrong_component,t))

    @cached_file_checks
    def test_missing_jobs_remain_unavailable_and_block_execution(self):
        # Move this real finite manifest to an empty execution directory. No output is fabricated.
        with tempfile.TemporaryDirectory() as td:
            p=Path(td)/'manifest.json';write_new(p,self.m)
            r=collect(p)
            self.assertEqual(r['status'],'incomplete')
            self.assertTrue(all(s['R_mask_model_kcal'] is None for s in r['scores'].values()))
            self.assertFalse(r['canonical_extension_permitted'])
            self.assertFalse(execution_gate(p,'core'))
            self.assertFalse(execution_gate(p,'numerical'))

    @cached_file_checks
    def test_corrupted_physical_state_rejected_even_with_recomputed_cache(self):
        m=copy.deepcopy(self.m);t=m['tasks'][0];t['charge']+=2
        payload={k:v for k,v in t.items() if k!='cache_key'}
        t['cache_key']=cache_key({'task':payload,'model':m['model'],'software':m['software'],'implementation':m['implementation']})
        with tempfile.TemporaryDirectory() as td:
            p=Path(td)/'manifest.json';write_new(p,m)
            with self.assertRaisesRegex(InvalidArtifact,'physical state'):validate(p)


if __name__=='__main__':unittest.main()
