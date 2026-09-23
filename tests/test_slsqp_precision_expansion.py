"""Fixed real source/method/denominator checks without molecular evaluations."""
from pathlib import Path
import sys
import unittest
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from affordable_common import read_json,verify
import slsqp_precision_expansion as expansion
import slsqp_precision as pilot
import adaptive_completion

class ExpansionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.path=ROOT/'workspaces/slsqp_precision_expansion_20260923/proposals_v1/manifest.json';cls.m=read_json(cls.path)

    def test_actual_all_source_preflight(self):
        v=expansion.validate(self.path)
        self.assertEqual((v['total_sources'],v['reused_pools'],v['new_searches'],v['GFN2_maximum']),(34,4,60,240))

    def test_canonical_only_calibration_membership(self):
        rows=self.m['all_sources'];self.assertEqual(len(rows),34)
        self.assertEqual(sum(r['role']=='calibration' for r in rows),25)
        self.assertEqual({r['case_id'] for r in rows if r['reused_precision_pilot']},set(pilot.CASES))
        self.assertEqual(len({t['task_id'] for t in self.m['tasks']}),60)

    def test_actual_old_iteration_failures_are_not_scored(self):
        for cid in expansion.FAILED:
            row=next(r for r in self.m['all_sources'] if r['case_id']==cid)
            self.assertIsNone(row['old_pool'])
            la=read_json(verify(row['old_receipts']['La']));self.assertEqual(la['status'],'unavailable')
            self.assertFalse(la['optimizer']['success']);self.assertEqual(la['optimizer']['iterations'],200)
            self.assertEqual(read_json(verify(row['old_receipts']['Ca']))['status'],'proposal_available')

    def test_only_numerical_tolerance_changes(self):
        old=dict(adaptive_completion.SETTINGS);engine,_=expansion.engine()
        self.assertEqual({k for k,v in old.items() if engine.SETTINGS[k]!=v},{'optimizer_ftol'})
        self.assertEqual(engine.SETTINGS['optimizer_ftol'],1e-8);self.assertEqual(adaptive_completion.SETTINGS,old)
        self.assertEqual(self.m['new_origin_calls'],0);self.assertEqual(self.m['new_DFT_calls'],0)

    def test_four_real_precision_pools_reused(self):
        reuse=read_json(verify(self.m['reuse']));self.assertEqual(reuse['available'],4)
        manifest=read_json(verify(reuse['manifest']));source=read_json(verify(manifest['source_manifest']))
        self.assertEqual(source['settings'],self.m['settings']);self.assertEqual(manifest['model'],self.m['model'])
        self.assertTrue(all(c['pool']['status']=='available' for c in reuse['cases']))

    def test_archived_qualified_endpoint_comparison(self):
        from slsqp_precision_expansion_report import endpoint
        original=read_json(ROOT/'workspaces/slsqp_precision_20260923/COMPARISON_v1.json')
        for row in original['rows']:
            for z,expected in row['endpoints'].items():
                r=endpoint(expected['old_receipt'],expected['new_receipt'])
                self.assertEqual(r['native_proposal_delta_kcal_mol'],expected['native_proposal_delta_kcal_mol'])
                self.assertEqual(r['maximum_coordinate_delta_A'],expected['maximum_coordinate_delta_A'])
                self.assertTrue(r['native_energy_pass'])

    def test_actual_failed_receipt_not_promoted_by_reporting(self):
        from slsqp_precision_expansion_report import endpoint
        for cid in expansion.FAILED:
            row=next(r for r in self.m['all_sources'] if r['case_id']==cid);pin=row['old_receipts']['La']
            r=endpoint(pin,pin)
            self.assertFalse(r['successful_old_comparison']);self.assertIsNone(r['native_energy_pass'])
            self.assertIsNone(r['native_proposal_delta_kcal_mol'])

    def test_completed34_and_canonical_reference(self):
        run=ROOT/'workspaces/slsqp_precision_expansion_20260923';r=read_json(run/'COMPARISON_v1.json');ref=read_json(run/'REFERENCE_v1.json')
        self.assertEqual((r['denominator'],r['counts']['all34']['available']),(34,34))
        self.assertEqual(r['counts']['all34']['new_old_reference']['operational'],{'correct':34})
        self.assertEqual(r['counts']['all34']['new_own_reference']['operational'],{'correct':34})
        for v in ref['variants'].values():
            self.assertEqual(len(v['rows']),25);self.assertTrue(all(q['role']=='calibration' for q in v['rows']))
            self.assertAlmostEqual(v['gap_model_kcal_mol'],v['bands']['La_min']-v['bands']['Ca_max'],places=8)
        self.assertFalse(ref['crystals_or_noncanonical_used_for_fit'])

    def test_failed_equivalence_gates_remain_visible(self):
        r=read_json(ROOT/'workspaces/slsqp_precision_expansion_20260923/COMPARISON_v1.json')
        self.assertEqual((r['numerical_passes'],r['numerical_comparison_denominator']),(27,32))
        failures={x['case_id'] for x in r['rows'] if x['numerical_pass'] is False}
        self.assertEqual(failures,{'c5axv8-pqq-la_model','i0jwn7-pqq-la_model','a8r3s4-pqq-la_model','p16027-pqq-la_model','q60ar6-pqq-la_model'})
        for row in r['rows']:
            if row['case_id'] in expansion.FAILED:
                self.assertIsNone(row['numerical_pass']);self.assertIsNone(row['old_scores']['operational'])
                self.assertEqual(row['status'],'available');self.assertIsNotNone(row['endpoints']['La']['failed_old_last_iterate_diagnostic'])

if __name__=='__main__':unittest.main()
