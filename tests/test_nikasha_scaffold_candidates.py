"""Actual completed collective candidates and explicitly corrupted copies."""
import copy
from pathlib import Path
import sys
import tempfile
import unittest
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from affordable_common import InvalidArtifact,read_json,record,write_new,xyz,verify
from nikasha_finite_candidates import sources,load_spec,protocol_for,SCAFFOLD_PROTOCOL
from nikasha_scaffold_candidates import physical_candidate,TARGET_IDS,context_coordinates
from nikasha_parallel_compare import scaffold_targets


class ScaffoldCandidates(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.inputs=read_json(ROOT/'diagnostics/nikasha_parallel_pilots_20260922/INPUTS.json')
        cls.path=ROOT/'workspaces/collective_scaffold_20260922/searches_v3/searches/1H4I__origin/receipt.json'
        if not cls.path.exists():raise unittest.SkipTest('actual collective calculation unavailable')
        cls.receipt=read_json(cls.path)
        if not cls.receipt['candidate_admitted']:raise unittest.SkipTest('actual first collective candidate not admitted')
        _,cls.base,_,cls.tasks=sources(cls.inputs,'1H4I')
        cls.candidate={'id':'scaffold_origin','target':'origin','mechanics_receipt':record(cls.path)}
        cls.limits={'maximum_heavy_displacement_A':.8}
        cls.temp=tempfile.TemporaryDirectory(prefix='scaffold_admission_',dir=ROOT/'workspaces')
        cls.root=Path(cls.temp.name)

    @classmethod
    def tearDownClass(cls):cls.temp.cleanup()

    def corrupt_receipt(self,receipt,name):
        path=self.root/(name+'.json');write_new(path,receipt)
        return {**self.candidate,'mechanics_receipt':record(path)}

    def test_actual_completed_parent_and_offset_caps_are_admitted(self):
        rows,checks=physical_candidate(self.candidate,self.tasks,self.limits)
        self.assertTrue(checks['parent_geometry']['pass'])
        self.assertFalse(checks['forcefield_energy_added_to_score'])
        self.assertEqual(rows['Ca'],xyz(verify(self.receipt['context_coordinates']['Ca'])))
        self.assertEqual(rows['Ca'][1:],rows['La'][1:])
        self.assertEqual(protocol_for('collective_scaffold'),SCAFFOLD_PROTOCOL)

    def test_unavailable_or_uphill_actual_receipt_is_rejected(self):
        for key,value in [('candidate_admitted',False),('parent_work_kcal_mol',1.),('forcefield_energy_added_to_score',True)]:
            receipt=copy.deepcopy(self.receipt);receipt[key]=value
            with self.assertRaises(InvalidArtifact):
                physical_candidate(self.corrupt_receipt(receipt,key),self.tasks,self.limits)

    def test_fixed_donor_corruption_is_rejected(self):
        receipt=copy.deepcopy(self.receipt);case=read_json(verify(receipt['parent_case']))
        positions=read_json(verify(receipt['final_parent_positions_A']))
        donor=read_json(verify(case['donors']))['atoms'][0]['parent_index']
        positions[donor][0]+=.01
        path=self.root/'explicitly_corrupted_donor.json';write_new(path,positions)
        receipt['final_parent_positions_A']=record(path)
        with self.assertRaises(InvalidArtifact):
            physical_candidate(self.corrupt_receipt(receipt,'donor_receipt'),self.tasks,self.limits)

    def test_real_source_context_replay_preserves_original_caps(self):
        case=read_json(verify(self.receipt['parent_case']))
        original=np.array(read_json(verify(case['parent']['positions_A'])))
        for metal in ('Ca','La'):
            np.testing.assert_allclose(context_coordinates(case,original,metal),
                [a[1:] for a in xyz(verify(case['origins'][metal]['xyz']))],atol=1e-12,rtol=0)

    def test_missing_target_does_not_become_successful_partial_source(self):
        result=scaffold_targets(self.base,self.base)
        self.assertEqual({r['status'] for r in result.values()},{'unavailable'})
        spec={'branch':'collective_scaffold','inputs':record(ROOT/'diagnostics/nikasha_parallel_pilots_20260922/INPUTS.json'),
              'agreement':record(ROOT/'diagnostics/scaffold_restart_round_20260922/PLAN.md'),
              'coordinate_limits':self.limits,'maximum_candidates_per_case':3,
              'cases':[{'case_id':r['case_id'],'status':'unavailable','candidates':[]} for r in self.inputs['cases']]}
        spec['cases'][0].update(status='prepared',candidates=[self.candidate])
        path=self.root/'explicit_missing_two_targets.json';write_new(path,spec)
        with self.assertRaises(InvalidArtifact):load_spec(path)

    def test_matched_target_identity_algebra_uses_real_archived_energies(self):
        # Parser/algebra identity only: reuse real old cells, not simulated FF outputs.
        case=copy.deepcopy(self.base)
        for target,new in TARGET_IDS.items():
            old={'origin':'origin','Ca_adaptive':'adaptive_Ca','La_adaptive':'adaptive_La'}[target]
            case['aliases'][new]={**case['aliases'][old], 'physical_checks':{'test_scope':'archived_energy_identity'}}
        result=scaffold_targets(case,self.base)
        for row in result.values():
            self.assertEqual(row['status'],'available')
            self.assertEqual(row['delta_R'],0.)
            self.assertEqual(row['delta_R_minus_origin_response'],0.)
            for metal in ('Ca','La'):self.assertEqual(row['endpoint_work'][metal]['composite_kcal_mol'],0.)


if __name__=='__main__':unittest.main()
