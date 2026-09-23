"""Real fixed source groups and existing physical maps, without molecular execution."""
from pathlib import Path
import copy
import json
import sys
import unittest
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from affordable_common import read_json,verify
import motion_envelope as envelope
import second_shell_context as shell
from mace_site_kinematics import Kinematics

class EnvelopePreparation(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.p=read_json(ROOT/'workspaces/motion_envelope_20260923/run_v2/PREPARATION.json')
        cls.m=envelope.validate(verify(cls.p['selection']));cls.d=read_json(verify(cls.m['discovery']))
    def test_actual_fixed100_and_canonical_choices(self):
        old=read_json(verify(self.m['prior_selection']))
        self.assertEqual([(t['triple_id'],t['members'],t['missing_members'])for t in old['triples']],[(t['triple_id'],t['members'],t['missing_members'])for t in self.p['triples']])
        self.assertEqual(self.p['counts'],{'triple_denominator':100,'prepared_triples':94,'prior_unavailable_triples':6,'source_selection_pairs':135,'prepared_pairs':135,'pilot_denominator':34,'prepared_pilot':34})
        canonical=read_json(verify(self.m['canonical_selection']))
        self.assertEqual([(x['case_id'],x['selected_triple']['triple_id'])for x in canonical['canonical_sources']],[(x['case_id'],x['declared_triple_id'])for x in self.p['pilot']if x['role']=='canonical_calibration'])
    def test_private_discovery_only_changes_cutoff(self):
        self.assertEqual(shell.POLICY['polar_neighbor_cutoff_A'],3.5)
        self.assertEqual({k:v for k,v in envelope.POLICY.items()if k!='polar_neighbor_cutoff_A'},{k:v for k,v in shell.POLICY.items()if k!='polar_neighbor_cutoff_A'})
        self.assertEqual(envelope.POLICY['polar_neighbor_cutoff_A'],4.3)
        raw=read_json(verify(self.m['source_inputs']));cid=self.d['cases'][0]['case_id'];row=raw['sources'][cid]
        policy=copy.deepcopy(shell.POLICY);actual=envelope.discover_source(row,self.m['config'])
        self.assertEqual(shell.POLICY,policy);self.assertEqual(json.loads(json.dumps(actual)),next(r for r in self.d['cases']if r['case_id']==cid))
        self.assertTrue(all(r['anchor_identity_replayed']for r in self.d['cases']))
    def test_actual135_paired_maps_reconstruct_origins(self):
        for row in self.p['cases']:
            a=read_json(verify(row['maps']['Ca']));b=read_json(verify(row['maps']['La']));self.assertEqual(a,b)
            kin=Kinematics(a['context']);coords=kin.evaluate(np.zeros(len(kin.modes)))[1]
            np.testing.assert_allclose(coords,np.asarray(a['context']['core_positions_A']),rtol=0,atol=1e-12)
            self.assertTrue(row['mapping_eligible']);self.assertEqual(row['status'],'prepared')
    def test_missing_defining_member_blocks_pilot_without_replacement(self):
        results={(r['selection_id'],r['case_id']):copy.deepcopy(r)for r in self.p['cases']}
        target=next(x for x in self.m['pilot']if x['role']=='canonical_calibration')
        group=next(g for g in self.m['groups']if g['selection_id']==target['selection_id'])
        missing=next(c for c in group['members']if c!=group['members'][0])
        results[group['selection_id'],missing].update(status='mapping_unavailable',mapping_eligible=False,reason='explicit corruption of real member mapping for test')
        _,triples,pilot=envelope.finalize(self.m,results)
        affected=[x for x in pilot if x.get('declared_triple_id')==target['declared_triple_id']]
        self.assertTrue(affected);self.assertTrue(all(x['status']=='declared_triple_unavailable'for x in affected))
        self.assertTrue(all(x['target_preparation_status']=='prepared'for x in affected))
        self.assertEqual(len(triples),100);self.assertEqual(len(pilot),34)
    def test_unavailable_anchor_cannot_qualify_other_members(self):
        results={(r['selection_id'],r['case_id']):copy.deepcopy(r)for r in self.p['cases']}
        group=self.m['groups'][0];anchor=group['members'][0]
        for r in results.values():r['status']='prepared_awaiting_group_check'
        results[group['selection_id'],anchor].update(status='mapping_unavailable',mapping_eligible=False)
        ordered,_,_=envelope.finalize(self.m,results)
        others=[r for r in ordered if r['selection_id']==group['selection_id'] and r['case_id']!=anchor]
        self.assertTrue(others);self.assertTrue(all(r['status']=='state_anchor_unavailable'for r in others))
    def test_actual_pilot_native_only_manifest_has68_real_inputs(self):
        a=read_json(ROOT/'workspaces/motion_envelope_20260923/INVENTORY_v1.json')
        self.assertEqual(a['counts']['origin_pairs_reused'],0);self.assertEqual(a['counts']['whole_pool_reuses'],0)
        r=read_json(ROOT/'workspaces/motion_envelope_20260923/pilot34_origins_v1/READY.json');m=read_json(verify(r['MACE_manifest']))
        self.assertEqual(len(m['tasks']),68);self.assertEqual(len({t['case_id']for t in m['tasks']}),34)
        self.assertEqual(r['new_GFN2_calls'],0);self.assertEqual(r['new_optimizer_starts'],0)
        for t in m['tasks']:
            row=next(x for x in self.p['cases']if (x['selection_id'],x['case_id'])==(t['selection_id'],t['source_case_id']))
            ep=row['representations']['context']['endpoints'][t['metal']]
            self.assertEqual((t['xyz'],t['charge'],t['spin_multiplicity']),(ep['xyz'],ep['charge'],ep['multiplicity']))
            verify(t['xyz'])

if __name__=='__main__':unittest.main()
