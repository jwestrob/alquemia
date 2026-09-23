"""Real ten-fold/crystal artifacts and explicitly corrupted copies only."""
from pathlib import Path
import copy
import sys
import tempfile
import unittest
import gemmi
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from affordable_common import InvalidArtifact,read_json,record,verify,write_new,xyz
import pqq_union_candidate as adapter
import consistent_context as context

ROOT=Path(__file__).resolve().parents[1]
RUN=ROOT/'workspaces/pqq_union_candidate_20260923'


class UnionCandidate(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if not (RUN/'two_groups_v3/RESULT.json').exists():
            raise unittest.SkipTest('actual declared-group replay artifact unavailable')
        cls.req=read_json(RUN/'TWO_GROUPS_v1.json')
        cls.prep=read_json(RUN/'two_groups_v3/PREPARATION.json')
        cls.result=read_json(RUN/'two_groups_v3/RESULT.json')

    def test_declared_real_sources_and_no_three_fold_substitution(self):
        self.assertEqual(len(adapter.check_request(self.req)),2)
        self.assertEqual(len(self.req['cases']),20)
        bad=copy.deepcopy(self.req);bad['groups'][0]['members']=bad['groups'][0]['members'][:3]
        with self.assertRaisesRegex(InvalidArtifact,'five distinct'):adapter.check_request(bad)
        bad=copy.deepcopy(self.req);bad['profile']='unqualified_cheap_profile'
        with self.assertRaisesRegex(InvalidArtifact,'unsupported profile'):adapter.check_request(bad)
        bad=copy.deepcopy(self.req);g=bad['groups'][0]
        original=next(x for x in bad['cases'] if x['case_id']==g['canonical_case_id'])
        other=next(x for x in bad['cases'] if x['root_case_id']==g['protein_id'] and x['source_conditioning_metal']=='La' and x is not original)
        original['canonical_coordinate_match']=False;other['canonical_coordinate_match']=True;g['canonical_case_id']=other['case_id']
        with self.assertRaisesRegex(InvalidArtifact,'actual archived calibration source'):adapter.check_request(bad)

    def test_actual_extra_protein_chain_is_not_hidden_by_declared_identity(self):
        bad=copy.deepcopy(self.req)
        with tempfile.TemporaryDirectory() as td:
            # Deliberately corrupted real source: duplicate its real chain A as Z.
            st=gemmi.read_structure(str(verify(bad['cases'][0]['source_structure'])))
            chain=st[0]['A'].clone();chain.name='Z';st[0].add_chain(chain)
            path=Path(td)/'corrupted_added_protein_chain.cif';st.make_mmcif_document().write_file(str(path))
            bad['cases'][0]['source_structure']=record(path)
            with self.assertRaisesRegex(InvalidArtifact,'sequence, numbering, roles or assembly'):
                adapter.check_request(bad)

    def test_union_identity_and_coordinates_reconstruct_real_archive(self):
        self.assertEqual((self.prep['supported'],self.prep['denominator']),(15,20))
        arc=read_json(verify(self.req['archive']))
        archived={r['case_id']:r for pin in arc['union_preparations'] for r in read_json(verify(pin))['cases']}
        for group in self.prep['groups']:
            rows=[r for r in self.prep['cases'] if r['group_id']==group['protein_id'] and r['status']=='prepared']
            self.assertEqual(len({r['state_key'] for r in rows}),1)
            for row in rows:
                old=archived[row['case_id']]
                self.assertEqual(row['state_key'],old['state_key'])
                self.assertEqual(read_json(verify(row['union']))['fragments'],read_json(verify(old['union']))['fragments'])
                for z in ('Ca','La'):
                    self.assertTrue(context.reusable_state(row['representations']['context']['endpoints'][z],
                        old['representations']['context']['endpoints'][z]))
            self.assertNotEqual(xyz(verify(rows[0]['representations']['context']['endpoints']['La']['xyz'])),
                xyz(verify(rows[1]['representations']['context']['endpoints']['La']['xyz'])))

    def test_actual_source_scores_and_strict_aggregates_match_completed_ledger(self):
        arc=read_json(verify(self.req['archive']));ledger=read_json(verify(arc['comparison']))
        ref=read_json(verify(arc['reference']));index={r['case_id']:r for r in ledger['rows']}
        canonical={r['actual_union_case_id']:r for r in ref['rows']}
        self.assertEqual((self.result['available'],self.result['denominator']),(15,20))
        self.assertTrue(all(r['source_domain_status']=='consumed_reference' for r in self.result['rows']))
        for row in self.result['rows']:
            for variant,method in (('operational','union_adaptive'),('mathematical','union_adaptive_mathematical')):
                if row['case_id'] in index:expected=index[row['case_id']]['methods'][method]['R']
                else:expected=canonical[row['case_id']]['pool'][variant]['composite_R_model_kcal_mol']
                self.assertEqual(row['candidate'][variant]['R'],expected)
        for group in self.result['aggregates']:
            prior=next(x for x in ledger['pools']+ledger['triples'] if x['root_case_id']==group['protein_id']
                and set(x['members'])==set(group['members']))
            for new,old in (('original_static','context_composite'),('union_static','context_union'),('union_adaptive','union_adaptive')):
                self.assertEqual(group['methods'][new]['R'],prior['methods'][old]['R'])
                expected_missing={cid for cid in group['members'] if index[cid]['methods'][old]['R'] is None}
                self.assertEqual(set(group['methods'][new]['missing_members']),expected_missing)

    def test_actual_missing_five_Q9_Ca_members_remain_missing(self):
        absent=[r for r in self.result['rows'] if r['status']!='available']
        self.assertEqual(len(absent),5)
        self.assertTrue(all(r['case_id'].startswith('q9z4j7') and r['source_conditioning_metal']=='Ca' for r in absent))
        for row in absent:
            self.assertIn('exact selected PQQ required',row['reason'])
            self.assertIsNone(row['candidate']['operational']['R'])
        for group in self.result['aggregates']:
            if group['protein_id'].startswith('q9z4j7') and group['descriptor'] in ('Ca5','balanced'):
                self.assertIsNone(group['methods']['union_adaptive']['R'])
                self.assertEqual(len(group['methods']['union_adaptive']['missing_members']),5)

    def test_corrupted_real_coordinate_mapping_rejected(self):
        data=copy.deepcopy(self.prep);data['cases'][0]['representations']=copy.deepcopy(data['cases'][1]['representations'])
        plan=read_json(RUN/'two_groups_v3/plan.json')
        with tempfile.TemporaryDirectory() as td:
            pp=Path(td)/'corrupted_swapped_fold_coordinates.json';write_new(pp,data)
            plan['preparation']=record(pp);mp=Path(td)/'corrupted_plan.json';write_new(mp,plan)
            with self.assertRaisesRegex(InvalidArtifact,'state or coordinates differ'):
                adapter.replay(mp,Path(td)/'must_not_exist.json')

    def test_explicit_corrupted_unavailable_candidate_keeps_real_static_separate(self):
        data=copy.deepcopy(self.prep);data['cases'][0].update(status='unsupported',reason='explicit_corrupted_fixture_failure')
        plan=read_json(RUN/'two_groups_v3/plan.json')
        with tempfile.TemporaryDirectory() as td:
            pp=Path(td)/'corrupted_failure_status.json';write_new(pp,data)
            plan['preparation']=record(pp);mp=Path(td)/'corrupted_plan.json';write_new(mp,plan)
            out=Path(td)/'corrupted_status_replay.json';adapter.replay(mp,out);r=read_json(out)['rows'][0]
            self.assertEqual(r['original_static'],self.result['rows'][0]['original_static'])
            self.assertIsNotNone(r['original_static']['R']);self.assertIsNone(r['candidate']['operational']['R'])
            self.assertEqual(r['reason'],'explicit_corrupted_fixture_failure')

    def test_actual_existing_crystal_singletons(self):
        path=RUN/'crystals_v2/RESULT.json'
        if not path.exists():self.skipTest('actual crystal replay not yet present')
        r=read_json(path);self.assertEqual((r['denominator'],r['available']),(3,3))
        self.assertEqual({x['case_id']:x['candidate']['operational']['decision'] for x in r['rows']},
            {'1H4I':'Ca-supported','4MAE':'La-supported','1KB0':'Ca-supported'})
        self.assertEqual(r['aggregates'],[]);self.assertEqual(r['new_molecular_calls'],0)


if __name__=='__main__':unittest.main()
