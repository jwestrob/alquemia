"""Actual canonical source preparation for deterministic three-La membership."""
from pathlib import Path
import sys
import unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from affordable_common import read_json,verify,paired
import union_triple_reference as ref

ROOT=Path(__file__).resolve().parents[1]
RUN=ROOT/'workspaces/union_triple_preparation_20260923/canonical28_v1'


class TripleReference(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if not (RUN/'REUSE.json').exists():raise unittest.SkipTest('actual canonical preparation/reuse not complete')
        cls.s=read_json(RUN/'SELECTION.json');cls.p=read_json(RUN/'PREPARATION.json');cls.r=read_json(RUN/'REUSE.json')
        cls.primary=read_json(verify(cls.s['primary_preparation']));cls.m=read_json(verify(cls.s['primary_selection']))

    def test_exact25_designated_canonical_and3_crystals(self):
        self.assertEqual((self.p['denominator'],self.p['supported']),(28,28))
        self.assertEqual(len(self.s['canonical_sources']),25)
        self.assertTrue(all(c['source']['canonical_coordinate_match'] for c in self.s['canonical_sources']))
        self.assertEqual([r['case_id'] for r in self.p['cases'][-3:]],['1H4I','4MAE','1KB0'])
        self.assertEqual(self.p['new_molecular_calls'],0);self.assertIsNone(self.p['new_reference'])

    def test_lexicographic_complete_choice_without_label_or_score(self):
        for c in self.s['canonical_sources']:
            options=sorted([t for t in self.m['triples'] if t['protein_id']==c['protein_id']],key=lambda t:tuple(sorted(t['members'])))
            chosen=next(t for t in options if not t['missing_members'])
            self.assertEqual(c['selected_triple'],chosen)
            skipped=[t for t in options[:options.index(chosen)] if t['missing_members']]
            self.assertEqual([t['triple_id'] for t in c['skipped_incomplete_before_selected']],[t['triple_id'] for t in skipped])
            self.assertTrue(all('__conditioned_La__' in x for x in chosen['members']))
            self.assertNotIn(c['case_id'],chosen['members'])

    def test_real_group_state_and_paired_coordinates(self):
        anchor={(r['selection_id'],r['case_id']):r for r in self.primary['cases']}
        selected={r['case_id']:r for r in self.s['canonical_sources']}
        for r in self.p['cases'][:25]:
            c=selected[r['case_id']];a=anchor[c['selection_id'],c['selected_triple']['members'][0]]
            self.assertEqual(r['state_key'],a['state_key']);self.assertEqual(r['protein_key'],a['protein_key'])
            ep=r['representations']['context']['endpoints']
            paired(verify(ep['La']['xyz']),verify(ep['Ca']['xyz']),ep['La']['charge'],ep['Ca']['charge'])
        self.assertEqual(sum(r['tenfold_comparison']['exact_reuse_eligible'] for r in self.p['cases']),19)

    def test_crystal_preparation_is_literal_archived_reuse(self):
        old=read_json(verify(self.s['crystals']))['cases']
        for r,c in zip(self.p['cases'][-3:],old):
            self.assertEqual(r['representations'],c['representations']);self.assertEqual(r['state_signature'],c['state_signature'])
            self.assertTrue(r['tenfold_comparison']['exact_reuse_eligible'])


if __name__=='__main__':unittest.main()
