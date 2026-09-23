"""Actual A8 components, source identities and saved geometry, zero molecular calls."""
from pathlib import Path
import importlib.util
import json
import statistics
import sys
import tempfile
import unittest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from affordable_common import read_json,verify
A=ROOT/'workspaces/union_triple_transfer_20260923/A8_AUDIT_v2.json'
G=ROOT/'workspaces/union_triple_transfer_20260923/A8_GEOMETRY_v1.json'

class ActualA8Audit(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.a=read_json(A);cls.g=read_json(G)

    def test_all_four_triples_and_seven_pairs_no_favorable_subset(self):
        a=self.a;self.assertEqual(len(a['rows']),7);self.assertEqual(len(a['all_four_triples']),4)
        self.assertEqual({t['triple_id'] for t in a['all_four_triples']},{'triple_008','triple_009','triple_010','triple_011'})
        self.assertEqual({t['expected_class'] for t in a['all_four_triples']},{'Ca'})
        for t in a['all_four_triples']:
            for m in t['methods'].values():
                self.assertEqual(m['required_members'],3);self.assertEqual(m['missing_members'],[])
                self.assertEqual(m['R'],statistics.median(m['member_R']))
        self.assertEqual(sum(t['methods']['triple_union']['decision']=='inconclusive' for t in a['all_four_triples']),2)
        self.assertTrue(all(t['methods']['union_precision']['decision']=='Ca-supported' for t in a['all_four_triples']))

    def test_actual_fragment_charge_and_component_accounting(self):
        for r in self.a['rows']:
            self.assertEqual((r['before_La_charge'],r['after_La_charge']),(-2,-2))
            self.assertEqual(r['added_fragments'],[])
            removed={tuple(x) for x in r['removed_fragments']}
            expected={('sidechain','A',352,'','SER')}
            if r['after_atoms']==149:expected.add(('sidechain','A',282,'','TRP'))
            else:self.assertEqual(r['after_atoms'],168)
            self.assertEqual(removed,expected)
            delta=r['R_component_changes']
            self.assertAlmostEqual(delta['composite_R_model_kcal_mol'],delta['native_R_model_kcal_mol']+delta['solvation_delta_R_kcal_mol'],places=7)
            self.assertAlmostEqual(delta['solvation_delta_R_kcal_mol'],delta['GFN2_ALPB_R_kcal_mol']-delta['GFN2_vacuum_R_kcal_mol'],places=7)
            for stage in ('before','after'):
                s=r[stage];work=s['selected_work'];dw=work['Ca']['composite_kcal_mol']-work['La']['composite_kcal_mol']
                self.assertAlmostEqual(s['pool']['operational']['composite_R_model_kcal_mol']-s['origin_components']['composite_R_model_kcal_mol'],dw,places=7)
            if r['after_atoms']==168 and r['case_id'].endswith(('sample-1','sample-3')):
                self.assertEqual(r['before']['selected_modes'],r['after']['selected_modes'])
                old=r['before']['selected_work'];new=r['after']['selected_work']
                self.assertLess(new['Ca']['composite_kcal_mol']-new['La']['composite_kcal_mol'],old['Ca']['composite_kcal_mol']-old['La']['composite_kcal_mol'])
                self.assertGreater(delta['solvation_delta_R_kcal_mol'],6)
                self.assertLess(abs(delta['native_R_model_kcal_mol']),1)

    def test_actual_density_diagnostics_seeds_and_different_continuation(self):
        old={r['case_id']:r['before'] for r in self.a['rows']}
        pools=list(old.values())+[r['after'] for r in self.a['rows']]+[self.a['changed_threefold_canonical']]
        cells=[c for p in pools for c in p['scalar_cells']]
        self.assertEqual(len(cells),144)
        for c in cells:
            self.assertTrue(c['normal_termination']);self.assertTrue(c['energy_check_stop'])
            self.assertTrue(c['SCF']['energy']['within_printed_tolerance'])
            self.assertFalse(c['SCF']['max_density']['within_printed_tolerance'])
            self.assertEqual(c['SCF']['charge'],c['charge']);self.assertEqual(c['SCF']['multiplicity'],c['multiplicity'])
            self.assertEqual(c['charge'],-3 if c['metal']=='Ca' else -2)
            for pin in c['retained_restart_seeds'].values():verify(pin)
        self.assertEqual(self.a['ongoing_A8_atom_counts'],[174]);self.assertEqual(self.a['changed_canonical_atom_count'],168)
        self.assertFalse(self.a['ongoing_case_is_same_as_any_changed_pool'])
        self.assertEqual(self.a['new_molecular_calls'],0)

    def test_saved_source_geometry_replays_exactly(self):
        path=ROOT/'diagnostics/union_triple_transfer_20260923/audit_a8_geometry.py'
        spec=importlib.util.spec_from_file_location('actual_a8_geometry',path);mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
        with tempfile.TemporaryDirectory() as td:
            p=Path(td)/'replay.json';mod.audit(A,p);self.assertEqual(read_json(p),self.g)
        self.assertEqual(len(self.g['rows']),7)
        for r in self.g['rows']:
            self.assertTrue(r['source_OG_is_fixed_under_both_actual_selected_subspaces'])
            for stage in r['stages'].values():
                for p in stage:
                    self.assertEqual((p['nearest_anchor']['resnum'],p['nearest_anchor']['atom']),(300,'OD1'))
                    self.assertGreater(p['distance_A'],3.5);self.assertLess(p['distance_A'],4.3)
        self.assertEqual(self.g['new_molecular_calls'],0)

if __name__=='__main__':unittest.main()
