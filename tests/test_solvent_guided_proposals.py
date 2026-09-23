"""Real archived pool/geometry checks; no invented scientific energies."""
import copy
from pathlib import Path
import sys
import importlib.util
import tempfile
import unittest
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
import solvent_guided_proposals as guided
from affordable_common import InvalidArtifact, read_json, verify, xyz
from mace_site_kinematics import Kinematics


class ArchivedRanking(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        p=ROOT/'workspaces/solvent_guided_20260922/archive_rankings_v2.json'
        if not p.exists():raise unittest.SkipTest('actual archived ranking result unavailable')
        cls.rank=read_json(p)
        cls.pools={}
        for pin in cls.rank['source_collections']:
            cls.pools.update({c['case_id']:c for c in read_json(verify(pin))['cases']})

    def test_all255_actual_rows_and_choices_replay(self):
        self.assertEqual(len(self.rank['rows']),255)
        for r in self.rank['rows']:
            for z in ('Ca','La'):
                self.assertEqual(guided.row_rank(self.pools[r['case_id']],z),r['metals'][z])

    def test_actual_regret_and_score_sign(self):
        for r in self.rank['rows']:
            if r['status']!='available':continue
            ca,la=(r['metals'][z] for z in ('Ca','La'))
            self.assertGreaterEqual(ca['native_selected_composite_regret_kcal_mol'],-1e-8)
            self.assertGreaterEqual(la['native_selected_composite_regret_kcal_mol'],-1e-8)
            self.assertAlmostEqual(r['delta_R_composite_vs_native_selection_kcal_mol'],
                la['selected_composite_energy_difference_kcal_mol']-ca['selected_composite_energy_difference_kcal_mol'],places=7)

    def test_actual_missing_sources_keep_null_scores(self):
        missing=[r for r in self.rank['rows'] if r['status']!='available']
        self.assertEqual(len(missing),21)
        self.assertTrue(all(r['composite_R_at_native_selected'] is None for r in missing))
        self.assertTrue(all(r['composite_R_at_composite_selected'] is None for r in missing))

    def test_labels_do_not_enter_actual_row_selection(self):
        c=copy.deepcopy(self.pools['1H4I']);expected=guided.row_rank(c,'Ca')
        c['old_result']['expected_class']='REMOVED_FOR_PARSER_TEST'
        self.assertEqual(guided.row_rank(c,'Ca'),expected)


class ActualProbes(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        p=ROOT/'workspaces/solvent_guided_20260922/probes_v1/specification.json'
        if not p.exists():raise unittest.SkipTest('actual prepared probes unavailable')
        cls.spec=read_json(p)

    def test_declared_probe_coverage_and_actual_mappings(self):
        self.assertEqual(len(self.spec['cases']),8)
        self.assertEqual(sum(len(c['probes']) for c in self.spec['cases']),32)
        for c in self.spec['cases']:
            t=c['source_tasks']['Ca'];kin=Kinematics(read_json(verify(t['mapping']))['context'])
            for q in c['candidates']:
                actual=np.asarray([a[1:] for a in xyz(verify(q['coordinate']))])
                np.testing.assert_allclose(kin.evaluate(q['full_q'])[1],actual,atol=1e-12,rtol=0)
                self.assertLessEqual(kin.displacement(q['full_q']),.8+1e-10)
                self.assertTrue(all(q['full_q'][i]==0 for i in set(range(t['mode_count']))-set(t['active_indices'])))

    def test_backoff_depends_only_on_geometry_and_preserves_failures(self):
        for c in self.spec['cases']:
            for p in c['probes']:
                if 'backoff_attempts' not in p:continue
                steps=[r['maximum_angle_step_radian'] for r in p['backoff_attempts']]
                self.assertEqual(steps,guided.PROBE_SETTINGS['geometry_only_step_sequence_radian'][:len(steps)])
                if p['status']=='available':self.assertEqual(p['step_radian'],steps[-1])
                else:self.assertIsNone(p['full_q'])

    def test_explicitly_corrupted_real_subspace_rejected(self):
        c=self.spec['cases'][0];p=c['probes'][0];t=c['source_tasks']['Ca']
        q=list(p['center']['full_q']);i=next(i for i in range(t['mode_count']) if i not in t['active_indices']);q[i]=.01
        with self.assertRaisesRegex(InvalidArtifact,'different physical subspace'):
            guided.probe_direction(q,p['native_winner']['full_q'],t,True)


class ActualFiniteReport(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        p=ROOT/'diagnostics/solvent_guided_20260922/report.py'
        loader=importlib.util.spec_from_file_location('solvent_guided_report_test',p)
        cls.module=importlib.util.module_from_spec(loader);loader.loader.exec_module(cls.module)
        cls.spec=ROOT/'workspaces/solvent_guided_20260922/probes_v1/specification.json'
        cls.pool=ROOT/'workspaces/solvent_guided_20260922/pool_v1'

    def test_actual_mace_only_phase_preserves_all_missing_solvent(self):
        p=self.pool/'after_MACE_1210100.json'
        if not p.exists():self.skipTest('actual completed MACE-only phase unavailable')
        with tempfile.TemporaryDirectory() as tmp:
            out=Path(tmp)/'report.json';self.module.report(p,self.spec,out);r=read_json(out)
            self.assertEqual(r['available'],0)
            self.assertEqual(r['denominator'],8)
            self.assertEqual(r['new_MACE_complete'],54)
            self.assertTrue(all(x['probe_R'] is None for x in r['rows']))
            self.assertTrue(all(x['prior_adaptive_R'] is not None for x in r['rows']))

    def test_actual_completed_composite_report_when_available(self):
        p=self.pool/'after_solvent_0_1210101.json'
        if not p.exists():self.skipTest('actual solvent stage has not finished')
        with tempfile.TemporaryDirectory() as tmp:
            out=Path(tmp)/'report.json';self.module.report(p,self.spec,out);r=read_json(out)
            self.assertEqual(r['denominator'],8)
            for x in r['rows']:
                self.assertIsNone(x['validated_new_protocol_decision'])
                if x['status']=='available':
                    w=x['endpoint_incremental_work_kcal_mol']
                    self.assertAlmostEqual(x['delta_R_vs_adaptive'],w['Ca']['composite_kcal_mol']-w['La']['composite_kcal_mol'],places=7)


if __name__=='__main__':unittest.main()
