"""Real six-source pilot inputs; no generated molecular test outputs."""
import copy
from pathlib import Path
import sys
import tempfile
import unittest
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
import standalone_response_proposals as a
MANIFEST=ROOT/'workspaces/standalone_response_20260923/run_v3/manifest.json'

class Response(unittest.TestCase):
    def test_actual_six_source_scope(self):
        r=a.validate(MANIFEST)
        self.assertEqual((r['sources'],r['searches'],r['baseline_reuses'],r['new_baseline_calls']),(6,12,56,16))
        self.assertEqual((r['maximum_new_standalone'],r['maximum_new_MACE']),(1000,480))

    def test_all_36_native_force_geometries_use_common_four_modes(self):
        m=a.read_json(MANIFEST)
        for c in m['cases']:
            ca,la=[t for t in m['tasks'] if t['case_id']==c['case_id']]
            self.assertEqual(ca['active_mode_ids'],la['active_mode_ids']);self.assertEqual(len(ca['active_indices']),4)
            for z in ('Ca','La'):
                t=ca if z=='Ca' else la;kin=a.Kinematics(a.read_json(a.verify(t['mapping']))['context'])
                for cell in c['cells'][z].values():
                    self.assertEqual(cell['charge'],t['charge']);self.assertEqual(cell['multiplicity'],1)
                    self.assertTrue(np.allclose(kin.evaluate(cell['full_q'])[1],[r[1:] for r in a.xyz(a.verify(cell['xyz']))],atol=1e-12,rtol=0))
                    self.assertFalse(any(cell['full_q'][i] for i in set(range(t['mode_count']))-set(t['active_indices'])))

    def test_missing_baseline_calls_are_exact_standalone_policy(self):
        m=a.read_json(MANIFEST)
        self.assertEqual({t['case_id'] for t in m['baseline_tasks']},set(a.CASES[-2:]))
        for t in m['baseline_tasks']:
            self.assertEqual(a.verify(t['input']).read_text(),a.xtb.control());cmd=a.xtb.command(t,a.verify(m['executable']))
            self.assertIn('--norestart',cmd);self.assertEqual(cmd[cmd.index('--acc')+1],'0.02');self.assertFalse(t['gradient_requested'])

    def test_corrupted_real_search_scope_rejected(self):
        m=a.read_json(MANIFEST);m['settings']['maxiter']=21
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'corrupted_actual_manifest.json';a.write_new(p,m)
            with self.assertRaises(a.InvalidArtifact):a.validate(p)

    def test_actual_standalone_transfer_algebra_and_gradient_projection(self):
        m=a.read_json(MANIFEST);c=next(c for c in m['cases'] if c['case_id']=='q88jh5-pqq-la_model')
        cell=c['cells']['Ca']['origin'];low={s:v['row'] for s,v in cell['low'].items()}
        comp=a.components(cell['MACE_eV'],low)
        self.assertEqual(comp['GFN2_ALPB_hartree'],low['alpb']['energy_hartree'])
        t=next(t for t in m['tasks'] if t['case_id']==c['case_id'] and t['metal']=='Ca')
        kin=a.Kinematics(a.read_json(a.verify(t['mapping']))['context']);jac=kin.evaluate(cell['full_q'])[3][t['active_indices']]
        f=np.load(a.verify(cell['forces']),allow_pickle=False)
        ng,sg,g=a.project_components(jac,f,low['vacuum']['gradient']['kcal_mol_A'],low['alpb']['gradient']['kcal_mol_A'])
        self.assertTrue(np.isfinite(g).all());self.assertTrue(np.allclose(g,ng+sg,atol=1e-12,rtol=0))
        prior=a.read_json(ROOT/'workspaces/standalone_xtb_20260923/run_v1/COMPARISON_v1.json')
        case=next(x for x in prior['derivative_cases'] if x['case_id']==c['case_id'] and x['accuracy']=='0.02')
        check=next(x for x in case['checks'] if x['quantity']=='Ca_solvent')
        self.assertAlmostEqual(sg[t['active_mode_ids'].index('A/221/chi3')],check['analytic_kcal_mol_radian'],places=10)

    def test_actual_executed_pilot_when_available(self):
        p=MANIFEST.parent/'collection.json'
        if not p.exists():self.skipTest('bounded solvent-aware search has not run')
        r=a.read_json(p);self.assertEqual(r['denominator'],6)
        searches=a.read_json(a.verify(r['searches']))['rows']
        self.assertLessEqual(sum(s['new_standalone_calls'] for s in searches),960)
        for s in searches:
            self.assertLessEqual(len(s['points']),40)
            if s['status']=='finite_proposal_available':
                self.assertTrue(s['candidate']['eligible']);self.assertFalse(s['stationary_minimum_claimed'])
        for row in r['rows']:
            if row['extension_status']=='available':self.assertEqual(a.choose_rows(row['extended_matrix'],row['candidate_ids']),row['extended_pool'])

    def test_corrupted_actual_result_reports_missing_baseline_without_fallback(self):
        p=MANIFEST.parent/'collection.json'
        if not p.exists():self.skipTest('actual pilot result required for corrupted-copy report test')
        r=copy.deepcopy(a.read_json(p));case=r['rows'][0]
        case['base_matrix']['La']['origin']={'status':'unavailable','reason':'explicitly corrupted test copy'}
        case['base_pool']=a.choose_rows(case['base_matrix'],a.BASE)
        case['extended_pool']=None;case['extension_status']='unavailable'
        case['extension_reason']='explicitly corrupted test copy'
        with tempfile.TemporaryDirectory() as d:
            src=Path(d)/'corrupted_actual_collection.json';dst=Path(d)/'report.json'
            a.write_new(src,r);reported=a.report(src,dst)
        self.assertEqual(reported['case_denominator'],6)
        self.assertEqual(len(reported['rows']),6)
        self.assertTrue(all(v is None for v in reported['rows'][0]['values'].values()))
        self.assertIsNone(reported['raw_class_gaps_model_kcal_mol']['solvent_response'])
        self.assertFalse(reported['production_changed'])

if __name__=='__main__':unittest.main()
