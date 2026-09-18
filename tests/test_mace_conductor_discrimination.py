"""Algebra and manifest checks on real archived calculations; no solver mock."""
import copy
import json
from pathlib import Path
import sys
import unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from affordable_common import cache_key,read_json,verify
from mace_conductor_discrimination import corrected_contrast,CASES

ROOT=Path(__file__).resolve().parents[1]/'workspaces/mace_omol_20260917'

class RealArchiveTests(unittest.TestCase):
    def test_identity_and_sign_from_actual_native_values(self):
        p=ROOT/'frozen_response_v2/collection_job_1201310.json'
        if not p.exists(): self.skipTest('real native functional archive unavailable')
        co=read_json(p)
        gk={s:co['rows']['GGR_2FW0_'+s]['GK_cross_kcal'] for s in ('Ca','La')}
        old=read_json(ROOT/'trial_gk_expansion_native_v1/collection_job_1201162.json')['variants']['primary']['cases']['GGR_2FW0']['R_kcal']
        self.assertEqual(corrected_contrast(old,gk,gk),old)
        # The native permanent distribution is an actual computed comparator,
        # used here solely to test the signed algebra, not claimed as a model.
        permanent={s:co['rows']['GGR_2FW0_'+s]['static_GK_kcal']['1'] for s in ('Ca','La')}
        expected=(permanent['Ca']-permanent['La'])-(gk['Ca']-gk['La'])
        self.assertAlmostEqual(corrected_contrast(old,permanent,gk)-old,expected,places=9)
        for s in ('Ca','La'):
            r=co['rows']['GGR_2FW0_'+s]
            self.assertAlmostEqual(r['static_GK_kcal']['2']-r['static_GK_kcal']['3'],gk[s],places=12)
            self.assertTrue(all(c['pass_'] for c in r['checks']))

    def test_real_prepared_scope_and_pair_cancellation(self):
        p=ROOT/'conductor_discrimination_v1/manifest.json'
        if not p.exists():self.skipTest('prepared real comparison unavailable')
        m=read_json(p)
        self.assertEqual(tuple(g['case_id'] for g in m['groups']),CASES)
        self.assertEqual(m['requested_forward_solves'],56)
        self.assertEqual(sum(t['reuse'] is not None for g in m['groups'] for t in g['endpoints'].values()),4)
        for g in m['groups']:
            old=g['old'];ca,la=(old['endpoints'][s] for s in ('Ca','La'))
            self.assertEqual(ca['native_environment_static_components'],la['native_environment_static_components'])
            self.assertEqual(ca['induction_environment_only_kcal'],la['induction_environment_only_kcal'])
            self.assertAlmostEqual(sum(old['components_R_kcal'].values()),old['R_kcal'],places=8)
            self.assertAlmostEqual(ca['total_kcal']-la['total_kcal'],old['R_kcal'],places=8)

    def test_real_cache_includes_method_source_cavity_and_grid(self):
        p=ROOT/'conductor_discrimination_v1/manifest.json'
        if not p.exists():self.skipTest('prepared real comparison unavailable')
        m=read_json(p);unsigned={k:v for k,v in m.items() if k!='cache_key'}
        self.assertEqual(cache_key(unsigned),m['cache_key'])
        for field in ('model','grids','groups','module','implementation'):
            changed=copy.deepcopy(unsigned);del changed[field]
            self.assertNotEqual(cache_key(changed),m['cache_key'])

    def test_actual_solver_outputs_and_complete_signed_comparison(self):
        p=ROOT/'conductor_discrimination_v1/collection_final.json'
        if not p.exists():self.skipTest('real completed continuum calculations unavailable')
        co=read_json(p);self.assertTrue(co['complete']);count=0
        for case,r in co['cases'].items():
            for grid,g in r['grids'].items():
                for metal in ('Ca','La'):
                    a=g['rows'][metal+'_average'];b=g['rows'][metal+'_difference']
                    for s in (a,b):
                        self.assertEqual(s['status'],'computed')
                        self.assertTrue(s['forward_solve_started'])
                        self.assertTrue(s['contraction_pass']);count+=1
                    self.assertEqual(g['C_CPCM_kcal'][metal],a['energy_kcal_mol']-b['energy_kcal_mol'])
                expected=(g['C_CPCM_kcal']['Ca']-g['C_CPCM_kcal']['La'])-(g['C_GK_kcal']['Ca']-g['C_GK_kcal']['La'])
                self.assertAlmostEqual(g['R_new_kcal']-r['old_R_kcal'],expected,places=9)
        self.assertEqual(count,56)
        for grid,rows in co['margins'].items():
            self.assertEqual(len(rows),12)
            for row in rows:
                a=co['cases'][row['La_case']]['grids'][grid]['R_new_kcal']
                b=co['cases'][row['Ca_case']]['grids'][grid]['R_new_kcal']
                self.assertEqual(row['new_margin_kcal'],a-b)
                self.assertEqual(row['new_direction_pass'],a-b>.02)

if __name__=='__main__':unittest.main()
