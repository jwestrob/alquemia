"""Actual pinned union contexts and adaptive selector, no synthetic science."""
from pathlib import Path
import sys
import tempfile
import unittest
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from affordable_common import read_json,verify,xyz,write_new,InvalidArtifact
import union_adaptive as u
from mace_site_kinematics import Kinematics
ROOT=Path(__file__).resolve().parents[1]
RUN=ROOT/'workspaces/union_adaptive_20260923/proposals_v1'

class Fixtures(unittest.TestCase):
    def test_actual_common8_preflight(self):
        v=u.validate(RUN/'manifest.json')
        self.assertEqual((v['cases'],v['optimizer_starts'],v['new_q0_calls']),(8,16,0))

    def test_actual_union_mapping_and_pairing(self):
        m=read_json(RUN/'manifest.json')
        for c in m['cases']:
            ca,la=[next(t for t in m['tasks'] if (t['case_id'],t['metal'])==(c['case_id'],z)) for z in ('Ca','La')]
            self.assertEqual(ca['charge']+1,la['charge'])
            self.assertEqual(ca['active_mode_ids'],la['active_mode_ids'])
            self.assertEqual(read_json(verify(ca['mapping'])),read_json(verify(la['mapping'])))
            self.assertTrue(np.array_equal([a[1:] for a in xyz(verify(ca['xyz']))],[a[1:] for a in xyz(verify(la['xyz']))]))
            k=Kinematics(read_json(verify(ca['mapping']))['context'])
            self.assertTrue(all(k.modes[i]['unit']=='radian' for i in ca['active_indices']))
            self.assertEqual(len(ca['active_indices']),4)

    def test_real_origins_and_units(self):
        m=read_json(RUN/'manifest.json')
        self.assertEqual(m['settings'],u.completion.SETTINGS)
        for t in m['tasks']:
            r=read_json(verify(t['origin_reuse']['point']['MACE']))
            self.assertEqual(r['energy_eV'],t['q0']['components']['MACE_eV'])
            self.assertEqual(r['forces'],t['origin_reuse']['point']['forces'])
            self.assertFalse(any(t['origin_reuse']['point']['full_q']))

    def test_explicit_corrupted_selector_rejected(self):
        m=read_json(RUN/'manifest.json');m['tasks'][0]['active_indices'][0]=0
        with tempfile.TemporaryDirectory() as td:
            p=Path(td)/'corrupted.json';write_new(p,m)
            with self.assertRaises(InvalidArtifact):u.validate(p)

    def test_no_invented_reference_or_default_change(self):
        m=read_json(RUN/'manifest.json')
        for t in m['tasks']:
            self.assertEqual(t['q0_status'],'available')
        self.assertIsNone(m['reference'])
        self.assertFalse(m['production_changed'])

    def test_shared_dispatch_real_minimal_pool(self):
        import nikasha_pool
        path=RUN.parent/'pool_v1/manifest.json'
        direct=u.validate_pool(path);shared=nikasha_pool.validate(path)
        self.assertEqual(direct,shared)
        self.assertEqual((shared['denominator'],shared['prepared'],shared['new_MACE_cells'],shared['new_GFN2_calls']),(8,8,16,64))
        m=read_json(path)
        for c in m['cases']:
            self.assertEqual([q['id'] for q in c['candidates']],['origin','adaptive_Ca','adaptive_La'])
            self.assertEqual(set(c['matrix']['Ca']),set(c['matrix']['La']))

    def test_corrupted_cross_state_rejected(self):
        m=read_json(RUN.parent/'pool_v1/manifest.json');m['tasks'][0]['charge']+=1
        with tempfile.TemporaryDirectory() as td:
            p=Path(td)/'corrupted_pool.json';write_new(p,m)
            with self.assertRaises(InvalidArtifact):u.validate_pool(p)

    def test_actual_final_work_sign_and_complete_denominator(self):
        j=read_json(RUN.parent/'COMPARISON_v1.json')
        final=read_json(verify(j['collection']))
        self.assertEqual((j['denominator'],len(j['rows']),len(final['cases'])),(8,8,8))
        self.assertTrue(all(c['pool']['status']=='available' for c in final['cases']))
        self.assertIsNone(j['new_calibration'])
        self.assertTrue(j['old_band_transfer_only'])
        for r in j['rows']:
            works=r['union_selected_endpoint_work'];scores=r['raw_R_model_kcal_mol']
            for z in ('Ca','La'):
                self.assertAlmostEqual(works[z]['composite_kcal_mol'],works[z]['native_MACE_kcal_mol']+works[z]['solvent_transfer_kcal_mol'],places=8)
            self.assertAlmostEqual(scores['union_adaptive_minimal']-scores['union_static'],works['Ca']['composite_kcal_mol']-works['La']['composite_kcal_mol'],places=6)
        self.assertEqual(j['raw_ordering']['union_adaptive_minimal']['correct_pairwise_directions'],15)
        self.assertEqual(j['counts']['adaptive']['union_adaptive_minimal'],{'correct':8})
        self.assertEqual(j['counts']['released']['union_adaptive_minimal'],{'correct':7,'inconclusive':1})

    def test_exact_canonical_selection_and_pilot_reuse(self):
        p=RUN.parent/'CANONICAL_INPUTS_v1.json';data=read_json(p)
        self.assertEqual(u.declared_population(p,verify(data['calibration'])),('canonical_remaining24',24))
        self.assertEqual(sum(r['role']=='calibration' for r in data['all_cases']),25)
        self.assertEqual({r['case_id'] for r in data['all_cases'] if r['role']!='calibration'},{'1H4I','4MAE','1KB0'})
        self.assertEqual(len(data['cases']),24)
        self.assertEqual(set(data['reuse_case_ids']),set(u.CANONICAL_REUSE))
        # A corrupted real source selector cannot substitute a favorable fold.
        data['cases'][0]['actual_union_case_id']='corrupted_actual_source_id'
        with tempfile.TemporaryDirectory() as td:
            q=Path(td)/'corrupted_inputs.json';write_new(q,data)
            with self.assertRaises(InvalidArtifact):u.declared_population(q,verify(data['calibration']))

    def test_actual_canonical_reference_and_missing_member(self):
        from nikasha_pool_compare import extrema_reference
        j=read_json(RUN.parent/'CANONICAL_REFERENCE_v1.json')
        self.assertEqual((j['denominator'],j['calibration_denominator'],j['crystal_transfer_denominator']),(28,25,3))
        self.assertFalse(j['crystals_or_noncanonical_used_for_fit'])
        for name,r in j['variants'].items():
            self.assertEqual(len(r['rows']),25)
            self.assertFalse(any(c['case_id'] in ('1H4I','4MAE','1KB0') for c in r['rows']))
            replay=extrema_reference(r['rows'],name,'Nikasha_union_adaptive_minimal_canonical25_v1')
            for key in ('status','bands','gap_model_kcal_mol','available_calibration'):
                self.assertEqual(replay[key],r[key])
            # Explicitly corrupted real result: a deleted energy cannot become zero.
            broken=[dict(c) for c in r['rows']];broken[0]['R_model_kcal_mol']=None
            absent=extrema_reference(broken,name)
            self.assertEqual(absent['status'],'unavailable_calibration_member')
            self.assertIsNone(absent['bands'])

if __name__=='__main__':unittest.main()
