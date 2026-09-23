"""Real-source geometry/input/parser checks; no molecular engine is run here."""
import copy
from pathlib import Path
import sys
import numpy as np
import unittest
import tempfile
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import native_solvent_force as f
from affordable_common import read_json,verify,xyz,InvalidArtifact,write_new
from mace_site_kinematics import Kinematics
ROOT=Path(__file__).resolve().parents[1]
RUN=ROOT/'workspaces/native_solvent_force_20260923/run_v1'


class RealFixtures(unittest.TestCase):
    def test_fixed_actual_displacement_manifest(self):
        result=f.validate(RUN/'stage0/manifest.json')
        assert result['tasks']==32 and result['blocked']==0
        d=read_json(RUN/'design.json')
        assert d['new_molecular_calls']==96 and d['external_q0_cells']==8
        assert {(c['case_id'],c['offset_radian'],c['metal'],c['medium']) for c in d['cells']}=={
            (case,step,z,s) for case in f.CASES for step in f.OFFSETS for z in ('Ca','La') for s in ('vacuum','alpb')}


    def test_exact_paired_coordinates_and_physical_caps(self):
        d=read_json(RUN/'design.json')
        for cid in f.CASES:
            for step in f.OFFSETS:
                ca,la=[next(c for c in d['cells'] if (c['case_id'],c['offset_radian'],c['metal'],c['medium'])==(cid,step,z,'vacuum')) for z in ('Ca','La')]
                assert ca['charge']+1==la['charge']
                assert np.array_equal([a[1:] for a in xyz(verify(ca['xyz']))],[a[1:] for a in xyz(verify(la['xyz']))])
                assert ca['geometry_checks']['pass']
                assert ca['geometry_checks']['bond_max_error_A']<1e-9
                assert ca['geometry_checks']['maximum_heavy_displacement_A']<.01
                k=Kinematics(read_json(verify(ca['physical_mapping']))['context'])
                assert np.allclose(k.evaluate(ca['full_q'])[1],[a[1:] for a in xyz(verify(ca['xyz']))],atol=1e-12,rtol=0)


    def test_terminal_role_and_rms_measure_are_declared(self):
        d=read_json(RUN/'design.json')
        for key,meta in d['maps'].items():
            assert meta['source_role']['resname']=='GLU'
            assert meta['mode_id'].endswith('/chi3')
            assert meta['RMS_moving_heavy_speed_A_radian']>meta['RMS_all_heavy_speed_A_radian']>0
            assert len(meta['moving_heavy_ids'])==len(meta['moving_heavy_indices'])
        assert not d['maps']['4MAE__Ca']['adaptive_selected']
        assert d['maps']['q88jh5-pqq-la_model__Ca']['adaptive_selected']


    def test_restart_recipe_preserves_native_solver(self):
        d=read_json(RUN/'design.json')
        for c in d['cells']:
            initial=f.recipe(c['charge'],1,c['medium'],0);later=f.recipe(c['charge'],1,c['medium'],2)
            assert initial.replace(' NoAutostart','')==later
            assert 'UseXTBMixer true' in later and 'MaxIter 500' in later and 'SmearTemp 300' in later
            assert 'Tight' not in later and 'EnGrad' not in later
            assert ('ALPB(Water)' in later)==(c['medium']=='alpb')


    def test_corrupted_real_manifest_rejected(self):
        # Explicitly corrupted copy of actual prepared inputs, not simulated chemistry.
        m=read_json(RUN/'stage0/manifest.json');m['tasks'][0]['charge']+=1
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'corrupted.json';write_new(path,m)
            with self.assertRaisesRegex(InvalidArtifact,'state or coordinate'):f.validate(path)


    def test_actual_old_gradient_matches_mapping_chain_rule(self):
        from accommodation_response import parse_gradient
        p=ROOT/'workspaces/accommodation_response_20260920/prepared_v2/followup/manifest.json'
        m=read_json(p)
        task=next(t for t in m['tasks'] if (t['case_id'],t['metal'],t['medium'],t['kind'])==('4MAE','La','alpb','center'))
        pin=f.completed(p,task['task_id']);assert pin is not None
        parsed=parse_gradient(task,pin);k=Kinematics(read_json(verify(task['physical_mapping']))['context'])
        j=k.evaluate(np.zeros(len(k.modes)))[3]
        manual=np.einsum('mij,ij->m',j,np.asarray(parsed['gradient_kcal_mol_A']))
        assert np.array_equal(manual,parsed['projected_gradient_kcal_per_unit'])
        assert parsed['quantity']=='gradient_not_force'

    def test_actual_continuations_and_reported_stage(self):
        rows=[read_json(RUN/('stage'+str(i))/'collection.json')['rows'] for i in range(3)]
        assert all(len(r)==32 and all(x['status']=='complete' for x in r) for r in rows)
        for stage in (1,2):
            previous={r['cell_id']:r for r in rows[stage-1]}
            m=read_json(RUN/('stage'+str(stage))/'manifest.json')
            for t in m['tasks']:
                source=previous[t['cell_id']]['seed_after']
                assert t['seed_source']=={'xtbw':source['preserved_after'],'gbw':source['preserved_gbw_after']}
        result=read_json(RUN/'COMPARISON_v1.json')
        assert result['qualified_cases']==0 and result['passed_derivative_quantities']==7
        assert all(c['numerical_pass'] and not c['force_pass'] for c in result['cases'])

    def test_actual_derivatives_sign_and_single_unit_conversion(self):
        result=read_json(RUN/'COMPARISON_v1.json')
        for case in result['cases']:
            q={r['quantity']:r for r in case['derivatives']}
            for z in ('Ca','La'):
                for medium in ('vacuum','alpb'):
                    r=q[z+'_'+medium];energies={float(k):v for k,v in r['raw']['energies_hartree'].items()}
                    direct=(energies[.0005]-energies[-.0005])/.001*f.HA_TO_KCAL
                    self.assertAlmostEqual(direct,r['centered_half_h_kcal_mol_radian'],places=10)
            expected=q['Ca_alpb']['analytic_kcal_mol_radian']-q['Ca_vacuum']['analytic_kcal_mol_radian']-q['La_alpb']['analytic_kcal_mol_radian']+q['La_vacuum']['analytic_kcal_mol_radian']
            self.assertAlmostEqual(expected,q['Ca_minus_La_solvent']['analytic_kcal_mol_radian'],places=12)

    def test_missing_actual_gradient_stays_unavailable(self):
        # Explicitly corrupted status in a copy of a real completed collection.
        centers=ROOT/'workspaces/native_pool_continuation_20260923/run_v1/stage2/collection.json'
        d=read_json(centers)
        r=next(r for r in d['rows'] if (r['case_id'],r['candidate'],r['metal'],r['medium'])==('4MAE','origin','Ca','vacuum'))
        r['status']='audit_failed';r['reason']='test corruption: gradient withheld'
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'corrupted_centers.json';out=Path(tmp)/'result.json';write_new(p,d)
            result=f.compare(RUN/'design.json',p,out)
            case=next(c for c in result['cases'] if c['case_id']=='4MAE')
            assert len(case['derivatives'])==7 and not case['numerical_pass']
            row=next(q for q in case['derivatives'] if q['quantity']=='Ca_minus_La_solvent')
            assert row['status']=='unavailable' and row['raw'] is None

if __name__=="__main__":unittest.main()
