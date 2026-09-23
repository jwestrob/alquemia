"""Real archived-H preparation and finite-manifest checks; no model calls."""
import json
from pathlib import Path
import sys
import unittest
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
import preparation_coverage as prep
import preparation_coverage_score as score
from affordable_common import read_json,verify,xyz
import adaptive_completion
import union_adaptive

class CoverageTests(unittest.TestCase):
    def setUp(self):
        self.path=ROOT/'workspaces/preparation_coverage_20260923/prepared_v1/PREPARATION.json'
        self.p=score.preparation(self.path)

    def test_archived_heavy_atoms_and_inventory(self):
        r=self.p['local']['archived_H_repair'];old=prep.atoms(verify(r['original_protonated']));new=prep.atoms(verify(r['continued_PDB']))
        self.assertEqual([a[:-1] for a in old],[a[:-1] for a in new]);self.assertEqual(len(old),9582)
        self.assertEqual([a for a in old if a[-2]!='H'],[a for a in new if a[-2]!='H'])
        self.assertTrue(any(a!=b for a,b in zip(old,new) if a[-2]=='H'))
        self.assertEqual(self.p['new_molecular_calls'],0)

    def test_fixed_union_state_and_coordinates(self):
        source=read_json(verify(self.p['existing_union_preparation']));canonical=next(c for c in source['cases'] if c['case_id']==self.p['frozen_group']['canonical_case_id'])
        for k in ('protein_key','state_key'):self.assertEqual(self.p['union'][k],canonical[k])
        for kind in ('local','union'):
            r=self.p[kind]['representations']['context'];a=xyz(verify(r['endpoints']['Ca']['xyz']));b=xyz(verify(r['endpoints']['La']['xyz']))
            self.assertEqual(a[1:],b[1:]);self.assertEqual(a[0][1:],b[0][1:]);self.assertTrue(self.p['original_failure_preserved'])

    def test_finite_actual_origin_inputs(self):
        path=ROOT/'workspaces/preparation_coverage_20260923/origins_v1/manifest.json';v=score.validate_origins(path)
        self.assertEqual((v['MACE_calls'],v['GFN2_calls']),(4,8))
        low=read_json(path.parent/'solvent/shard_0/manifest.json');self.assertEqual(len(low['tasks']),8)
        from compact_solvation import input_text
        for t in low['tasks']:
            self.assertTrue(verify(t['xyz']).is_relative_to(path.parent/'solvent/shard_0'))
            self.assertEqual(verify(t['input']).read_text(),input_text(t['charge'],t['multiplicity'],t['medium'],'native').replace('%scf\n','%scf\n MaxIter 500\n'))
        self.assertEqual(read_json(path.parent/'solvent/shard_0/PREFLIGHT.json')['tasks'],8)

    def test_original_optimizer_isolation(self):
        before=json.dumps(adaptive_completion.SETTINGS,sort_keys=True);u,p=score.engine()
        self.assertEqual(u.SETTINGS,adaptive_completion.SETTINGS)
        self.assertAlmostEqual(u.SETTINGS['optimizer_ftol'],3.6749322179347735e-11,places=20)
        self.assertNotEqual(u.PROTOCOL,union_adaptive.PROTOCOL)
        self.assertEqual(json.dumps(adaptive_completion.SETTINGS,sort_keys=True),before)

    def test_real_union_kinematics_before_force_selection(self):
        from second_shell_context import parent_state
        from coordination_preparation_context import geometry
        from mace_site_kinematics import Kinematics
        u=self.p['union'];cfg=read_json(verify(self.p['existing_union_preparation']))['config']
        state=parent_state(u['original_core'],cfg['topology'],require_endpoint_receipts=False)
        rep=u['representations']['context'];context=read_json(verify(rep['preparation']));maps=[]
        for z in ('Ca','La'):
            atom=xyz(verify(rep['endpoints'][z]['xyz']));core=xyz(verify(u['original_core']['endpoints'][z]['xyz']))
            g=geometry(state,context,atom,core);maps.append(g);k=Kinematics(g['context'])
            self.assertTrue(np.allclose(k.evaluate(np.zeros(len(k.modes)))[1],[a[1:] for a in atom],atol=1e-12,rtol=0))
            self.assertGreaterEqual(sum(v['unit']=='radian' for v in k.modes),4)
        self.assertEqual(maps[0],maps[1])

    def test_actual_scored_origin_and_proposal_reuse(self):
        run=ROOT/'workspaces/preparation_coverage_20260923'
        self.assertEqual(score.validate(run/'proposals_v1/manifest.json')['optimizer_starts'],2)
        self.assertEqual(score.validate_pool(run/'pool_v1/manifest.json')['new_GFN2_calls'],8)
        for stage in ('origins_v1','pool_v1'):
            col=read_json(run/stage/'collection_final.json');self.assertEqual(col['GFN2_complete'],8)
            self.assertTrue(all(c['pool']['status']=='available' for c in col['cases']))

    def test_actual_score_sign_reference_and_failure_history(self):
        run=ROOT/'workspaces/preparation_coverage_20260923';r=read_json(run/'COMPARISON_v1.json')
        self.assertFalse(r['historical225_modified']);self.assertTrue(r['historical_preparation_failure_retained'])
        self.assertEqual(r['independent_biological_observations_added'],0)
        for row in r['rows'].values():
            self.assertEqual(row['decision'],'La-supported')
            self.assertAlmostEqual(row['R_model_kcal_mol'],row['components']['native_R_model_kcal_mol']+row['components']['solvation_delta_R_kcal_mol'],places=8)
        work=r['rows']['union_adaptive']['selected_work']
        self.assertAlmostEqual(r['rows']['union_adaptive']['R_model_kcal_mol']-r['rows']['union_static']['R_model_kcal_mol'],work['Ca']['composite_kcal_mol']-work['La']['composite_kcal_mol'],places=8)
        self.assertEqual(r['optimizers']['Ca']['optimizer']['ftol_hartree_equivalent'],adaptive_completion.SETTINGS['optimizer_ftol'])

if __name__=='__main__':unittest.main()
