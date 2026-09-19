"""Pinned PQQ sources, including two actual vicinal disulfides."""
from pathlib import Path
import sys
import unittest
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from affordable_common import read_json,verify,xyz,InvalidArtifact
from second_shell_context import parent_state,expansion
from environment_context_chemistry import complete_expansion
from environment_pqq_context import validate
from environment_pqq_report import summarize


class ContextChemistry(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.inv=read_json(ROOT/'workspaces/environment_pqq_20260919/inventory_v3.json')
        cls.cfg=read_json(verify(cls.inv['static_config']))

    def test_old_23_supported_contexts_are_exactly_unchanged(self):
        count=0
        for row in self.inv['cases']:
            if row['disulfide_closures']:continue
            state=parent_state(row['source'],self.cfg['topology'])
            self.assertEqual(expansion(state),complete_expansion(state));count+=1
        self.assertEqual(count,23)

    def test_complete_source_disulfides_retained_without_sulfur_caps(self):
        count=0
        for row in self.inv['cases']:
            if not row['disulfide_closures']:continue
            state=parent_state(row['source'],self.cfg['topology']);rows,audit=complete_expansion(state);count+=1
            ss=[b for b in audit['mapping']['retained_bonds'] if b['a']['atom']=='SG' and b['b']['atom']=='SG']
            self.assertEqual(len(ss),1)
            self.assertFalse(any(a['retained']['element']=='S' for a in audit['mapping']['cut_bonds_and_caps']))
            self.assertEqual(rows['Ca'][1:],rows['La'][1:])
            for metal in ('Ca','La'):
                for old,new in audit['core_to_context'].items():self.assertEqual(rows[metal][new],state['original'][metal][old])
        self.assertEqual(count,2)

    def test_corrupted_real_disulfide_distance_is_rejected(self):
        row=next(r for r in self.inv['cases'] if r['disulfide_closures']);state=parent_state(row['source'],self.cfg['topology'])
        source=row['disulfide_closures'][0]['atoms'][0];key=(source['chain_index'],source['residue_index'],source['atom'])
        # Explicitly corrupted in-memory copy of a real source; no energy/label.
        state['graph'].atoms[key].pos.x+=5
        with self.assertRaises(InvalidArtifact):complete_expansion(state)

    def test_corrupted_real_disulfide_h_inventory_is_rejected(self):
        row=next(r for r in self.inv['cases'] if r['disulfide_closures']);state=parent_state(row['source'],self.cfg['topology'])
        source=row['disulfide_closures'][0]['atoms'][0];rkey=(source['chain_index'],source['residue_index'])
        h=next(h for h,p in state['graph'].hparents.items() if p[:2]==rkey and p[2]=='CB')
        del state['graph'].hparents[h]
        with self.assertRaises(InvalidArtifact):complete_expansion(state)

    def test_complete_finite_panel_and_exact_transfer_reuse(self):
        path=ROOT/'workspaces/environment_pqq_20260919/prepared_v1/manifest.json';m=read_json(path)
        self.assertEqual(validate(path)['new_native_MACE_endpoints'],52)
        self.assertEqual(len([r for r in m['cases'] if r['source']['evaluation_role']=='calibration']),25)
        self.assertEqual({r['source']['case'] for r in m['cases'] if r['source']['evaluation_role']!='calibration'}, {'1H4I','4MAE','1KB0'})
        for r in m['reused']:
            old=read_json(verify(r['accepted']['manifest']));task=next(t for t in old['tasks'] if t['task_id']==r['accepted']['task_id'])
            self.assertEqual(task['xyz']['sha256'],r['task']['xyz']['sha256']);self.assertEqual(task['charge'],r['task']['charge'])

    def test_actual_endpoint_algebra_and_frozen_calibration_bands(self):
        result=summarize(ROOT/'workspaces/environment_pqq_20260919/prepared_v1/result_1202474.json')
        self.assertEqual(result['bands']['Ca_max_model_kcal_mol'],-405347.73205159395)
        self.assertEqual(result['bands']['La_min_model_kcal_mol'],-405345.416592485)
        self.assertAlmostEqual(result['context_gap_model_kcal_mol'],2.3154591089696623,places=8)
        self.assertEqual(sum(r['expected_region_pass'] for r in result['rows'] if r['role']=='calibration'),25)
        self.assertEqual(sum(r['expected_region_pass'] for r in result['rows'] if r['role']!='calibration'),3)
        self.assertTrue(all(r['retained_original_coordinates_exact'] for r in result['rows']))
        self.assertEqual(result['cost']['new_endpoints'],52)
        self.assertEqual(result['cost']['reused_context_endpoints'],4)


if __name__=='__main__':unittest.main()
