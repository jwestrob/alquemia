"""Real source-state, role and physical-displacement tests; no mock energies."""
from pathlib import Path
import sys
import unittest
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import accommodation_torsion_profiles as at
from affordable_common import read_json,verify,xyz
from mace_site_kinematics import Kinematics
ROOT=Path(__file__).resolve().parents[1]
WORK=ROOT/'workspaces/accommodation_torsion_20260920/prepared_v3'

class TorsionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.d=read_json(WORK/'design.json')

    def test_frozen_scope_and_native_validation_selection(self):
        self.assertEqual(len(self.d['points']),64)
        self.assertEqual(sum(len(c['modes']) for c in self.d['cases']),7)
        dft=read_json(WORK/'dft/manifest.json')
        self.assertEqual(len(dft['all_tasks']),16)
        for t in dft['all_tasks']:
            p=next(p for p in self.d['points'] if p['task_id']==t['point_id'])
            self.assertTrue((p['role']=='extra_acidic_ligand_homolog' and abs(p['angle_radian'])==.2) or
                (p['case_id'] in at.CASES[2:] and p['point']=='origin'))
            self.assertEqual(xyz(verify(t['xyz'])),xyz(verify(p['xyz'])))
        self.assertEqual(len(self.d['reused_DFT_centers']),2)

    def test_whole_group_rotation_preserves_source_bonds_and_frozen_atoms(self):
        for c in self.d['cases']:
            for metal in ('Ca','La'):
                kin=Kinematics(read_json(verify(c['maps'][metal]))['context'])
                for role,mode in c['modes'].items():
                    index=next(i for i,m in enumerate(kin.modes) if m['id']==mode)
                    q=np.zeros(len(kin.modes));q[index]=.4
                    p=kin.evaluate(q)[0];moving=set(kin.modes[index]['moving_indices'])
                    fixed=[i for i in range(len(p)) if i not in moving]
                    np.testing.assert_allclose(p[fixed],kin.positions[fixed],atol=1e-12,rtol=0)
                    self.assertTrue(kin.check(q)['pass'])
                    meta=kin.data['source_atom_metadata']
                    names={meta[i]['atom'] for i in moving if meta[i] is not None}
                    self.assertTrue({'OE1','OE2'}<=names if role=='anchor_glutamate' else {'OD1','OD2'}<=names)

    def test_paired_coordinates_charge_and_no_unrecorded_overlap(self):
        for c in self.d['cases']:
            for p in (p for p in self.d['points'] if p['case_id']==c['case_id'] and p['metal']=='Ca'):
                l=next(x for x in self.d['points'] if x['case_id']==c['case_id'] and x['point']==p['point'] and x['metal']=='La')
                a,b=xyz(verify(p['xyz'])),xyz(verify(l['xyz']))
                self.assertEqual(a[1:],b[1:]);self.assertEqual(a[0][1:],b[0][1:]);self.assertEqual(l['charge']-p['charge'],1)
                self.assertEqual(p['geometry_checks']['pass'],not p['geometry_checks']['new_severe_overlaps'])

    def test_actual_warm_native_MACE_and_archived_control_repetition(self):
        from mace_hybrid import EV_TO_KCAL
        result=read_json(WORK/'mace_1203745/result.json')
        self.assertEqual(result['new_MACE_calls'],64)
        rows={r['task_id']:r for r in result['rows']}
        self.assertEqual(len(rows),64)
        for c in self.d['cases'][:2]:
            for metal,e in c['archived_control_endpoints'].items():
                point=rows[c['case_id']+'__origin__'+metal]
                self.assertEqual(point['status'],'complete')
                self.assertLessEqual(abs(point['energy_eV']-e['native_MACE_energy_eV'])*EV_TO_KCAL,.01)

    def test_original_control_origins_and_plm_unknown_labels(self):
        for c in self.d['cases']:
            if c['case_id'] in at.CASES[:2]:
                for metal,e in c['archived_control_endpoints'].items():
                    self.assertEqual(xyz(verify(e['xyz'])),xyz(verify(c['origins'][metal]['xyz'])))
                    self.assertEqual(e['charge'],c['origins'][metal]['charge'])
            else:self.assertEqual(c['label_scope'],'unknown_PLM_prediction')
        r=at.validate(WORK/'manifest.json');self.assertEqual(r['low']['tasks'],128)
        self.assertEqual(r['dft']['tasks'],16)

if __name__=='__main__':unittest.main()
