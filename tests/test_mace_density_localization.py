"""Component localization checked against actual independent native collections."""
from pathlib import Path
import sys
import unittest
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from affordable_common import read_json,verify
from mace_density_localization import reconstruct
BASE=ROOT/'workspaces/mace_omol_20260917'

class Localization(unittest.TestCase):
    def test_actual_atom_residue_and_rigid_closure(self):
        p=BASE/'ggr_density_localization_v1/result.json'
        if not p.exists():self.skipTest('actual localization not executed')
        r=read_json(p)
        self.assertTrue(r['checks_pass']);self.assertEqual(len(r['cases']),3)
        self.assertLessEqual(max(abs(x['error_kcal']) for x in r['checks']),1e-7)
        self.assertIsNone(r['new_score']);self.assertIsNone(r['calibrated_class'])
        for c in r['cases'].values():
            for name,key in [('direct','direct_density_kcal'),('induction','environment_induction_transfer_kcal')]:
                self.assertAlmostEqual(sum(x['R_components_kcal'][name] for x in c['residues'].values()),c['reported_components'][key],places=7)

    def test_reconstruct_actual_failed_structure_without_solver(self):
        p=BASE/'trial_gk_expansion_native_v1/collection_job_1201162.json'
        if not p.exists():self.skipTest('actual native collection unavailable')
        r=reconstruct(p,'GGR_2FW0','primary')
        archived=read_json(BASE/'ggr_density_localization_v1/result.json')['cases']['GGR_2FW0']
        self.assertEqual(r,archived)

if __name__=='__main__':unittest.main()
