"""Real saved endpoints and explicitly constructed distant-ion diagnostics."""
from pathlib import Path
import copy
import sys
import unittest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from affordable_common import InvalidArtifact, read_json, verify, xyz
from mace_file_checks import cached_file_checks
from mace_hybrid import check_atoms
from mace_omol_locality import endpoint
from mace_omol_spectator import expected, SETTINGS

REPORT=ROOT/'workspaces/mace_omol_20260917/intact_report_v1/result.json'


@unittest.skipUnless(REPORT.exists(),'requires actual five-case OMOL report and receipts')
class SpectatorTests(unittest.TestCase):
    @classmethod
    @cached_file_checks
    def setUpClass(cls):
        cls.source,cls.tasks,cls.coords,cls.reused,_=expected(REPORT)

    def test_all_twenty_real_endpoints_retained_with_explicit_sodium(self):
        self.assertEqual(len(self.tasks),20)
        self.assertEqual(len(self.reused),20)
        for t in self.tasks:
            original=xyz(verify(t['source_xyz']));actual=self.coords[t['task_id']]
            self.assertEqual(actual[:-1],original)
            self.assertEqual(actual[-1][0],'Na')
            self.assertEqual(t['spectator_index'],len(actual)-1)
            self.assertGreaterEqual(t['spectator']['minimum_distance_A'],10000.)
            self.assertLessEqual(t['spectator']['direct_Coulomb_score_bound_kcal_mol'],
                                 SETTINGS['maximum_direct_Coulomb_difference_kcal_mol'])
            self.assertEqual(t['state'],check_atoms(actual,t['charge']))
            self.assertEqual(t['state']['all_electron_count']-check_atoms(original,t['charge']-1)['all_electron_count'],10)
            self.assertIn(t['charge'],range(-10,11))

    def test_identical_spectator_for_both_metals_and_both_positions(self):
        for name in self.source['scores']:
            tasks=[t for t in self.tasks if t['case_id']==name]
            self.assertEqual(len(tasks),4)
            self.assertEqual(len({tuple(t['spectator']['xyz_A']) for t in tasks}),1)
            bykey={(t['metal'],t['position']):t for t in tasks}
            for position in ('bound','detached'):
                la,ca=(bykey[(metal,position)] for metal in ('La','Ca'))
                self.assertEqual(la['charge']-ca['charge'],1)
                i=la['metal_index'];a=self.coords[la['task_id']];b=self.coords[ca['task_id']]
                self.assertEqual(a[:i]+a[i+1:],b[:i]+b[i+1:])

    def test_corrupted_copy_of_actual_energy_rejected(self):
        # Corruption regression, not a fabricated successful scientific fixture.
        result=copy.deepcopy(self.source['scores']['GGR_1GLG']['endpoints']['Ca']['bound'])
        result['energy_eV']+=1.
        with self.assertRaisesRegex(InvalidArtifact,'actual execution'):
            endpoint(result)


if __name__=='__main__':unittest.main()
