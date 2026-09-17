"""Real native preparations/results; damaged copies only for rejection tests."""
from pathlib import Path
import copy
import sys
import unittest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from affordable_common import InvalidArtifact,read_json,verify
from mace_tinker_framework_solver import validate,parse_parameters,check_parameters,parse_energy

W=ROOT/'workspaces/mace_omol_20260917/tinker_framework_solver_v1'


class NativePreparation(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if not (W/'manifest.json').exists():
            raise unittest.SkipTest('real native preparation unavailable')
        cls.manifest=validate(W/'manifest.json')

    def test_real_parameters_geometry_and_masks(self):
        for task in self.manifest['tasks']:
            receipt=read_json(verify(task['preflight_receipt']))
            params=parse_parameters(verify(receipt['log']).read_text())
            self.assertTrue(check_parameters(params,task,read_json(verify(task['mapping']))))

    def test_damaged_real_parameter_inventory_rejected(self):
        task=self.manifest['tasks'][0]
        text=verify(read_json(verify(task['preflight_receipt']))['log']).read_text()
        lines=text.splitlines()
        lines=[line for line in lines if not (line.split()[:2]==['PREP','1'])]
        with self.assertRaisesRegex(InvalidArtifact,'incomplete native parameter'):
            parse_parameters('\n'.join(lines))

    def test_altered_real_source_mask_rejected(self):
        task=next(t for t in self.manifest['tasks'] if t['frozen_indices'])
        params=copy.deepcopy(read_json(verify(task['parameters'])))
        params['atoms'][task['frozen_indices'][0]-1]['response_allowed']=True
        with self.assertRaisesRegex(InvalidArtifact,'source mask'):
            check_parameters(params,task,read_json(verify(task['mapping'])))

    def test_parameter_initialization_is_not_energy_evidence(self):
        task=self.manifest['tasks'][0]
        text=verify(read_json(verify(task['preflight_receipt']))['log']).read_text()
        with self.assertRaisesRegex(InvalidArtifact,'incomplete/nonconverged'):
            parse_energy(text,task)


class NativeEnergy(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if not (W/'manifest.json').exists():
            raise unittest.SkipTest('real native preparation unavailable')
        cls.manifest=validate(W/'manifest.json')
        cls.finished=[t for t in cls.manifest['tasks'] if (W/'tasks'/t['task_id']/'result.json').exists()]
        if len(cls.finished)!=12:
            raise unittest.SkipTest('twelve actual native energy attempts have not completed')

    def test_actual_native_results_and_nonconvergence(self):
        for task in self.finished:
            result=read_json(W/'tasks'/task['task_id']/'result.json')
            text=verify(result['receipt']['log']).read_text()
            if result['status']=='computed_framework_control':
                parsed=parse_energy(text,task)
                self.assertEqual(parsed['energies_kcal_mol'],result['energies_kcal_mol'])
                self.assertLessEqual(parsed['max_frozen_dipole_eA'],1e-12)
                self.assertIsNone(parsed['numerical_score'])
                with self.assertRaisesRegex(InvalidArtifact,'incomplete/nonconverged'):
                    parse_energy(text.replace('ALQUEMIA_ENERGY_COMPLETE','DAMAGED_COMPLETION'),task)
            else:
                with self.assertRaises(InvalidArtifact):parse_energy(text,task)


if __name__=='__main__':unittest.main()
