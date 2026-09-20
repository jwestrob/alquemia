"""Real archived input invariants; scientific outputs are never synthesized."""
from pathlib import Path
import sys,tempfile,unittest
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
import compact_solvation as c
from affordable_common import InvalidArtifact,read_json,write_new,verify,xyz

class CompactSolvationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.manifest=ROOT/'workspaces/compact_solvation_20260920/pilot_v1/manifest.json'
    def test_real_fixed_eight_task_scope(self):
        r=c.validate(self.manifest);self.assertEqual((r['tasks'],r['reused']),(8,0))
        m=read_json(self.manifest)
        self.assertEqual(len({t['scientific_key'] for t in m['tasks']}),8)
        for t in m['tasks']:
            self.assertEqual(xyz(verify(t['xyz'])),xyz(verify(t['source_endpoint']['xyz'])))
            self.assertEqual(t['charge'],t['source_endpoint']['charge'])
            self.assertIn('SmearTemp 300',verify(t['input']).read_text())
    def test_changed_electronic_state_rejected(self):
        m=read_json(self.manifest);m['all_tasks'][0]['charge']+=1
        with tempfile.TemporaryDirectory() as td:
            p=Path(td)/'explicitly_corrupted_real_charge.json';write_new(p,m)
            with self.assertRaisesRegex(InvalidArtifact,'electronic state'):c.validate(p)
    def test_cache_contains_solvent_and_numerical_solver(self):
        m=read_json(self.manifest);m['all_tasks'][0]['medium']='alpb'
        with tempfile.TemporaryDirectory() as td:
            p=Path(td)/'explicitly_corrupted_real_solvent.json';write_new(p,m)
            with self.assertRaises(InvalidArtifact):c.validate(p)
    def test_native_and_numerical_check_are_distinct(self):
        t=read_json(self.manifest)['all_tasks'][0]
        primary=c.input_text(t['charge'],1,t['medium'],'native')
        check=c.input_text(t['charge'],1,t['medium'],'ordinary_tight')
        self.assertIn('UseXTBMixer true',primary);self.assertIn('UseXTBMixer false',check)
        self.assertIn('TightSCF',check);self.assertNotIn('TightSCF',primary)
        self.assertTrue(all('SmearTemp 300' in x and 'ReadXTBParam false' in x for x in (primary,check)))
    def test_actual_native_parameters_and_charge_closure(self):
        m=read_json(self.manifest)
        for t in m['tasks']:
            pin=c.completed(self.manifest,t['task_id']);self.assertIsNotNone(pin)
            audit=c.diagnostics(pin,t)
            self.assertEqual(audit['charge_sanity_status'],'pass')
            self.assertAlmostEqual(audit['charge_sum_e'],t['charge'],delta=5e-4)
            self.assertIsNone(audit['ALPB_component_hartree'])
            self.assertEqual(sum(audit['metal_reference_occupation']),2 if t['metal']=='Ca' else 3)
    def test_actual_failed_numerical_endpoint_has_no_energy(self):
        p=ROOT/'workspaces/compact_solvation_20260920/numerical_v1/manifest.json'
        for t in read_json(p)['tasks']:self.assertIsNone(c.completed(p,t['task_id']))
if __name__=='__main__':unittest.main()
