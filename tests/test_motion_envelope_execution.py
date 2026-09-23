"""Actual68 native tasks and136 strict scalar inputs; no fabricated energies."""
from pathlib import Path
import copy
import sys
import tempfile
import unittest
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from affordable_common import InvalidArtifact,read_json,verify,write_new
import motion_envelope_run as run
import motion_envelope_scalar as scalar
BASE=ROOT/'workspaces/motion_envelope_20260923'
class EnvelopeExecution(unittest.TestCase):
    def test_actual_native_only_searches_validate_without_solvent(self):
        path=BASE/'pilot34_searches_v1/manifest.json';v=run.validate(path);m=read_json(path)
        self.assertEqual(v['bounded_searches'],68);self.assertFalse(v['composite_available'])
        self.assertEqual(m['settings']['optimizer_ftol'],1e-8)
        self.assertTrue(all(set(t['q0']['components'])=={'MACE_eV'} and t['q0']['low'] is None for t in m['tasks']))
        self.assertEqual(len({t['source_case_id']for t in m['tasks']}),34)
    def test_strict136_real_inputs_validate(self):
        path=BASE/'pilot34_scalar_origins_v1/manifest.json';scalar.validate(path);m=read_json(path)
        self.assertEqual(len(m['tasks']),136);self.assertEqual(m['reused'],{})
        for t in m['tasks']:
            self.assertEqual(verify(t['input']).read_text(),scalar.recipe(t['charge'],t['medium'],'fresh'))
            self.assertEqual(verify(t['source_xyz']).read_bytes(),verify(t['xyz']).read_bytes())
    def test_real_qualified_strict_output_is_accepted(self):
        q=read_json(ROOT/'workspaces/strict_native_pool_20260923/run_v1/COMPARISON.json')
        pin=q['cells'][0]['branches']['fresh']['actual'];m=read_json(verify(pin['manifest']));t=next(t for t in m['tasks']if t['task_id']==pin['task_id'])
        result=scalar.endpoint(verify(pin['manifest']),t)
        self.assertEqual(result['status'],'complete');self.assertEqual(result['energy_hartree'],pin['energy_hartree'])
        self.assertEqual(result['observed_TolE_hartree'],1e-10)
    def test_corrupted_actual_native_task_is_rejected(self):
        m=read_json(BASE/'pilot34_searches_v1/manifest.json');m['tasks'][0]['charge']+=1
        with tempfile.TemporaryDirectory()as td:
            p=Path(td)/'corrupted_charge.json';write_new(p,m)
            with self.assertRaisesRegex(InvalidArtifact,'actual native force task differs'):run.validate(p)
if __name__=='__main__':unittest.main()
