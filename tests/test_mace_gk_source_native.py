"""Real source-only native inputs and executed component agreement."""
import copy
import json
from pathlib import Path
import sys
import tempfile
import unittest
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from affordable_common import InvalidArtifact,read_json,verify
from mace_gk_source_native import validate,parse
W=ROOT/'workspaces/mace_omol_20260917/gk_source_native_v2'

class NativeSource(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if not (W/'manifest.json').exists():raise unittest.SkipTest('real component preparation unavailable')
        cls.m=validate(W/'manifest.json')

    def test_exact_archived_inputs_and_all_five_geometries(self):
        self.assertEqual(len(self.m['tasks']),15)
        for t in self.m['tasks']:
            r=read_json(verify(t['parent_result']));cmd=r['receipt']['command']
            self.assertEqual(verify(t['xyz']).read_bytes(),Path(cmd[1]).read_bytes())
            self.assertEqual(verify(t['mask']).read_bytes(),Path(cmd[2]).read_bytes())
            self.assertEqual(verify(t['overrides']).read_bytes(),Path(cmd[3]).read_bytes())

    def test_missing_real_task_rejected(self):
        m=copy.deepcopy(self.m);m['tasks'].pop()
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'corrupted_real_manifest.json';p.write_text(json.dumps(m))
            with self.assertRaises(InvalidArtifact):validate(p)

    def test_real_native_agreement_and_corrupted_completion(self):
        paths=list(W.glob('collection_job_*.json'))
        if not paths:self.skipTest('native source-only calculations have not run')
        r=read_json(paths[-1]);self.assertTrue(r['complete']);self.assertTrue(r['checks_pass'])
        self.assertEqual(r['native_energy_calls'],15)
        for t in self.m['tasks']:
            out=read_json(verify(r['tasks'][t['task_id']]['result']));text=verify(out['receipt']['log']).read_text()
            self.assertEqual(parse(text,t)['energies'],out['energies'])
            with self.assertRaisesRegex(InvalidArtifact,'incomplete'):
                parse(text.replace('ALQUEMIA_DENSITY_GK_COMPLETE','CORRUPTED_COMPLETION'),t)

if __name__=='__main__':unittest.main()
