"""Large checkpoint admission uses real pinned medium/large artifacts."""
import json
from pathlib import Path
import sys
import tempfile
import unittest
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from affordable_common import InvalidArtifact,read_json,verify
from mace_analytic_pilot import validate,core_gate,LARGE_PROTOCOL
from mace_hybrid import collect
MANIFEST=ROOT/'workspaces/mace_large_20260916/pilot_v2/manifest.json'

@unittest.skipUnless(MANIFEST.exists(),'requires prepared real large pilot')
class LargePilotTests(unittest.TestCase):
    def test_large_variant_accepts_only_its_own_verified_checkpoint(self):
        m=read_json(MANIFEST);validate(m)
        medium=read_json(verify(read_json(verify(m['finite_reference']))['manifest']))
        m['model']['checkpoint']=medium['model']['checkpoint']
        with self.assertRaisesRegex(InvalidArtifact,'checkpoint'):validate(m)
    def test_unknown_model_variant_rejected(self):
        m=read_json(MANIFEST);m['model_variant']='unsupported'
        with self.assertRaisesRegex(InvalidArtifact,'protocol'):validate(m)
    def test_empty_large_execution_remains_unavailable_and_correctly_labeled(self):
        with tempfile.TemporaryDirectory() as folder:
            p=Path(folder)/'manifest.json';p.write_bytes(MANIFEST.read_bytes())
            self.assertEqual(core_gate(p)['status'],'failed')
            c=collect(p);self.assertEqual(c['protocol_id'],LARGE_PROTOCOL)
            self.assertIsNone(c['direct']);self.assertIsNone(c['S_kcal_mol'])
            self.assertEqual(c['analytic_comparison']['comparison_label'],'analytic_large_minus_analytic_medium')

if __name__=='__main__':unittest.main()
