"""Admission/recovery tests on copies of the real analytic-pilot manifest."""
import json
from pathlib import Path
import sys
import tempfile
import unittest
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from affordable_common import InvalidArtifact,read_json
from mace_analytic_pilot import validate,core_gate
from mace_hybrid import collect
MANIFEST=ROOT/'workspaces/mace_analytic_20260916/pilot_v1/manifest.json'

@unittest.skipUnless(MANIFEST.exists(),'requires prepared real analytic pilot')
class AnalyticPilotTests(unittest.TestCase):
    def test_real_manifest_and_changed_model_rejection(self):
        m=read_json(MANIFEST);validate(m);m['model']['solvent']='water'
        with self.assertRaisesRegex(InvalidArtifact,'physical model'):validate(m)
    def test_stale_kernel_test_receipt_rejected(self):
        m=read_json(MANIFEST);m['implementation']['mace_analytic.py']['sha256']='0'*64
        with self.assertRaisesRegex(InvalidArtifact,'stale'):validate(m)
    def test_absent_execution_cannot_pass_gate_or_become_zero_score(self):
        with tempfile.TemporaryDirectory() as folder:
            p=Path(folder)/'manifest.json';p.write_bytes(MANIFEST.read_bytes())
            self.assertEqual(core_gate(p)['status'],'failed')
            c=collect(p)
            self.assertEqual(c['status'],'incomplete')
            self.assertIsNone(c['S_kcal_mol']);self.assertIsNone(c['direct'])
            self.assertFalse(c['analytic_comparison']['core_gate'])

if __name__=='__main__':unittest.main()
