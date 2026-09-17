"""Integrity checks use pinned real rotation manifests, never synthetic science."""
import copy
from pathlib import Path
import sys
import unittest
import tempfile
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from affordable_common import InvalidArtifact, read_json
from mace_rotation import validate, collect_rotation
MANIFEST=ROOT/'workspaces/mace_rotation_20260916/audit_v1/original/manifest.json'

@unittest.skipUnless(MANIFEST.exists(),'requires prepared real rotation investigation')
class RotationManifestTests(unittest.TestCase):
    def setUp(self):self.m=read_json(MANIFEST)
    def test_real_manifest_accepted(self):validate(self.m)
    def test_unapproved_rotation_rejected(self):
        self.m['tasks'][0]['rotation_matrix'][0][0]+=0.01
        with self.assertRaisesRegex(InvalidArtifact,'transformation'):validate(self.m)
    def test_charge_change_rejected(self):
        self.m['tasks'][0]['charge']+=2
        with self.assertRaisesRegex(InvalidArtifact,'scientific state'):validate(self.m)
    def test_changed_dipole_model_rejected(self):
        self.m['model']['dtype']='float32'
        with self.assertRaisesRegex(InvalidArtifact,'physical model'):validate(self.m)
    def test_report_uses_actual_rotation_denominator(self):
        from mace_hybrid import report
        collection=MANIFEST.parent/'collection_job_1200396.json'
        if not collection.exists():self.skipTest('requires actual completed rotation collection')
        with tempfile.TemporaryDirectory() as folder:
            output=Path(folder)/'report.md';report(collection,output)
            text=output.read_text()
            self.assertIn('Successful core calls: 4/4.',text)
            self.assertIn('New DFT calls: 0.',text)
            self.assertNotIn('4/12',text)

if __name__=='__main__':unittest.main()
