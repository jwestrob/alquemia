"""Real canonical/crystal preparation and exact actual standalone reuse."""
from pathlib import Path
import sys
import tempfile
import unittest
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
import standalone_xtb_reference as s
M=ROOT/'workspaces/standalone_xtb_reference_20260923/run_v1/manifest.json'
class Reference(unittest.TestCase):
 def test_actual_28_336_36_300_scope(self):
  v=s.validate(M);self.assertEqual((v['sources'],v['logical_cells'],v['actual_reuses'],v['new_calls']),(28,336,36,300))
 def test_canonical_only_roles_and_identical_backend(self):
  m=s.read_json(M);self.assertEqual(sum(c['role']=='calibration' for c in m['cases']),25)
  self.assertEqual({c['case_id'] for c in m['cases'] if c['role']!='calibration'},set(s.CRYSTALS))
  for t in m['cells']:
   self.assertEqual(t['accuracy'],'0.02');self.assertEqual(t['multiplicity'],1)
   if t['reuse']:
    old=t['reuse']['source_task'];self.assertEqual(s.xyz(s.verify(old['xyz'])),s.xyz(s.verify(t['xyz'])))
    self.assertEqual((old['charge'],old['multiplicity']),(t['charge'],t['multiplicity']))
 def test_malformed_actual_label_cannot_recalibrate(self):
  m=s.read_json(M);m['cases'][0]['expected_class_for_report_only']='corrupted'
  with tempfile.TemporaryDirectory() as d:
   p=Path(d)/'corrupted_real_manifest.json';s.write_new(p,m)
   with self.assertRaises(s.InvalidArtifact):s.validate(p)
 def test_actual_calibration_when_complete(self):
  p=M.parent/'REFERENCE_v1.json'
  if not p.exists():self.skipTest('300 new standalone reference calls/calibration not yet executed')
  r=s.read_json(p);self.assertEqual(len(r['rows']),28);self.assertFalse(r['noncanonical_folds_used_for_calibration'])
  for name,v in r['variants'].items():
   self.assertEqual(len(v['rows']),25);self.assertEqual(s.extrema_reference(v['rows'],name,'standalone_xtb671_canonical25_v1'),v)
if __name__=='__main__':unittest.main()
