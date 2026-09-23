"""Real225source preparation, including the older native-solvent failure."""
from pathlib import Path
import sys
import tempfile
import unittest
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
import standalone_xtb_transfer as s
M=ROOT/'workspaces/standalone_xtb_transfer_20260923/run_v1/manifest.json'
class Transfer(unittest.TestCase):
 def test_actual_scope_and_new_backend_does_not_inherit_native_failure(self):
  v=s.validate(M);self.assertEqual((v['cases'],v['physical_sources'],v['logical_cells'],v['actual_reuses'],v['new_calls']),(225,208,832,4,828))
  m=s.read_json(M);a=next(c for c in m['cases'] if c['case_id']=='a8r3s4-pqq-la_model__conditioned_Ca__seed-1_sample-3');self.assertEqual(a['status'],'prepared')
  self.assertIsNone(a['prior_row']['methods']['released']['R'])
 def test_static_only_with_17_missing_and_exact_reuse(self):
  m=s.read_json(M);self.assertEqual(sum(c['status']=='unavailable' for c in m['cases']),17)
  for t in m['cells']:
   self.assertEqual((t['candidate'],t['accuracy']),('origin','0.02'))
   if t['reuse']:self.assertEqual(s.xyz(s.verify(t['xyz'])),s.xyz(s.verify(t['reuse']['source_task']['xyz'])))
 def test_reference_is_frozen_static_not_failed_adaptive(self):
  r=s.reference_check(s.verify(s.read_json(M)['reference']));self.assertEqual(r['variants']['static']['status'],'available');self.assertIsNone(r['variants']['minimal']['bands'])
  r['variants']['static']['bands']['Ca_max']+=.1
  with tempfile.TemporaryDirectory() as d:
   p=Path(d)/'corrupted_real_reference.json';s.write_new(p,r)
   with self.assertRaises(s.InvalidArtifact):s.reference_check(p)
 def test_actual_final_denominators_when_executed(self):
  p=M.parent/'COMPARISON_v1.json'
  if not p.exists():self.skipTest('828 new molecular calls/transfer comparison not yet executed')
  d=s.read_json(p);self.assertEqual((len(d['rows']),len(d['pools']),len(d['triples'])),(225,75,100));self.assertFalse(d['calibration_changed'])
  self.assertEqual(d['counts']['all225'],s.summarize_rows(d['rows'],d['bands']))
  ix={r['case_id']:r for r in d['rows']}
  for g in d['triples']:
   if any(ix[c]['methods']['standalone_static']['R'] is None for c in g['members']):self.assertIsNone(g['methods']['standalone_static']['R'])
if __name__=='__main__':unittest.main()
