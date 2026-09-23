"""Pinned PQQ geometry/control tests; no fabricated successful solver output."""
from pathlib import Path
import sys
import tempfile
import unittest
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
import standalone_xtb as x
M=ROOT/'workspaces/standalone_xtb_20260923/run_v1/manifest.json'
class Standalone(unittest.TestCase):
 def test_actual_scope_and_fresh_inputs(self):
  v=x.validate(M);self.assertEqual((v['tasks'],v['gradient_calls_within_total']),(224,16))
 def test_accuracy_only_pair_and_no_native_seed(self):
  m=x.read_json(M);a={t['cell_id']:t for t in m['tasks'] if t['accuracy']=='0.2'};b={t['cell_id']:t for t in m['tasks'] if t['accuracy']=='0.02'}
  self.assertEqual(set(a),set(b))
  for k in a:
   self.assertEqual(x.verify(a[k]['xyz']).read_bytes(),x.verify(b[k]['xyz']).read_bytes())
   ca=x.command(a[k],m['executable']['path']);cb=x.command(b[k],m['executable']['path']);ca[ca.index('--acc')+1]='0.02';self.assertEqual(ca,cb);self.assertIn('--norestart',ca)
 def test_source_physical_jacobian_on_real_displacements(self):
  fd=x.read_json(x.verify(x.read_json(M)['force_design']))
  for key,meta in fd['maps'].items():
   kin=x.Kinematics(x.read_json(x.verify(meta['physical_mapping']))['context']);q=np.zeros(len(kin.modes));i=meta['mode_index'];j=kin.evaluate(q)[3][i]
   q[i]=1e-6;a=kin.evaluate(q)[1];q[i]=-1e-6;b=kin.evaluate(q)[1]
   self.assertLess(np.max(np.abs((a-b)/2e-6-j)),1e-7)
 def test_corrupted_real_state_is_rejected(self):
  m=x.read_json(M);m['tasks'][0]['charge']+=1
  with tempfile.TemporaryDirectory() as td:
   p=Path(td)/'corrupted_real_manifest.json';x.write_new(p,m)
   with self.assertRaises(x.InvalidArtifact):x.validate(p)
 def test_actual_outputs_when_executed(self):
  m=x.read_json(M)
  if not all((Path(t['directory'])/'execution.json').exists() for t in m['tasks']):self.skipTest('224 molecular calls not yet executed')
  rows=[x.parse_task(M,t) for t in m['tasks']];self.assertEqual(len(rows),224)
  self.assertEqual(sum(r['status']=='complete' for r in rows),224)
  for t,r in zip(m['tasks'],rows):
   if r['status']=='complete':
    self.assertIsNotNone(r['energy_hartree'])
    if t['gradient_requested']:self.assertIsNotNone(r['gradient'])
   else:self.assertIsNotNone(r['reason'])
 def test_actual_reported_tight_pool_and_derivative_gates(self):
  cp=M.parent/'COMPARISON_v1.json'
  if not cp.exists():self.skipTest('scientific comparison not executed yet')
  from nikasha_pool import choose_rows
  d=x.read_json(cp);self.assertEqual(d['reported_accuracy'],'0.02');self.assertIsNone(d['own_reference'])
  for c in d['pool_cases']:
   for level in c['levels'].values():self.assertEqual(choose_rows(level['matrix'],x.CANDIDATES),level['pool'])
  self.assertEqual(sum(c['energy_qualification_pass'] for c in d['pool_cases']),4)
  self.assertEqual(sum(y['pass'] for c in d['derivative_cases'] if c['accuracy']=='0.02' for y in c['checks']),14)
if __name__=='__main__':unittest.main()
