"""Actual completed envelope searches and finite paired common-pool tasks."""
from pathlib import Path
import copy
import sys
import tempfile
import unittest
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from affordable_common import InvalidArtifact,read_json,verify,write_new,xyz
import motion_envelope_pool as pool
import motion_envelope_scalar as scalar
BASE=ROOT/'workspaces/motion_envelope_20260923/pilot34_pool_v1'
class EnvelopePool(unittest.TestCase):
    def test_actual_finite_matrix_and_strict_scalar_scope(self):
        v=pool.validate_pool(BASE/'manifest.json');m=read_json(BASE/'manifest.json')
        self.assertEqual(v['denominator'],34);self.assertLessEqual(v['new_MACE_cells'],68);self.assertLessEqual(v['new_GFN2_calls'],272)
        self.assertEqual(m['reference'],None);self.assertEqual(m['numerical_policy_id'],scalar.PROFILE)
        for c in m['cases']:
            if c['status']!='prepared':continue
            self.assertEqual(set(c['matrix']['Ca']),set(c['matrix']['La']))
            self.assertEqual(set(c['aliases']),{'adaptive_Ca','adaptive_La'})
        for t in m['tasks']:
            self.assertEqual(xyz(verify(t['xyz']))[0][0],t['metal'])
        lm=read_json(BASE/'solvent/shard_0/manifest.json');scalar.validate(BASE/'solvent/shard_0/manifest.json')
        self.assertEqual(len(lm['tasks']),2*len(m['tasks']));self.assertEqual(lm['reused'],{})
        self.assertTrue(all('NoAutostart'in verify(t['input']).read_text()for t in lm['tasks']))
    def test_real_shared_selection_requires_every_candidate_cell(self):
        actual=read_json(BASE/'manifest.json')['cases'][0]['matrix']
        self.assertEqual(pool.pool.choose_rows(actual,list(actual['Ca']))['status'],'unavailable')
        # The pending actual cells cannot be substituted with the already complete origin.
        origin_only=pool.pool.choose_rows(actual,['origin'])
        self.assertEqual(origin_only['status'],'available')
if __name__=='__main__':unittest.main()
