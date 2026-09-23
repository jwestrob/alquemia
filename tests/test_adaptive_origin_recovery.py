"""Real inherited-failure origins; no artificial energies or protein fixtures."""
import json
from pathlib import Path
import sys
import tempfile
import unittest
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
import adaptive_origin_recovery as a
MANIFEST=ROOT/'workspaces/adaptive_origin_recovery_20260923/run_v2/manifest.json'

class Recovery(unittest.TestCase):
    def test_actual_four_origins_without_old_proposal_gate(self):
        r=a.validate(MANIFEST);self.assertEqual((r['prepared'],r['optimizer_starts']),(2,4))
        m=a.read_json(MANIFEST)
        for c in m['cases']:
            self.assertEqual(c['old_row']['pool']['status'],'unavailable')
            self.assertEqual(c['status'],'prepared')
            self.assertTrue(all(x['status']=='complete' for x in c['origins'].values()))

    def test_same_common_modes_and_unmodified_optimizer(self):
        m=a.read_json(MANIFEST);self.assertEqual(m['settings'],a.completion.SETTINGS)
        for cid in a.CASES:
            ca,la=[next(t for t in m['tasks'] if t['case_id']==cid and t['metal']==z) for z in ('Ca','La')]
            self.assertEqual(ca['active_mode_ids'],la['active_mode_ids']);self.assertEqual(len(ca['active_indices']),4)
            self.assertFalse(any(ca['origin_reuse']['point']['full_q']))
            self.assertFalse(any(la['origin_reuse']['point']['full_q']))

    def test_corrupted_actual_manifest_settings_rejected(self):
        m=a.read_json(MANIFEST);m['settings']['maxiter']=201
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'corrupted_real_manifest.json';p.write_text(json.dumps(m))
            with self.assertRaises(a.InvalidArtifact):a.validate(p)

    def test_actual_reference_is_frozen_canonical_only(self):
        m=a.read_json(MANIFEST);ref=a.check_reference(a.verify(m['reference']),m)
        self.assertFalse(ref['noncanonical_folds_used_for_calibration'])
        self.assertEqual(ref['candidate_ids'],list(a.CANDIDATES))
        self.assertEqual(len(ref['variants']['operational']['rows']),25)

    def test_actual_final_pool_when_executed(self):
        path=MANIFEST.parent/'pool/final_collection.json'
        if not path.exists():self.skipTest('new adaptive recovery not executed yet')
        m=a.read_json(MANIFEST.parent/'pool/manifest.json');a.validate_pool(MANIFEST.parent/'pool/manifest.json')
        result=a.read_json(path);self.assertEqual(result['denominator'],2)
        self.assertLessEqual(result['required_new_MACE_cells'],4);self.assertLessEqual(result['required_new_GFN2_calls'],16)
        for c in result['cases']:
            if c['pool']['status']=='available':
                self.assertEqual([q['id'] for q in c['candidates']],list(a.CANDIDATES))
                self.assertEqual(a.pool.choose_rows(c['matrix'],a.CANDIDATES),c['pool'])
            else:self.assertIsNone(c['pool']['operational'])

if __name__=='__main__':unittest.main()
