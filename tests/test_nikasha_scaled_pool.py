"""Real prior pools exercise the full-calibration and missing-fold adapters."""
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
from affordable_common import InvalidArtifact, read_json, verify
from nikasha_pool import expansion_bases, choose_rows


class ScaledPoolSources(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.paths = [ROOT / 'workspaces/nikasha_shared_pool_20260922/pilot_v2/after_solvent_0_1209845.json',
                     ROOT / 'workspaces/nikasha_shared_pool_20260922/remaining26_v1/final_collection.json']
        cls.fold = ROOT / 'workspaces/nikasha_shared_pool_20260922/primary225_v1/final_collection.json'
        if not all(p.exists() for p in [*cls.paths, cls.fold]):
            raise unittest.SkipTest('real completed base-pool fixtures unavailable')

    def test_original30_union_preserves_every_actual_matrix_and_source_order(self):
        data, manifest, origins = expansion_bases(self.paths)
        source = read_json(verify(manifest['source_manifest']))
        self.assertEqual(manifest['population'], 'original30')
        self.assertEqual([c['case_id'] for c in data['cases']], [c['case_id'] for c in source['cases']])
        self.assertEqual(len(origins), 30)
        for c in data['cases']:
            archived = read_json(verify(origins[c['case_id']]))
            self.assertEqual(c, next(r for r in archived['cases'] if r['case_id'] == c['case_id']))
            self.assertEqual(choose_rows(c['matrix'], [q['id'] for q in c['candidates']]), c['pool'])

    def test_duplicate_collection_cannot_be_counted_as_original30(self):
        with self.assertRaisesRegex(InvalidArtifact, 'duplicate source'):
            expansion_bases([self.paths[0], self.paths[0]])

    def test_full225_base_retains_twenty_actual_unavailable_sources(self):
        data, manifest, origins = expansion_bases(self.fold)
        self.assertEqual(data, read_json(self.fold))
        self.assertEqual(len(data['cases']), 225)
        self.assertEqual(set(origins), set(manifest['declared_case_ids']))
        unavailable = [c for c in data['cases'] if c['pool']['status'] == 'unavailable']
        self.assertEqual(len(unavailable), 20)
        self.assertTrue(all(c['pool']['mathematical'] is None and c['pool']['operational'] is None
                            for c in unavailable))


if __name__ == '__main__':
    unittest.main()
