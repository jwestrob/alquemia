"""Real pinned archives/preparations only; no invented scientific outputs."""
import json
from pathlib import Path
import sys
import tempfile
import unittest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
import accommodation_fold_dft as dft

class FrozenFoldDFT(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.release=json.loads((ROOT/'diagnostics/pqq_pmdh_fixed_core_calibration_20260914/result.json').read_text())
        cls.run_dir=ROOT/'workspaces/accommodation_fold_DFT_20260920/run_v1'
        if not (cls.run_dir/'manifest.json').exists():raise unittest.SkipTest('real prepared fold manifest is unavailable')
        cls.manifest=json.loads((cls.run_dir/'manifest.json').read_text())

    def test_archived_sign_units_and_published_bands(self):
        for row in self.release['scores']:
            value=(row['energies_hartree']['Ca']-row['energies_hartree']['La'])*dft.HA_TO_KCAL
            self.assertEqual(value,row['R_kcal_mol'])
            self.assertEqual(dft.decision(value,self.manifest['bands']),row['class']+'-supported')
        self.assertEqual(dft.decision(None,self.manifest['bands']),'unavailable')

    def test_finite_nonoverlapping_preserved_inputs(self):
        result=dft.dry_run(self.run_dir/'manifest.json')
        self.assertEqual(result['new_endpoints'],416)
        self.assertEqual(result['endpoints_per_shard'],104)
        self.assertEqual(len(self.manifest['reused_canonical_pairs']),25)
        for t in self.manifest['new_tasks']:
            self.assertEqual(t['input']['sha256'],t['source_endpoint']['input']['sha256'])
            self.assertEqual(t['xyz']['sha256'],t['source_endpoint']['xyz']['sha256'])

    def test_real_partial_collection_does_not_invent_new_scores(self):
        p=self.run_dir/'prelaunch_collection.json'
        if not p.exists():self.skipTest('real prelaunch collection unavailable')
        c=json.loads(p.read_text())
        self.assertEqual(c['new_endpoint_statuses'],{'pending':416})
        self.assertEqual(c['counts']['single_sources']['canonical25']['correct'],25)
        self.assertEqual(c['counts']['single_sources']['primary225']['unavailable'],225)
        self.assertTrue(all(r['R_kcal_mol'] is None for r in c['rows'] if r['primary_evaluation_pool']))

    def test_corrupted_real_input_cannot_change_method_silently(self):
        t=self.manifest['new_tasks'][0]
        text=Path(t['input']['path']).read_text()
        with tempfile.TemporaryDirectory() as tmp:
            bad=Path(tmp)/'explicitly_corrupted_real_input.inp'
            bad.write_text(text.replace('DefGrid3','DefGrid2'))
            with self.assertRaises(dft.InvalidArtifact):dft.recipe(bad,t['xyz']['path'])

if __name__=='__main__':unittest.main()
