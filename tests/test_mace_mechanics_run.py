"""Real finite mechanics manifests, cache separation and partial recovery."""
import copy
import json
from pathlib import Path
import sys
import tempfile
import unittest
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from affordable_common import InvalidArtifact,read_json
from mace_mechanics_run import validate,collect_mechanics
CORE=ROOT/'workspaces/mace_mechanics_20260916/core_v1/manifest.json'
SHORT=ROOT/'workspaces/mace_mechanics_20260916/short_v1/manifest.json'


@unittest.skipUnless(CORE.exists() and SHORT.exists(),'requires real mechanics task manifests')
class MechanicsRunnerTests(unittest.TestCase):
    def test_finite_counts_and_exact_cached_GGR_only(self):
        self.assertEqual(validate(CORE)['tasks'],116)
        self.assertEqual(validate(SHORT)['tasks'],106)
        with tempfile.TemporaryDirectory() as folder:
            p=Path(folder)/'manifest.json';p.write_bytes(CORE.read_bytes());c=collect_mechanics(p)
            self.assertEqual(c['status'],'incomplete');self.assertEqual(len(c['reused_rows']),20)
            self.assertTrue(all(r['status']=='unavailable' for r in c['rows'].values()))
            self.assertTrue(all(r['result']['status']=='computed' for r in c['reused_rows'].values()))
            self.assertIsNone(c['relaxation_correction_kcal_mol'])

    def test_changed_model_or_point_is_rejected(self):
        for name in ('model','point'):
            m=copy.deepcopy(read_json(CORE))
            if name=='model':m['model']['solvent']='water'
            else:m['tasks'][0]['point']='explicitly_corrupted_unlisted_point'
            with tempfile.TemporaryDirectory() as folder:
                p=Path(folder)/'manifest.json';p.write_text(json.dumps(m))
                with self.assertRaises(InvalidArtifact):validate(p)


if __name__=='__main__':unittest.main()
