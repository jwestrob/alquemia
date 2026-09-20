"""Source/state checks on real native GFN2 data; no scientific output fabricated."""
import json
from pathlib import Path
import sys
import tempfile
import unittest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from affordable_common import InvalidArtifact, read_json, verify
from compact_solvation_qualification import prepare, validate, recipe, restart_diagnostics


class QualificationSources(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp=tempfile.TemporaryDirectory()
        cls.folder=Path(cls.tmp.name)/'prepared'
        cls.result=prepare(ROOT/'workspaces/compact_solvation_20260920/full_v1/collection_1203165.json',
                           ROOT/'diagnostics/compact_qualification_20260920/PLAN.md',cls.folder)
        cls.path=cls.folder/'manifest.json';cls.manifest=read_json(cls.path)

    @classmethod
    def tearDownClass(cls):cls.tmp.cleanup()

    def test_real_source_coverage_and_seed_identity(self):
        m=self.manifest;self.assertEqual(self.result['tasks'],32)
        self.assertEqual(sum(t['restart'] is not None for t in m['tasks']),8)
        self.assertEqual(sum(t['variant']=='native_tight_fresh' for t in m['tasks']),24)
        for t in m['tasks']:
            self.assertEqual(verify(t['xyz']).read_bytes(),verify(t['primary_task']['xyz']).read_bytes())
            self.assertEqual(t['charge'],t['primary_task']['charge'])
            self.assertEqual(t['multiplicity'],1)
            if t['restart']:
                self.assertEqual(verify(t['restart']['source']).read_bytes(),verify(t['restart']['immutable_copy']).read_bytes())
                self.assertFalse(Path(t['restart']['runtime_path']).exists())

    def test_identified_corruption_cannot_change_electronic_state(self):
        m=read_json(self.path);m['tasks'][0]['charge']+=1
        p=self.folder/'corrupted_manifest.json';p.write_text(json.dumps(m))
        with self.assertRaises(InvalidArtifact):validate(p)

    def test_scf_only_changes_no_chemical_or_temperature_change(self):
        for t in self.manifest['tasks']:
            s=verify(t['input']).read_text()
            self.assertIn('Convergence Tight',s);self.assertIn('SmearTemp 300',s)
            self.assertEqual('ALPB(Water)' in s,t['medium']=='alpb')
            self.assertEqual('UseXTBMixer true' in s,t['variant']=='native_tight_fresh')
            self.assertNotIn('Opt',s);self.assertNotIn('NumGrad',s)
        with self.assertRaises(InvalidArtifact):recipe(-2,'vacuum','undeclared_solver')

    def test_orbital_restart_fix_uses_actual_gbw_and_only_eight_tasks(self):
        folder=Path(self.tmp.name)/'orbital_fix'
        result=prepare(ROOT/'workspaces/compact_solvation_20260920/full_v1/collection_1203165.json',
                       ROOT/'diagnostics/compact_qualification_20260920/RESTART_FIX.md',folder,
                       stage='orbital_restart_fix')
        self.assertEqual(result['tasks'],8)
        for t in read_json(folder/'manifest.json')['tasks']:
            self.assertEqual(t['restart']['kind'],'gbw')
            self.assertIsNone(t['restart']['runtime_path'])
            self.assertIn('%moinp "seed_source.gbw"',verify(t['input']).read_text())
            self.assertGreater(verify(t['restart']['immutable_copy']).stat().st_size,100000)
            self.assertFalse((Path(t['input']['path']).parent/'endpoint.runtime.xtbw').exists())

    def test_actual_ignored_keyword_and_actual_explicit_restart_are_distinguished(self):
        old=(ROOT/'workspaces/compact_qualification_20260920/gbw_v1/tasks/1H4I__Ca__vacuum__ordinary_tight_gbw/endpoint.out').read_text()
        new=(ROOT/'workspaces/compact_qualification_20260920/explicit_v1/tasks/1H4I__Ca__alpb__ordinary_explicit_gbw/endpoint.out').read_text()
        self.assertEqual(restart_diagnostics(old,'ordinary_tight_gbw')['restart_status'],'requested_restart_not_confirmed')
        self.assertTrue(restart_diagnostics(old,'ordinary_tight_gbw')['ignored_moinp_warning'])
        self.assertEqual(restart_diagnostics(new,'ordinary_explicit_gbw')['restart_status'],'orbital_restart_confirmed')
        explicit=recipe(-2,'vacuum','ordinary_explicit_gbw')
        self.assertIn('\n Guess MORead\n MOInp "seed_source.gbw"\n',explicit)
        self.assertNotIn('NoAutostart',explicit)


if __name__=='__main__':unittest.main()
