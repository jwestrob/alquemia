"""Real archived native seeds and exact prepared eight-task manifest only."""
import copy
import json
from pathlib import Path
import sys
import tempfile
import unittest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
import native_xtb_restart as n

SOURCES=ROOT/'diagnostics/structure_informed_starts_20260922/RESTART_SOURCES.json'
MANIFEST=ROOT/'workspaces/native_xtb_restart_20260922/run_v2/manifest.json'
RECOVERY=ROOT/'workspaces/native_xtb_restart_20260922/recovery_v1/manifest.json'


class NativeRestart(unittest.TestCase):
    def test_real_seed_sources_and_only_guess_policy_change(self):
        rows=n.source_rows(SOURCES)
        self.assertEqual(len(rows),4)
        for row in rows.values():
            self.assertEqual(n.recipe(row['medium'],False).replace(' NoAutostart',''),n.recipe(row['medium']))
            self.assertNotIn(n.SETTINGS['restart_marker'],n.verify(row['output']).read_text())

    def test_corrupted_copy_of_actual_seed_pin_rejected(self):
        source=n.read_json(SOURCES);source['rows'][0]['xtbw']['sha256']='0'*64
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'corrupted_real_sources.json';p.write_text(json.dumps(source))
            with self.assertRaises(n.InvalidArtifact):n.source_rows(p)

    def test_exact_eight_prepared_destinations_and_seeds(self):
        self.assertEqual(n.validate(MANIFEST)['tasks'],8)
        m=n.read_json(MANIFEST)
        for t in m['tasks']:
            same=t['source']['geometry']==t['seed_source']['geometry']
            self.assertEqual(same,t['seed_kind']=='self')
            self.assertEqual(t['source']['medium'],t['seed_source']['medium'])
            self.assertEqual(Path(t['active_seed_path']).name,'endpoint.runtime.xtbw')
            self.assertEqual(n.verify(t['immutable_seed']).read_bytes(),n.verify(t['seed_source']['xtbw']).read_bytes())

    def test_collection_preserves_actual_incomplete_or_complete_status(self):
        with tempfile.TemporaryDirectory() as d:
            result=n.collect(MANIFEST,Path(d)/'collection.json')
            self.assertEqual(result['denominator'],8)
            self.assertIsNone(result['full_Ca_La_score'])
            for r in result['rows']:
                if r['status']!='confirmed_restart':self.assertIsNone(r['energy_hartree'])
                else:
                    self.assertTrue(r['restart_marker_observed'])
                    self.assertTrue(r['details']['native_mixer_observed'])
                    self.assertAlmostEqual(r['change_from_unseeded_kcal_mol'],
                        (r['energy_hartree']-r['source']['energy_hartree'])*n.HA_TO_KCAL)

    def test_actual_restart_outputs_when_available(self):
        collections=list(MANIFEST.parent.glob('collection_*.json'))
        if not collections:self.skipTest('eight scientific restart calls not yet executed')
        c=n.read_json(collections[0]);self.assertEqual(len(c['rows']),8)
        for r in c['rows']:
            self.assertIsNotNone(r['observed_energy_hartree'])
            text=n.verify(r['actual']['output']).read_text()
            self.assertTrue(r['details']['native_mixer_observed'])
            self.assertEqual(r['seed_before']['active_seed']['sha256'],r['seed_source']['xtbw']['sha256'])
            n.verify(r['seed_after']['preserved_after'])
            if r['status']=='confirmed_restart':
                self.assertIn(n.SETTINGS['restart_marker'],text)
            else:
                self.assertEqual(r['status'],'restart_not_confirmed')
                self.assertNotIn(n.SETTINGS['restart_marker'],text)
                self.assertIn('INITIAL GUESS: SAD',text)
                self.assertIsNone(r['energy_hartree'])

    def test_recovery_pairs_actual_same_source_gbw_and_xtbw(self):
        self.assertEqual(n.validate(RECOVERY)['tasks'],8)
        m=n.read_json(RECOVERY)
        self.assertEqual(m['activation'],'matched_gbw_auto')
        for t in m['tasks']:
            self.assertEqual(n.verify(t['immutable_gbw']).read_bytes(),n.verify(t['seed_source']['gbw']).read_bytes())
            self.assertEqual(Path(t['seed_source']['gbw']['path']).parent,Path(t['seed_source']['xtbw']['path']).parent)
            self.assertEqual(n.verify(t['input']).read_text(),n.recipe(t['medium']))

    def test_actual_recovery_outputs_when_available(self):
        collections=list(RECOVERY.parent.glob('collection_*.json'))
        if not collections:self.skipTest('documented activation recovery not yet executed')
        c=n.read_json(collections[0]);self.assertEqual(len(c['rows']),8)
        for r in c['rows']:
            if r['status']=='confirmed_restart':
                self.assertIn(n.SETTINGS['restart_marker'],n.verify(r['actual']['output']).read_text())
                self.assertEqual(r['seed_before']['active_gbw']['sha256'],r['seed_source']['gbw']['sha256'])
                n.verify(r['seed_after']['preserved_after'])
                n.verify(r['seed_after']['preserved_gbw_after'])
            else:self.assertIsNone(r['energy_hartree'])
        for g in n.GEOMETRIES:
            for k in n.SEEDS:
                rows={r['medium']:r['energy_hartree'] for r in c['rows'] if r['geometry']==g and r['seed_kind']==k}
                if all(e is not None for e in rows.values()):
                    self.assertAlmostEqual(c['matched_La_transfer_kcal_mol'][g][k],
                        (rows['alpb']-rows['vacuum'])*n.HA_TO_KCAL)


if __name__=='__main__':unittest.main()
