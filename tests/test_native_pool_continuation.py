"""Actual prepared/archive fixtures; no fabricated scientific energies."""
import copy
import json
from pathlib import Path
import sys
import tempfile
import unittest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
import native_pool_continuation as n

INV=ROOT/'diagnostics/native_xtb_restart_20260922/SEED_AVAILABILITY.json'
RUN=ROOT/'workspaces/native_pool_continuation_20260923/run_v1'
MANIFEST=RUN/'stage1/manifest.json'


class Continuation(unittest.TestCase):
    def test_actual_manifest_scope_and_native_recipe(self):
        result=n.validate(MANIFEST)
        self.assertEqual(result['tasks'],80);self.assertEqual(result['gradient_requests'],0)
        m=n.read_json(MANIFEST)
        for t in m['tasks']:
            self.assertEqual(n.verify(t['xyz']).read_bytes(),n.verify(t['source']['xyz']).read_bytes())
            body=n.verify(t['input']).read_text()
            self.assertNotIn('NoAutostart',body);self.assertNotIn('TightSCF',body)
            self.assertIn('MaxIter 500',body);self.assertIn('SmearTemp 300',body)
            self.assertIn('UseXTBMixer true',body)

    def test_gradient_requests_exact_eight_origins(self):
        rows=n.read_json(INV)['rows']
        selected=[r for r in rows if n.want_gradient(2,r)]
        self.assertEqual(len(selected),8)
        self.assertEqual({r['case_id'] for r in selected},{'4MAE','q88jh5-pqq-la_model'})
        for r in selected:
            self.assertEqual(r['candidate'],'origin')
            plain=n.recipe(r['charge'],r['multiplicity'],r['medium'])
            grad=n.recipe(r['charge'],r['multiplicity'],r['medium'],True)
            self.assertEqual(grad.replace(' EnGrad',''),plain)

    def test_corrupted_real_source_membership_rejected(self):
        inv=n.read_json(INV);inv['rows'].pop()
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'corrupted_inventory.json';p.write_text(json.dumps(inv))
            with self.assertRaises(n.InvalidArtifact):n.source_rows(p)

    def test_corrupted_actual_recipe_rejected(self):
        m=n.read_json(MANIFEST);m['settings']['reported_stage']=1
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'corrupted_manifest.json';p.write_text(json.dumps(m))
            with self.assertRaises(n.InvalidArtifact):n.validate(p)

    def test_actual_collection_partial_status_and_algebra(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'collection.json';n.collect(MANIFEST,p);c=n.read_json(p)
            self.assertEqual(c['denominator'],80);self.assertEqual(len(c['rows']),80)
            for r in c['rows']:
                if r['status']!='confirmed_restart':self.assertIsNone(r['energy_hartree'])
                else:
                    self.assertAlmostEqual(r['change_from_archive_kcal_mol'],
                        (r['energy_hartree']-r['source']['source_energy_hartree'])*n.HA_TO_KCAL)
                    self.assertIn('INITIAL GUESS: XTBRESTART',n.verify(r['actual']['output']).read_text())

    def test_completed_actual_stages_when_available(self):
        result=RUN/'result.json'
        if not result.exists():self.skipTest('two-stage scientific execution not yet complete')
        with tempfile.TemporaryDirectory() as directory:
            replay=Path(directory)/'comparison.json'
            n.compare(RUN/'stage1/collection.json',RUN/'stage2/collection.json',replay)
            d=n.read_json(replay)
        self.assertEqual(d['case_denominator'],4)
        second=n.read_json(n.verify(d['stage2']))
        for r in second['rows']:
            if r['gradient']:
                self.assertEqual(r['candidate'],'origin')
                self.assertFalse(r['gradient']['numerical_gradient_used'])
                n.verify(r['gradient']['engrad'])
        for r in d['cases']:
            for cells in r['matrix_stage2'].values():
                for cell in cells.values():
                    self.assertIn('archived_low',cell)
                    if cell['status']=='complete':
                        self.assertEqual(cell['components']['GFN2_vacuum_hartree'],cell['low']['vacuum']['energy_hartree'])
                        self.assertEqual(cell['components']['GFN2_ALPB_hartree'],cell['low']['alpb']['energy_hartree'])
            if r['status']!='qualified_within_declared_test':self.assertIsNone(r['qualified_R'])
            else:
                self.assertTrue(all(c['pass'] for c in r['cell_settling']))
                self.assertTrue(all(c['pass'] for c in r['same_geometry_contrasts']))
            if r['endpoint_work_from_old_pool']:
                w=r['endpoint_work_from_old_pool']
                self.assertAlmostEqual(w['Ca']['composite_kcal_mol']-w['La']['composite_kcal_mol'],r['stage2_minus_old_R'],places=6)


if __name__=='__main__':unittest.main()
