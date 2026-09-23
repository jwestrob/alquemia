"""Real archived source/state/seed and completed two-pass checks."""
import sys,tempfile,unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
import precision_geometry_continuation as p
from affordable_common import read_json,verify,HA_TO_KCAL,InvalidArtifact,write_new

D=ROOT/'workspaces/precision_geometry_continuation_20260923'
M=D/'run_v2/stage1/manifest.json'

class ContinuationTests(unittest.TestCase):
    def test_actual_four_source_states_and_seed_pairs(self):
        rows=p.source_rows(D/'SEEDS.json')
        self.assertEqual(len(rows),4)
        for r in rows:
            receipt=read_json(verify(r['receipt']))
            self.assertEqual(receipt['parallelism']['nprocs'],1)
            self.assertEqual((r['charge'],r['multiplicity'],r['metal'],r['medium']),(-1,1,'La','vacuum'))
            self.assertGreater(verify(r['gbw']).stat().st_size,2_000_000)
            self.assertGreater(verify(r['xtbw']).stat().st_size,12000)
            self.assertEqual(verify(r['input']).read_text().replace(' NoAutostart',''),p.recipe(-1,1,'vacuum'))

    def test_actual_manifest_keeps_self_seed_and_exact_xyz(self):
        result=p.validate(M);self.assertEqual(result['tasks'],4)
        m=read_json(M);self.assertEqual(m['execution_resources'],{'mpi_ranks':1,'concurrent_tasks':4})
        for t in m['tasks']:
            self.assertEqual(verify(t['xyz']).read_bytes(),verify(t['source']['xyz']).read_bytes())
            for field,source in [('immutable_seed','xtbw'),('immutable_gbw','gbw')]:
                self.assertEqual(verify(t[field]).read_bytes(),verify(t['source'][source]).read_bytes())
            self.assertFalse(t['gradient_requested'])
            self.assertNotIn('NoAutostart',verify(t['input']).read_text())
            self.assertIn('MaxIter 500',verify(t['input']).read_text())

    def test_cross_geometry_seed_mapping_is_rejected(self):
        # Explicitly corrupted metadata over real pinned seeds, no fake energy.
        m=read_json(M);m['tasks'][0]['seed_source']=m['tasks'][1]['seed_source']
        import json
        with tempfile.NamedTemporaryFile(mode='w',suffix='.json',dir=M.parent) as f:
            json.dump(m,f);f.flush()
            with self.assertRaisesRegex(InvalidArtifact,'source or exact self-seed changed'):p.validate(f.name)

    def test_completed_passes_and_uniform_stage2_algebra(self):
        path=D/'run_v2/COMPARISON.json'
        if not path.exists():self.skipTest('scientific two-pass outputs not yet executed')
        c=read_json(path);self.assertEqual(c['denominator'],4);self.assertEqual(c['reported_stage'],2)
        self.assertEqual(c['complete'],4)
        for r in c['rows']:
            for actual in r['actual']:
                self.assertEqual(actual['status'],'confirmed_restart')
                text=verify(actual['actual']['output']).read_text()
                self.assertIn('INITIAL GUESS: XTBRESTART',text)
                self.assertTrue(actual['details']['native_mixer_observed'])
            delta=(r['stage2_energy_hartree']-r['stage1_energy_hartree'])*HA_TO_KCAL
            self.assertAlmostEqual(r['stage2_minus_stage1_kcal_mol'],delta,places=10)
            self.assertEqual(r['settled'],abs(delta)<=.1)
        for t in c['targets']:
            rows={r['geometry']:r for r in c['rows'] if r['case_id']==t['case_id']}
            for stage in ('cold','stage1','stage2'):
                delta=(rows['new'][stage+'_energy_hartree']-rows['old'][stage+'_energy_hartree'])*HA_TO_KCAL
                self.assertAlmostEqual(t['geometry_new_minus_old_kcal_mol'][stage],delta,places=10)
        self.assertIsNone(c['new_pool_score'])

    def test_explicitly_corrupted_failed_stage_keeps_four_missing(self):
        path=D/'run_v2/stage1/collection.json'
        if not path.exists():self.skipTest('real stage1 collection unavailable')
        # Malformed-input fixture only: mark a copy of the actual results
        # unavailable; never manufacture a successful scientific calculation.
        c=read_json(path)
        for r in c['rows']:r.update(status='unavailable',energy_hartree=None,reason='explicitly_corrupted_test_copy')
        with tempfile.TemporaryDirectory() as d:
            d=Path(d);prev=d/'corrupted_stage1.json';write_new(prev,c)
            out=d/'stage2'
            p.prepare(D/'SEEDS.json',ROOT/'diagnostics/precision_geometry_continuation_20260923/PLAN.md',out,prev)
            manifest=out/'manifest.json';m=read_json(manifest)
            self.assertEqual(m['tasks'],[]);self.assertEqual(len(m['missing']),4)
            result=out/'collection.json';p.collect(manifest,result)
            collected=read_json(result)
            self.assertEqual(collected['denominator'],4);self.assertEqual(collected['confirmed'],0)
            self.assertTrue(all(r['energy_hartree'] is None for r in collected['rows']))

if __name__=='__main__':unittest.main()
