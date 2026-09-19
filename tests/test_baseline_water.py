"""Promoted routing and score accounting on real archived calculations only."""
from pathlib import Path
import copy
import shutil
import sys
import tempfile
import unittest
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
import baseline_water as b
from affordable_common import InvalidArtifact,HA_TO_KCAL,read_json,record,verify,write_new,cache_key,xyz

class BaselineWaterTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.base=ROOT/'workspaces/water_promotion_20260919'
        cls.wet=cls.base/'dft_v2/manifest.json'
        cls.dry=cls.base/'dry_dft_v3/manifest.json'

    def test_default_and_actual_wet_scores_preserved(self):
        self.assertEqual(b.release(b.DEFAULT_RELEASE)['default_water_policy'],'contextual_if_supported')
        check=b.dry_run(self.wet)
        self.assertEqual((check['new_DFT_endpoints'],check['reused_DFT_endpoints']),(0,8))
        with tempfile.TemporaryDirectory() as td:
            r=b.collect(self.wet,Path(td)/'result.json')
        self.assertEqual(r,read_json(self.base/'result_v2.json'))
        self.assertEqual(r['status'],'complete')
        for c in r['cases']:
            self.assertEqual(c['selected_score'],c['water_prepared'])
            self.assertNotEqual(c['original']['R_kcal_mol'],c['selected_score']['R_kcal_mol'])
            for a in ('original','water_prepared'):
                s=c[a];e=s['endpoints']
                self.assertEqual(s['R_kcal_mol'],(e['Ca']['energy_hartree']-e['La']['energy_hartree'])*HA_TO_KCAL)
                self.assertIsNone(s['S_kcal_mol'])
                self.assertEqual(s['decision'],'uncalibrated_protocol')
        self.assertIsNone(r['occupancy_probabilities']);self.assertIsNone(r['entropy_correction'])

    def test_dry_pqq_keeps_real_released_classifications(self):
        self.assertEqual(b.dry_run(self.dry)['reused_DFT_endpoints'],4)
        with tempfile.TemporaryDirectory() as td:r=b.collect(self.dry,Path(td)/'result.json')
        self.assertEqual(r['status'],'complete')
        calls={c['case_id']:c['selected_score'] for c in r['cases']}
        self.assertAlmostEqual(calls['PQQ_1H4I']['S_kcal_mol'],7.52,delta=.01)
        self.assertAlmostEqual(calls['PQQ_4MAE']['S_kcal_mol'],38.09,delta=.01)
        self.assertNotEqual(calls['PQQ_1H4I']['decision'],calls['PQQ_4MAE']['decision'])
        for c in r['cases']:
            self.assertEqual(c['selected_arm'],'original');self.assertIsNone(c['water_prepared'])
            self.assertEqual(c['selected_score']['protocol_id'],b.CANONICAL)
        for t in read_json(self.dry)['all_tasks']:
            self.assertEqual(xyz(verify(t['xyz'])),xyz(verify(t['original_xyz'])))

    def test_replay_executes_no_new_science(self):
        self.assertEqual(b.water_execute(self.base/'prepared_v2/plan.json')['new_model_calls'],0)
        self.assertEqual(b.execute(self.wet)['new_DFT_endpoints'],0)

    def test_withheld_real_endpoint_remains_missing_and_recovers(self):
        m=read_json(self.wet);t=next(t for t in m['all_tasks'] if t['arm']=='water_prepared')
        m['reused'].pop(t['task_id']);m['tasks']=[t]
        with tempfile.TemporaryDirectory() as td:
            p=Path(td)
            for field,name in [('input','endpoint.inp'),('xyz','core.xyz')]:
                shutil.copyfile(verify(t[field]),p/name);t[field]=record(p/name)
            t['output_path']=str(p/'endpoint.out')
            t['cache_key']=cache_key({'task':{k:v for k,v in t.items() if k!='cache_key'},'release':m['release']})
            write_new(p/'withheld_real_receipt.json',m)
            r=b.collect(p/'withheld_real_receipt.json',p/'result.json')
        c=next(c for c in r['cases'] if c['case_id']==t['case'])
        self.assertEqual(r['status'],'incomplete');self.assertEqual(c['original']['status'],'complete')
        self.assertEqual(c['selected_score']['status'],'unavailable');self.assertIsNone(c['selected_score']['R_kcal_mol'])
        self.assertEqual(b.dry_run(self.wet)['reused_DFT_endpoints'],8)

    def test_corrupted_real_charge_input_cannot_reuse(self):
        m=read_json(self.wet);t=m['all_tasks'][0]
        with tempfile.TemporaryDirectory() as td:
            p=Path(td);ip=p/'explicitly_corrupted_real.inp'
            body=verify(t['input']).read_text().replace(f'xyzfile {t["charge"]} 1',f'xyzfile {t["charge"]+1} 1')
            ip.write_text(body);t['input']=record(ip)
            t['cache_key']=cache_key({'task':{k:v for k,v in t.items() if k!='cache_key'},'release':m['release']})
            write_new(p/'bad.json',m)
            with self.assertRaisesRegex(InvalidArtifact,'input state'):b.dry_run(p/'bad.json')

    def test_changed_route_cannot_hide_original(self):
        m=read_json(self.wet);m['cases'][0]['selected_arm']='original'
        with tempfile.TemporaryDirectory() as td:
            p=Path(td)/'corrupted_route.json';write_new(p,m)
            with self.assertRaisesRegex(InvalidArtifact,'routing'):b.dry_run(p)

    def test_real_receipts_cannot_be_swapped_between_metals(self):
        m=read_json(self.wet);a,c=m['all_tasks'][:2]
        m['reused'][a['task_id']]=m['reused'][c['task_id']]
        with tempfile.TemporaryDirectory() as td:
            p=Path(td)/'swapped_real_receipt.json';write_new(p,m)
            with self.assertRaisesRegex(InvalidArtifact,'quantum reuse'):b.dry_run(p)

    def test_explicit_original_request(self):
        req=read_json(self.base/'request_v1.json');req['water_policy']='original'
        with tempfile.TemporaryDirectory() as td:
            p=Path(td);write_new(p/'request.json',req)
            b.prepare(p/'request.json',p/'prepared')
            plan=b.checked_plan(p/'prepared/plan.json')
            self.assertIsNone(plan['water_manifest'])
            self.assertTrue(all(c['operation']=='original_reference' for c in plan['cases']))

    def test_fresh_retry_reuses_completed_science(self):
        with tempfile.TemporaryDirectory() as td:
            r=b.prepare_retry(self.wet,Path(td)/'retry')
        self.assertEqual((r['new_DFT_endpoints'],r['reused_DFT_endpoints']),(0,8))

    def test_unsupported_wet_core_never_falls_back(self):
        from contextual_water_context import source_parent
        c=read_json(self.wet)['cases'][0];whole=read_json(verify(c['preparation']))
        core=read_json(verify(whole['source_preparation']))
        # Explicitly corrupted copy of real alpha preparation, not scientific data.
        core['protocol_id']='unsupported_cofactor_core'
        with tempfile.TemporaryDirectory() as td:
            p=Path(td)/'corrupted_core.json';write_new(p,core);whole['source_preparation']=record(p)
            with self.assertRaisesRegex(InvalidArtifact,'unsupported'):source_parent(whole)

if __name__=='__main__':unittest.main()
