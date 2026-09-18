"""Pinned real-source tests; scientific outputs are never synthesized."""
from pathlib import Path
import sys
import tempfile
import unittest
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from affordable_common import InvalidArtifact,read_json,record,verify,write_new
from mace_frozen_response import frontend,parse,source_arrays,task_key,validate
W=ROOT/'workspaces/mace_omol_20260917'


class FrozenResponse(unittest.TestCase):
    def test_real_response_and_native_frontend_preserve_two_states(self):
        m=validate(W/'frozen_response_v2/manifest.json')
        self.assertEqual(m['new_response_solves'],0)
        sw=read_json(verify(m['software'])); parent=read_json(verify(sw['parent']))
        code=frontend(verify(parent['frontend']).read_text())
        self.assertEqual(code,verify(sw['frontend']).read_text())
        self.assertNotIn('call alquemia_induce0c',code)
        self.assertIn('call ufield0d(vd,vp,sd,sp)',code)
        self.assertLess(code.index('call switch(switch_mode)'),code.index('call ufield0d(vd,vp,sd,sp)'))
        for t in m['tasks']:
            s=read_json(verify(t['static'])); r=read_json(verify(t['response']))
            born,md,mp,fd,fp=source_arrays(s,r)
            self.assertTrue(np.any(md!=mp))
            self.assertTrue(np.any(fd!=fp))
            frozen=np.array(t['frozen_indices'])-1
            self.assertTrue(np.all(md[frozen]==0))
            self.assertTrue(np.all(mp[frozen]==0))

    def test_corrupted_real_dipole_rejected_with_updated_file_hash(self):
        m=read_json(W/'frozen_response_v2/manifest.json')
        with tempfile.TemporaryDirectory() as d:
            d=Path(d); t=m['tasks'][0]
            rows=verify(t['input']).read_text().splitlines()
            fields=rows[1].split(); fields[2]=str(float(fields[2])+0.01); rows[1]=' '.join(fields)
            p=d/'corrupted_real_response.dat'; p.write_text('\n'.join(rows)+'\n')
            t['input']=record(p); t['cache_key']=task_key(t,m)
            p=d/'corrupted_real_manifest.json'; write_new(p,m)
            with self.assertRaises(InvalidArtifact): validate(p)

    def test_actual_native_operator_and_energy_receipts(self):
        m=read_json(W/'frozen_response_v2/manifest.json')
        paths=list((W/'frozen_response_v2').glob('collection_job_*.json'))
        if not paths: self.skipTest('actual native operator integration not yet complete')
        co=read_json(paths[-1])
        for t in m['tasks']:
            r=co['rows'][t['task_id']]
            if r['status']!='computed':
                self.assertIn('failure_reason',r)
                continue
            actual=parse(verify(r['receipt']['log']).read_text(),t)
            self.assertEqual(actual['checks'],r['checks'])
            self.assertEqual(actual['GK_cross_kcal'],r['GK_cross_kcal'])
            self.assertIsNone(actual['new_score'])
        self.assertEqual(co['checks_pass'],co['complete'] and all(c['pass_'] for r in co['rows'].values() for c in r.get('checks',[])))


if __name__=='__main__': unittest.main()
