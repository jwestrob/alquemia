"""File-only recovery/admission checks. Temporary fixtures are never scientific output."""
import contextlib, copy, importlib.util, io, json, tempfile, unittest
from pathlib import Path
from unittest.mock import patch
SPEC=importlib.util.spec_from_file_location('batch',Path(__file__).with_name('CPU_BATCH.py'))
b=importlib.util.module_from_spec(SPEC); SPEC.loader.exec_module(b)
MANIFEST=b.ROOT/'workspaces/orca_parameter_audit_20261012/cpu_proposal_batch_preflight_v1/manifest.json'

class NoCompute:
    def evaluate(self,*args): raise AssertionError('test attempted a molecular call')
    def optimize(self,*args): raise AssertionError('test attempted an optimizer restart')

class Checks(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory(prefix='cpu_batch_test_'); self.root=Path(self.tmp.name)
        self.m=b.read(MANIFEST)
    def tearDown(self): self.tmp.cleanup()
    def test_real_pinned_input_and_wrong_state_rejected(self):
        t=self.m['tasks'][0]
        req={'manifest':b.pin(MANIFEST),'task_id':t['task_id'],'xyz':t['xyz'],'charge':t['charge'],'multiplicity':t['multiplicity']}
        task,rows=b.authorize_request(MANIFEST,self.m,req)
        self.assertEqual(task,t); self.assertEqual(rows[0][0],'Ca')
        with self.assertRaisesRegex(ValueError,'state differs'): b.authorize_request(MANIFEST,self.m,{**req,'charge':t['charge']+1})
        with self.assertRaisesRegex(ValueError,'root request differs'): b.authorize_request(MANIFEST,self.m,{**req,'xyz':self.m['tasks'][1]['xyz']})
    def test_started_incomplete_model_call_not_repeated(self):
        t=self.m['tasks'][0]; req={'manifest':b.pin(MANIFEST),'task_id':t['task_id'],'xyz':t['xyz'],'charge':t['charge'],'multiplicity':t['multiplicity']}
        b.write(self.root/'DISPATCH_STARTED.json',{'test_only':True})
        r,p=b.cell(self.root,req,NoCompute(),self.m['model'])
        self.assertEqual(r['status'],'unavailable'); self.assertIn('not_repeated',r['reason'])
        self.assertEqual(b.cell(self.root,req,NoCompute(),self.m['model'])[1],p)
    def test_actual_archived_force_receipt_reuse_without_new_call(self):
        # Alter only request provenance inside this explicitly temporary test fixture.
        t=self.m['tasks'][-2]; fixture=self.m['fixtures'][0]
        rp=self.root/'mace_request.json'; req={'manifest':b.pin(MANIFEST),'task_id':t['task_id'],'xyz':fixture['xyz'],'charge':t['charge'],'multiplicity':t['multiplicity']}
        b.write(rp,req); native=copy.deepcopy(b.read(b.verify(fixture['native_receipt'])))
        native.update(request=b.pin(rp),model=self.m['model'],status='complete',test_fixture_only=True)
        b.write(self.root/'mace_result.json',native)
        r,p=b.cell(self.root,req,NoCompute(),self.m['model'])
        self.assertEqual(r['energy_eV'],fixture['energy_eV']); b.verify(r['forces'])
    def test_interrupted_optimizer_no_fresh_start(self):
        mp=self.root/'manifest.json'; b.write(mp,{'test_only':True}); t={'task_id':'test_Ca'}
        d=self.root/'proposals'/'test_Ca'; b.write(d/'OPTIMIZER_STARTED.json',{'test_only':True})
        p=b.proposal(mp,t,NoCompute(),NoCompute()); r=b.read(b.verify(p))
        self.assertEqual(r['status'],'unavailable'); self.assertIn('no_q0_restart',r['reason'])
        self.assertEqual(b.proposal(mp,t,NoCompute(),NoCompute()),p)
    def test_all_terminal_rebuilds_missing_and_stale_index_without_models(self):
        mp=self.root/'manifest.json'; b.write(mp,self.m)
        for c in self.m['cases']:
            b.write(self.root/'sources'/c['case_id']/'SOURCE_RESULT.json',{'case_id':c['case_id'],'status':'unavailable','test_fixture_only':True})
        with patch.object(b,'validate',return_value={}), patch.object(b.subprocess,'Popen',side_effect=AssertionError('no process may launch')), contextlib.redirect_stdout(io.StringIO()):
            b.execute(mp)
            first=b.pin(self.root/'attempts/001/FINAL_COLLECTION.json')
            r=b.read(self.root/'FINAL_COLLECTION.json'); self.assertTrue(r['terminal']); self.assertEqual(r['terminal_sources'],106)
            (self.root/'FINAL_COLLECTION.json').write_text('{"stale_test_index":true}\n')
            b.execute(mp)
        b.verify(first); self.assertTrue(b.read(self.root/'FINAL_COLLECTION.json')['terminal'])
        self.assertEqual(b.read(self.root/'attempts/002/RUN_RESULT.json')['new_molecular_calls'],0)

if __name__=='__main__': unittest.main(verbosity=2)
