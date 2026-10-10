"""Administrative mocks only: no scientific calculations or claimed Slurm runs."""
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

spec=importlib.util.spec_from_file_location('monitor',Path(__file__).with_name('monitor.py'))
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)


def reply(stdout='',stderr='',code=0):
    return dict(stdout=stdout,stderr=stderr,returncode=code)


class MonitorTests(unittest.TestCase):
    def test_ignore_steps_and_other_jobs(self):
        q=reply('123.batch|COMPLETED|0:0|1|1|n\n124|COMPLETED|0:0|1|1|n\n123|RUNNING|0:0|1|1|n\n')
        self.assertEqual(m.record_for('123',q)['state'],'RUNNING')
        self.assertIsNone(m.record_for('123',reply(code=1)))

    def test_readiness_requires_all_three_conditions(self):
        for status,code,marker in [('RUNNING','0:0',True),('FAILED','7:0',True),('COMPLETED','0:0',False),('COMPLETED','1:0',True)]:
            self.assertFalse(m.successful(dict(state=status,exit_code=code),marker))
        self.assertTrue(m.successful(dict(state='COMPLETED',exit_code='0:0'),True))

    def test_submission_uncertainty_is_not_retryable(self):
        self.assertEqual(m.submission_outcome(reply('123;biotite')),('watching','123'))
        self.assertEqual(m.submission_outcome(reply(stderr='Batch job submission failed: Invalid account',code=1)),('waiting',None))
        self.assertEqual(m.submission_outcome(reply(stderr='RPC timed out',code=None)),('submission_uncertain',None))
        self.assertEqual(m.submission_outcome(reply()),('submission_uncertain',None))

    def test_queue_intent_persisted_before_call(self):
        with tempfile.TemporaryDirectory() as d:
            path=Path(d)/'state.json'; state=dict(codex='/codex',thread='thread')
            def queue(cmd):
                self.assertEqual(json.loads(path.read_text())['status'],'enqueue_pending')
                return reply('accepted')
            with patch.object(m,'run',side_effect=queue):m.notify(state,path,'message')
            self.assertEqual(json.loads(path.read_text())['status'],'queued')

    def test_rejection_then_one_pending_probe_then_success(self):
        with tempfile.TemporaryDirectory(prefix='readiness_') as d:
            base=Path(d);instructions=base/'ON_READY.md';instructions.write_text('test')
            statepath=base/'state.json'; clock=[100.0];counts={'sbatch':0,'sacct':0,'queue':0};submitted=[]
            def fake_run(cmd):
                if cmd[0]=='sbatch':
                    counts['sbatch']+=1;submitted.append(clock[0])
                    if counts['sbatch']==1:return reply(stderr='Batch job submission failed: Invalid account',code=1)
                    return reply('123')
                if cmd[0]=='sacct':
                    counts['sacct']+=1
                    if counts['sacct']==1:return reply('123|PENDING|0:0|0|1|None\n')
                    state=json.loads(statepath.read_text())
                    Path(state['output_pattern'].replace('%j','123')).write_text(state['marker']+'\n')
                    return reply('123|COMPLETED|0:0|1|1|testnode\n')
                counts['queue']+=1
                self.assertIn('AUTOMATED SLURM READINESS EVENT',cmd[-1])
                return reply('accepted')
            args=['monitor','--directory',d,'--thread','test-thread','--codex','/bin/true','--instructions',str(instructions)]
            with patch('sys.argv',args),patch.object(m,'run',side_effect=fake_run),patch.object(m.time,'time',side_effect=lambda:clock[0]),patch.object(m.time,'sleep',side_effect=lambda n:clock.__setitem__(0,clock[0]+n)):
                m.main()
                with self.assertRaises(SystemExit):m.main()
            self.assertEqual(counts,dict(sbatch=2,sacct=2,queue=1))
            self.assertGreaterEqual(submitted[1]-submitted[0],1800)
            self.assertEqual(json.loads(statepath.read_text())['status'],'queued')


if __name__=='__main__':unittest.main()
