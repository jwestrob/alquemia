"""Administrative replay from actual completion receipts; no invented chemistry."""
import datetime as dt
import json
import sys
import unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from slurm_native_health_watch import inspect_text,select_events,snapshot
W=ROOT/'workspaces/metal_environment_response_20260926/dy_frozen_f_direction_v1'

def timestamp(s):return dt.datetime.fromisoformat(s).timestamp()

def receipt_snapshot(clock):
    # Replay only terminal-success facts backed by actual receipt timestamps.
    # Live pre-completion SCF health is unknown here and intentionally not modeled.
    rows={}
    for p in W.glob('*/endpoint.out.execution.json'):
        r=json.loads(p.read_text());done=timestamp(r['finished_at_utc'])<=clock
        rows[r['task_id']]={'terminal':done,'normal_termination':done and r['normal_termination'],
            'scf_converged':done and r['scf_converged'],'explicit_failure':False,
            'execution_returncode':r['returncode'] if done else None,'alerts':[]}
    return {'endpoints':rows,'terminal_endpoints':sum(r['terminal'] for r in rows.values()),'total_endpoints':len(rows)}

class CompletionDebounce(unittest.TestCase):
    def test_actual_1220300_chronology_coalesces_to_final_only(self):
        health=json.loads((W/'health_receipt.json').read_text())
        moments=sorted({timestamp(e['utc']) for e in health['notifications'].values()})
        state={'notifications':{}};events=[]
        for clock in moments:
            snap=receipt_snapshot(clock)
            scheduler='COMPLETED' if snap['terminal_endpoints']==snap['total_endpoints'] else 'RUNNING'
            new=select_events(snap,scheduler,state,clock)
            events.extend(k for k,_ in new)
            state['notifications'].update({k:{'status':'queued'} for k,_ in new})
        self.assertEqual(events,['job_terminal'])
        self.assertNotIn('successful_partial_since_unix',state)
        self.assertEqual(select_events(receipt_snapshot(moments[-1]),'COMPLETED',state,moments[-1]),[])

    def test_long_partial_success_alerts_once_after_120_seconds(self):
        # Retained real four-cell scout had two finishedLa endpoints and two
        # unfinishedDy endpoints; replay elapsed monitor timing against that state.
        old=ROOT/'workspaces/metal_environment_response_20260926/lanm_ef3_hans_scout_v1'
        snap=snapshot(old/'manifest.json',224)
        snap['endpoints']={k:{**v,'alerts':[]} for k,v in snap['endpoints'].items()}
        clock=max(timestamp(v['finished_at_utc']) for v in snap['endpoints'].values() if v.get('finished_at_utc'))
        state={'notifications':{}}
        self.assertEqual(select_events(snap,'RUNNING',state,clock),[])
        self.assertEqual(select_events(snap,'RUNNING',state,clock+119),[])
        events=select_events(snap,'RUNNING',state,clock+120)
        self.assertEqual([k for k,_ in events],['successful_partial_persistent'])
        state['notifications']={k:{'status':'queued'} for k,_ in events}
        self.assertEqual(select_events(snap,'RUNNING',state,clock+240),[])

    def test_actual_failure_bypasses_success_debounce(self):
        p=ROOT/'workspaces/metal_environment_response_20260926/dy_small_guess_v3/PModel/endpoint.out'
        row=inspect_text(p.read_text());row['terminal']=row['explicit_failure']
        self.assertTrue(row['terminal'])
        snap={'endpoints':{'PModel':row},'terminal_endpoints':1,'total_endpoints':1}
        events=select_events(snap,'RUNNING',{'notifications':{}},p.stat().st_mtime)
        self.assertIn('PModel:explicit_native_failure',[k for k,_ in events])
        self.assertIn('PModel:near_root_stagnation_review_only',[k for k,_ in events])
        self.assertIn('PModel:terminal',[k for k,_ in events])

if __name__=='__main__':unittest.main()
