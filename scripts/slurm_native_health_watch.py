#!/usr/bin/env python3
"""Alert-only SCF health monitor for explicitly owned native Slurm manifests.

Uses supported codex queue; never cancels, submits, or changes jobs. Alerts are
numerical-health observations, not scientific endpoint acceptance decisions.
"""
import argparse
import datetime as dt
import fcntl
import getpass
import hashlib
import json
from pathlib import Path
import re
import subprocess
import time

POLICY = {'id':'native_scf_health_v2', 'trah_macro_window':8,
          'trah_error_floor':1.0, 'minimum_relative_error_improvement':0.5,
          'extreme_micro_step_hartree':1e6, 'extreme_distinct_macro_count':3,
          'silent_output_seconds':900, 'automatic_cancellation':False,
          'diis_window':12, 'diis_error_floor':1e-2, 'diis_rms_density_floor':1e-1,
          'diis_max_density_floor':1.0, 'diis_energy_span_hartree':1.0}
TERMINAL={'COMPLETED','FAILED','CANCELLED','TIMEOUT','NODE_FAIL','OUT_OF_MEMORY','PREEMPTED','BOOT_FAIL','DEADLINE','REVOKED'}
F=r'[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eEdD][-+]?\d+)?'

def now():return dt.datetime.now(dt.timezone.utc).isoformat()
def number(v):return float(v.replace('D','E').replace('d','e'))
def save(path,data):
    p=path.with_suffix(path.suffix+'.tmp');p.write_text(json.dumps(data,indent=2,allow_nan=False)+'\n');p.replace(path)

def inspect_text(text):
    normal='ORCA TERMINATED NORMALLY' in text
    converged=bool(re.search(r'SCF CONVERGED AFTER',text))
    failure=bool(re.search(r'ORCA finished by error termination|SCF NOT CONVERGED|SCF DID NOT CONVERGE',text,re.I))
    macros=[{'iteration':int(i),'energy_hartree':number(e),'error_norm':number(g)} for i,e,g in
            re.findall(r'^\s*(\d+)\s+('+F+r')\s+('+F+r').*\(TRAH MAcro\)',text,re.M)]
    micro=[{'iteration':int(i),'step_hartree':number(e)} for i,e in
           re.findall(r'^\s*(\d+)\s+dE\s+('+F+r').*\(TRAH MIcro\)',text,re.M)]
    diis_matches=list(re.finditer(r'^\s*(\d+)\s+('+F+r')\s+('+F+r')\s+('+F+r')\s+('+F+r')\s+('+F+r')\s+('+F+r')\s+('+F+r')\s*$',text,re.M))
    diis=[dict(zip(('iteration','energy_hartree','delta_energy_hartree','rms_density','max_density','diis_error','damping','iteration_seconds'),[int(m[1])]+[number(v) for v in m.groups()[1:]])) for m in diis_matches]
    phase_transition=max(text.rfind('(TRAH MAcro)'),text.rfind('Leaving SCF to start the TRAH'),text.rfind('Initializing SOSCF'))
    diis_active=bool(diis_matches and diis_matches[-1].end()>phase_transition)
    alerts=[]
    if failure:alerts.append('explicit_native_failure')
    window=macros[-POLICY['trah_macro_window']:]
    if not normal and not converged:
        d=diis[-POLICY['diis_window']:]
        if diis_active and len(d)==POLICY['diis_window'] and all(r['diis_error']>POLICY['diis_error_floor'] and r['rms_density']>POLICY['diis_rms_density_floor'] and r['max_density']>POLICY['diis_max_density_floor'] for r in d) and d[-1]['diis_error']>=d[0]['diis_error']*POLICY['minimum_relative_error_improvement'] and max(r['energy_hartree'] for r in d)-min(r['energy_hartree'] for r in d)>POLICY['diis_energy_span_hartree']:
            alerts.append('gross_diis_nonprogress')
        if len(window)==POLICY['trah_macro_window'] and all(r['error_norm']>POLICY['trah_error_floor'] for r in window) and window[-1]['error_norm']>=window[0]['error_norm']*POLICY['minimum_relative_error_improvement']:
            alerts.append('persistent_large_trah_error_without_twofold_improvement')
        extreme={r['iteration'] for r in micro if abs(r['step_hartree'])>POLICY['extreme_micro_step_hartree']}
        if len(extreme)>=POLICY['extreme_distinct_macro_count']:
            alerts.append('repeated_extreme_trah_micro_steps')
    return {'normal_termination':normal,'scf_converged':converged,'explicit_failure':failure,
            'last_trah_macros':window,'trah_macro_count':len(macros),
            'last_diis_iteration':diis[-1]['iteration'] if diis else None,
            'last_diis_rows':diis[-POLICY['diis_window']:], 'diis_current_phase':diis_active,
            'extreme_micro_distinct_macros':sorted({r['iteration'] for r in micro if abs(r['step_hartree'])>POLICY['extreme_micro_step_hartree']}),
            'negative_gap_warning_count':len(re.findall('negative HOMO - LUMO gap',text)),
            'alerts':alerts}

def snapshot(manifest,allocated_cpus,clock=None):
    clock=time.time() if clock is None else clock
    m=json.loads(Path(manifest).read_text()); rows={}
    for task in m['tasks']:
        out=Path(task['output_path']);receipt=Path(str(out)+'.execution.json')
        row=inspect_text(out.read_text(errors='replace')) if out.exists() else inspect_text('')
        row.update(output=str(out),output_bytes=out.stat().st_size if out.exists() else 0,
                   output_age_seconds=max(0,clock-out.stat().st_mtime) if out.exists() else None,
                   receipt_exists=receipt.exists())
        if receipt.exists():
            r=json.loads(receipt.read_text());row['execution_returncode']=r.get('returncode')
            row['receipt_job_id']=r.get('allocation',{}).get('slurm_job_id')
            row['finished_at_utc']=r.get('finished_at_utc')
            if r.get('returncode') not in (None,0) and 'explicit_native_failure' not in row['alerts']:
                row['alerts'].append('execution_receipt_failure')
        row['terminal']=row['normal_termination'] or row['explicit_failure'] or row['receipt_exists']
        if not row['terminal'] and row['output_age_seconds'] is not None and row['output_age_seconds']>=POLICY['silent_output_seconds']:
            row['alerts'].append('output_silent_review_needed_not_proven_stall')
        rows[task['task_id']]=row
    resources=m['execution_resources'];remaining=sum(not r['terminal'] for r in rows.values())
    occupied=min(remaining,resources['concurrent_tasks'])*resources['mpi_ranks']
    return {'endpoints':rows,'terminal_endpoints':len(rows)-remaining,'total_endpoints':len(rows),
            'maximum_remaining_endpoint_rank_demand':occupied,
            'allocated_cpus':allocated_cpus,
            'allocated_cpus_without_remaining_endpoint_work':max(0,allocated_cpus-occupied),
            'resource_note':'scheduler/manifest accounting; not measured CPU utilization'}

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for x in ('manifest','job','thread','codex','receipt'):p.add_argument('--'+x,required=True)
    p.add_argument('--interval',type=float,default=30)
    p.add_argument('--inspect-only',action='store_true');a=p.parse_args()
    if not a.job.isdecimal() or a.interval<1:p.error('numeric owned job and interval>=1 required')
    manifest=Path(a.manifest).resolve();receipt=Path(a.receipt).resolve();receipt.parent.mkdir(parents=True,exist_ok=True)
    if a.inspect_only:
        print(json.dumps(snapshot(manifest,0),indent=2));return
    identity={'job':a.job,'thread':a.thread,'manifest':str(manifest),'manifest_sha256':hashlib.sha256(manifest.read_bytes()).hexdigest(),'policy':POLICY}
    with receipt.with_suffix('.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        state=json.loads(receipt.read_text()) if receipt.exists() else dict(identity=identity,started_utc=now(),notifications={},status='watching')
        if state['identity']!=identity:raise SystemExit('Watcher receipt identity differs')
        if any(n['status'] != 'queued' for n in state['notifications'].values()):
            raise SystemExit('Prior wake delivery uncertain/failed; inspect queue and receipt before explicit recovery')
        if state['status']=='terminal':return
        while True:
            if hashlib.sha256(manifest.read_bytes()).hexdigest()!=identity['manifest_sha256']:raise SystemExit('Manifest changed during monitoring')
            q=subprocess.run(['sacct','-n','-P','-j',a.job,'--format=JobIDRaw,User,State,AllocCPUS'],capture_output=True,text=True,timeout=30)
            record=next((line.split('|') for line in q.stdout.splitlines() if line.split('|')[0]==a.job),None)
            if q.returncode or not record:
                state.update(last_scheduler_error=q.stderr or 'accounting row not yet available',last_checked_utc=now());save(receipt,state);time.sleep(a.interval);continue
            if record[1]!=getpass.getuser():raise SystemExit('Refusing to monitor a job owned by another user')
            job_state=record[2].split()[0].rstrip('+')
            try:
                snap=snapshot(manifest,int(record[3]))
            except (OSError,ValueError) as error:
                state.update(last_snapshot_error=str(error),last_checked_utc=now());save(receipt,state);time.sleep(a.interval);continue
            events=[]
            for tid,row in snap['endpoints'].items():
                events += [(tid+':'+reason,reason) for reason in row['alerts']]
                if row['terminal']:events.append((tid+':terminal','endpoint terminal'))
            if snap['terminal_endpoints'] and snap['terminal_endpoints']<snap['total_endpoints']:
                events.append(('partial_terminal:'+str(snap['terminal_endpoints']),'partial completion leaves allocated CPUs without remaining endpoint work'))
            if job_state in TERMINAL:events.append(('job_terminal','scheduler terminal '+job_state))
            unseen=[(key,why) for key,why in events if key not in state['notifications']]
            state.update(snapshot=snap,job_state=job_state,last_checked_utc=now());save(receipt,state)
            if unseen:
                keys=[k for k,_ in unseen]
                message='[AUTOMATED OWNED SLURM HEALTH EVENT; not a human message]\nJob '+a.job+' '+job_state+'\n'+ '\n'.join(k+': '+why for k,why in unseen)+'\nReceipt: '+str(receipt)+'\nReview actual residuals/endpoint progress and decide recovery. Alert-only watcher has not cancelled or changed any job. Idle-rank numbers are manifest estimates, not measured utilization.'
                for k in keys:state['notifications'][k]={'status':'enqueue_pending','utc':now()}
                save(receipt,state)
                delivery=subprocess.run([a.codex,'queue','--thread',a.thread,'--message',message],capture_output=True,text=True,timeout=60)
                for k in keys:state['notifications'][k]={'status':'queued' if delivery.returncode==0 else 'enqueue_failed','utc':now(),'stdout':delivery.stdout,'stderr':delivery.stderr}
                save(receipt,state)
                if delivery.returncode:raise SystemExit('Wake enqueue failed; receipt retained for explicit recovery')
            if job_state in TERMINAL:
                state['status']='terminal';save(receipt,state);return
            time.sleep(a.interval)

if __name__=='__main__':main()
