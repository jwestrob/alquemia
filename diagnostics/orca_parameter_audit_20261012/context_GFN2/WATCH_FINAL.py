"""Durable completion notification for the explicitly owned finite local worker."""
import argparse, json, os, subprocess, time
from pathlib import Path

def main():
    p=argparse.ArgumentParser();p.add_argument('--worker-pid',type=int,required=True)
    p.add_argument('--worker-start',required=True);p.add_argument('--directory',required=True)
    p.add_argument('--thread',required=True);p.add_argument('--codex',required=True);a=p.parse_args()
    base=Path(a.directory);receipt=base/'completion_event.json'
    if receipt.exists():raise RuntimeError('existing completion event: inspect before rearming')
    while True:
        result=base/'FINAL_COLLECTION.json'
        try:data=json.loads(result.read_text());state='completed';break
        except (FileNotFoundError,json.JSONDecodeError):pass
        try:start=Path(f'/proc/{a.worker_pid}/stat').read_text().split()[21]
        except FileNotFoundError:start=None
        if start!=a.worker_start:data=None;state='worker_disappeared_without_final_collection';break
        time.sleep(5)
    event={'event':'context-factorial-solvent-v2','state':state,'result':str(result),'created_epoch':time.time(),
           'worker_pid':a.worker_pid,'thread':a.thread,'scientific_jobs_modified':False}
    if data:event.update(executed_inputs=data['executed_inputs'],qualified_components=data['accepted_components'],declared_components=data['declared_components'])
    message='[AUTOMATED OWNED COMPLETION EVENT; not a human message]\nEvent: context-factorial-solvent-v2\nReceipt: '+str(receipt)+'\nState: '+state+'\nResult: '+str(result)+'\nRead actual results; preserve failed absolute v1 gate. Integrate only the already-authorized context factorial. No new optimization is authorized by this event.'
    event['queue_attempts']=[]
    tmp=receipt.with_suffix('.tmp');tmp.write_text(json.dumps(event,indent=2)+'\n');tmp.replace(receipt)
    for attempt in range(3):
        try:
            r=subprocess.run([a.codex,'queue','--thread',a.thread,'--message',message],text=True,capture_output=True,timeout=30)
            q={'returncode':r.returncode,'stdout':r.stdout,'stderr':r.stderr}
        except Exception as exc:q={'error':str(exc),'returncode':None}
        event['queue_attempts'].append(q)
        tmp=receipt.with_suffix('.tmp');tmp.write_text(json.dumps(event,indent=2)+'\n');tmp.replace(receipt)
        if q['returncode']==0:break
        time.sleep(5)
    print(json.dumps(event),flush=True)

if __name__=='__main__':main()
