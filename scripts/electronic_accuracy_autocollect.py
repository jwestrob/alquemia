#!/usr/bin/env python3
"""One-job, read-only completion monitor; no scientific execution or resubmission."""
from __future__ import annotations
import argparse
from datetime import datetime,timezone
from email.message import EmailMessage
import fcntl
import json
import os
from pathlib import Path
import re
import subprocess
import time

from affordable_common import record,verify,read_json,write_new
from electronic_accuracy_collect import collect

TERMINAL={'COMPLETED','FAILED','CANCELLED','TIMEOUT','OUT_OF_MEMORY','NODE_FAIL',
          'BOOT_FAIL','DEADLINE','PREEMPTED','REVOKED'}


def utc():return datetime.now(timezone.utc).isoformat()


def scheduler(job_id):
    """Neither an empty queue nor a failed scheduler query implies completion."""
    fields=['JobIDRaw','State','ElapsedRaw','AllocCPUS','TotalCPU','MaxRSS','NodeList']
    command=['sacct','-n','-P','-j',str(job_id),'--format='+','.join(fields)]
    p=subprocess.run(command,capture_output=True,text=True,timeout=45)
    if p.returncode:raise RuntimeError('sacct query failed: '+p.stderr.strip())
    rows=[dict(zip(fields,line.split('|'))) for line in p.stdout.splitlines() if line.strip()]
    main=next((r for r in rows if r['JobIDRaw']==str(job_id)),None)
    state=None if main is None else main['State'].split()[0].rstrip('+')
    q=subprocess.run(['squeue','-h','-j',str(job_id),'-o','%T'],capture_output=True,text=True,timeout=45)
    if q.returncode and 'Invalid job id' not in q.stderr:raise RuntimeError('squeue query failed: '+q.stderr.strip())
    return {'utc':utc(),'job_id':str(job_id),'command':command,'rows':rows,'state':state,
        'queue_states':q.stdout.splitlines(),'terminal':state in TERMINAL and not q.stdout.strip(),
        'allocated_core_seconds':None if main is None or not main['ElapsedRaw'].isdigit() or not main['AllocCPUS'].isdigit()
            else int(main['ElapsedRaw'])*int(main['AllocCPUS'])}


def report(result,accounting,manifest):
    n=result['complete_endpoints'];total=result['total_endpoints']
    text=['# Independent electronic reference: '+('complete' if n==total else 'incomplete or unsupported'),'',
        f'Collected {utc()}. Validated correlated endpoints: **{n}/{total}**.',
        f'Scheduler state: `{accounting.get("state")}`. No scientific reruns or added endpoints.','',
        '## Endpoint status','', '| Case | Metal | Status | CC total (hartree) |',
        '|---|---|---|---:|']
    for r in result['rows']:
        v='unavailable' if r['result'] is None else format(r['result']['total_hartree'],'.12f')
        text.append(f"| {r['case']} | {r['metal']} | {r['status']} | {v} |")
    text+=['','## Fixed intersite comparisons','',
        'D = (R_alpha − R_GGR) × 627.509474 kcal/mol; R = E_Ca − E_La.',
        'Positive D has the expected relative direction. No old aquo reference or bands are inherited.','',
        '| Alpha preparation | Native baseline D | Correlated D |', '|---|---:|---:|']
    for case in ('ALPHA_1F6S_original','ALPHA_1F6S_water_prepared'):
        b=result['baseline']['alpha_minus_GGR_kcal_mol'].get(case)
        c=result['correlated']['alpha_minus_GGR_kcal_mol'].get(case)
        text.append(f'| {case} | '+('unavailable' if b is None else f'{b:+.6f}')+' | '+('unavailable' if c is None else f'{c:+.6f}')+' |')
    shift=result['correlated']['water_preparation_delta_R_kcal_mol']
    if shift is not None:
        text+=['',f'Water preparation changes the correlated contrast by **{shift:+.6f} kcal/mol**.',
            ('It moves alpha in the more La-like relative direction.' if shift>0 else
             'It does not move alpha in the more La-like relative direction.')]
    if n<total:
        text+=['','**The full electronic-method comparison is unavailable.** Missing, failed or unparsed',
               'endpoints were not filled with baseline energies. Completed partial pairs remain explicit.']
    errors=[f"- {r['task_id']}: {r.get('reason',r['status'])}" for r in result['rows'] if r['status']!='complete']
    if errors:text+=['','## Unavailable endpoints','']+errors
    text+=['','## Interpretation and cost','',
        'One consumed, condition-qualified alpha/GGR biological comparison; the two alpha preparations',
        'are not independent biological samples. This is a finite-basis electronic-recipe diagnostic,',
        'not a binding free energy or broad affinity validation. Basis, correlation and solvent response',
        'change together, so a shift does not identify a unique source of baseline error.',
        'PTES/core accounting and CC diagnostics are retained in collection.json. No numerical basis/PNO',
        'refinement has yet established a high-accuracy limit. The baseline Hamiltonian is unchanged;',
        "Jacob's separately authorized water-preparation promotion is not a promotion of this CC method.",'',
        f"Measured job allocation: {accounting.get('allocated_core_seconds')} CPU-core-seconds (null means unavailable).",
        'See accounting.json for elapsed time, allocated CPUs, CPU time, memory records and node.',
        'GPU allocations: none. Preparation, fixture tests and collection CPU cost were not separately metered.',
        'This expensive reference is not an affordable routine scorer. Larger PQQ and replica checks remain unsubmitted.','',
        'Recommendation: retain the native Hamiltonian. Use the complete reference result, when available,',
        'to decide whether electronic refinement adds useful information before another expensive panel.',
        'Next informative candidate is the narrow-margin alpha/GGR2FW0 comparison, with this same method;',
        'it is a proposal, not a submitted calculation. The declared PQQ guardrail also remains unrun.','',
        f'Manifest: `{manifest}`. Full unrounded components and execution receipts: `collection.json`.','']
    return '\n'.join(text)


def cleanup_successful(result,manifest):
    """Only generated matrices for strictly validated successful owned endpoints."""
    mp=Path(manifest).resolve();m=read_json(mp)
    rows={r['task_id']:r for r in result['rows']};deleted=[];errors=[]
    for t in m['tasks']:
        if rows[t['task_id']]['status']!='complete':continue
        td=Path(t['output_path']).parent.resolve()
        if td.parent!=mp.parent:raise RuntimeError('cleanup directory is outside this pilot')
        for p in sorted(td.iterdir()):
            if p.is_symlink() or not p.is_file() or not re.fullmatch(r'endpoint\.runtime\..+\.tmp(?:\.\d+)?',p.name):continue
            item={'path':str(p),'bytes':p.stat().st_size}
            try:p.unlink();deleted.append(item)
            except OSError as exc:errors.append({**item,'error':str(exc)})
    return {'utc':utc(),'policy':'validated successful own temporary matrices only; scientific artifacts retained',
        'removed':deleted,'removed_bytes':sum(x['bytes'] for x in deleted),'errors':errors}


def notify(cfg,out,summary):
    if (out/'email_receipt.json').exists():return
    msg=EmailMessage();msg['To']=cfg['recipient'];msg['Subject']='Alquemia: independent electronic reference status'
    msg.set_content(summary+'\n\nFull report: '+str(out/'REPORT.md')+'\n')
    eml=out/'message.eml';eml.write_bytes(msg.as_bytes())
    try:
        p=subprocess.run(['/usr/sbin/sendmail','-t'],input=msg.as_bytes(),capture_output=True,timeout=45)
        r={'returncode':p.returncode,'stdout':p.stdout.decode(errors='replace'),'stderr':p.stderr.decode(errors='replace'),
            'delivery_status':'local relay accepted; mailbox delivery not independently verified' if p.returncode==0 else 'relay_failed'}
    except (OSError,subprocess.TimeoutExpired) as exc:r={'returncode':None,'delivery_status':'relay_error_or_uncertain','error':str(exc)}
    write_new(out/'email_receipt.json',{'utc':utc(),'recipient':cfg['recipient'],'message':record(eml),**r})


def finalize(cfg,accounting):
    if not accounting['terminal']:raise RuntimeError('finalization requires a terminal job absent from queue')
    manifest=verify(cfg['manifest']);out=Path(cfg['output']);out.mkdir(parents=True,exist_ok=True)
    if (out/'DONE.json').exists():return
    # Never act while an executor still owns this manifest, even if accounting
    # briefly reports a terminal state. This locks rather than alters its file.
    with (manifest.parent/'execute.lock').open('a+') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        if (out/'accounting.json').exists():accounting=read_json(out/'accounting.json')
        else:write_new(out/'accounting.json',accounting)
        try:
            if (out/'collection.json').exists():result=read_json(out/'collection.json')
            else:
                result=collect(manifest);write_new(out/'collection.json',result)
            body=report(result,accounting,manifest)
            (out/'REPORT.md').write_text(body)
            Path(cfg['report_path']).write_text(body)
            if not (out/'scratch_cleanup.json').exists():
                cleanup=cleanup_successful(result,manifest);write_new(out/'scratch_cleanup.json',cleanup)
            n=result['complete_endpoints'];summary=f'Hi Jacob,\n\nThe six-endpoint independent electronic-reference job ended ({accounting["state"]}); {n}/6 endpoints passed collection.\n'
            if n==6:
                ds=result['correlated']['alpha_minus_GGR_kcal_mol'];shift=result['correlated']['water_preparation_delta_R_kcal_mol']
                summary+=f'Alpha-minus-GGR contrasts: original {ds["ALPHA_1F6S_original"]:+.3f}, prepared-water {ds["ALPHA_1F6S_water_prepared"]:+.3f} kcal/mol. The corresponding native DFT values are -14.559 and +10.924. Water preparation shifts the CC contrast by {shift:+.3f} kcal/mol.\n'
                summary+='This supports the same relative ordering after water preparation.\n' if ds['ALPHA_1F6S_water_prepared']>0 else 'The correlated recipe does not reproduce the expected prepared-alpha/GGR ordering; this needs review before attributing a cause.\n'
            else:
                summary+='A complete Ca/La comparison is unavailable. Missing/failed/unparsed endpoints remain explicit, with no substituted energies.\n'
                summary+='; '.join(f"{r['task_id']}: {r.get('reason',r['status'])}" for r in result['rows'] if r['status']!='complete')+'\n'
            summary+=f'Job allocation: {accounting["allocated_core_seconds"]} CPU-core-seconds; no GPUs. This is a costly finite-basis diagnostic on one consumed biological comparison, not broad validation or a new production score. No extra calculations were launched.\n'
        except Exception as exc:
            err={'utc':utc(),'error':f'{type(exc).__name__}: {exc}','scientific_reruns':0}
            if not (out/'collection_error.json').exists():write_new(out/'collection_error.json',err)
            body='# Electronic reference: collection requires review\n\n'+err['error']+'\n\nNo missing energy was inferred. Raw scientific outputs are retained.\n'
            (out/'REPORT.md').write_text(body);Path(cfg['report_path']).write_text(body)
            summary=f'Hi Jacob,\n\nElectronic-reference job {cfg["job_id"]} ended ({accounting["state"]}), but its final collection requires review: {err["error"]}. No completed discrimination result is claimed and no additional calculation was launched. Raw outputs remain preserved.\n'
        Path(cfg['vault_path']).write_text(body+'\n\nPrimary artifacts: '+str(out)+'\n')
        notify(cfg,out,summary)
        write_new(out/'DONE.json',{'utc':utc(),'report':record(out/'REPORT.md'),'email_receipt':record(out/'email_receipt.json'),
            'scientific_reruns':0,'monitor':record(__file__)})


def run(config,preview=False):
    cfg=read_json(config);verify(cfg['manifest'])
    if preview:
        a=scheduler(cfg['job_id']);r=collect(verify(cfg['manifest']))
        return {'scheduler':a,'complete_endpoints':r['complete_endpoints'],'statuses':{x['task_id']:x['status'] for x in r['rows']},
            'correlated':r['correlated'],'would_finalize':a['terminal'],'email_sent':False,'cleanup_executed':False}
    directory=Path(config).resolve().parent
    with (directory/'monitor.lock').open('a+') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        write_new(directory/f'monitor_started_{os.getpid()}.json',{'utc':utc(),'pid':os.getpid(),'configuration':record(config),'implementation':record(__file__)})
        while True:
            try:
                a=scheduler(cfg['job_id'])
                if a['terminal']:finalize(cfg,a);return {'status':'finished'}
            except Exception as exc:
                with (directory/'monitor_errors.jsonl').open('a') as f:f.write(json.dumps({'utc':utc(),'error':str(exc)})+'\n')
            time.sleep(cfg['poll_seconds'])


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--config',required=True);p.add_argument('--preview',action='store_true')
    a=p.parse_args();print(json.dumps(run(a.config,a.preview),indent=2))
