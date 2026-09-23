"""Finite read-only completion monitor for the two declared solvent jobs."""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import subprocess
import sys
import time


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--primary',required=True); p.add_argument('--matched',required=True)
    p.add_argument('--jobs',nargs=2,required=True); p.add_argument('--report',required=True)
    p.add_argument('--vault',required=True); p.add_argument('--once',action='store_true')
    args=p.parse_args(); matched=Path(args.matched).resolve(); primary=Path(args.primary).resolve()
    while True:
        running=subprocess.check_output(['squeue','--noheader','--jobs',','.join(args.jobs),'--format','%i %T'],text=True).strip()
        if running:
            print(json.dumps({'status':'waiting','actual_queue':running}),flush=True)
            if args.once:return
            time.sleep(60); continue
        break
    out=matched.parent/'final_collection.json'
    if not out.exists():
        subprocess.run([sys.executable,str(matched.parent/'implementation/solvent_cpcm_matched.py'),
                        'collect','--manifest',str(matched),'--output',str(out)],check=True)
    result=json.loads(out.read_text())
    acct=subprocess.check_output(['sacct','-j',','.join(args.jobs),'--parsable2',
        '--format=JobID,State,ElapsedRaw,AllocCPUS,CPUTimeRAW,MaxRSS,ExitCode,NodeList'],text=True)
    accounting=matched.parent/'scheduler_accounting_final.txt'
    if not accounting.exists():accounting.write_text(acct)
    events=[]
    for mp in (primary,matched):
        ep=mp.parent/'budget_events.jsonl'
        if ep.exists():events.extend(json.loads(s) for s in ep.read_text().splitlines() if s.strip())
    alloc=sum(x['allocated_core_seconds'] for x in events if 'allocated_core_seconds' in x)
    costs=[]
    for mp in (primary,matched):
        m=json.loads(mp.read_text())
        for task in m['tasks']:
            rp=Path(task['output_path']+'.execution.json')
            if not rp.exists():
                costs.append({'task_id':task['task_id'],'receipt':None});continue
            r=json.loads(rp.read_text())
            secs=(datetime.fromisoformat(r['finished_at_utc'])-datetime.fromisoformat(r['started_at_utc'])).total_seconds()
            costs.append({'task_id':task['task_id'],'receipt':str(rp),'wall_seconds':secs,
                          'rank_seconds':secs*r['parallelism']['nprocs'],'returncode':r['returncode'],
                          'scf_converged':r['scf_converged'],'normal_termination':r['normal_termination']})
    cost_path=matched.parent/'final_costs.json'
    if not cost_path.exists():cost_path.write_text(json.dumps({'allocated_core_seconds':alloc,'GPU_seconds':0,'attempts':costs},indent=2)+'\n')
    available=result['available_pairs']; n=result['complete_endpoints']
    lines=['# Matched CPCM pilot result','',
        ('Four complete matched pairs are available for development comparison; no new reference or prediction band has been calibrated.' if available==4 else
         f'The CPCM replacement is not qualified for benchmark expansion: only {available}/4 matched metal pairs are available. Failed or unsupported endpoints remain missing.'),'',
        'ORCA actually switched native GFN2 to ordinary ORCA SCF for CPCM. The second submission used the same eight physical states for matched ordinary-SCF vacuum controls, with explicit archived-GBW MORead. No converged CPCM cell was rerun.','',
        f'Actual audited endpoint coverage: {n}/16, comprising '+str(result['reused_CPCM_complete'])+'/8 original CPCM endpoints and '+str(n-result['reused_CPCM_complete'])+'/8 matched vacuum endpoints. Sixteen total scientific attempts were declared; all receipt statuses are retained.','',
        '| Case | Matched status | Released ALPB R | New CPCM R | Change |',
        '|---|---|---:|---:|---:|']
    fmt=lambda x:'unavailable' if x is None else f'{x:.8f}'
    for row in result['rows']:
        lines.append('| '+row['case_id']+' | '+row['status']+' | '+fmt(row['released_ALPB_R'])+' | '+fmt(row['CPCM_R'])+' | '+fmt(row['CPCM_minus_ALPB_R'])+' |')
    lines+=['','R is E(Ca)−E(La), in model kcal/mol. The complete matched CPCM-minus-vacuum transfer replaces ALPB; the printed dielectric component is not itself the correction. Missing values are not zero. These consumed development cases cannot establish broad affinity discrimination.','',
            '## Numerical and physical qualification','']
    for e in result['endpoints']:
        if e['status']!='complete':lines.append(f'- {e["case_id"]} {e["metal"]} {e["medium"]}: {e.get("reason","unavailable")}.')
    lines+=['','Each successful endpoint retains its native parameter export, actual solver/temperature/state/charge audit and receipt. Successful CPCM endpoints additionally retain installed dielectric/radii, Gaussian-vdW surface, surface density and a nonzero printed reaction-field contribution. Ordinary-minus-native vacuum energy differences are reported separately in the collection.','',
            '## Cost and next action','',f'Actual executor allocation: **{alloc:.3f} core-seconds; 0 GPU-seconds**. Endpoint wall/rank-seconds, failures and scheduler accounting are in the linked records. This includes idle resources within the two allocated jobs; it is not a matched-hardware speed comparison.','',
            ('Review the four raw matched contrasts and numerical differences before preparing the complete original30 canonical/crystal/unknown population. Do not calibrate on this four-case pilot.' if available==4 else
             'Review the explicit SCF/state failures before any retry. No full30 calibration or225-fold transfer is launched, and no new decision is available. Preserve the released ALPB default.'),'',
            f'- Final collection: `{out}`',f'- Costs: `{cost_path}`',f'- Scheduler: `{accounting}`','']
    report=Path(args.report); report.parent.mkdir(parents=True,exist_ok=True)
    if not report.exists():report.write_text('\n'.join(lines))
    vault=Path(args.vault); vault.parent.mkdir(parents=True,exist_ok=True)
    note='\n'.join(['# Nikasha matched CPCM solvent pilot','',datetime.now(timezone.utc).isoformat(),'',lines[2],'',
        f'ORCA6.1.1 forced ordinary SCF. Matched vacuum controls were actually run; audited coverage {n}/16 endpoints, {available}/4 pairs. No baseline/default change or new calibration. Cost {alloc:.3f} allocated core-seconds, zero GPU.','',
        f'Report: {report}',f'Collection: {out}',f'Costs: {cost_path}',''])
    if not vault.exists():vault.write_text(note)
    print(json.dumps({'status':'collected','available_pairs':available,'endpoints':n,'report':str(report),'vault':str(vault)}),flush=True)


if __name__=='__main__':main()
