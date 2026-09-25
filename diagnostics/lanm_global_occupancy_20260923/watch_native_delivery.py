"""Job-specific, chat-independent collection and email; never submits work."""
import argparse
import datetime as dt
from email.message import EmailMessage
from email.utils import formatdate,make_msgid
import fcntl
import json
from pathlib import Path
import re
import subprocess
import time

ROOT=Path(__file__).resolve().parents[2]
DIAG=Path(__file__).resolve().parent
PY='/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python'
TERMINAL={'COMPLETED','FAILED','CANCELLED','TIMEOUT','OUT_OF_MEMORY','NODE_FAIL',
          'BOOT_FAIL','DEADLINE','PREEMPTED','REVOKED','SPECIAL_EXIT'}

def terminal_record(job):
    p=subprocess.run(['sacct','-n','-X','-j',job,'--format=JobIDRaw,State,ExitCode,ElapsedRaw,AllocCPUS,NodeList','-P'],
                     capture_output=True,text=True)
    if p.returncode:return None
    for line in p.stdout.splitlines():
        f=line.split('|')
        if len(f)>=6 and f[0]==job and f[1].split()[0].rstrip('+') in TERMINAL:
            return dict(job_id=f[0],state=f[1],exit_code=f[2],elapsed_seconds=int(f[3]),
                        allocated_cpus=int(f[4]),node=f[5])
    return None

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--job',required=True)
    p.add_argument('--manifest',type=Path,required=True)
    p.add_argument('--seeded-vacuum',action='store_true')
    a=p.parse_args()
    if not re.fullmatch(r'\d+',a.job):raise ValueError('numeric job ID required')
    lock=(DIAG/f'delivery_{a.job}.lock').open('a+')
    fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    receipt=DIAG/f'DELIVERY_{a.job}.json'
    if receipt.exists():return
    while True:
        result=terminal_record(a.job)
        if result is not None:break
        time.sleep(60)
    cp=a.manifest.parent/'COLLECTION.json'
    if a.manifest.exists() and not cp.exists():
        command=([PY,str(DIAG/'seed_vacuum_recovery.py'),'collect'] if a.seeded_vacuum else
                 [PY,str(DIAG/'collect_native_recovery.py')])
        r=subprocess.run(command+['--manifest',str(a.manifest)],
                         cwd=ROOT,capture_output=True,text=True)
        result['collector_recovery']={'exit_code':r.returncode,'stdout':r.stdout,'stderr':r.stderr}
    collection=json.loads(cp.read_text()) if cp.exists() else None
    count=collection['complete_cells'] if collection else 0
    denominator=2 if a.seeded_vacuum else 4
    lines=['# Whole-protein LanM native feasibility', '',
           f"Job {a.job}: {result['state']}; {count}/{denominator} native origin cells passed collection.",
           f"Node {result['node']}; {result['allocated_cpus']} allocated CPUs; {result['elapsed_seconds']} elapsed seconds.", '',
           ('These are two same-metal ALPB-seeded vacuum recovery endpoints; solvent energies are reused. Single-seed convergence is not a ground-state qualification.' if a.seeded_vacuum else
            'These are the same Hans8DQ2 EF12 La/Dy vacuum/ALPB endpoints.')+' No new MACE or DFT ran.',
           'Both earlier relaxed geometries remain rejected for covalent distortion. This result tests native whole-protein execution; it is not a validated within-series preference or accommodated score.', '']
    if a.manifest.exists():
        manifest=json.loads(a.manifest.read_text())
        lines.append('Resources: '+json.dumps(manifest['execution_resources'],sort_keys=True))
        for task in manifest['tasks']:
            row=next((r for r in (collection or {}).get('rows',[]) if r['task_id']==task['task_id']),{})
            lines.append(f"- {task['task_id']}: {row.get('status','unavailable')}; energy Eh={row.get('energy_hartree')}; reason={row.get('reason','none')}")
            op=Path(task['output_path'])
            if op.exists() and row.get('status')!='complete':
                errors=[s.strip() for s in op.read_text(errors='replace').splitlines()
                        if re.search(r'Not enough memory|Error\s*\(ORCA|SCF NOT CONVERGED|SCF failed to converge|ORCA finished by error',s,re.I)]
                lines.extend('  '+s for s in errors[:6])
    lines += ['',f'Collection: {cp}',
              'No automatic resubmission or expansion is attached. PQQ production is unchanged.', '']
    body='\n'.join(lines)
    report=DIAG/f'NATIVE_RESULT_{a.job}.md'
    with report.open('x') as f:f.write(body)
    vault=Path('/home/jwestrob/jwestrob/obsidian-vault/agent-captures')/f'2026-09-24_Nikasha-LanM-native-{a.job}.md'
    try:
        with vault.open('x') as f:f.write(body)
        result['vault']=str(vault)
    except OSError as exc:result['vault_error']=str(exc)
    msg=EmailMessage();msg['To']='jacobwestroberts@gmail.com';msg['From']='jwestrob@biotite.berkeley.edu'
    msg['Subject']=f'Nikasha LanM: native run {a.job} {result["state"]}, {count}/{denominator} cells'
    msg['Date']=formatdate(localtime=False);msg['Message-ID']=make_msgid();msg.set_content(body)
    eml=DIAG/f'native_result_{a.job}.eml'
    with eml.open('xb') as f:f.write(msg.as_bytes())
    mail=subprocess.run(['/usr/sbin/sendmail','-t'],input=msg.as_bytes(),capture_output=True)
    result.update(complete_cells=count,collection=str(cp) if collection else None,report=str(report),
                  completed_utc=dt.datetime.now(dt.timezone.utc).isoformat(),
                  email_status='accepted_by_local_mail_relay' if mail.returncode==0 else 'relay_failed',
                  email_returncode=mail.returncode,email_stderr=mail.stderr.decode(errors='replace'),
                  execution_scope='collection/notification only; no job submissions')
    with receipt.open('x') as f:json.dump(result,f,indent=2)
    print(json.dumps(result),flush=True)

if __name__=='__main__':main()
