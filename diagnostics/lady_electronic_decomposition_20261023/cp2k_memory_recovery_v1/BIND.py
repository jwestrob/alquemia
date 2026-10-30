#!/usr/bin/env python3
"""Pin one resource-only recovery without changing any electronic input."""
import hashlib,json,shlex
from pathlib import Path
B=Path(__file__).resolve().parent
D=B.parent
R=B.parents[2]
W=R/'workspaces/lady_electronic_decomposition_20261023'
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
 return h.hexdigest()
def once(p,s):
 p.parent.mkdir(parents=True,exist_ok=True)
 if p.exists():assert p.read_text()==s,'Refuse changed '+str(p)
 else:p.write_text(s)
def main():
 old=D/'cp2k_inputs_v1/CONFIG.json'
 assert sha(old)=='075aaf39b904651e0cd406fc127aa23523986c344b294fda1b61ebeb554db423'
 c=json.loads(old.read_text());taskroot=W/'cp2k_quartet_memory_v1'
 c['version']=2;c['native_cap_seconds']=7073;c['memory_MiB']=262144
 c['recovery_of']={'configuration':str(old),'sha256':sha(old),'job_id':'21618','reason':'initialization_OUT_OF_MEMORY','charged_prior_seconds':127,'pilot_cumulative_native_budget_seconds':7200}
 c['diagnostic_results']=str(D/'cp2k_results_memory_v1')
 c['pilot_admission']=str(taskroot/'PILOT_ADMISSION.json')
 for gid,s in c['systems'].items():s['task']=str(taskroot/'tasks'/gid)
 for p in c['pins']:assert sha(Path(p['path']))==p['sha256'],p
 for p in [old,B/'RUN.py',B/'BIND.py',B/'RECOVERY_SPEC.md',B/'RECOVERY_REVIEW.md',B/'RUNNER_REVIEW.md',B/'RESOURCE_EVIDENCE.json',B/'FAILED_ATTEMPT.json']:
  c['pins'].append({'path':str(p),'sha256':sha(p)})
 path=B/'CONFIG.json';once(path,json.dumps(c,indent=2)+'\n')
 for gid,s in c['systems'].items():
  task=Path(s['task']);task.mkdir(parents=True,exist_ok=True)
  argv=['/home/jwestrob/.pyenv/versions/3.12.0/bin/python',str(B/'RUN.py'),str(path),str(task),gid]
  script='\n'.join(['#!/bin/bash','#SBATCH --job-name=lady-almo-memory-'+gid,'#SBATCH --account=banfield','#SBATCH --partition=standard','#SBATCH --qos=standard_expand','#SBATCH --nodes=1','#SBATCH --ntasks=1','#SBATCH --cpus-per-task=24','#SBATCH --mem=256G','#SBATCH --no-requeue','set -euo pipefail','exec '+shlex.join(argv),''])
  once(task/'job.sbatch',script)
 print(json.dumps({'configuration':str(path),'sha256':sha(path),'pins':len(c['pins']),'submitted':False},indent=2))
if __name__=='__main__':main()
