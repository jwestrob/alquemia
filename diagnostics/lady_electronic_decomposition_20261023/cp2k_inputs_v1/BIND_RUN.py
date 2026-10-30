#!/usr/bin/env python3
"""Freeze run configuration and finite Slurm scripts; no submission."""
import hashlib,json,shlex
from pathlib import Path
B=Path(__file__).resolve().parent
R=B.parents[2]
W=R/'workspaces/lady_electronic_decomposition_20261023'
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1024*1024),b''): h.update(b)
 return h.hexdigest()
def once(p,s):
 p.parent.mkdir(parents=True,exist_ok=True)
 if p.exists():
  assert p.read_text()==s, 'Refuse changed existing '+str(p)
 else: p.write_text(s)
def main():
 systems=json.loads((B/'DECK_PROVENANCE.json').read_text())['systems']
 taskroot=W/'cp2k_quartet_v1'
 for gid,s in systems.items(): s['task']=str(taskroot/'tasks'/gid)
 image=W/'cp2k_build_v1/cp2k-2026.2-mpich-generic.sif'
 binary=W/'cp2k_build_v2/build/bin/cp2k.psmp'
 files=[B/n for n in ['RUN.py','PARSE.py','PREPARE.py','BIND_RUN.py','RUN_SPEC.md','ENERGY_LEDGER.json','DECK_PROVENANCE.json','PILOT_SETTINGS_REVIEW.md','DECK_REVIEW.md','RUNNER_REVIEW.md','assets/BASIS_LADY','assets/POTENTIAL_LADY','assets/PROVENANCE.json']]
 files += [image,binary]
 files += [Path(s['deck']) for s in systems.values()]
 files += [Path(s['xyz']) for s in systems.values()]
 config={'version':1,'pilot':'aquo8_central_La','systems':systems,'cpus':24,'mpi_ranks':12,'omp_threads':2,'native_cap_seconds':7200,'assets':str(B/'assets'),'parser':str(B/'PARSE.py'),'pilot_admission':str(taskroot/'PILOT_ADMISSION.json'),'diagnostic_results':str(B.parent/'cp2k_results_v1'),'runtime':{'image':str(image),'build':str(W/'cp2k_build_v2'),'source':str(W/'cp2k_build_v1/cp2k-2026.2')},'pins':[{'path':str(p),'sha256':sha(p)} for p in files]}
 assert sha(image)=='90c3a7a7368e8025689a3ea13ba773ecf832a94ea5378c984de1b861c9d8deab'
 assert sha(binary)=='19e4c45c2920299863d0d74543b55ded033ab1e145cbba1cd5cfc3d660e23b32'
 path=B/'CONFIG.json'
 once(path,json.dumps(config,indent=2)+'\n')
 for gid,s in systems.items():
  task=Path(s['task']);task.mkdir(parents=True,exist_ok=True)
  argv=['/home/jwestrob/.pyenv/versions/3.12.0/bin/python',str(B/'RUN.py'),str(path),str(task),gid]
  script='\n'.join(['#!/bin/bash','#SBATCH --job-name=lady-almo-'+gid,'#SBATCH --account=banfield','#SBATCH --partition=standard','#SBATCH --qos=standard_expand','#SBATCH --nodes=1','#SBATCH --ntasks=1','#SBATCH --cpus-per-task=24','#SBATCH --no-requeue','set -euo pipefail','exec '+shlex.join(argv),''])
  once(task/'job.sbatch',script)
 print(json.dumps({'configuration':str(path),'sha256':sha(path),'taskroot':str(taskroot),'pins':len(config['pins']),'submitted':False},indent=2))
if __name__=='__main__':main()
