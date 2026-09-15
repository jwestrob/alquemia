"""Watch only a named task-owned job and write a terminal accounting receipt."""
import argparse
import json
from pathlib import Path
import subprocess
import time

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--job',required=True)
p.add_argument('--output',type=Path,required=True)
a=p.parse_args()
while True:
    r=subprocess.run(['sacct','-j',a.job,'--format=JobID,State,AllocCPUS,ElapsedRaw,CPUTimeRAW,MaxRSS,NodeList','-P'],capture_output=True,text=True)
    rows=[x.split('|') for x in r.stdout.splitlines()[1:] if x.startswith(a.job+'|')]
    if rows and rows[0][1].split()[0] in ('COMPLETED','FAILED','CANCELLED','TIMEOUT','OUT_OF_MEMORY','NODE_FAIL','PREEMPTED','BOOT_FAIL'):
        with a.output.open('x') as f: json.dump({'job_id':a.job,'sacct':r.stdout,'finished_observed_unix':time.time()},f,indent=2)
        print(r.stdout,flush=True)
        break
    time.sleep(30)
