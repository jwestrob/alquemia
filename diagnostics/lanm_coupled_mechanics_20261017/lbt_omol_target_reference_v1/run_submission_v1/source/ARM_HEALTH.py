from pathlib import Path
import json,subprocess,time,sys,datetime,os
wd=Path(sys.argv[1]);job=sys.argv[2];manifest=wd/'execution'/sys.argv[6]/'MANIFEST.json'
root=Path(sys.argv[3]);thread=sys.argv[4];codex=sys.argv[5]
while not manifest.is_file():
 p=subprocess.run(['sacct','-n','-P','-j',job,'--format=JobIDRaw,State'],capture_output=True,text=True)
 states=[x.split('|')[1].split()[0].rstrip('+') for x in p.stdout.splitlines() if x.startswith(job+'|')]
 if states and states[0] in {'COMPLETED','FAILED','CANCELLED','TIMEOUT','OUT_OF_MEMORY','NODE_FAIL','BOOT_FAIL','DEADLINE','PREEMPTED','REVOKED'}:
  (wd/'HEALTH_ARM_TERMINAL_NO_MANIFEST.json').write_text(json.dumps({'job':job,'state':states[0],'utc':datetime.datetime.now(datetime.timezone.utc).isoformat()}));raise SystemExit(0)
 time.sleep(15)
cmd=['/usr/bin/python3',str(root/'scripts/slurm_native_health_watch.py'),'--manifest',str(manifest),'--job',job,'--thread',thread,'--codex',codex,'--receipt',str(wd/'health_receipt.json'),'--interval','30']
(wd/'HEALTH_ARMED.json').write_text(json.dumps({'job':job,'pid':os.getpid(),'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'command':cmd},indent=2))
os.execv(cmd[0],cmd)
