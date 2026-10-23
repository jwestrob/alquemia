#!/usr/bin/env python3
"""Only the exclusive reference-save correction and exact preparation delta."""
import argparse,importlib.util,json,time
from pathlib import Path
HERE=Path(__file__).resolve().parent
s=importlib.util.spec_from_file_location('save_delta_binding',HERE/'BIND.py');B=importlib.util.module_from_spec(s);s.loader.exec_module(B)
p=argparse.ArgumentParser();p.add_argument('--old',required=True);p.add_argument('--new',required=True);p.add_argument('--out',required=True);a=p.parse_args();tick,cpu=time.monotonic(),time.process_time()
out=Path(a.out);out.mkdir(parents=True,exist_ok=False);old=B.read(B.pin(a.old));new=B.read(B.pin(a.new));checks=[]
def check(ok,name):B.require(ok,name);checks.append(dict(name=name,passed=True))
check({k:v for k,v in old.items() if k!='sources'}=={k:v for k,v in new.items() if k!='sources'},'Prepared science/input metadata unchanged')
check(all(old['sources'][k]==new['sources'][k] for k in ('BIND','MODEL','COLLECT')),'Only reference source changed')
r=B.module(new['sources']['REFERENCE'],'save_delta_reference');f=out/'EXCLUSIVE.json';r.save(f,dict(first=True));before=B.pin(f)
try:r.save(f,dict(first=False))
except FileExistsError:check(B.pin(f)==before,'Exclusive save rejects overwrite and preserves evidence')
else:raise AssertionError('Exclusive save overwrote existing evidence')
B.save(out/'RESULT.json',dict(status='passed_exclusive_save_delta',implementation=B.pin(__file__),old_manifest=B.pin(a.old),manifest=B.pin(a.new),sources=new['sources'],checks=checks,new_native_calls=0,fresh_outcomes_read=False,wall_seconds=time.monotonic()-tick,CPU_seconds=time.process_time()-cpu,local_threads=1))
