#!/usr/bin/env python3
"""Exercise root admission helpers on immutable initial-ten evidence only."""
import argparse
import copy
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import time

ROOT = Path(__file__).resolve().parents[6]
REL = Path('lanm_selectivity_free_energy_20261021/explicit_model_v1/paired_water_coordination_joint_v1')
D, W = ROOT/'diagnostics'/REL, ROOT/'workspaces'/REL

def pin(path):
    path=Path(path).resolve()
    return {'path':str(path),'sha256':hashlib.sha256(path.read_bytes()).hexdigest()}

def main(out):
    start,cpu=time.monotonic(),time.process_time()
    assert all(os.environ.get(k)=='1' for k in ('OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'))
    p=W/'preparation_v1/RESULT.json'
    assert pin(p)['sha256']=='24d4d3e68f38f402b49e9bccce05853384a3a6ddfdf929176527f63dc3693fa6'
    prep=json.loads(p.read_text());config=json.loads(Path(prep['config']['path']).read_text())
    assert pin(config['training_result']['path'])==config['training_result']
    data=json.loads(Path(config['training_result']['path']).read_text())
    source=D/'selected_review_v1/REVIEW.py'
    spec=importlib.util.spec_from_file_location('selected_review_saved_admission',source)
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
    checks=[]
    for gid,value in prep['initial']['model'].items():
        for kind,fn,args in [('native',mod.native_authority,(gid,value)),('reference',mod.reference_replay,(gid,data['rows'][gid]))]:
            try:
                evidence=fn(*args)
                checks.append({'cell':gid,'kind':kind,'passed':True})
            except Exception as e:
                checks.append({'cell':gid,'kind':kind,'passed':False,'error':repr(e)})
    # These controls change only detached in-memory request/reference values.
    for gid in ('water_central_La','aquo9_central_Dy'):
        changed=copy.deepcopy(prep['initial']['model'][gid]);changed['energy']+=1
        try:
            mod.native_authority(gid,changed)
            checks.append({'cell':gid,'kind':'reject_mismatched_native_energy','passed':False})
        except ValueError as e:
            checks.append({'cell':gid,'kind':'reject_mismatched_native_energy','passed':'Producer energy differs' in str(e),'error':str(e)})
        changed=copy.deepcopy(data['rows'][gid]);changed['reference']['forces'][0][0]+=1
        try:
            mod.reference_replay(gid,changed)
            checks.append({'cell':gid,'kind':'reject_mismatched_reference_force','passed':False})
        except Exception as e:
            checks.append({'cell':gid,'kind':'reject_mismatched_reference_force','passed':'Reference arrays differ' in str(e),'error':repr(e)})
    result={'status':'passed' if all(c['passed'] for c in checks) else 'failed',
            'checks':checks,'checks_passed':sum(c['passed'] for c in checks),'checks_total':len(checks),
            'source':pin(source),'control':pin(__file__),'initial_fixture':pin(p),'training':config['training_result'],
            'new_native_calls':0,'live_fit_result_read':False,'numerical_threads':1,
            'wall_seconds':time.monotonic()-start,'CPU_seconds':time.process_time()-cpu}
    out=Path(out);out.parent.mkdir(parents=True,exist_ok=True)
    with out.open('x') as handle:json.dump(result,handle,indent=2,allow_nan=False);handle.write('\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('checks',)}))
    print(json.dumps([c for c in checks if not c['passed']]))
    print(json.dumps(pin(out)))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--out',required=True);args=p.parse_args();main(args.out)
