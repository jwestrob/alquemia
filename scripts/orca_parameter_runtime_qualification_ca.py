#!/usr/bin/env python3
"""Declared two-cell full corrected-file calcium import control; no retries."""
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import argparse, copy, ctypes, json, os, resource, shutil, time
import orca_parameter_runtime_qualification as base

W=base.W

def prepare():
    old=base.read(W/'manifest_La8.json')
    parameter=W/'inputs/La_vacuum_corrected/control.xtb.json'
    if base.sha(parameter)!='a4f83ca72d58d6e80f634e61bb8e99b7af0f26c38cc4a9295e85215ecf72d05c':
        raise ValueError('corrected full-file identity changed')
    plan={'created_utc':base.now(),'question':'Does importing the exact corrected full La parameter file preserve both unchanged Ca control endpoints?',
          'authorization':'Root instruction: add a declared two-cell Ca corrected-full-file vacuum/water diagnostic; increase finite maximum58 to60; <=4 aggregate numerical threads.',
          'original_plan':base.pin(W/'PLAN.json'),'previous_maximum_native_calls':58,'new_maximum_native_calls':60,
          'new_jobs':['Ca_vacuum_corrected_file','Ca_water_corrected_file'],'settings':'Identical existing 222-atom Ca geometry/state and default input except explicit full corrected parameter import; no retry or tolerance change.',
          'admission':'All original process/state/printed SCF/property gates, <=1e-8 Eh and <=1e-7 e max-atom difference against accepted same-runtime default and archived default.',
          'parameter':base.pin(parameter)}
    base.write(W/'PLAN_EXTENSION_Ca2.json',plan)
    jobs=[]
    for medium in ['vacuum','water']:
        original=next(j for j in old['jobs'] if j['id']==f'Ca_{medium}_default')
        default=base.read(W/'results'/original['id']/'receipt.json')
        if not default['passed'] or not default['comparison_to_archived']['accepted_parity']:
            raise ValueError('accepted same-runtime Ca default missing')
        j=copy.deepcopy(original);j['id']=f'Ca_{medium}_corrected_file'
        j['extension']='Ca2_corrected_full_file';j['variant']='corrected_file'
        j['baseline_job_id']=original['id'];j['input_dir']='inputs_Ca2/'+j['id']
        d=W/j['input_dir'];d.mkdir(parents=True,exist_ok=False)
        shutil.copyfile(W/original['input_dir']/'core.xyz',d/'core.xyz')
        source=(W/original['input_dir']/'endpoint.inp').read_text()
        prior=' WriteXTBParam true\n ReadXTBParam false'
        if source.count(prior)!=1:raise ValueError('unexpected source method block')
        (d/'endpoint.inp').write_text(source.replace(prior,' XTBParamFile "control.xtb.json"\n ReadXTBParam true\n WriteXTBParam false'))
        shutil.copyfile(parameter,d/'control.xtb.json')
        j['files']={p.name:base.pin(p) for p in d.iterdir()}
        jobs.append(j)
    base.write(W/'manifest_Ca2.json',{'schema':'native_runtime_Ca_full_import.v1','jobs':jobs,'plan':base.pin(W/'PLAN_EXTENSION_Ca2.json'),'helper':base.pin(__file__),'base_helper':base.pin(base.__file__),'runtime_binding':base.pin(W/'RUNTIME_BINDING.json')})
    print(json.dumps({'prepared':2,'new_total_maximum':60}),flush=True)

def run():
    base.verify_binding();libc=ctypes.CDLL(None,use_errno=True)
    if libc.prctl(36,1,0,0,0)!=0:raise OSError('subreaper unavailable')
    mp=W/'manifest_Ca2.json';m=base.read(mp);usage=resource.getrusage(resource.RUSAGE_CHILDREN)
    start=time.monotonic()
    base.status('running',[],stage='Ca2_corrected_full_file',runner_pid=os.getpid(),requested_jobs=[j['id'] for j in m['jobs']])
    with ThreadPoolExecutor(max_workers=2) as ex:
        rows=list(ex.map(lambda j:base.run_job(j,mp),m['jobs']))
    comparisons=[]
    for r in rows:
        default=base.read(W/'results'/r['job']['baseline_job_id']/'receipt.json')
        comparisons.append({'job':r['job']['id'],'receipt':base.pin(W/'results'/r['job']['id']/'receipt.json'),'against_same_runtime_default':base.compare(r['output'],default['output'],1e-8,1e-7),'against_archived_default':r['comparison_to_archived']})
    end=resource.getrusage(resource.RUSAGE_CHILDREN)
    result={'stage':'Ca2_corrected_full_file','completed':len(rows),'accepted':sum(r['passed'] for r in rows),'comparisons':comparisons,
            'all_pass':all(x[k]['accepted_parity'] for x in comparisons for k in ['against_same_runtime_default','against_archived_default']),
            'elapsed_wall_seconds':time.monotonic()-start,'summed_endpoint_wall_seconds':sum(r['process']['elapsed_seconds'] for r in rows),
            'reaped_child_user_seconds':end.ru_utime-usage.ru_utime,'reaped_child_system_seconds':end.ru_stime-usage.ru_stime}
    base.write(W/'EXECUTION_Ca2.json',result)
    base.status('stage_complete',[],stage='Ca2_corrected_full_file',completed=2,accepted=result['accepted'])
    print(json.dumps(result),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('op',choices=['prepare','run']);a=p.parse_args()
    prepare() if a.op=='prepare' else run()
