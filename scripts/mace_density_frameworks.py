"""Declared source-backed AMOEBA/GK framework preparation without energies."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import shutil
import time
import numpy as np
from affordable_common import InvalidArtifact,cache_key,read_json,record,verify,write_new
from mace_density_inputs import audit
from mace_amoeba_capability import prepare_case
from mace_tinker_capability import prepare as export_native
from mace_tinker_framework_solver import native_call,parse_parameters,check_parameters

PROTOCOL='source_graph_AMOEBA2018_GK_framework_preparation_v1'


def key(t,m):
    return cache_key(dict(task={k:v for k,v in t.items() if k!='cache_key'},
        **{k:m[k] for k in ('protocol','preparation','reference_solver','software','native_parameters','implementation')}))


def prepare(preparation,reference_solver,output):
    start=time.monotonic();cpu=time.process_time();audit(preparation);prep=read_json(preparation)
    ref=read_json(reference_solver);old=read_json(verify(ref['parent']));ar=read_json(verify(old['parent']))
    am=read_json(verify(ar['manifest']));ff=[verify(p) for p in am['forcefields']]
    if len(ff)!=2:raise InvalidArtifact('qualified paired AMOEBA2018 forcefields required')
    sw=read_json(verify(ref['software']))
    if sw['returncode']:raise InvalidArtifact('qualified native framework executable unavailable')
    root=Path(output).resolve();root.mkdir(parents=True,exist_ok=False);impl=root/'implementation';impl.mkdir()
    paths={k:verify(v) for k,v in prep['implementation'].items()}
    for n in ('mace_density_inputs.py','mace_density_frameworks.py','mace_amoeba_capability.py','mace_tinker_capability.py','mace_tinker_framework_solver.py','mace_responsive_charges.py'):
        paths[n]=Path(__file__).with_name(n)
    for n,p in paths.items():shutil.copyfile(p,impl/n)
    pins={p.name:record(p) for p in impl.glob('*.py')}
    a=root/'amoeba';a.mkdir();names=tuple(prep['cases']);states={};rows=[]
    write_new(a/'manifest.json',dict(protocol=PROTOCOL,source_preparation=record(preparation),
        reference=record(verify(ar['manifest'])),forcefields=am['forcefields'],implementation=pins))
    for name,pin in prep['cases'].items():
        case=read_json(verify(pin));states[name]=case['state']
        row=prepare_case(case['normalized_global_preparation'],ff,a/name,names);rows.append(row)
        if row['status']!='framework_parameterized_no_energy':
            write_new(root/'failure.json',dict(case_id=name,result=record(a/name/'result.json'),energy_calls=0))
            raise InvalidArtifact(row['failure_reason'])
    write_new(a/'result.json',dict(protocol=PROTOCOL,manifest=record(a/'manifest.json'),cases=rows,
        numerical_score=None,new_energy_calls=0,new_force_calls=0))
    native=root/'native_export'
    ex=export_native(a/'result.json',verify(old['parameters']),verify(old['build_receipt']),None,
        verify(prep['agreement']),native,case_states=states)
    m=dict(protocol=PROTOCOL,preparation=record(preparation),reference_solver=record(reference_solver),
        parent=record(native/'manifest.json'),software=ref['software'],executable=sw['executable'],
        native_parameters=old['parameters'],plan=prep['agreement'],implementation=pins,tasks=[],
        new_energy_calls=0,new_response_solves=0,numerical_score=None,baseline_changed=False)
    for source in ex['tasks']:
        name=source['case_id'];d=root/'tasks'/name;d.mkdir(parents=True)
        shutil.copyfile(verify(source['xyz']),d/'framework.xyz');(d/'framework.freeze').write_text('0\n')
        opts=['parameters '+str(verify(old['parameters'])),'multipoleterm ONLY','polarizeterm','solvateterm',
              'solvate GK','dielectric 1.0','gk-radius SOLUTE','gkc 2.455','polarization MUTUAL','polar-iter 100','polar-eps 1e-5']
        (d/'framework.key').write_text('\n'.join(opts)+'\n')
        task=dict(task_id=name+'_all_standard',case_id=name,variant='all_standard',poleps=1e-5,
            rotation=np.eye(3).tolist(),translation_A=[0.,0.,0.],frozen_indices=[],mapping=source['mapping'],
            xyz=record(d/'framework.xyz'),mask=record(d/'framework.freeze'),key=record(d/'framework.key'))
        receipt=native_call(sw['executable'],task,'preflight',d/'preflight.log',1)
        write_new(d/'preflight_receipt.json',receipt);text=verify(receipt['log']).read_text()
        if receipt['returncode'] or 'ALQUEMIA_PREFLIGHT_COMPLETE' not in text:raise InvalidArtifact('native parameter read failed')
        params=parse_parameters(text);check_parameters(params,task,read_json(verify(task['mapping'])))
        write_new(d/'parameters.json',params);task.update(parameters=record(d/'parameters.json'),preflight_receipt=record(d/'preflight_receipt.json'))
        task['cache_key']=key(task,m);m['tasks'].append(task)
    m['preparation_receipt']=dict(wall_seconds=time.monotonic()-start,CPU_seconds=time.process_time()-cpu)
    write_new(root/'manifest.json',m);validate(root/'manifest.json');return m


def validate(path):
    m=read_json(path)
    if m['protocol'] not in (PROTOCOL,'source_multisite_AMOEBA2018_GK_verified_framework_import_v1') or m['new_energy_calls'] or m['new_response_solves']:raise InvalidArtifact('framework preparation scope differs')
    for p in [m['preparation'],m['reference_solver'],m['parent'],m['software'],m['executable'],m['native_parameters'],m['plan'],*m['implementation'].values()]:verify(p)
    prep=read_json(verify(m['preparation']))
    if len(m['tasks'])!=len(prep['cases']) or {t['case_id'] for t in m['tasks']}!=set(prep['cases']):raise InvalidArtifact('framework inventory differs')
    for t in m['tasks']:
        for k in ('mapping','xyz','mask','key','parameters','preflight_receipt'):verify(t[k])
        if t['cache_key']!=key(t,m):raise InvalidArtifact('framework cache differs')
        receipt=read_json(verify(t['preflight_receipt']));text=verify(receipt['log']).read_text()
        if receipt['returncode'] or 'ALQUEMIA_PREFLIGHT_COMPLETE' not in text:raise InvalidArtifact('native parameter receipt failed')
        params=parse_parameters(text)
        if params!=read_json(verify(t['parameters'])):raise InvalidArtifact('native parameters differ from output')
        check_parameters(params,t,read_json(verify(t['mapping'])))
    return m


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest='command',required=True)
    s=sub.add_parser('prepare')
    for name in ('preparation','reference-solver','output'):s.add_argument('--'+name,required=True)
    s=sub.add_parser('audit');s.add_argument('--manifest',required=True)
    a=p.parse_args()
    if a.command=='prepare':r=prepare(a.preparation,a.reference_solver,a.output)
    else:r=validate(a.manifest)
    print(json.dumps(dict(protocol=r['protocol'],cases=[t['case_id'] for t in r['tasks']],new_energy_calls=0,new_response_solves=0)))
