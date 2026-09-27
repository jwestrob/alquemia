"""Eight-endpoint numerical diagnostic; original failed rigid gate stays failed."""
from __future__ import annotations
import argparse
import json
import os
from pathlib import Path
import shutil
import numpy as np
from affordable_common import (BOHR_TO_A, HA_TO_KCAL, InvalidArtifact, cache_key,
                               read_json, record, verify, write_new)
from affordable_workflow import dry_run, execute
from mace_omol_vacuum import METHOD, REFINED_EMBEDDED_GRID, embedded_input, parse_endpoint
from metal_environment_reference import read_pcgrad
from metal_environment_force_checks import (check_pins, xyz_data, pc_data, rotation,
                                             centers, gradients)

PROTOCOL='nikasha_embedded_grid_diagnostic_v1'
PROFILE=REFINED_EMBEDDED_GRID
TOLERANCES=dict(rigid_energy_Eh=1e-5,rigid_gradient_Eh_bohr=1e-4,response_change_kcal_mol=.05)


def inventory():
    return [(metal,'translation',None) for metal in ('Ca','La')]+[(metal,env,PROFILE) for metal in ('Ca','La') for env in ('A','B','rigid')]


def task_id(metal,env,profile):return metal+'_'+env+'_'+('refined' if profile else 'original')


def rendered_geometry(config,metal,env):
    xl,symbols,x=xyz_data(verify(config['endpoints'][metal]['xyz']))
    pl,q,y=pc_data(verify(config['environments']['B' if env=='B' else 'A']['pointcharges']))
    if env in ('A','B'):return ''.join(xl),''.join(pl)
    if env=='rigid':xx,yy=x@rotation().T+(.173,.117,.231),y@rotation().T+(.173,.117,.231)
    elif env=='translation':xx,yy=x+(.173,.117,.231),y+(.173,.117,.231)
    else:raise InvalidArtifact('unsupported declared geometry')
    for i in range(len(x)):xl[i+2]=symbols[i]+' '+' '.join(f'{v:.16g}' for v in xx[i])+'\n'
    for i in range(len(y)):pl[i+1]=' '.join(f'{v:.17g}' for v in (q[i],*yy[i]))+'\n'
    return ''.join(xl),''.join(pl)


def cache_identity(m,t):
    return cache_key(dict(task={k:v for k,v in t.items() if k!='cache_key'},protocol=PROTOCOL,
                         method=METHOD,orca=m['orca'],tolerances=m['tolerances'],inputs=m['inputs'],
                         source_collection=m['source_collection'],force_collection=m['force_collection'],
                         implementation=m['implementation']))


def reference_sources(source_collection,force_collection,inputs):
    original,cm=centers(source_collection,inputs)
    collection=read_json(source_collection)
    from ggr_sensitivity import executed
    _,actual=executed(verify(collection['manifest']))
    for metal in ('Ca','La'):
        row=collection['rows'][metal+'_B'];task=next(t for t in cm['tasks'] if t['task_id']==metal+'_B')
        if row['status']!='complete' or actual[metal+'_B']['status']!='complete':raise InvalidArtifact('source B unavailable')
        if any(row[k]!=actual[metal+'_B'][k] for k in ('output','receipt','energy_hartree')):raise InvalidArtifact('source B receipt differs')
        parsed=parse_endpoint(task,verify(row['output']),verify(row['engrad']),permanent_field=True)
        if parsed['energy_hartree']!=row['energy_hartree']:raise InvalidArtifact('source B energy differs')
    force=read_json(force_collection)
    fm,fr=executed(verify(force['manifest']))
    if fm['inputs']!=cm['inputs'] or fm['method']!=METHOD:raise InvalidArtifact('force source/method differs')
    if force['status']!='not_qualified':raise InvalidArtifact('diagnostic requires recorded original failed gate')
    for task in fm['tasks']:
        row=force['rows'][task['task_id']];observed=fr[task['task_id']]
        if row['status']!='complete' or observed['status']!='complete' or any(row[k]!=observed[k] for k in ('output','receipt','energy_hartree')):
            raise InvalidArtifact('original force matrix incomplete or receipt differs')
    return collection,force,cm


def prepare(source_manifest,source_collection,force_collection,agreement,output,*,workers,mpi_ranks):
    if workers<1 or workers>8 or mpi_ranks<1:raise InvalidArtifact('one to eight workers and positive ranks required')
    if record(source_manifest)!=read_json(source_collection)['manifest']:raise InvalidArtifact('source manifest differs')
    inputs=verify(read_json(source_manifest)['inputs']);config=read_json(inputs);check_pins(config)
    _,_,cm=reference_sources(source_collection,force_collection,inputs)
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False);impl=out/'implementation';impl.mkdir();pins={}
    for src in sorted(Path(__file__).parent.glob('*.py')):
        shutil.copyfile(src,impl/src.name);pins[src.name]=record(impl/src.name)
    shutil.copyfile(impl/'render_orca_runtime_input.py',impl/'_base_render_orca_runtime_input.py')
    shutil.copyfile(impl/'metal_environment_runtime.py',impl/'render_orca_runtime_input.py')
    for name in ('_base_render_orca_runtime_input.py','render_orca_runtime_input.py'):pins[name]=record(impl/name)
    tasks=[]
    for metal,env,profile in inventory():
        td=out/task_id(metal,env,profile);td.mkdir();x,p=rendered_geometry(config,metal,env)
        (td/'core.xyz').write_text(x);(td/'environment.pc').write_text(p)
        charge=config['endpoints'][metal]['charge'];(td/'endpoint.inp').write_text(embedded_input(charge,grid_profile=profile))
        tasks.append(dict(task_id=td.name,metal=metal,environment=env,grid_profile=profile,charge=charge,multiplicity=1,
                          input=record(td/'endpoint.inp'),xyz=record(td/'core.xyz'),pointcharges=record(td/'environment.pc'),
                          output_path=str(td/'endpoint.out'),engrad_path=str(td/'endpoint.engrad'),task_type='analytic_gradient'))
    m=dict(protocol_id=PROTOCOL,method=METHOD,inputs=record(inputs),source_collection=record(source_collection),
           force_collection=record(force_collection),agreement=record(agreement),tasks=tasks,tolerances=TOLERANCES,
           orca=cm['orca'],implementation=pins,
           execution_policy={'task_runner':pins['run_orca_task_manifest.py'],'runtime_renderer':pins['render_orca_runtime_input.py']},
           execution_resources={'mpi_ranks':mpi_ranks,'concurrent_tasks':workers},
           energy_scope='embedded_electronic_component_only',compute_budget=None,wall_time_limit=None)
    for t in tasks:t['cache_key']=cache_identity(m,t)
    write_new(out/'manifest.json',m)
    return dict(status='prepared',manifest=record(out/'manifest.json'),task_count=8)


def validate(manifest):
    m=read_json(manifest)
    if m['protocol_id']!=PROTOCOL or m['method']!=METHOD or m['tolerances']!=TOLERANCES:raise InvalidArtifact('grid diagnostic protocol/method/tolerances changed')
    check_pins(m['implementation']);verify(m['agreement']);verify(m['orca'])
    config=read_json(verify(m['inputs']));check_pins(config)
    reference_sources(verify(m['source_collection']),verify(m['force_collection']),verify(m['inputs']))
    if len(m['tasks'])!=8 or {(t['metal'],t['environment'],t['grid_profile']) for t in m['tasks']}!=set(inventory()):raise InvalidArtifact('eight-cell grid inventory changed')
    for t in m['tasks']:
        x,p=rendered_geometry(config,t['metal'],t['environment'])
        if verify(t['xyz']).read_text()!=x or verify(t['pointcharges']).read_text()!=p:raise InvalidArtifact('grid diagnostic geometry changed')
        if t['charge']!=config['endpoints'][t['metal']]['charge'] or t['multiplicity']!=1 or verify(t['input']).read_text()!=embedded_input(t['charge'],grid_profile=t['grid_profile']):raise InvalidArtifact('exact named numerical input required')
        if t['cache_key']!=cache_identity(m,t):raise InvalidArtifact('grid diagnostic cache differs')
    return dry_run(manifest)


def rigid_check(metal,mode,base,row,count):
    if row['status']!='complete' or base['status']!='complete':return dict(metal=metal,mode=mode,status='unavailable')
    x,y=gradients(base,count);xx,yy=gradients(row,count);r=rotation() if mode=='refined_rigid' else np.eye(3)
    er=abs(row['energy_hartree']-base['energy_hartree']);qr=float(abs(xx-x@r.T).max()*BOHR_TO_A/HA_TO_KCAL);mr=float(abs(yy-y@r.T).max()*BOHR_TO_A/HA_TO_KCAL)
    ok=er<=TOLERANCES['rigid_energy_Eh'] and max(qr,mr)<=TOLERANCES['rigid_gradient_Eh_bohr']
    return dict(metal=metal,mode=mode,energy_signed_shift_Eh=row['energy_hartree']-base['energy_hartree'],energy_residual_Eh=er,
                core_gradient_residual_Eh_bohr=qr,mm_gradient_residual_Eh_bohr=mr,status='passed' if ok else 'failed')


def collect(manifest):
    from ggr_sensitivity import executed
    validate(manifest);m,rows=executed(manifest);config=read_json(verify(m['inputs']))
    source=read_json(verify(m['source_collection']));force=read_json(verify(m['force_collection']));count=config['environments']['A']['count']
    for t in m['tasks']:
        row=rows[t['task_id']]
        if row['status']!='complete':continue
        try:
            row.update(parse_endpoint(t,verify(row['output']),t['engrad_path'],permanent_field=True,grid_profile=t['grid_profile']))
            pc=Path(t['output_path']).with_name('endpoint.runtime.pcgrad');read_pcgrad(pc,count);row['pointcharge_gradient']=record(pc)
        except (ValueError,OSError) as exc:row.update(status='invalid',reason=str(exc),energy_hartree=None)
    checks=[];response={};changes={}
    for metal in ('Ca','La'):
        a,b,rigid=[rows[task_id(metal,k,PROFILE)] for k in ('A','B','rigid')]
        checks.append(rigid_check(metal,'original_translation',source['rows'][metal+'_A'],rows[task_id(metal,'translation',None)],count))
        checks.append(rigid_check(metal,'refined_rigid',a,rigid,count))
        response[metal]=(b['energy_hartree']-a['energy_hartree'])*HA_TO_KCAL if a['status']==b['status']=='complete' else None
        old=(source['rows'][metal+'_B']['energy_hartree']-source['rows'][metal+'_A']['energy_hartree'])*HA_TO_KCAL
        change=None if response[metal] is None else response[metal]-old;changes[metal]=change
        checks.append(dict(metal=metal,mode='response_grid_change',original_response_kcal_mol=old,refined_response_kcal_mol=response[metal],signed_change_kcal_mol=change,tolerance_kcal_mol=.05,
                           status='unavailable' if change is None else ('passed' if abs(change)<=.05 else 'failed')))
    delta=None if any(v is None for v in response.values()) else response['La']-response['Ca']
    old_delta=source['delta_env_el_kcal_mol'];change=None if delta is None else delta-old_delta
    checks.append(dict(mode='double_difference_grid_change',signed_change_kcal_mol=change,tolerance_kcal_mol=.05,status='unavailable' if change is None else ('passed' if abs(change)<=.05 else 'failed')))
    return dict(protocol_id=PROTOCOL,manifest=record(manifest),rows=rows,checks=checks,original_force_qualification_status=force['status'],
                original_rigid_checks=[c for c in force['checks'] if c['mode']=='rigid'],
                refined_per_metal_response_kcal_mol=response,refined_delta_env_el_kcal_mol=delta,
                refined_grid_status='passed_tested_numerical_checks' if all(c['status']=='passed' for c in checks if c['mode']!='original_translation') else 'not_qualified',
                refined_finite_difference_status='unavailable_not_executed',full_hybrid_status='unsupported',classification=None,
                measured_execution_events=str(Path(manifest).parent/'budget_events.jsonl'))


def main():
    p=argparse.ArgumentParser(description=__doc__);s=p.add_subparsers(dest='op',required=True)
    a=s.add_parser('prepare')
    for k in ('source-manifest','source-collection','force-collection','agreement','output'):a.add_argument('--'+k,required=True)
    for k in ('workers','mpi-ranks'):a.add_argument('--'+k,type=int,required=True)
    for op in ('dry-run','execute','collect'):
        a=s.add_parser(op);a.add_argument('--manifest',required=True);a.add_argument('--output')
    a=p.parse_args()
    if a.op=='prepare':result=prepare(a.source_manifest,a.source_collection,a.force_collection,a.agreement,a.output,workers=a.workers,mpi_ranks=a.mpi_ranks)
    elif a.op=='dry-run':result=validate(a.manifest)
    elif a.op=='collect':result=collect(a.manifest)
    else:
        validate(a.manifest);r=read_json(a.manifest)['execution_resources']
        if int(os.environ.get('SLURM_NTASKS','0'))<r['concurrent_tasks']*r['mpi_ranks']:raise InvalidArtifact('MPI layout exceeds actual slots')
        os.environ['METAL_ENV_WORKERS']=str(r['concurrent_tasks']);result=execute(a.manifest)
    if a.op!='prepare' and a.output:write_new(a.output,result)
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
