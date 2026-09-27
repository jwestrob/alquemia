"""Four fixed-grid, stricter-SCF endpoints; no relaxation or retry expansion."""
from __future__ import annotations
import argparse
import json
import os
from pathlib import Path
import shutil
from affordable_common import InvalidArtifact,cache_key,read_json,record,verify,write_new
from affordable_workflow import dry_run,execute
from mace_omol_vacuum import (METHOD,REFINED_EMBEDDED_GRID,STRICT_EMBEDDED_SCF,
                              embedded_input,parse_endpoint,terminal_scf_residuals)
from metal_environment_reference import read_pcgrad
from metal_environment_force_checks import check_pins
from metal_environment_grid_check import rigid_check,TOLERANCES

PROTOCOL='nikasha_embedded_fixed_grid_strict_scf_v1'
GRID=REFINED_EMBEDDED_GRID
SCF=STRICT_EMBEDDED_SCF


def inventory():return [(metal,env) for metal in ('Ca','La') for env in ('A','rigid')]


def source_cells(source_collection):
    from ggr_sensitivity import executed
    c=read_json(source_collection);m,actual=executed(verify(c['manifest']));selected={}
    for metal,env in inventory():
        key=f'{metal}_{env}_refined';row=c['rows'][key]
        if row['status']!='complete' or actual[key]['status']!='complete' or any(row[k]!=actual[key][k] for k in ('output','receipt','energy_hartree')):
            raise InvalidArtifact('completed refined source receipt unavailable/differs')
        task=next(t for t in m['tasks'] if t['task_id']==key)
        parsed=parse_endpoint(task,verify(row['output']),verify(row['engrad']),permanent_field=True,grid_profile=GRID)
        if parsed['energy_hartree']!=row['energy_hartree']:raise InvalidArtifact('refined source energy differs')
        selected[(metal,env)]=task
    return c,m,selected


def identity(m,t):
    return cache_key(dict(task={k:v for k,v in t.items() if k!='cache_key'},protocol=PROTOCOL,method=METHOD,
                         grid_profile=GRID,scf_profile=SCF,source_collection=m['source_collection'],
                         orca=m['orca'],implementation=m['implementation'],tolerances=m['tolerances']))


def prepare(source_collection,agreement,output,*,workers,mpi_ranks):
    if not 1<=workers<=4 or mpi_ranks<1:raise InvalidArtifact('one to four workers and positive ranks required')
    _,sm,source=source_cells(source_collection)
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False);impl=out/'implementation';impl.mkdir();pins={}
    for src in sorted(Path(__file__).parent.glob('*.py')):
        shutil.copyfile(src,impl/src.name);pins[src.name]=record(impl/src.name)
    shutil.copyfile(impl/'render_orca_runtime_input.py',impl/'_base_render_orca_runtime_input.py')
    shutil.copyfile(impl/'metal_environment_runtime.py',impl/'render_orca_runtime_input.py')
    for name in ('_base_render_orca_runtime_input.py','render_orca_runtime_input.py'):pins[name]=record(impl/name)
    tasks=[]
    for metal,env in inventory():
        src=source[(metal,env)];td=out/f'{metal}_{env}_strict';td.mkdir()
        for key,name in [('xyz','core.xyz'),('pointcharges','environment.pc')]:shutil.copyfile(verify(src[key]),td/name)
        (td/'endpoint.inp').write_text(embedded_input(src['charge'],grid_profile=GRID,scf_profile=SCF))
        tasks.append(dict(task_id=td.name,metal=metal,environment=env,grid_profile=GRID,scf_profile=SCF,
                          charge=src['charge'],multiplicity=1,input=record(td/'endpoint.inp'),xyz=record(td/'core.xyz'),
                          pointcharges=record(td/'environment.pc'),output_path=str(td/'endpoint.out'),
                          engrad_path=str(td/'endpoint.engrad'),task_type='analytic_gradient'))
    m=dict(protocol_id=PROTOCOL,method=METHOD,grid_profile=GRID,scf_profile=SCF,source_collection=record(source_collection),
           agreement=record(agreement),tasks=tasks,orca=sm['orca'],implementation=pins,tolerances=TOLERANCES,
           execution_policy={'task_runner':pins['run_orca_task_manifest.py'],'runtime_renderer':pins['render_orca_runtime_input.py']},
           execution_resources={'mpi_ranks':mpi_ranks,'concurrent_tasks':workers},energy_scope='embedded_electronic_component_only',
           compute_budget=None,wall_time_limit=None)
    for t in tasks:t['cache_key']=identity(m,t)
    write_new(out/'manifest.json',m)
    return dict(status='prepared',manifest=record(out/'manifest.json'),task_count=4)


def validate(manifest):
    m=read_json(manifest)
    if any(m[k]!=v for k,v in dict(protocol_id=PROTOCOL,method=METHOD,grid_profile=GRID,scf_profile=SCF,tolerances=TOLERANCES).items()):raise InvalidArtifact('strict SCF protocol changed')
    check_pins(m['implementation']);verify(m['agreement']);verify(m['orca'])
    _,_,source=source_cells(verify(m['source_collection']))
    if len(m['tasks'])!=4 or {(t['metal'],t['environment']) for t in m['tasks']}!=set(inventory()):raise InvalidArtifact('four-cell inventory changed')
    for t in m['tasks']:
        src=source[(t['metal'],t['environment'])]
        for key in ('xyz','pointcharges'):
            verify(t[key])
            if t[key]['sha256']!=src[key]['sha256']:raise InvalidArtifact('strict SCF source geometry/field changed')
        if t['charge']!=src['charge'] or t['multiplicity']!=1 or t['grid_profile']!=GRID or t['scf_profile']!=SCF or verify(t['input']).read_text()!=embedded_input(src['charge'],grid_profile=GRID,scf_profile=SCF):raise InvalidArtifact('strict SCF exact state/input changed')
        if t['cache_key']!=identity(m,t):raise InvalidArtifact('strict SCF cache changed')
    return dry_run(manifest)


def collect(manifest):
    from ggr_sensitivity import executed
    validate(manifest);m,rows=executed(manifest);old=read_json(verify(m['source_collection']));checks=[]
    for t in m['tasks']:
        row=rows[t['task_id']];output=Path(t['output_path'])
        # Preserve printed residuals even if SCF status/strict admission fails.
        if output.exists():
            try:row['observed_terminal_scf_residuals']=terminal_scf_residuals(output.read_text())
            except (ValueError,OSError) as exc:row['scf_residual_observation_status']=str(exc)
        if row['status']!='complete':continue
        try:
            row.update(parse_endpoint(t,verify(row['output']),t['engrad_path'],permanent_field=True,grid_profile=GRID,scf_profile=SCF))
            count=int(verify(t['pointcharges']).read_text().splitlines()[0]);pc=output.with_name('endpoint.runtime.pcgrad');read_pcgrad(pc,count);row['pointcharge_gradient']=record(pc)
            row['energy_change_from_looser_scf_Eh']=row['energy_hartree']-old['rows'][f"{t['metal']}_{t['environment']}_refined"]['energy_hartree']
        except (ValueError,OSError) as exc:row.update(status='invalid',reason=str(exc),energy_hartree=None)
    count=int(verify(m['tasks'][0]['pointcharges']).read_text().splitlines()[0])
    for metal in ('Ca','La'):checks.append(rigid_check(metal,'refined_rigid',rows[f'{metal}_A_strict'],rows[f'{metal}_rigid_strict'],count))
    return dict(protocol_id=PROTOCOL,manifest=record(manifest),rows=rows,checks=checks,
                status='passed_strict_scf_and_rigid_checks' if all(c['status']=='passed' for c in checks) else 'not_qualified',
                original_grid_status=old['refined_grid_status'],original_rigid_checks=[c for c in old['checks'] if c['mode']=='refined_rigid'],
                strict_scf_finite_difference_status='not_executed',strict_scf_environment_response_status='not_executed_B_absent',
                full_hybrid_status='unsupported',classification=None,measured_execution_events=str(Path(manifest).parent/'budget_events.jsonl'))


def main():
    p=argparse.ArgumentParser(description=__doc__);s=p.add_subparsers(dest='op',required=True)
    a=s.add_parser('prepare')
    for k in ('source-collection','agreement','output'):a.add_argument('--'+k,required=True)
    for k in ('workers','mpi-ranks'):a.add_argument('--'+k,type=int,required=True)
    for op in ('dry-run','execute','collect'):
        a=s.add_parser(op);a.add_argument('--manifest',required=True);a.add_argument('--output')
    a=p.parse_args()
    if a.op=='prepare':result=prepare(a.source_collection,a.agreement,a.output,workers=a.workers,mpi_ranks=a.mpi_ranks)
    elif a.op=='dry-run':result=validate(a.manifest)
    elif a.op=='collect':result=collect(a.manifest)
    else:
        validate(a.manifest);r=read_json(a.manifest)['execution_resources']
        if int(os.environ.get('SLURM_NTASKS','0'))<r['concurrent_tasks']*r['mpi_ranks']:raise InvalidArtifact('MPI layout exceeds actual slots')
        os.environ['METAL_ENV_WORKERS']=str(r['concurrent_tasks']);result=execute(a.manifest)
    if a.op!='prepare' and a.output:write_new(a.output,result)
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
