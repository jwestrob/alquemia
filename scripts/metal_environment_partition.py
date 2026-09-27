"""Four real expanded-region endpoints: coarse partition sensitivity, not qualification."""
from __future__ import annotations
import argparse,json,os,shutil
from pathlib import Path
import numpy as np
from affordable_common import BOHR_TO_A,HA_TO_KCAL,InvalidArtifact,cache_key,read_json,record,verify,write_new
from affordable_workflow import dry_run,execute
from mace_omol_vacuum import METHOD,embedded_input,parse_endpoint
from metal_environment_reference import read_pcgrad
from metal_environment_force_checks import check_pins,xyz_data,pc_data
PROTOCOL='nikasha_expanded_Thr159_electronic_partition_diagnostic_v1'
TASK_IDS=('Ca_A','Ca_B','La_A','La_B')


def validate_inputs(config):
    check_pins(config)
    if set(config['endpoints'])!=set(TASK_IDS):raise InvalidArtifact('four expanded endpoint identities required')
    if config.get('complete_hybrid_status')!='unsupported_missing_cross_repulsion_dispersion_and_boundary_mechanics':raise InvalidArtifact('hybrid limitations changed')
    geometry={}
    for tid in TASK_IDS:
        t=config['endpoints'][tid];metal,env=tid.split('_');_,symbols,coords=xyz_data(verify(t['xyz']))
        if t['metal']!=metal or t['environment']!=env or t['charge']!=(-3 if metal=='Ca' else -2) or t['multiplicity']!=1:raise InvalidArtifact('expanded electronic state differs')
        if len(coords)!=63 or symbols[0]!=metal or t['pointcharges']!=config['environments'][env]['pointcharges']:raise InvalidArtifact('expanded geometry/field pairing differs')
        numbers={'H':1,'C':6,'N':7,'O':8,'S':16,'Ca':20,'La':57}
        physical=sum(numbers[e] for e in symbols)-t['charge'];ecp=46 if metal=='La' else 0
        if t['physical_electrons']!=physical or t['native_r2scan3c_ecp_core_electrons']!=ecp or t['native_r2scan3c_explicit_electrons']!=physical-ecp or (physical-ecp)%2:raise InvalidArtifact('expanded electron/ECP metadata differs')
        mapping=read_json(verify(t['core_mapping']))
        if len(mapping)!=63 or not np.array_equal(coords,np.array([r['xyz_A'] for r in mapping])):raise InvalidArtifact('expanded source mapping coordinates differ')
        if symbols[1:]!=[r['element'] for r in mapping[1:]]:raise InvalidArtifact('expanded source mapping elements differ')
        geometry[tid]=coords
    for env in ('A','B'):
        if not np.array_equal(geometry['Ca_'+env],geometry['La_'+env]):raise InvalidArtifact('paired expanded coordinates differ')
    for metal in ('Ca','La'):
        changed=np.where(np.any(geometry[metal+'_A']!=geometry[metal+'_B'],axis=1))[0].tolist()
        if changed!=[57]:raise InvalidArtifact('B must use its real changed QM HG1 coordinates')
    if config['environments']['A']['pointcharges']['sha256']!=config['environments']['B']['pointcharges']['sha256']:raise InvalidArtifact('expanded external A/B fields differ')
    for e in config['environments'].values():
        _,q,y=pc_data(verify(e['pointcharges']))
        if len(q)!=9078 or abs(float(q.sum())+6)>1e-10:raise InvalidArtifact('expanded external charge closure differs')
    return geometry


def responses(rows):
    per={}
    for metal in ('Ca','La'):
        a,b=[rows.get(metal+'_'+env,{}) for env in ('A','B')]
        per[metal]=(b['energy_hartree']-a['energy_hartree'])*HA_TO_KCAL if a.get('status')==b.get('status')=='complete' else None
    delta=per['La']-per['Ca'] if all(v is not None for v in per.values()) else None
    return per,delta


def source_reference(path,config):
    from ggr_sensitivity import executed
    source=read_json(path);sm,actual=executed(verify(source['manifest']))
    if sm['inputs']!=config['base_inputs'] or sm['method']!=METHOD:raise InvalidArtifact('small-region source/preparation/method differs')
    for tid in TASK_IDS:
        row=source['rows'][tid];observed=actual[tid];task=next(t for t in sm['tasks'] if t['task_id']==tid)
        if row['status']!='complete' or observed['status']!='complete':raise InvalidArtifact('small-region reference incomplete')
        if any(row[k]!=observed[k] for k in ('output','receipt','energy_hartree')):raise InvalidArtifact('small-region reference receipt differs')
        parsed=parse_endpoint(task,verify(row['output']),verify(row['engrad']),permanent_field=True)
        if parsed['energy_hartree']!=row['energy_hartree']:raise InvalidArtifact('small-region parsed energy differs')
    per,delta=responses(source['rows'])
    if abs(delta-source['delta_env_el_kcal_mol'])>1e-10:raise InvalidArtifact('small-region response algebra differs')
    return source,per,delta


def identity(m,t):return cache_key({'task':{k:v for k,v in t.items() if k!='cache_key'},'protocol':PROTOCOL,'method':METHOD,'inputs':m['inputs'],'source_collection':m['source_collection'],'orca':m['orca'],'implementation':m['implementation']})


def prepare(inputs,source_collection,agreement,output,*,workers,mpi_ranks):
    if workers!=4 or mpi_ranks!=86:raise InvalidArtifact('declared layout is four workers times 86 MPI ranks')
    config=read_json(inputs);validate_inputs(config);source_reference(source_collection,config)
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False);impl=out/'implementation';impl.mkdir();pins={}
    for src in sorted(Path(__file__).parent.glob('*.py')):shutil.copyfile(src,impl/src.name);pins[src.name]=record(impl/src.name)
    shutil.copyfile(impl/'render_orca_runtime_input.py',impl/'_base_render_orca_runtime_input.py');shutil.copyfile(impl/'metal_environment_runtime.py',impl/'render_orca_runtime_input.py')
    for name in ('_base_render_orca_runtime_input.py','render_orca_runtime_input.py'):pins[name]=record(impl/name)
    tasks=[]
    for tid in TASK_IDS:
        endpoint=config['endpoints'][tid];td=out/tid;td.mkdir()
        shutil.copyfile(verify(endpoint['xyz']),td/'core.xyz');shutil.copyfile(verify(endpoint['pointcharges']),td/'environment.pc');(td/'endpoint.inp').write_text(embedded_input(endpoint['charge']))
        tasks.append(dict(task_id=tid,metal=endpoint['metal'],environment=endpoint['environment'],charge=endpoint['charge'],multiplicity=1,input=record(td/'endpoint.inp'),xyz=record(td/'core.xyz'),pointcharges=record(td/'environment.pc'),core_mapping=endpoint['core_mapping'],output_path=str(td/'endpoint.out'),engrad_path=str(td/'endpoint.engrad'),task_type='analytic_gradient'))
    m=dict(protocol_id=PROTOCOL,method=METHOD,inputs=record(inputs),source_collection=record(source_collection),agreement=record(agreement),tasks=tasks,orca=config['orca'],implementation=pins,execution_policy={'task_runner':pins['run_orca_task_manifest.py'],'runtime_renderer':pins['render_orca_runtime_input.py']},execution_resources={'mpi_ranks':mpi_ranks,'concurrent_tasks':workers},energy_scope='embedded_electronic_component_only',compute_budget=None,wall_time_limit=None,numerical_qualification_status='not_established_original_rigid_gate_failed',full_hybrid_status='unsupported',classification=None)
    for t in tasks:t['cache_key']=identity(m,t)
    write_new(out/'manifest.json',m);return {'status':'prepared','manifest':record(out/'manifest.json'),'task_count':4}


def validate(manifest):
    m=read_json(manifest)
    if m['protocol_id']!=PROTOCOL or m['method']!=METHOD:raise InvalidArtifact('partition protocol/method differs')
    check_pins(m['implementation']);verify(m['agreement']);verify(m['orca']);config=read_json(verify(m['inputs']));validate_inputs(config);source_reference(verify(m['source_collection']),config)
    if len(m['tasks'])!=4 or {t['task_id'] for t in m['tasks']}!=set(TASK_IDS):raise InvalidArtifact('four-cell inventory differs')
    if m['execution_resources']!={'mpi_ranks':86,'concurrent_tasks':4}:raise InvalidArtifact('four-by-86 MPI layout differs')
    for t in m['tasks']:
        e=config['endpoints'][t['task_id']]
        if verify(t['xyz']).read_bytes()!=verify(e['xyz']).read_bytes() or verify(t['pointcharges']).read_bytes()!=verify(e['pointcharges']).read_bytes():raise InvalidArtifact('exact endpoint-specific expanded coordinates/field required')
        if any(t[k]!=e[k] for k in ('metal','environment','charge','multiplicity','core_mapping')):raise InvalidArtifact('endpoint state/source map differs')
        if verify(t['input']).read_text()!=embedded_input(e['charge']):raise InvalidArtifact('original DefGrid3 TightSCF EnGrad profile required')
        if t['cache_key']!=identity(m,t):raise InvalidArtifact('partition cache identity differs')
    return dry_run(manifest)


def hydroxyl_projection(config,tid,gradient):
    mapping=read_json(verify(config['endpoints'][tid]['core_mapping']))
    byid={v['id']:v for v in mapping};carbon,oxygen,hydrogen=[np.array(byid['A/159/ /'+name]['xyz_A']) for name in ('CB','OG1','HG1')]
    axis=oxygen-carbon;axis/=np.linalg.norm(axis);tangent=np.cross(axis,hydrogen-oxygen);radius=float(np.linalg.norm(tangent))
    if radius<=1e-12:raise InvalidArtifact('degenerate physical HG1 rotation')
    arc=float(np.dot(gradient,tangent/radius))
    return dict(gradient_per_arc_kcal_mol_A=arc,force_per_arc_kcal_mol_A=-arc,gradient_per_radian_kcal_mol=float(np.dot(gradient,tangent)),arc_radius_A=radius,quantity='electronic-component derivative only; different QM/MM representations')


def collect(manifest):
    from ggr_sensitivity import executed
    validate(manifest);m,rows=executed(manifest);config=read_json(verify(m['inputs']));source,oldper,olddelta=source_reference(verify(m['source_collection']),config)
    for t in m['tasks']:
        row=rows[t['task_id']]
        if row['status']!='complete':continue
        try:
            row.update(parse_endpoint(t,verify(row['output']),t['engrad_path'],permanent_field=True))
            pc=Path(t['output_path']).with_name('endpoint.runtime.pcgrad');read_pcgrad(pc,9078);row['pointcharge_gradient']=record(pc)
        except (ValueError,OSError) as exc:row.update(status='invalid',reason=str(exc),energy_hartree=None)
    per,delta=responses(rows)
    projections={}
    for tid in TASK_IDS:
        small=source['rows'][tid];pcg=read_pcgrad(verify(small['pointcharge_gradient']),9087)*HA_TO_KCAL/BOHR_TO_A
        projections[tid]={'small_region_MM_HG1':hydroxyl_projection(config,tid,pcg[2412]),'expanded_region_QM_HG1':None}
        if rows[tid]['status']=='complete':projections[tid]['expanded_region_QM_HG1']=hydroxyl_projection(config,tid,np.array(rows[tid]['gradient_kcal_mol_per_A'])[57])
    return dict(protocol_id=PROTOCOL,manifest=record(manifest),rows=rows,status='component_matrix_complete' if delta is not None else 'component_matrix_unavailable',per_metal_response_kcal_mol=per,delta_env_el_kcal_mol=delta,small_region_per_metal_response_kcal_mol=oldper,small_region_delta_env_el_kcal_mol=olddelta,HG1_physical_arc_projections=projections,partition_response_change_kcal_mol=None if delta is None else delta-olddelta,numerical_qualification_status='not_established_original_rigid_gate_failed',interpretation='Coarse electronic-component partition sensitivity only; small differences remain numerically unresolved. No biological pass/fail or absolute cross-region comparison.',full_hybrid_status='unsupported',classification=None,measured_execution_events=str(Path(manifest).parent/'budget_events.jsonl'))


def main():
    p=argparse.ArgumentParser(description=__doc__);s=p.add_subparsers(dest='op',required=True);a=s.add_parser('prepare')
    for k in ('inputs','source-collection','agreement','output'):a.add_argument('--'+k,required=True)
    for k in ('workers','mpi-ranks'):a.add_argument('--'+k,type=int,required=True)
    for op in ('dry-run','execute','collect'):
        a=s.add_parser(op);a.add_argument('--manifest',required=True);a.add_argument('--output')
    a=p.parse_args()
    if a.op=='prepare':r=prepare(a.inputs,a.source_collection,a.agreement,a.output,workers=a.workers,mpi_ranks=a.mpi_ranks)
    elif a.op=='dry-run':r=validate(a.manifest)
    elif a.op=='collect':r=collect(a.manifest)
    else:
        validate(a.manifest);m=read_json(a.manifest);res=m['execution_resources']
        if int(os.environ.get('SLURM_NTASKS','0'))<res['concurrent_tasks']*res['mpi_ranks']:raise InvalidArtifact('MPI layout exceeds allocation slots')
        os.environ['METAL_ENV_WORKERS']=str(res['concurrent_tasks']);r=execute(a.manifest)
    if a.op!='prepare' and a.output:write_new(a.output,r)
    print(json.dumps(r,indent=2))
if __name__=='__main__':main()
