"""Finite opt-in development operations using the existing manifested ORCA runner."""
from __future__ import annotations

import argparse
import fcntl
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import time

from affordable_common import InvalidArtifact, read_json, record, verify, write_new, digest, energy, xyz, paired, snapshot_implementation


def resolve_record(rec, parent):
    p=Path(rec['path'])
    if not p.is_absolute(): p=Path(parent)/p
    return verify(dict(rec,path=str(p)))


def prepare_pilot(root, output):
    root,output=Path(root).resolve(),Path(output).resolve()
    if output.exists(): raise InvalidArtifact(f'refusing existing pilot {output}')
    base=root/'diagnostics/pqq_pmdh_fixed_core_calibration_20260914'
    embedding=read_json(root/'diagnostics/pqq_balanced_embedding_20260914/preparation.json')
    sources=[('1h4i_fixed_core',base/'reserved_crystal_holdout/prepared/01_1H4I/pmdh_fc_holdout_1h4i_carve_manifest.json')]
    for radius in ('3.3','3.6'):
        sources.append(('1h4i_qm'+radius.replace('.',''),verify(embedding['environments']['mxaf'][radius]['anchor_manifest'])))
    output.mkdir(parents=True)
    tasks,cases=[],[]
    for name,manifest_path in sources:
        m=read_json(manifest_path); task_ids=[]; charges={}
        for metal in ('La','Ca'):
            nested=m['outputs'].get(metal)
            ir=nested.get('orca_input',nested.get('input')) if nested else m['outputs'][metal+'_input']
            xr=nested['xyz'] if nested else m['outputs'][metal+'_xyz']
            source_input=resolve_record(ir,manifest_path.parent); source_xyz=resolve_record(xr,manifest_path.parent)
            text=source_input.read_text()
            simple=[l.strip() for l in text.splitlines() if l.strip().startswith('!')]
            if len(simple)!=1 or set(simple[0].lower().split()[1:])!={'r2scan-3c','noautostart','cpcm(water)','defgrid3'}:
                raise InvalidArtifact(f'unsupported endpoint Hamiltonian: {source_input}')
            if re.search(r'%pointcharges|%basis|%pal|\bOpt\b|\bNumGrad\b',text,re.I):
                raise InvalidArtifact('cannot alter Hamiltonian or resource policy')
            coord=re.search(r'^\s*\*\s+xyzfile\s+(-?\d+)\s+(\d+)\s+\S+\s*$',text,re.M|re.I)
            if not coord or coord[2]!='1': raise InvalidArtifact('missing closed-shell coordinate declaration')
            charges[metal]=int(coord[1])
            d=output/name/metal; d.mkdir(parents=True)
            xp=d/'core.xyz'; xp.write_bytes(source_xyz.read_bytes())
            inp=d/'endpoint.inp'
            inp.write_text(f'! r2SCAN-3c NoAutostart CPCM(Water) DefGrid3 MBIS\n%method\n MBIS_LARGEPRINT true\nend\n* xyzfile {coord[1]} 1 core.xyz\n')
            tid=name+'_'+metal; task_ids.append(tid)
            tasks.append({'task_id':tid,'input':record(inp),'xyz':record(xp),
                          'output_path':str(d/'endpoint.out'),'metal':metal,'charge':charges[metal],
                          'case':name,'source_input':record(source_input),'source_xyz':record(source_xyz)})
        invariant=paired(output/name/'La/core.xyz',output/name/'Ca/core.xyz',charges['La'],charges['Ca'])
        cases.append({'case':name,'source_manifest':record(manifest_path),'task_ids':task_ids,
                      'paired_invariants':invariant,'use':'method_development_already_consumed',
                      'baseline_protocol_id':m['protocol_id']})
    result={'schema_version':'alquemia.affordable_pilot.v1','protocol_id':'native_r2scan3c_cpcm_mbis_endpoints_v1',
            'status':'approved_prepared','cases':cases,'tasks':tasks,
            'agreement':record(root/'diagnostics/affordable_challenger_20260915/CONTINUATION_AGREEMENT.md'),
            'execution_policy':{'task_runner':record(root/'scripts/run_orca_task_manifest.py'),
                                'runtime_renderer':record(root/'scripts/render_orca_runtime_input.py')},
            'implementation':record(__file__),
            'orca':record('/groups/banfield/users/jwestrob/bin/ORCA/orca_6_1_1_linux_x86-64_shared_openmpi418_nodmrg/orca'),
            'cost_tracking':{'compute_budget':None,'wall_time_limit':None,
                      'baseline_cost_evidence':record(root/'diagnostics/affordable_challenger_20260915/live_audit.json'),
                      'basis':'SLURM 1198050: 154 s x 344 CPUs = 52976 core-s for four endpoints; reference measurement only',
                      'runtime_policy':'execute agreed task list; record costs; no compute-budget stopping rule or timeout'},
            'numerical_checks':{'identity_kcal_mol':.01,'grid_box_rotation_kcal_mol':.5,'partition_score_kcal_mol':2.0},
            'environment_model':{'dielectric_inside':1.0,'dielectric_outside':78.54,'salt_molar':0.0,'temperature_K':298.15,
                                 'radii_A':{'H':1.2,'C':1.7,'N':1.55,'O':1.52,'S':1.8,'Ca':1.8,'La':1.8},
                                 'radii_policy':'Bondi_CHNOS_common_1p8A_metal_v1','probe_radius_A':1.4,
                                 'primary_grid_spacing_A':.5,'refined_grid_spacing_A':.4,'box_padding_A':20.,'expanded_box_padding_A':30.,
                                 'charge_model':'ORCA_MBIS_monopoles','charge_quality':'actual_ESP_relative_RMS_at_most_0p10_or_absolute_RMS_0p005_au',
                                 'physical_assembly':'deposited_catalytic_chain_A_only_not_biological_oligomer'}}
    write_new(output/'pilot_manifest.json',result)
    return result


def dry_run(manifest):
    from run_orca_task_manifest import load_manifest_tasks
    m,tasks=load_manifest_tasks(Path(manifest))
    verify(m['agreement']); verify(m['orca'])
    for r in m['execution_policy'].values(): verify(r)
    for task in tasks:
        verify({'path':str(task['input']),'sha256':task['input_sha256']})
        verify({'path':str(task['xyz']),'sha256':task['xyz_sha256']})
    return {'status':'dry_run_pass','tasks':len(tasks),'historical_budget':m.get('budget'),
            'budget_enforced':False,'manifest':record(manifest)}


def execute(manifest):
    from run_orca_task_manifest import run_manifest, load_manifest_tasks, _completed_attempt_is_valid
    dry_run(manifest)
    mp=Path(manifest).resolve(); m,tasks=load_manifest_tasks(mp)
    if not os.environ.get('SLURM_JOB_ID'): raise InvalidArtifact('execution requires SLURM')
    cpus=int(os.environ['SLURM_CPUS_ON_NODE'])
    events=mp.parent/'budget_events.jsonl'
    lock=(mp.parent/'execute.lock').open('a+')
    fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    previous=[json.loads(l) for l in events.read_text().splitlines()] if events.exists() else []
    admitted_jobs={e['slurm_job_id'] for e in previous if e.get('admitted_endpoint_count')}
    accounted_jobs={e['slurm_job_id'] for e in previous if 'allocated_core_seconds' in e}
    unaccounted_jobs=sorted(admitted_jobs-accounted_jobs)
    unfinished=[]
    for t in tasks:
        op=t['output']; ep=Path(str(op)+'.execution.json')
        if op.exists() or ep.exists():
            if op.exists() and ep.exists() and _completed_attempt_is_valid(ep,op,manifest_sha256=digest(mp),task=t,
                    runner_identity=m['execution_policy']['task_runner'],runtime_renderer_identity=m['execution_policy']['runtime_renderer']):
                continue
            raise InvalidArtifact(f"partial attempt retained: {t['task_id']}; prepare explicit fresh retry")
        unfinished.append(t['task_id'])
    if not unfinished: return {'status':'already_complete'}
    resources=m.get('execution_resources',{})
    ranks=int(resources.get('mpi_ranks',min(16,cpus)))
    if ranks < 1 or ranks > cpus: raise InvalidArtifact('invalid explicit MPI rank request')
    workers=int(resources.get('concurrent_tasks',min(len(unfinished),max(1,cpus//ranks))))
    if workers < 1 or workers*ranks > cpus: raise InvalidArtifact('invalid explicit task concurrency')
    workers=min(workers,len(unfinished))
    def append(event):
        with events.open('a') as f: f.write(json.dumps(event)+'\n'); f.flush(); os.fsync(f.fileno())
    implementation=snapshot_implementation(mp.parent/f"implementation_{os.environ['SLURM_JOB_ID']}")
    append({'admitted_endpoint_count':len(unfinished),'tasks':unfinished,'compute_budget':None,'unaccounted_prior_jobs':unaccounted_jobs,'actual_implementation':implementation,
            'slurm_job_id':os.environ['SLURM_JOB_ID'],'manifest_sha256':digest(mp),'time_unix':time.time()})
    start=time.monotonic(); error=None
    try:
        result=run_manifest(mp,orca_path=verify(m['orca']),workers=workers,nprocs=ranks,
                            selected_tasks=set(unfinished),expected_runner_sha256=m['execution_policy']['task_runner']['sha256'])
    except Exception as exc:
        error=str(exc); result=[]
    finally:
        elapsed=time.monotonic()-start
        append({'allocated_core_seconds':elapsed*cpus,'elapsed_seconds':elapsed,'allocated_cpus':cpus,
                'slurm_job_id':os.environ['SLURM_JOB_ID'],'error':error,'time_unix':time.time()})
    if error: raise InvalidArtifact(error)
    return {'status':'completed','tasks':result,'elapsed_seconds':elapsed,'allocated_core_seconds':elapsed*cpus}


def collect(manifest):
    from affordable_environment import mbis_charges
    m=read_json(manifest); rows=[]
    for t in m['tasks']:
        op=Path(t['output_path'])
        row={'task_id':t['task_id'],'case':t['case'],'metal':t['metal'],'status':'not_run','energy_hartree':None,'charges':None}
        if op.exists():
            try:
                row['energy_hartree']=energy(op)
                row['charges']=mbis_charges(op,verify(t['xyz']),t['charge'])
                row['status']='endpoint_and_mbis_available'
            except ValueError as exc:
                row['status']='unscorable'; row['reason']=str(exc)
            row['output']=record(op)
            ep=Path(str(op)+'.execution.json')
            row['execution_receipt']=record(ep) if ep.exists() else None
        rows.append(row)
    return {'schema_version':'alquemia.affordable_collection.v1','manifest':record(manifest),'rows':rows,
            'environmental_score_status':'requires_ESP_and_APBS_physical_checks','environmental_score':None}


def prepare_retry(manifest, task_id):
    """Fresh immutable attempt of exactly one failed endpoint; retain cost history."""
    mp=Path(manifest).resolve();m=read_json(mp)
    if m.get('retry_of'):
        raise InvalidArtifact('prepare retries from the primary pilot manifest')
    matches=[t for t in m['tasks'] if t['task_id']==task_id]
    if len(matches)!=1: raise InvalidArtifact('unknown retry task')
    t=matches[0];op=Path(t['output_path'])
    if not op.exists(): raise InvalidArtifact('no attempt exists to retry')
    try: energy(op)
    except ValueError: pass
    else: raise InvalidArtifact('successful energy is not rerun as a retry')
    events=mp.parent/'budget_events.jsonl'
    used=sum(json.loads(l).get('admitted_endpoint_count',0) for l in events.read_text().splitlines())
    count=len(list(mp.parent.glob('retry_*_manifest.json')))+1
    d=mp.parent/'retries'/f'{count}_{task_id}';d.mkdir(parents=True,exist_ok=False)
    xp=d/'core.xyz';xp.write_bytes(verify(t['xyz']).read_bytes())
    inp=d/'endpoint.inp';inp.write_bytes(verify(t['input']).read_bytes())
    revised=dict(t,input=record(inp),xyz=record(xp),output_path=str(d/'endpoint.out'),retry_source_output=record(op))
    result=dict(m,tasks=[revised],retry_of=record(mp),retry_task_id=task_id)
    path=mp.parent/f'retry_{count}_manifest.json';write_new(path,result)
    return {'status':'retry_prepared_unchanged','manifest':record(path),'prior_admissions':used}


def main():
    p=argparse.ArgumentParser(description=__doc__); sub=p.add_subparsers(dest='op',required=True)
    q=sub.add_parser('prepare'); q.add_argument('--root',type=Path,required=True);q.add_argument('--output',type=Path,required=True)
    for name in ('dry-run','execute','collect'):
        q=sub.add_parser(name);q.add_argument('--manifest',type=Path,required=True);q.add_argument('--output',type=Path)
    q=sub.add_parser('prepare-retry');q.add_argument('--manifest',type=Path,required=True);q.add_argument('--task',required=True)
    a=p.parse_args()
    if a.op=='prepare': r=prepare_pilot(a.root,a.output)
    elif a.op=='prepare-retry': r=prepare_retry(a.manifest,a.task)
    else:
        r={'dry-run':dry_run,'execute':execute,'collect':collect}[a.op](a.manifest)
        if a.output: write_new(a.output,r)
    print(json.dumps({k:v for k,v in r.items() if k in ('status','elapsed_seconds','allocated_core_seconds')}) or 'prepared')


if __name__=='__main__':
    import sys
    if len(sys.argv)>1 and sys.argv[1]=='baseline':
        from baseline_water import main as baseline_main
        baseline_main(sys.argv[2:])
    else: main()
