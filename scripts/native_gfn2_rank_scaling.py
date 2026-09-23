"""Three fixed MPI rank settings for eight real, immutable native GFN2 cells."""
from __future__ import annotations
import argparse
import copy
import datetime as dt
import json
import os
from pathlib import Path
import re
import shutil
import statistics
from affordable_common import InvalidArtifact,HA_TO_KCAL,read_json,record,verify,write_new,xyz
from affordable_workflow import dry_run,execute as execute_existing
from compact_solvation import completed,diagnostics
from adaptive_origin_recovery import snapshot
from mace_hybrid import EV_TO_KCAL

ROOT=Path(__file__).resolve().parents[1]
PROTOCOL='Nikasha_native_GFN2_matched_MPI_ranks_1_4_8_v1'
CASES=('q9z4j7-pqq-la_model','q89gy2-pqq-la_model')
RANKS=(1,4,8)
HOST='node-224-2t-8gpu-1'
SETTINGS={'ranks':[1,4,8],'workers':8,'cell_tolerance_kcal_mol':.1,'paired_R_tolerance_kcal_mol':.2,
          'reference_ranks':8,'host':HOST,'new_calls':24,'restart':'NoAutostart','rank_order':[1,4,8]}


def sources(source):
    m=read_json(source);ts=[]
    pool_collection=verify(m['pool_manifest']).parent/'collection_final.json';pool=read_json(pool_collection)
    if pool['manifest']!=m['pool_manifest']:raise InvalidArtifact('actual parent pool manifest differs')
    for cid in CASES:
        for z in ('Ca','La'):
            for medium in ('vacuum','alpb'):
                tid=f'{cid}__{z}__at_adaptive_Ca__{medium}'
                t=copy.deepcopy(next(t for t in m['tasks'] if t['task_id']==tid))
                pin=completed(source,tid)
                if pin is None:raise InvalidArtifact('actual source endpoint unavailable')
                body=verify(t['input']).read_text()
                if not all(v in body for v in ('Native-GFN2-xTB NoAutostart','%maxcore 2000','MaxIter 500','SmearTemp 300','UseXTBMixer true')):raise InvalidArtifact('source electronic recipe differs')
                if (t['charge'],t['multiplicity'])!=(-3 if z=='Ca' else -2,1):raise InvalidArtifact('source state differs')
                if len(xyz(verify(t['xyz'])))!=(145 if cid==CASES[0] else 184):raise InvalidArtifact('source composition differs')
                audit=diagnostics(pin,t)
                native=next(c for c in pool['cases'] if c['case_id']==cid)['matrix'][z]['adaptive_Ca']['MACE']
                nr=read_json(verify(native));rq=read_json(verify(nr['request']))
                if nr['status']!='complete' or xyz(verify(rq['xyz']))!=xyz(verify(t['xyz'])) or (rq['charge'],rq['multiplicity'])!=(t['charge'],1):raise InvalidArtifact('matched native energy/state differs')
                t.update(archived=pin,archived_audit=audit,native_energy_eV=nr['energy_eV'],native_source=native,native_collection=record(pool_collection));ts.append(t)
    for cid in CASES:
        ca,la=[xyz(verify(next(t for t in ts if t['case_id']==cid and t['metal']==z and t['medium']=='vacuum')['xyz'])) for z in ('Ca','La')]
        if ca[1:]!=la[1:] or ca[0][1:]!=la[0][1:]:raise InvalidArtifact('source paired coordinates differ')
    return m,ts


def prepare(source,agreement,output):
    old,ts=sources(source);root=Path(output).resolve();root.mkdir(parents=True,exist_ok=False)
    impl=snapshot(Path(__file__).parent,root/'implementation');manifests=[]
    for n in RANKS:
        out=root/f'rank_{n}';out.mkdir();tasks=[]
        for oldtask in ts:
            t=copy.deepcopy(oldtask);d=out/'tasks'/t['task_id'];d.mkdir(parents=True)
            for field,name in [('input','endpoint.inp'),('xyz','core.xyz')]:
                shutil.copyfile(verify(t[field]),d/name);t[field]=record(d/name)
            t['output_path']=str(d/'endpoint.out');tasks.append(t)
        m={'protocol_id':PROTOCOL,'settings':SETTINGS,'ranks':n,'source':record(source),'agreement':record(agreement),
           'orca':old['orca'],'implementation':impl,'execution_resources':{'mpi_ranks':n,'concurrent_tasks':8},
           'execution_policy':{k:impl[v] for k,v in [('task_runner','run_orca_task_manifest.py'),('runtime_renderer','render_orca_runtime_input.py')]},
           'tasks':tasks,'all_tasks':tasks,'production_changed':False}
        mp=out/'manifest.json';write_new(mp,m);validate(mp);manifests.append(record(mp))
    design={'protocol_id':PROTOCOL,'settings':SETTINGS,'source':record(source),'agreement':record(agreement),'implementation':impl,'manifests':manifests,'new_calls':24,'production_changed':False}
    dp=root/'design.json';write_new(dp,design);write_new(root/'PREFLIGHT.json',{'design':record(dp),'rank_manifests':manifests,'tasks_per_rank':8,'new_calls':24,'status':'validated','new_calls_in_preflight':0});return record(dp)


def validate(manifest):
    m=read_json(manifest)
    if m['protocol_id']!=PROTOCOL or m['settings']!=SETTINGS or m['ranks'] not in RANKS:raise InvalidArtifact('frozen numerical scope differs')
    if m['execution_resources']!={'mpi_ranks':m['ranks'],'concurrent_tasks':8}:raise InvalidArtifact('runtime resource policy differs')
    if m['tasks']!=m['all_tasks'] or len(m['tasks'])!=8:raise InvalidArtifact('eight fresh tasks required')
    for p in m['implementation'].values():verify(p)
    old=read_json(verify(m['source']));by={t['task_id']:t for t in old['tasks']}
    expected={f'{c}__{z}__at_adaptive_Ca__{s}' for c in CASES for z in ('Ca','La') for s in ('vacuum','alpb')}
    if {t['task_id'] for t in m['tasks']}!=expected or m['orca']!=old['orca']:raise InvalidArtifact('source membership/backend differs')
    for t in m['tasks']:
        original=by[t['task_id']]
        if any(t[k]!=original[k] for k in ('case_id','candidate','metal','medium','charge','multiplicity')) or t.get('native_reuse')!=original.get('native_reuse'):raise InvalidArtifact('source state changed')
        native=read_json(verify(t['native_source']));rq=read_json(verify(native['request']));verify(t['native_collection'])
        if native['energy_eV']!=t['native_energy_eV'] or xyz(verify(rq['xyz']))!=xyz(verify(t['xyz'])) or (rq['charge'],rq['multiplicity'])!=(t['charge'],1):raise InvalidArtifact('native energy reuse changed')
        if any(verify(t[k]).read_bytes()!=verify(original[k]).read_bytes() for k in ('input','xyz')):raise InvalidArtifact('source geometry or electronic recipe changed')
        out=Path(t['output_path']).resolve()
        if out.parent!=verify(t['input']).parent or not out.is_relative_to(Path(manifest).resolve().parent):raise InvalidArtifact('task escaped rank directory')
    return dry_run(manifest)


def execute(manifest):
    validate(manifest);m=read_json(manifest)
    if os.environ.get('SLURM_JOB_NODELIST')!=HOST or int(os.environ.get('SLURM_CPUS_ON_NODE','0'))!=8*m['ranks']:raise InvalidArtifact('matched host/allocation differs')
    return execute_existing(manifest)


def collect(manifest,output):
    validate(manifest);m=read_json(manifest);rows=[]
    for t in m['tasks']:
        row={k:t[k] for k in ('task_id','case_id','metal','medium','charge','multiplicity','native_energy_eV','archived')}
        op=Path(t['output_path']);rp=Path(str(op)+'.execution.json');row.update(status='unavailable',energy_hartree=None,reason=None,available_artifacts=[record(p) for p in (op,rp) if p.exists()])
        try:
            pin=completed(manifest,t['task_id'])
            if pin is None:raise InvalidArtifact('native endpoint incomplete or failed')
            rec=read_json(verify(pin['receipt']));audit=diagnostics(pin,t)
            if rec['parallelism']!={'nprocs':m['ranks'],'omp_threads_per_rank':1,'binding':'none'} or rec['allocation']['host']!=HOST:raise InvalidArtifact('actual rank or host differs')
            runtime=verify(rec['artifacts']['runtime_input']).read_text()
            pal=re.search(r'%pal\s+nprocs\s+(\d+)\s+end',runtime,re.I)
            if pal is None or int(pal[1])!=m['ranks']:raise InvalidArtifact('actual runtime PAL differs')
            if read_json(verify(audit['parameter_export']))!=read_json(verify(t['archived_audit']['parameter_export'])):raise InvalidArtifact('native parameters changed')
            text=verify(pin['output']).read_text();iters=re.findall(r'SCF CONVERGED AFTER\s+(\d+)\s+CYCLES',text)
            seconds=(dt.datetime.fromisoformat(rec['finished_at_utc'].replace('Z','+00:00'))-dt.datetime.fromisoformat(rec['started_at_utc'].replace('Z','+00:00'))).total_seconds()
            row.update(status='complete',**pin,audit=audit,iterations=int(iters[-1]) if iters else None,wall_seconds=seconds,task_rank_seconds=seconds*m['ranks'])
        except Exception as exc:row['reason']=str(exc)
        rows.append(row)
    result={'protocol_id':PROTOCOL,'manifest':record(manifest),'ranks':m['ranks'],'denominator':8,'complete':sum(x['status']=='complete' for x in rows),'rows':rows,'new_calls_in_collection':0}
    write_new(output,result);return result


def compare(design,output):
    d=read_json(design)
    if d['protocol_id']!=PROTOCOL or d['settings']!=SETTINGS or len(d['manifests'])!=3:raise InvalidArtifact('fixed rank design differs')
    collections={};rows=[]
    for pin in d['manifests']:
        mp=verify(pin);m=read_json(mp);cp=mp.parent/'collection.json'
        if not cp.exists():collections[m['ranks']]={'status':'unavailable','rows':[]};continue
        c=read_json(cp)
        if c['manifest']!=pin:raise InvalidArtifact('collection belongs to different rank task')
        collections[m['ranks']]={**c,'collection':record(cp)}
    matrices={n:{r['task_id']:r for r in collections[n]['rows']} for n in RANKS}
    for n in RANKS:
        cells=[];pairs=[]
        for t in read_json(verify(d['manifests'][0]))['tasks']:
            a=matrices[n].get(t['task_id']);b=matrices[8].get(t['task_id']);ok=bool(a and b and a['status']==b['status']=='complete')
            delta=(a['energy_hartree']-b['energy_hartree'])*HA_TO_KCAL if ok else None
            old=(a['energy_hartree']-t['archived']['energy_hartree'])*HA_TO_KCAL if a and a['status']=='complete' else None
            cells.append({'task_id':t['task_id'],'status':'available' if ok else 'unavailable','delta_vs_fresh8_kcal_mol':delta,'delta_vs_old8_kcal_mol':old,'passes_component_gate':abs(delta)<=.1 if delta is not None else None})
        for cid in CASES:
            scores={}
            for rank in (n,8):
                values={}
                for z in ('Ca','La'):
                    low={s:matrices[rank].get(f'{cid}__{z}__at_adaptive_Ca__{s}') for s in ('vacuum','alpb')}
                    if not all(v and v['status']=='complete' for v in low.values()):break
                    values[z]={'native_eV':low['vacuum']['native_energy_eV'],'vacuum_Ha':low['vacuum']['energy_hartree'],'alpb_Ha':low['alpb']['energy_hartree']}
                if len(values)!=2:scores[rank]=None;continue
                ca,la=values['Ca'],values['La'];native=(ca['native_eV']-la['native_eV'])*EV_TO_KCAL;vac=(ca['vacuum_Ha']-la['vacuum_Ha'])*HA_TO_KCAL;alpb=(ca['alpb_Ha']-la['alpb_Ha'])*HA_TO_KCAL
                scores[rank]={'native_R_model_kcal_mol':native,'vacuum_R_kcal_mol':vac,'ALPB_R_kcal_mol':alpb,'solvent_R_kcal_mol':alpb-vac,'composite_R_model_kcal_mol':native+alpb-vac}
            delta=scores[n]['composite_R_model_kcal_mol']-scores[8]['composite_R_model_kcal_mol'] if scores[n] is not None and scores[8] is not None else None
            pairs.append({'case_id':cid,'components':scores[n],'fresh8_components':scores[8],'delta_R_vs_fresh8_kcal_mol':delta,'passes_R_gate':abs(delta)<=.2 if delta is not None else None})
        times=[x['wall_seconds'] for x in collections[n]['rows'] if x['status']=='complete']
        rows.append({'ranks':n,'cells':cells,'pairs':pairs,'all_numerical_gates_pass':all(x['passes_component_gate'] is True for x in cells) and all(x['passes_R_gate'] is True for x in pairs),
          'successful_wall_seconds':{'count':len(times),'min':min(times) if times else None,'median':statistics.median(times) if times else None,'max':max(times) if times else None,'sum':sum(times)},'successful_task_rank_seconds':sum(times)*n})
    result={'protocol_id':PROTOCOL,'design':record(design),'collections':{str(n):c.get('collection') for n,c in collections.items()},'rank_rows':rows,'new_calls_in_comparison':0,'biological_accuracy_claim':False,'production_changed':False}
    write_new(output,result);return result


def main():
    p=argparse.ArgumentParser(description=__doc__);s=p.add_subparsers(dest='op',required=True)
    for op,fields in {'prepare':('source','agreement','output'),'validate':('manifest',),'execute':('manifest',),'collect':('manifest','output'),'compare':('design','output')}.items():
        q=s.add_parser(op)
        for f in fields:q.add_argument('--'+f,required=True,type=Path)
    a=vars(p.parse_args());print(json.dumps(globals()[a.pop('op')](**a),indent=2))
if __name__=='__main__':main()
