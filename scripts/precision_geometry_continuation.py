"""Two uniform self-continuations of four exact rank-one anomaly cells."""
from __future__ import annotations
import argparse,fcntl,json,os,shutil,time
from pathlib import Path
from affordable_common import InvalidArtifact,HA_TO_KCAL,read_json,record,verify,write_new
from compact_solvation import completed,diagnostics
from native_pool_continuation import recipe,collect_task
from run_orca_task_manifest import load_manifest_tasks,run_manifest
from structure_informed_starts import scf_details
from union_adaptive import snapshot

PROTOCOL='native_GFN2_geometry_sensitive_rank1_two_self_continuations_v1'
CASES={'q4w6g0-pqq-la_model__conditioned_Ca__seed-1_sample-2':'adaptive_Ca',
       'p38539-pqq-la_model__conditioned_La__seed-1_sample-2':'adaptive_La'}
SETTINGS={'stages':2,'reported_stage':2,'mpi_ranks':1,'workers':4,'MaxIter':500,
          'electronic_temperature_K':300,'component_tolerance_kcal_mol':.1,
          'maximum_new_scalar_calls':8,'seed_policy':'same_geometry_same_cell_only'}


def inventory(collection,output):
    c=read_json(collection);rows=[]
    if len(c['targets'])!=2 or {r['case_id']:r['candidate'] for r in c['targets']}!=CASES:
        raise InvalidArtifact('exact two anomaly targets required')
    for target in c['targets']:
        for geometry in ('old','new'):
            r=target['matrix'][geometry]['1'];receipt=read_json(verify(r['receipt']))
            mp=verify(receipt['manifest']);m=read_json(mp)
            t=next(t for t in m['tasks'] if t['task_id']==receipt['task_id'])
            pin=completed(mp,t['task_id'])
            if not pin or pin['energy_hartree']!=r['energy_hartree'] or receipt['parallelism']['nprocs']!=1:
                raise InvalidArtifact('source success/energy/rank mismatch')
            detail=scf_details(verify(r['output']).read_text());audit=diagnostics(pin,t)
            if (t['metal'],t['medium'],t['charge'],t['multiplicity'])!=('La','vacuum',-1,1):
                raise InvalidArtifact('source electronic state differs')
            if verify(t['input']).read_text().replace(' NoAutostart','')!=recipe(-1,1,'vacuum'):
                raise InvalidArtifact('source recipe differs from declared self-continuation')
            runtime=verify(receipt['artifacts']['runtime_input']);base=runtime.with_suffix('')
            gbw=Path(str(base)+'.gbw')
            xtbw=Path(str(base)+'.xtbw')
            row={'case_id':target['case_id'],'geometry':geometry,'candidate':target['candidate'],
                 'metal':'La','medium':'vacuum','charge':-1,'multiplicity':1,'electron_count':detail['electrons'],
                 'xyz':t['xyz'],'input':t['input'],'output':r['output'],'receipt':r['receipt'],
                 'source_manifest':record(mp),'source_task_id':t['task_id'],'orca':receipt['orca_executable'],
                 'parameter_export':audit['parameter_export'],'source_energy_hartree':r['energy_hartree'],
                 'source_SCF_cycles':detail['cycles'],'gbw':record(gbw) if gbw.exists() else None,
                 'xtbw':record(xtbw) if xtbw.exists() else None}
            row['status']='compatible_seed_pair_available' if row['gbw'] and row['xtbw'] else 'seed_pair_unavailable'
            rows.append(row)
    result={'protocol_id':PROTOCOL,'collection':record(collection),'rows':rows,'denominator':4,
            'available':sum(r['status']=='compatible_seed_pair_available' for r in rows),'new_calls':0}
    write_new(output,result);return result


def source_rows(path):
    data=read_json(path);verify(data['collection']);rows=data['rows']
    if data['protocol_id']!=PROTOCOL or len(rows)!=4 or data['available']!=4 or {
        (r['case_id'],r['geometry']) for r in rows}!={(c,g) for c in CASES for g in ('old','new')}:
        raise InvalidArtifact('all four exact source seeds required')
    for r in rows:
        for k in ('gbw','xtbw','xyz','input','output','receipt','parameter_export','source_manifest','orca'):verify(r[k])
        if r['status']!='compatible_seed_pair_available' or (r['charge'],r['multiplicity'],r['metal'],r['medium'])!=(-1,1,'La','vacuum'):
            raise InvalidArtifact('source seed/state unsupported')
        if verify(r['gbw']).parent!=verify(r['xtbw']).parent:raise InvalidArtifact('unmatched seed pair')
        if verify(r['input']).read_text().replace(' NoAutostart','')!=recipe(-1,1,'vacuum'):
            raise InvalidArtifact('source recipe differs')
    return rows


def tid(r):return r['case_id']+'__'+r['geometry']


def prepare(sources,agreement,output,previous=None):
    rows=source_rows(sources);stage=2 if previous else 1;prior={};pm=None
    if previous:
        p=read_json(previous);pm=read_json(verify(p['manifest']))
        if p['stage']!=1 or pm['sources']!=record(sources) or pm['agreement']!=record(agreement):
            raise InvalidArtifact('stage2 requires matching stage1')
        prior={r['task_id']:r for r in p['rows']}
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False);tasks=[];missing=[]
    for r in rows:
        pr=prior.get(tid(r));seed={k:r[k] for k in ('gbw','xtbw')}
        if stage==2:
            if pr is None or pr['status']!='confirmed_restart':
                missing.append({'task_id':tid(r),'case_id':r['case_id'],'geometry':r['geometry'],
                                'status':'unavailable','reason':'stage1_restart_unavailable','energy_hartree':None});continue
            seed={'gbw':pr['seed_after']['preserved_gbw_after'],'xtbw':pr['seed_after']['preserved_after']}
        d=out/'tasks'/tid(r);d.mkdir(parents=True)
        shutil.copyfile(verify(r['xyz']),d/'core.xyz')
        for k in ('gbw','xtbw'):shutil.copyfile(verify(seed[k]),d/('seed.immutable.'+k))
        (d/'endpoint.inp').write_text(recipe(-1,1,'vacuum'))
        tasks.append({**{k:r[k] for k in ('case_id','geometry','candidate','metal','medium','charge','multiplicity')},
            'task_id':tid(r),'source':r,'seed_source':seed,'gradient_requested':False,
            'xyz':record(d/'core.xyz'),'input':record(d/'endpoint.inp'),'output_path':str(d/'endpoint.out'),
            'immutable_seed':record(d/'seed.immutable.xtbw'),'immutable_gbw':record(d/'seed.immutable.gbw'),
            'active_seed_path':str(d/'endpoint.runtime.xtbw'),'active_gbw_path':str(d/'endpoint.runtime.gbw')})
    if pm:
        impl={};source=verify(pm['implementation']['precision_geometry_continuation.py']).parent
        dest=out/'implementation';dest.mkdir()
        for p in source.glob('*.py'):shutil.copyfile(p,dest/p.name);impl[p.name]=record(dest/p.name)
    else:impl=snapshot(out/'implementation')
    m={'protocol_id':PROTOCOL,'stage':stage,'settings':SETTINGS,'tasks':tasks,'missing':missing,
       'sources':record(sources),'agreement':record(agreement),'previous':record(previous) if previous else None,
       'orca':rows[0]['orca'],'implementation':impl,'denominator':4,
       'execution_resources':{'mpi_ranks':1,'concurrent_tasks':4},
       'execution_policy':{'task_runner':impl['run_orca_task_manifest.py'],'runtime_renderer':impl['render_orca_runtime_input.py']},
       'baseline_changed':False,'new_MACE_DFT_optimization_calls':0}
    mp=out/'manifest.json';write_new(mp,m);v=validate(mp,True);write_new(out/'PREFLIGHT.json',v);return v


def validate(manifest,fresh=False):
    mp=Path(manifest).resolve();m=read_json(mp);rows={tid(r):r for r in source_rows(verify(m['sources']))}
    if m['protocol_id']!=PROTOCOL or m['settings']!=SETTINGS or m['stage'] not in (1,2) or m['execution_resources']!={'mpi_ranks':1,'concurrent_tasks':4}:
        raise InvalidArtifact('fixed method/stage/parallelism differs')
    verify(m['agreement']);verify(m['orca'])
    for pin in m['implementation'].values():verify(pin)
    if len(m['tasks'])+len(m['missing'])!=4 or {t['task_id'] for t in m['tasks']+m['missing']}!=set(rows):
        raise InvalidArtifact('fixed four-cell denominator differs')
    prior={}
    if m['stage']==1:
        if m['previous'] or m['missing']:raise InvalidArtifact('first stage requires all seeds')
    else:
        p=read_json(verify(m['previous']));pm=read_json(verify(p['manifest']))
        if p['stage']!=1 or pm['sources']!=m['sources'] or pm['agreement']!=m['agreement']:raise InvalidArtifact('continuation linkage differs')
        prior={r['task_id']:r for r in p['rows']}
        if {r['task_id'] for r in m['missing']}!={k for k,r in prior.items() if r['status']!='confirmed_restart'}:
            raise InvalidArtifact('failed previous cell lost or successful cell dropped')
    if m['tasks']:load_manifest_tasks(mp)
    for t in m['tasks']:
        r=rows[t['task_id']]
        seed={k:r[k] for k in ('gbw','xtbw')} if m['stage']==1 else {
            'gbw':prior[t['task_id']]['seed_after']['preserved_gbw_after'],'xtbw':prior[t['task_id']]['seed_after']['preserved_after']}
        if t['source']!=r or t['seed_source']!=seed or any(t[k]!=r[k] for k in ('case_id','geometry','candidate','metal','medium','charge','multiplicity')):
            raise InvalidArtifact('source or exact self-seed changed')
        if t['gradient_requested'] or verify(t['xyz']).read_bytes()!=verify(r['xyz']).read_bytes() or verify(t['input']).read_text()!=recipe(-1,1,'vacuum'):
            raise InvalidArtifact('geometry/method changed')
        for field,k in (('immutable_seed','xtbw'),('immutable_gbw','gbw')):
            if verify(t[field]).read_bytes()!=verify(seed[k]).read_bytes():raise InvalidArtifact('immutable seed differs')
        d=Path(t['output_path']).parent
        if Path(t['active_seed_path'])!=d/'endpoint.runtime.xtbw' or Path(t['active_gbw_path'])!=d/'endpoint.runtime.gbw':raise InvalidArtifact('runtime seed basename differs')
        if fresh and {p.name for p in d.iterdir()}!={'core.xyz','endpoint.inp','seed.immutable.gbw','seed.immutable.xtbw'}:
            raise InvalidArtifact('nonfresh task directory')
    return {'status':'validated','manifest':record(mp),'stage':m['stage'],'tasks':len(m['tasks']),'missing':len(m['missing']),'new_calls_in_validation':0}


def execute(manifest):
    mp=Path(manifest).resolve();m=read_json(mp)
    if not os.environ.get('SLURM_JOB_ID') or int(os.environ.get('SLURM_CPUS_ON_NODE','0'))<4:raise InvalidArtifact('four-CPU allocation required')
    with (mp.parent/'execute.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB);validate(mp,True);before=[]
        for t in m['tasks']:
            shutil.copyfile(verify(t['immutable_seed']),t['active_seed_path']);shutil.copyfile(verify(t['immutable_gbw']),t['active_gbw_path'])
            before.append({'task_id':t['task_id'],'active_seed':record(t['active_seed_path']),'active_gbw':record(t['active_gbw_path'])})
        write_new(mp.parent/'SEEDS_BEFORE.json',{'manifest':record(mp),'rows':before})
        start=time.monotonic();error=None;result=None
        try:result=run_manifest(mp,orca_path=verify(m['orca']),workers=4,nprocs=1) if m['tasks'] else []
        except Exception as exc:error=repr(exc)
        finally:
            after=[]
            for t in m['tasks']:
                d=Path(t['output_path']).parent;row={'task_id':t['task_id']}
                for field,name,outfield in (('active_seed_path','xtbw','preserved_after'),('active_gbw_path','gbw','preserved_gbw_after')):
                    p=Path(t[field]);dest=d/('seed.after.'+name)
                    if p.exists():shutil.copyfile(p,dest)
                    row[outfield]=record(dest) if dest.exists() else None
                after.append(row)
            write_new(mp.parent/'SEEDS_AFTER.json',{'manifest':record(mp),'rows':after})
            write_new(mp.parent/'EXECUTION.json',{'manifest':record(mp),'results':result,'error':error,
                'wall_seconds':time.monotonic()-start,'allocated_cpus':int(os.environ['SLURM_CPUS_ON_NODE']),'slurm_job_id':os.environ['SLURM_JOB_ID']})
        if error:raise InvalidArtifact(error)
    return {'status':'executed','tasks':len(m['tasks'])}


def collect(manifest,output):
    validate(manifest);mp=Path(manifest).resolve();m=read_json(mp)
    def idx(name):
        p=mp.parent/name;return {r['task_id']:r for r in read_json(p)['rows']} if p.exists() else {}
    before=idx('SEEDS_BEFORE.json');after=idx('SEEDS_AFTER.json');rows=[]
    for t in m['tasks']:
        r=collect_task(mp,t,before.get(t['task_id']),after.get(t['task_id']));r['geometry']=t['geometry']
        if r['status']=='confirmed_restart' and read_json(verify(r['actual']['receipt']))['parallelism']['nprocs']!=1:
            r.update(status='audit_failed',reason='actual rank differs',energy_hartree=None)
        rows.append(r)
    rows.extend(m['missing'])
    result={'protocol_id':PROTOCOL,'manifest':record(mp),'stage':m['stage'],'rows':rows,'denominator':4,
            'confirmed':sum(r['status']=='confirmed_restart' for r in rows),'new_calls_in_collection':0}
    write_new(output,result);return {k:v for k,v in result.items() if k!='rows'}


def compare(stage1,stage2,output):
    cols=[read_json(p) for p in (stage1,stage2)];ms=[read_json(verify(c['manifest'])) for c in cols]
    if cols[0]['stage']!=1 or cols[1]['stage']!=2 or ms[1]['previous']!=record(stage1) or ms[0]['sources']!=ms[1]['sources']:
        raise InvalidArtifact('exact two-pass linkage required')
    source=source_rows(verify(ms[0]['sources']));indices=[{r['task_id']:r for r in c['rows']} for c in cols];rows=[]
    for s in source:
        rr=[i[tid(s)] for i in indices];ok=all(r['status']=='confirmed_restart' for r in rr)
        vals=[r['energy_hartree'] for r in rr];step=(vals[1]-vals[0])*HA_TO_KCAL if ok else None
        rows.append({'case_id':s['case_id'],'geometry':s['geometry'],'candidate':s['candidate'],
            'status':'complete' if ok else 'unavailable','cold_energy_hartree':s['source_energy_hartree'],
            'stage1_energy_hartree':vals[0],'stage2_energy_hartree':vals[1],
            'stage2_minus_stage1_kcal_mol':step,'settled':abs(step)<=.1 if ok else None,
            'stage2_minus_cold_kcal_mol':(vals[1]-s['source_energy_hartree'])*HA_TO_KCAL if ok else None,
            'actual':rr})
    targets=[]
    for cid in CASES:
        rr={r['geometry']:r for r in rows if r['case_id']==cid};works={}
        for tag in ('cold','stage1','stage2'):
            vals=[rr[g][tag+'_energy_hartree'] for g in ('old','new')]
            works[tag]=(vals[1]-vals[0])*HA_TO_KCAL if all(x is not None for x in vals) else None
        targets.append({'case_id':cid,'geometry_new_minus_old_kcal_mol':works,
            'stage2_geometry_difference_within_component_gate':abs(works['stage2'])<=.1 if works['stage2'] is not None else None})
    result={'protocol_id':PROTOCOL,'stage1':record(stage1),'stage2':record(stage2),'reported_stage':2,
            'settings':SETTINGS,'rows':rows,'targets':targets,'denominator':4,
            'complete':sum(r['status']=='complete' for r in rows),'settled':sum(r['settled'] is True for r in rows),
            'new_pool_score':None,'new_molecular_calls_in_comparison':0,'production_changed':False}
    write_new(output,result);return {k:v for k,v in result.items() if k!='rows'}


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest='op',required=True)
    for op,fields in {'inventory':('collection','output'),'prepare':('sources','agreement','output'),
                      'validate':('manifest',),'execute':('manifest',),'collect':('manifest','output'),
                      'compare':('stage1','stage2','output')}.items():
        s=sub.add_parser(op)
        for field in fields:s.add_argument('--'+field,required=True,type=Path)
        if op=='prepare':s.add_argument('--previous',type=Path)
    args=vars(p.parse_args());print(json.dumps(globals()[args.pop('op')](**args),indent=2))
