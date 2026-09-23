"""Finite stages for the declared32 precision-pool continuation qualification."""
from __future__ import annotations
import argparse,fcntl,json,os,shutil,time
from pathlib import Path
from affordable_common import InvalidArtifact,read_json,record,verify,write_new
from native_pool_continuation import recipe,collect_task
from run_orca_task_manifest import load_manifest_tasks,run_manifest
from union_adaptive import snapshot
from precision_pool_continuation import PROTOCOL,CANDIDATES,FOLDS,CRYSTALS

SETTINGS={'stages':2,'reported_stage':2,'mpi_ranks':1,'workers':32,'MaxIter':500,
          'electronic_temperature_K':300,'component_tolerance_kcal_mol':.1,
          'contrast_tolerance_kcal_mol':.2,'pool_tolerance_kcal_mol':.2,
          'origin_selection_tolerance_kcal_mol':.1,'logical_call_denominator':768}


def tid(r):return '__'.join(r[k] for k in ('case_id','candidate','metal','medium'))


def source_rows(path):
    inv=read_json(path);verify(inv['agreement'])
    rows=inv['rows'];cases=inv['cases'];ids={c['case_id'] for c in cases}
    if inv['protocol_id']!=PROTOCOL or len(cases)!=32 or len(ids)!=32 or len(rows)!=384:
        raise InvalidArtifact('declared32/384 population differs')
    if sum(c['role']=='calibration' for c in cases)!=25 or not set(FOLDS+CRYSTALS)<=ids:
        raise InvalidArtifact('calibration/transfer membership differs')
    if {(r['case_id'],r['candidate'],r['metal'],r['medium']) for r in rows}!={
        (c,q,z,s) for c in ids for q in CANDIDATES for z in ('Ca','La') for s in ('vacuum','alpb')}:
        raise InvalidArtifact('cell membership differs')
    return inv,rows


def prepare(inventory,output,previous=None):
    inv,rows=source_rows(inventory);stage=2 if previous else 1;prior={};pm=None
    if previous:
        pc=read_json(previous);pm=read_json(verify(pc['manifest']))
        if pc['stage']!=1 or pm['inventory']!=record(inventory):raise InvalidArtifact('exact stage1 required')
        prior={r['task_id']:r for r in pc['rows']}
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False);tasks=[];missing=[];reused=[]
    for r in rows:
        identity={k:r[k] for k in ('case_id','candidate','metal','medium')};identity['task_id']=tid(r)
        if r['status']!='compatible_seed_pair_available':
            missing.append({**identity,'status':'unavailable','reason':r['reason'],'energy_hartree':None});continue
        if r['continued_reuse'] is not None:
            reuse=r['continued_reuse'];entry=reuse['stage'+str(stage)]
            verify(reuse['stage'+str(stage)+'_collection']);verify(entry['actual']['output']);verify(entry['actual']['receipt'])
            reused.append({**identity,'source':r,'actual_row':entry,'actual_collection':reuse['stage'+str(stage)+'_collection']});continue
        seed={k:r[k] for k in ('gbw','xtbw')}
        if stage==2:
            pr=prior[tid(r)]
            if pr['status']!='confirmed_restart':
                missing.append({**identity,'status':'unavailable','reason':'stage1_restart_unavailable','energy_hartree':None});continue
            seed={'gbw':pr['seed_after']['preserved_gbw_after'],'xtbw':pr['seed_after']['preserved_after']}
        d=out/'tasks'/tid(r);d.mkdir(parents=True)
        shutil.copyfile(verify(r['xyz']),d/'core.xyz')
        for k in ('gbw','xtbw'):shutil.copyfile(verify(seed[k]),d/('seed.immutable.'+k))
        (d/'endpoint.inp').write_text(recipe(r['charge'],r['multiplicity'],r['medium']))
        tasks.append({**identity,'source':r,'seed_source':seed,'charge':r['charge'],'multiplicity':r['multiplicity'],
            'gradient_requested':False,'xyz':record(d/'core.xyz'),'input':record(d/'endpoint.inp'),'output_path':str(d/'endpoint.out'),
            'immutable_seed':record(d/'seed.immutable.xtbw'),'immutable_gbw':record(d/'seed.immutable.gbw'),
            'active_seed_path':str(d/'endpoint.runtime.xtbw'),'active_gbw_path':str(d/'endpoint.runtime.gbw')})
    if pm:
        impl={};source=verify(pm['implementation']['precision_pool_continue_run.py']).parent;dest=out/'implementation';dest.mkdir()
        for p in source.glob('*.py'):shutil.copyfile(p,dest/p.name);impl[p.name]=record(dest/p.name)
    else:impl=snapshot(out/'implementation')
    m={'protocol_id':PROTOCOL,'stage':stage,'settings':SETTINGS,'tasks':tasks,'missing':missing,'reused':reused,
       'inventory':record(inventory),'agreement':inv['agreement'],'previous':record(previous) if previous else None,
       'orca':next(r['orca'] for r in rows if r['status']=='compatible_seed_pair_available'),'implementation':impl,
       'case_denominator':32,'cell_denominator':384,'execution_resources':{'mpi_ranks':1,'concurrent_tasks':32},
       'execution_policy':{'task_runner':impl['run_orca_task_manifest.py'],'runtime_renderer':impl['render_orca_runtime_input.py']},
       'production_changed':False,'new_MACE_DFT_optimization_calls':0}
    mp=out/'manifest.json';write_new(mp,m);v=validate(mp,True);write_new(out/'PREFLIGHT.json',v);return v


def validate(manifest,fresh=False):
    mp=Path(manifest).resolve();m=read_json(mp);inv,rows=source_rows(verify(m['inventory']));source={tid(r):r for r in rows}
    if m['protocol_id']!=PROTOCOL or m['settings']!=SETTINGS or m['stage'] not in (1,2) or m['agreement']!=inv['agreement']:
        raise InvalidArtifact('fixed method/stage/agreement differs')
    if m['execution_resources']!={'mpi_ranks':1,'concurrent_tasks':32}:raise InvalidArtifact('execution profile differs')
    verify(m['orca'])
    for pin in m['implementation'].values():verify(pin)
    allrows=m['tasks']+m['missing']+m['reused']
    if len(allrows)!=384 or {t['task_id'] for t in allrows}!=set(source):raise InvalidArtifact('denominator differs')
    expected_reuse={tid(r) for r in rows if r['continued_reuse'] is not None}
    if {r['task_id'] for r in m['reused']}!=expected_reuse:raise InvalidArtifact('exact reuse membership differs')
    prior={}
    if m['stage']==1:
        if m['previous']:raise InvalidArtifact('first stage has previous')
        expected_missing={tid(r) for r in rows if r['status']!='compatible_seed_pair_available'}
    else:
        pc=read_json(verify(m['previous']));pm=read_json(verify(pc['manifest']))
        if pc['stage']!=1 or pm['inventory']!=m['inventory']:raise InvalidArtifact('stage linkage differs')
        prior={r['task_id']:r for r in pc['rows']}
        expected_missing={k for k,r in prior.items() if r['status']!='confirmed_restart'}
    if {r['task_id'] for r in m['missing']}!=expected_missing:raise InvalidArtifact('missing cells removed or successes dropped')
    if m['tasks']:load_manifest_tasks(mp)
    for t in m['reused']:
        r=source[t['task_id']];stage='stage'+str(m['stage']);reuse=r['continued_reuse']
        if t['source']!=r or t['actual_row']!=reuse[stage] or t['actual_collection']!=reuse[stage+'_collection']:raise InvalidArtifact('reused output changed')
        verify(t['actual_collection']);verify(t['actual_row']['actual']['receipt']);verify(t['actual_row']['actual']['output'])
    for t in m['tasks']:
        r=source[t['task_id']];seed={k:r[k] for k in ('gbw','xtbw')} if m['stage']==1 else {
            'gbw':prior[t['task_id']]['seed_after']['preserved_gbw_after'],'xtbw':prior[t['task_id']]['seed_after']['preserved_after']}
        if t['source']!=r or t['seed_source']!=seed or any(t[k]!=r[k] for k in ('case_id','candidate','metal','medium','charge','multiplicity')):
            raise InvalidArtifact('source/self-seed changed')
        if t['gradient_requested'] or verify(t['xyz']).read_bytes()!=verify(r['xyz']).read_bytes() or verify(t['input']).read_text()!=recipe(r['charge'],r['multiplicity'],r['medium']):
            raise InvalidArtifact('geometry/recipe differs')
        for field,k in (('immutable_seed','xtbw'),('immutable_gbw','gbw')):
            if verify(t[field]).read_bytes()!=verify(seed[k]).read_bytes():raise InvalidArtifact('seed bytes differ')
        d=Path(t['output_path']).parent
        if Path(t['active_seed_path'])!=d/'endpoint.runtime.xtbw' or Path(t['active_gbw_path'])!=d/'endpoint.runtime.gbw':raise InvalidArtifact('runtime basename differs')
        if fresh and {p.name for p in d.iterdir()}!={'core.xyz','endpoint.inp','seed.immutable.gbw','seed.immutable.xtbw'}:
            raise InvalidArtifact('nonfresh task directory')
    return {'status':'validated','manifest':record(mp),'stage':m['stage'],'fresh_tasks':len(m['tasks']),
            'reused':len(m['reused']),'missing':len(m['missing']),'cell_denominator':384,'new_calls_in_validation':0}


def execute(manifest):
    mp=Path(manifest).resolve();m=read_json(mp)
    if not os.environ.get('SLURM_JOB_ID') or int(os.environ.get('SLURM_CPUS_ON_NODE','0'))<32:raise InvalidArtifact('32-CPU allocation required')
    with (mp.parent/'execute.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB);validate(mp,True);before=[]
        for t in m['tasks']:
            shutil.copyfile(verify(t['immutable_seed']),t['active_seed_path']);shutil.copyfile(verify(t['immutable_gbw']),t['active_gbw_path'])
            before.append({'task_id':t['task_id'],'active_seed':record(t['active_seed_path']),'active_gbw':record(t['active_gbw_path'])})
        write_new(mp.parent/'SEEDS_BEFORE.json',{'manifest':record(mp),'rows':before})
        start=time.monotonic();error=None;result=None
        try:result=run_manifest(mp,orca_path=verify(m['orca']),workers=32,nprocs=1) if m['tasks'] else []
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
        r=collect_task(mp,t,before.get(t['task_id']),after.get(t['task_id']));r['reused']=False
        if r['status']=='confirmed_restart' and read_json(verify(r['actual']['receipt']))['parallelism']['nprocs']!=1:
            r.update(status='audit_failed',reason='actual rank differs',energy_hartree=None)
        rows.append(r)
    rows.extend({**t['actual_row'],'task_id':t['task_id'],'case_id':t['case_id'],'candidate':t['candidate'],
                 'reused':True,'reused_collection':t['actual_collection']} for t in m['reused'])
    rows.extend(m['missing'])
    result={'protocol_id':PROTOCOL,'manifest':record(mp),'stage':m['stage'],'rows':rows,'denominator':384,
            'confirmed':sum(r['status']=='confirmed_restart' for r in rows),'fresh_tasks':len(m['tasks']),
            'reused':len(m['reused']),'new_calls_in_collection':0}
    write_new(output,result);return {k:v for k,v in result.items() if k!='rows'}


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest='op',required=True)
    for op,fields in {'prepare':('inventory','output'),'validate':('manifest',),'execute':('manifest',),'collect':('manifest','output')}.items():
        s=sub.add_parser(op)
        for field in fields:s.add_argument('--'+field,required=True,type=Path)
        if op=='prepare':s.add_argument('--previous',type=Path)
    args=vars(p.parse_args());print(json.dumps(globals()[args.pop('op')](**args),indent=2))
