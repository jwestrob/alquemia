"""Uniform32 native TolE1e-10 qualification: fresh versus original-cold restart."""
from __future__ import annotations
import argparse,fcntl,json,os,re,shutil,time
from pathlib import Path
from affordable_common import InvalidArtifact,read_json,record,verify,write_new
from precision_pool_continue_run import source_rows,tid
from native_pool_continuation import collect_task
from strict_native_stopping import recipe as strict_recipe, key
from compact_solvation import completed,diagnostics
from structure_informed_starts import scf_details
from union_adaptive import snapshot

PROTOCOL='native_GFN2_precision32_TolE1e10_fresh_vs_cold_seed_v1'
BRANCHES=('fresh','cold_seed')
SETTINGS={'TolE_hartree':1e-10,'MaxIter':500,'electronic_temperature_K':300,'native_mixer':True,
          'mpi_ranks':1,'workers':32,'component_tolerance_kcal_mol':.1,'contrast_tolerance_kcal_mol':.2,
          'pool_tolerance_kcal_mol':.2,'origin_selection_tolerance_kcal_mol':.1,'logical_call_denominator':768}


def recipe(charge,medium,branch):
    s=strict_recipe(charge,medium)
    return s.replace('! Native-GFN2-xTB','! Native-GFN2-xTB NoAutostart') if branch=='fresh' else s


def reuse_rows(path,sources):
    c=read_json(path);lookup={tid(r):r for r in sources};out={}
    for r in c['rows']:
        if r['seed_kind']!='cold':continue
        src=lookup[tid(r)]
        if (r['status']!='confirmed_restart' or r['observed_TolE_hartree']!=1e-10 or
            r['source']!=src or r['seed_source']!={k:src[k] for k in ('gbw','xtbw')}):
            raise InvalidArtifact('exact strict-from-original-cold reuse unavailable')
        for name in ('receipt','output'):verify(r['actual'][name])
        body=read_json(verify(r['actual']['manifest']))
        task=next(t for t in body['tasks'] if t['task_id']==r['actual']['task_id'])
        if verify(task['input']).read_text()!=recipe(src['charge'],src['medium'],'cold_seed') or verify(task['xyz']).read_bytes()!=verify(src['xyz']).read_bytes():
            raise InvalidArtifact('reuse recipe/geometry mismatch')
        out[tid(src)]=r
    if len(out)!=10:raise InvalidArtifact('ten exact strict20 cold-seed reuses required')
    return out


def prepare(inventory,strict20,agreement,branch,output):
    if branch not in BRANCHES:raise InvalidArtifact('unknown branch')
    inv,rows=source_rows(inventory);reuse=reuse_rows(strict20,rows) if branch=='cold_seed' else {}
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False);tasks=[];missing=[];reused=[]
    for r in rows:
        identity={k:r[k] for k in ('case_id','candidate','metal','medium')};identity['task_id']=tid(r)
        if r['status']!='compatible_seed_pair_available':
            missing.append({**identity,'status':'unavailable','reason':r['reason'],'energy_hartree':None});continue
        if tid(r) in reuse:
            reused.append({**identity,'source':r,'actual_row':reuse[tid(r)],'actual_collection':record(strict20)});continue
        d=out/'tasks'/tid(r);d.mkdir(parents=True);shutil.copyfile(verify(r['xyz']),d/'core.xyz')
        (d/'endpoint.inp').write_text(recipe(r['charge'],r['medium'],branch))
        t={**identity,'source':r,'charge':r['charge'],'multiplicity':r['multiplicity'],'gradient_requested':False,
           'xyz':record(d/'core.xyz'),'input':record(d/'endpoint.inp'),'output_path':str(d/'endpoint.out'),
           'active_seed_path':str(d/'endpoint.runtime.xtbw'),'active_gbw_path':str(d/'endpoint.runtime.gbw'),
           'seed_source':None}
        if branch=='cold_seed':
            t['seed_source']={k:r[k] for k in ('gbw','xtbw')}
            for k in ('gbw','xtbw'):shutil.copyfile(verify(r[k]),d/('seed.immutable.'+k))
            t.update(immutable_seed=record(d/'seed.immutable.xtbw'),immutable_gbw=record(d/'seed.immutable.gbw'))
        tasks.append(t)
    impl=snapshot(out/'implementation');m={'protocol_id':PROTOCOL,'branch':branch,'settings':SETTINGS,
        'inventory':record(inventory),'strict20':record(strict20),'agreement':record(agreement),'tasks':tasks,'missing':missing,'reused':reused,
        'orca':rows[0]['orca'],'implementation':impl,'case_denominator':32,'cell_denominator':384,
        'execution_resources':{'mpi_ranks':1,'concurrent_tasks':32},
        'execution_policy':{'task_runner':impl['run_orca_task_manifest.py'],'runtime_renderer':impl['render_orca_runtime_input.py']},
        'production_changed':False,'new_MACE_DFT_optimization_calls':0}
    mp=out/'manifest.json';write_new(mp,m);v=validate(mp,True);write_new(out/'PREFLIGHT.json',v);return v


def validate(manifest,fresh=False):
    mp=Path(manifest).resolve();m=read_json(mp);inv,rows=source_rows(verify(m['inventory']));src={tid(r):r for r in rows}
    if m['protocol_id']!=PROTOCOL or m['settings']!=SETTINGS or m['branch'] not in BRANCHES or m['execution_resources']!={'mpi_ranks':1,'concurrent_tasks':32}:
        raise InvalidArtifact('protocol/settings differ')
    verify(m['agreement']);verify(m['orca']);verify(m['strict20'])
    for pin in m['implementation'].values():verify(pin)
    reuse=reuse_rows(verify(m['strict20']),rows) if m['branch']=='cold_seed' else {}
    allrows=m['tasks']+m['missing']+m['reused']
    if len(allrows)!=384 or {r['task_id'] for r in allrows}!=set(src):raise InvalidArtifact('384-cell denominator differs')
    if {r['task_id'] for r in m['reused']}!=set(reuse):raise InvalidArtifact('reuse membership differs')
    if {r['task_id'] for r in m['missing']}!={tid(r) for r in rows if r['status']!='compatible_seed_pair_available'}:raise InvalidArtifact('unavailable denominator differs')
    if m['tasks']:load=__import__('run_orca_task_manifest').load_manifest_tasks(mp)
    for r in m['reused']:
        if r['source']!=src[r['task_id']] or r['actual_row']!=reuse[r['task_id']] or r['actual_collection']!=m['strict20']:raise InvalidArtifact('reuse changed')
    for t in m['tasks']:
        r=src[t['task_id']];branch=m['branch'];d=Path(t['output_path']).parent
        if t['source']!=r or any(t[k]!=r[k] for k in ('case_id','candidate','metal','medium','charge','multiplicity')) or t['multiplicity']!=1:
            raise InvalidArtifact('physical/electronic state changed')
        if t['gradient_requested'] or verify(t['xyz']).read_bytes()!=verify(r['xyz']).read_bytes() or verify(t['input']).read_text()!=recipe(r['charge'],r['medium'],branch):raise InvalidArtifact('recipe/coordinate changed')
        if Path(t['active_seed_path'])!=d/'endpoint.runtime.xtbw' or Path(t['active_gbw_path'])!=d/'endpoint.runtime.gbw':raise InvalidArtifact('runtime basename differs')
        expected={'core.xyz','endpoint.inp'}
        if branch=='cold_seed':
            if t['seed_source']!={k:r[k] for k in ('gbw','xtbw')}:raise InvalidArtifact('wrong original seed')
            for field,k in (('immutable_seed','xtbw'),('immutable_gbw','gbw')):
                if verify(t[field]).read_bytes()!=verify(r[k]).read_bytes():raise InvalidArtifact('seed bytes differ')
            expected|={'seed.immutable.gbw','seed.immutable.xtbw'}
        elif t['seed_source'] is not None or 'immutable_seed' in t or 'immutable_gbw' in t:raise InvalidArtifact('fresh branch has a seed')
        if fresh and {p.name for p in d.iterdir()}!=expected:raise InvalidArtifact('nonfresh task directory')
    return {'status':'validated','manifest':record(mp),'branch':m['branch'],'fresh_tasks':len(m['tasks']),
            'reused':len(m['reused']),'missing':len(m['missing']),'cell_denominator':384,'new_calls_in_validation':0}


def execute(manifest):
    from run_orca_task_manifest import run_manifest
    mp=Path(manifest).resolve();m=read_json(mp);seeded=m['branch']=='cold_seed'
    if not os.environ.get('SLURM_JOB_ID') or int(os.environ.get('SLURM_CPUS_ON_NODE','0'))<32:raise InvalidArtifact('32-CPU allocation required')
    with (mp.parent/'execute.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB);validate(mp,True);before=[]
        for t in m['tasks']:
            if seeded:
                shutil.copyfile(verify(t['immutable_seed']),t['active_seed_path']);shutil.copyfile(verify(t['immutable_gbw']),t['active_gbw_path'])
            before.append({'task_id':t['task_id'],'active_seed':record(t['active_seed_path']) if seeded else None,
                'active_gbw':record(t['active_gbw_path']) if seeded else None,'seeded':seeded})
        write_new(mp.parent/'SEEDS_BEFORE.json',{'manifest':record(mp),'rows':before})
        start=time.monotonic();error=None;results=None
        try:results=run_manifest(mp,orca_path=verify(m['orca']),workers=32,nprocs=1) if m['tasks'] else []
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
            write_new(mp.parent/'EXECUTION.json',{'manifest':record(mp),'results':results,'error':error,
                'wall_seconds':time.monotonic()-start,'allocated_cpus':int(os.environ['SLURM_CPUS_ON_NODE']),'slurm_job_id':os.environ['SLURM_JOB_ID']})
        if error:raise InvalidArtifact(error)
    return {'status':'executed','tasks':len(m['tasks'])}


def collect_fresh(mp,t,before,after):
    r={k:t[k] for k in ('task_id','case_id','candidate','metal','medium')};r.update(status='unavailable',energy_hartree=None,source=t['source'],seed_source=None)
    try:
        pin=completed(mp,t['task_id'])
        if not pin:
            r['reason']='execution_failed_or_not_complete';r['artifacts']=[record(p) for p in (Path(t['output_path']),Path(t['output_path']+'.execution.json')) if p.exists()];return r
        text=verify(pin['output']).read_text();detail=scf_details(text);audit=diagnostics(pin,t)
        r.update(actual=pin,observed_energy_hartree=pin['energy_hartree'],details=detail,audit=audit,seed_before=before,seed_after=after)
        guess=[v.strip() for v in re.findall(r'INITIAL GUESS:\s*([^\r\n]+)',text)];r['initial_guess']=guess
        if not before or before['active_seed'] is not None or before['active_gbw'] is not None or before['seeded'] or 'NoAutostart' not in verify(t['input']).read_text() or guess!=['SAD']:raise InvalidArtifact('fresh SAD initialization not confirmed')
        source=t['source']
        if not detail['native_mixer_observed'] or (detail['charge'],detail['multiplicity'],detail['electrons'])!=(t['charge'],t['multiplicity'],source['electron_count']) or audit['parameter_export']['sha256']!=source['parameter_export']['sha256'] or audit['charge_sanity_status']!='pass':raise InvalidArtifact('native method/state/parameter/charge audit failed')
        r.update(status='complete',energy_hartree=pin['energy_hartree'])
    except Exception as exc:r.update(status='audit_failed',reason=str(exc))
    return r


def collect(manifest,output):
    validate(manifest);mp=Path(manifest).resolve();m=read_json(mp)
    def idx(name):
        p=mp.parent/name;return {r['task_id']:r for r in read_json(p)['rows']} if p.exists() else {}
    before=idx('SEEDS_BEFORE.json');after=idx('SEEDS_AFTER.json');rows=[]
    for t in m['tasks']:
        r=(collect_task if m['branch']=='cold_seed' else collect_fresh)(mp,t,before.get(t['task_id']),after.get(t['task_id']))
        r['reused']=False;r['observed_TolE_hartree']=None
        if r['status'] in ('confirmed_restart','complete'):
            text=verify(r['actual']['output']).read_text();tol=re.findall(r'Energy Change\s+TolE\s+\.{4}\s+([-+0-9.eE]+)',text)
            r['observed_TolE_hartree']=float(tol[0]) if len(tol)==1 else None
            if r['observed_TolE_hartree']!=1e-10 or read_json(verify(r['actual']['receipt']))['parallelism']['nprocs']!=1:r.update(status='audit_failed',reason='effective tolerance/rank not confirmed',energy_hartree=None)
            else:r['status']='complete'
        rows.append(r)
    rows.extend({**t['actual_row'],'task_id':t['task_id'],'status':'complete','reused':True,'reused_collection':t['actual_collection']} for t in m['reused'])
    rows.extend(m['missing'])
    result={'protocol_id':PROTOCOL,'manifest':record(mp),'branch':m['branch'],'rows':rows,'denominator':384,
        'complete':sum(r['status']=='complete' for r in rows),'fresh_tasks':len(m['tasks']),'reused':len(m['reused']),'new_calls_in_collection':0}
    write_new(output,result);return {k:v for k,v in result.items() if k!='rows'}


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest='op',required=True)
    for op,fields in {'prepare':('inventory','strict20','agreement','output'),'validate':('manifest',),'execute':('manifest',),'collect':('manifest','output')}.items():
        s=sub.add_parser(op)
        for k in fields:s.add_argument('--'+k,required=True,type=Path)
        if op=='prepare':s.add_argument('--branch',required=True,choices=BRANCHES)
    a=vars(p.parse_args());print(json.dumps(globals()[a.pop('op')](**a),indent=2))
