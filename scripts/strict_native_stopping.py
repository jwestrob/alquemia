"""Exactly20 native TolE1e-10 calls from two fixed seeds for ten real cells."""
from __future__ import annotations
import argparse,fcntl,json,os,re,shutil,time
from pathlib import Path
import numpy as np
from affordable_common import InvalidArtifact,HA_TO_KCAL,read_json,record,verify,write_new,xyz
from native_pool_continuation import recipe as original_recipe,collect_task
from run_orca_task_manifest import load_manifest_tasks,run_manifest
from union_adaptive import snapshot

PROTOCOL='native_GFN2_five_failures_Ca_partners_TolE1e10_two_seeds_v1'
CELLS=(('c5b120-pqq-la_model','adaptive_La','alpb'),('bbl57595.1-pqq-la_model','origin','vacuum'),
       ('p15279-pqq-la_model','adaptive_Ca','vacuum'),('p38539-pqq-la_model','adaptive_Ca','vacuum'),
       ('q4w6g0-pqq-la_model','origin','vacuum'))
SETTINGS={'TolE_hartree':1e-10,'MaxIter':500,'electronic_temperature_K':300,'native_mixer':True,
          'mpi_ranks':1,'workers':20,'component_tolerance_kcal_mol':.1,'maximum_calls':20}


def key(r):return tuple(r[k] for k in ('case_id','candidate','metal','medium'))
def tid(r):return '__'.join(key(r)+(r['seed_kind'],))
def recipe(charge,medium):return original_recipe(charge,1,medium).replace(' MaxIter 500\n',' MaxIter 500\n TolE 1e-10\n')


def sources(comparison):
    c=read_json(comparison);inv=read_json(verify(c['inventory']));stage2=read_json(verify(c['stage2']))
    failed={(r['case_id'],r['candidate'],r['medium']) for r in c['cells'] if r['energy_settling_pass'] is False}
    if failed!=set(CELLS) or any(r['metal']!='La' for r in c['cells'] if r['energy_settling_pass'] is False):
        raise InvalidArtifact('exact prior five failures required')
    originals={key(r):r for r in inv['rows']};last={key(r):r for r in stage2['rows']};rows=[]
    for cid,q,medium in CELLS:
        pair=[]
        for metal in ('Ca','La'):
            r=originals[cid,q,metal,medium];p=last[cid,q,metal,medium]
            if r['status']!='compatible_seed_pair_available' or p['status']!='confirmed_restart':raise InvalidArtifact('required seed unavailable')
            for pin in (r['gbw'],r['xtbw'],r['xyz'],p['seed_after']['preserved_gbw_after'],p['seed_after']['preserved_after']):verify(pin)
            pair.append(r)
            for seed_kind in ('cold','pass2'):
                seed={k:r[k] for k in ('gbw','xtbw')} if seed_kind=='cold' else {
                    'gbw':p['seed_after']['preserved_gbw_after'],'xtbw':p['seed_after']['preserved_after']}
                rows.append({'case_id':cid,'candidate':q,'metal':metal,'medium':medium,'seed_kind':seed_kind,
                    'source':r,'seed_source':seed,'old_pass2':p,'old_pass2_collection':c['stage2']})
        ca,la=pair;xc=xyz(verify(ca['xyz']));xl=xyz(verify(la['xyz']))
        if xc[0][0]!='Ca' or xl[0][0]!='La' or [r[0] for r in xc[1:]]!=[r[0] for r in xl[1:]] or la['charge']-ca['charge']!=1 or ca['multiplicity']!=la['multiplicity']:
            raise InvalidArtifact('Ca partner state/order mismatch')
        if np.max(np.abs(np.array([r[1:] for r in xc])-np.array([r[1:] for r in xl])))>1e-12:
            raise InvalidArtifact('Ca partner geometry differs')
    return rows


def prepare(comparison,agreement,output):
    rows=sources(comparison);out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False);tasks=[]
    for r in rows:
        src=r['source'];d=out/'tasks'/tid(r);d.mkdir(parents=True)
        shutil.copyfile(verify(src['xyz']),d/'core.xyz')
        for k in ('gbw','xtbw'):shutil.copyfile(verify(r['seed_source'][k]),d/('seed.immutable.'+k))
        (d/'endpoint.inp').write_text(recipe(src['charge'],r['medium']))
        tasks.append({**r,'task_id':tid(r),'charge':src['charge'],'multiplicity':src['multiplicity'],'gradient_requested':False,
            'xyz':record(d/'core.xyz'),'input':record(d/'endpoint.inp'),'output_path':str(d/'endpoint.out'),
            'immutable_seed':record(d/'seed.immutable.xtbw'),'immutable_gbw':record(d/'seed.immutable.gbw'),
            'active_seed_path':str(d/'endpoint.runtime.xtbw'),'active_gbw_path':str(d/'endpoint.runtime.gbw')})
    impl=snapshot(out/'implementation');m={'protocol_id':PROTOCOL,'settings':SETTINGS,'tasks':tasks,'comparison':record(comparison),
        'agreement':record(agreement),'orca':rows[0]['source']['orca'],'implementation':impl,
        'execution_resources':{'mpi_ranks':1,'concurrent_tasks':20},
        'execution_policy':{'task_runner':impl['run_orca_task_manifest.py'],'runtime_renderer':impl['render_orca_runtime_input.py']},
        'new_MACE_DFT_optimization_calls':0,'production_changed':False}
    mp=out/'manifest.json';write_new(mp,m);v=validate(mp,True);write_new(out/'PREFLIGHT.json',v);return v


def validate(manifest,fresh=False):
    mp=Path(manifest).resolve();m=read_json(mp);rows={tid(r):r for r in sources(verify(m['comparison']))}
    if m['protocol_id']!=PROTOCOL or m['settings']!=SETTINGS or m['execution_resources']!={'mpi_ranks':1,'concurrent_tasks':20}:
        raise InvalidArtifact('fixed protocol/profile differs')
    verify(m['agreement']);verify(m['orca'])
    for pin in m['implementation'].values():verify(pin)
    if len(m['tasks'])!=20 or {t['task_id'] for t in m['tasks']}!=set(rows):raise InvalidArtifact('twenty-task denominator differs')
    load_manifest_tasks(mp)
    for t in m['tasks']:
        r=rows[t['task_id']]
        if any(t[k]!=v for k,v in r.items()):raise InvalidArtifact('source/seed assignment differs')
        if t['charge']!=r['source']['charge'] or t['multiplicity']!=1 or t['gradient_requested']:raise InvalidArtifact('electronic state differs')
        if verify(t['xyz']).read_bytes()!=verify(r['source']['xyz']).read_bytes() or verify(t['input']).read_text()!=recipe(t['charge'],t['medium']):
            raise InvalidArtifact('coordinates or recipe differs')
        for field,k in (('immutable_seed','xtbw'),('immutable_gbw','gbw')):
            if verify(t[field]).read_bytes()!=verify(r['seed_source'][k]).read_bytes():raise InvalidArtifact('seed differs')
        d=Path(t['output_path']).parent
        if Path(t['active_seed_path'])!=d/'endpoint.runtime.xtbw' or Path(t['active_gbw_path'])!=d/'endpoint.runtime.gbw':raise InvalidArtifact('basename differs')
        if fresh and {p.name for p in d.iterdir()}!={'core.xyz','endpoint.inp','seed.immutable.gbw','seed.immutable.xtbw'}:raise InvalidArtifact('nonfresh task directory')
    return {'status':'validated','manifest':record(mp),'tasks':20,'new_calls_in_validation':0}


def execute(manifest):
    mp=Path(manifest).resolve();m=read_json(mp)
    if not os.environ.get('SLURM_JOB_ID') or int(os.environ.get('SLURM_CPUS_ON_NODE','0'))<20:raise InvalidArtifact('20-CPU allocation required')
    with (mp.parent/'execute.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB);validate(mp,True);before=[]
        for t in m['tasks']:
            shutil.copyfile(verify(t['immutable_seed']),t['active_seed_path']);shutil.copyfile(verify(t['immutable_gbw']),t['active_gbw_path'])
            before.append({'task_id':t['task_id'],'active_seed':record(t['active_seed_path']),'active_gbw':record(t['active_gbw_path'])})
        write_new(mp.parent/'SEEDS_BEFORE.json',{'manifest':record(mp),'rows':before})
        start=time.monotonic();error=None;result=None
        try:result=run_manifest(mp,orca_path=verify(m['orca']),workers=20,nprocs=1)
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
    return {'status':'executed','tasks':20}


def collect(manifest,output):
    validate(manifest);mp=Path(manifest).resolve();m=read_json(mp)
    def idx(name):
        p=mp.parent/name;return {r['task_id']:r for r in read_json(p)['rows']} if p.exists() else {}
    before=idx('SEEDS_BEFORE.json');after=idx('SEEDS_AFTER.json');rows=[]
    for t in m['tasks']:
        r=collect_task(mp,t,before.get(t['task_id']),after.get(t['task_id']));r['seed_kind']=t['seed_kind'];r['old_pass2']=t['old_pass2']
        r['observed_TolE_hartree']=None
        if r['status']=='confirmed_restart':
            text=verify(r['actual']['output']).read_text();tol=re.findall(r'Energy Change\s+TolE\s+\.{4}\s+([-+0-9.eE]+)',text)
            r['observed_TolE_hartree']=float(tol[0]) if len(tol)==1 else None
            if r['observed_TolE_hartree']!=1e-10 or read_json(verify(r['actual']['receipt']))['parallelism']['nprocs']!=1:
                r.update(status='audit_failed',reason='requested native tolerance/rank not confirmed',energy_hartree=None)
        rows.append(r)
    index={(key(r),r['seed_kind']):r for r in rows};cells=[];paired=[]
    for cid,q,s in CELLS:
        for z in ('Ca','La'):
            a=index[(cid,q,z,s),'cold'];b=index[(cid,q,z,s),'pass2'];ok=a['status']==b['status']=='confirmed_restart'
            d=(b['energy_hartree']-a['energy_hartree'])*HA_TO_KCAL if ok else None
            cells.append({'case_id':cid,'candidate':q,'metal':z,'medium':s,'status':'complete' if ok else 'unavailable',
                'from_cold_hartree':a['energy_hartree'],'from_pass2_hartree':b['energy_hartree'],
                'seed_difference_kcal_mol':d,'seed_agreement_pass':abs(d)<=.1 if ok else None,
                'cold_hartree':a['source']['source_energy_hartree'],'old_pass2_hartree':a['old_pass2']['energy_hartree']})
        values={}
        for start in ('cold','pass2'):
            a=index[(cid,q,'Ca',s),start];b=index[(cid,q,'La',s),start]
            values[start]=(a['energy_hartree']-b['energy_hartree'])*HA_TO_KCAL if a['status']==b['status']=='confirmed_restart' else None
        paired.append({'case_id':cid,'candidate':q,'medium':s,'Ca_minus_La_component_kcal_mol':values,
                       'complete_solvent_or_classifier_score':None})
    result={'protocol_id':PROTOCOL,'manifest':record(mp),'rows':rows,'cells':cells,'paired_components':paired,
        'denominator':20,'complete':sum(r['status']=='confirmed_restart' for r in rows),'two_start_cell_denominator':10,
        'two_start_cells_complete':sum(c['status']=='complete' for c in cells),'seed_agreement_passes':sum(c['seed_agreement_pass'] is True for c in cells),
        'new_molecular_calls_in_collection':0,'production_changed':False,'classifier_or_reference':None}
    write_new(output,result);return {k:v for k,v in result.items() if k!='rows'}


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest='op',required=True)
    for op,fields in {'prepare':('comparison','agreement','output'),'validate':('manifest',),'execute':('manifest',),'collect':('manifest','output')}.items():
        s=sub.add_parser(op)
        for k in fields:s.add_argument('--'+k,required=True,type=Path)
    a=vars(p.parse_args());print(json.dumps(globals()[a.pop('op')](**a),indent=2))
