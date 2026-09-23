"""Exactly two same-state seeded attempts for the declared failed envelope cell."""
from __future__ import annotations
import argparse,fcntl,json,os,re,shutil,time
from pathlib import Path
from affordable_common import InvalidArtifact,HA_TO_KCAL,read_json,record,verify,write_new,xyz
from compact_solvation import completed,diagnostics
from native_pool_continuation import collect_task
from strict_native_pool import recipe
from structure_informed_starts import scf_details
from run_orca_task_manifest import load_manifest_tasks,run_manifest
from union_adaptive import snapshot

PROTOCOL='Nikasha_single_failed_native_cell_two_seed_recovery_sensitivity_v1'
CASE='a0acd6b9f2-pqq-la_model__conditioned_La__seed-1_sample-0__envelope_2eda778295ffed58'
SEEDS=('origin','adaptive_La')
SETTINGS={'TolE_hartree':1e-10,'MaxIter':500,'temperature_K':300,'mpi_ranks':1,'workers':2,
          'maximum_calls':2,'component_agreement_kcal_mol':.1,'reported_seed':'origin'}
BASE=Path('workspaces/motion_envelope_transfer_20260923')


def archived_task(path,candidate):
    m=read_json(path)
    ts=[t for t in m['tasks'] if (t['case_id'],t['candidate'],t['metal'],t['medium'])==(CASE,candidate,'La','vacuum')]
    if len(ts)!=1:raise InvalidArtifact('exact archived cell missing/nonunique')
    t=ts[0];rp=Path(t['output_path']+'.execution.json');r=read_json(rp);op=Path(t['output_path'])
    if (t['charge'],t['multiplicity'])!=(-3,1) or verify(t['input']).read_text()!=recipe(-3,'vacuum','fresh'):
        raise InvalidArtifact('actual strict electronic recipe differs')
    for key,actual in [('template_input',t['input']),('xyz',t['xyz']),('output',record(op))]:
        if r['artifacts'][key]['sha256']!=actual['sha256']:raise InvalidArtifact('actual source receipt/input mismatch')
    if r['orca_version']!='6.1.1' or r['parallelism']['nprocs']!=1:raise InvalidArtifact('actual backend/rank mismatch')
    # Origin receipt retains the documented byte-identical runner-path recovery;
    # require its original runner/renderer hashes, not a newly invented execution.
    for key,policy in [('task_runner','task_runner'),('runtime_renderer','runtime_renderer')]:
        if r[key]['sha256']!=m['execution_policy'][policy]['sha256']:raise InvalidArtifact('runner bytes differ')
        verify(m['execution_policy'][policy])
    return m,t,r,{'manifest':record(path),'task_id':t['task_id'],'output':record(op),'receipt':record(rp)}


def inventory(output):
    cp=BASE/'run_v2/shard_0/pool/solvent/shard_0/manifest.json'
    om=BASE/'run_v1/shard_0/origins/scalar/manifest.json'
    m,t,r,pin=archived_task(cp,'adaptive_Ca');text=verify(pin['output']).read_text()
    if r['normal_termination'] or r['scf_converged'] or completed(cp,t['task_id']) is not None:
        raise InvalidArtifact('declared target is not the real failed cell')
    if '500' not in text or 'LEANSCF' not in text:raise InvalidArtifact('declared SCF failure not found')
    atoms=xyz(verify(t['xyz']));order=[a[0] for a in atoms]
    target={'task':t,'actual':pin,'atom_count':len(atoms),'element_order':order,'orca':m['orca'],
            'failure_receipt':r,'parameter_export':record(Path(t['output_path']).parent/'endpoint.runtime.xtb.json')}
    rows=[]
    for seed,mp in [('origin',om),('adaptive_La',cp)]:
        sm,st,sr,sp=archived_task(mp,seed)
        if not(sr['returncode']==0 and sr['normal_termination'] and sr['scf_converged']):raise InvalidArtifact('seed calculation failed')
        txt=verify(sp['output']).read_text();detail=scf_details(txt)
        es=re.findall(r'FINAL SINGLE POINT ENERGY\s+([-+0-9.]+)',txt)
        if len(es)!=1:raise InvalidArtifact('seed energy unavailable')
        sp['energy_hartree']=float(es[0]);audit=diagnostics(sp,st)
        if ([a[0]for a in xyz(verify(st['xyz']))]!=order or sm['orca']!=m['orca'] or
            audit['parameter_export']['sha256']!=target['parameter_export']['sha256'] or
            (detail['charge'],detail['multiplicity'],detail['electrons'])!=(-3,1,636) or
            detail['energy']['tolerance']!=1e-10 or not detail['native_mixer_observed']):
            raise InvalidArtifact('seed state/ordering/parameters mismatch')
        basis=set(map(int,re.findall(r'Number of basis functions\s+\.\.\.\s+(\d+)',txt)))
        shells=set(map(int,re.findall(r'Number of shells\s+\.\.\.\s+(\d+)',txt)))
        if basis!={573} or shells!={325} or len(atoms)!=204:raise InvalidArtifact('actual basis shape differs')
        d=Path(st['output_path']).parent;seeds={k:record(d/('endpoint.runtime.'+k))for k in ('gbw','xtbw')}
        sizes={k:verify(v).stat().st_size for k,v in seeds.items()}
        if sizes!={'gbw':3855516,'xtbw':17368}:raise InvalidArtifact('seed file shape differs')
        rows.append({'seed_kind':seed,'task':st,'actual':sp,'seed_source':seeds,'bytes':sizes,
            'atom_count':204,'basis_functions':573,'shells':325,'electron_count':636,
            'parameter_export':audit['parameter_export'],'SCF':{k:v for k,v in detail.items()if k!='printed_orbital_occupations'},
            'context_mapping':t['source_mapping'],'context_preparation':t['source_preparation']})
    verify(t['source_preparation']);verify(t['source_mapping'])
    if rows[0]['task']['source_xyz']['sha256']!=record(Path(t['source_preparation']['path']).parent/'La_context.xyz')['sha256']:
        raise InvalidArtifact('origin is not exact same prepared context')
    if rows[1]['task']['source_mapping']!=t['source_mapping'] or rows[1]['task']['source_preparation']!=t['source_preparation']:
        raise InvalidArtifact('candidate context identity differs')
    d={'protocol_id':PROTOCOL,'target':target,'seeds':rows,'seed_count':2,'new_molecular_calls':0}
    write_new(output,d);return {'inventory':record(output),'compatible_seeds':2,'atoms':204}


def prepare(inventory,agreement,output):
    inv=read_json(inventory);target=inv['target'];out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False);tasks=[]
    for row in inv['seeds']:
        kind=row['seed_kind'];d=out/'tasks'/kind;d.mkdir(parents=True)
        shutil.copyfile(verify(target['task']['xyz']),d/'core.xyz')
        for k in ('gbw','xtbw'):shutil.copyfile(verify(row['seed_source'][k]),d/('seed.immutable.'+k))
        (d/'endpoint.inp').write_text(recipe(-3,'vacuum','cold_seed'))
        src={'electron_count':636,'parameter_export':target['parameter_export'],
             'source_energy_hartree':row['actual']['energy_hartree']}
        tasks.append({'task_id':kind,'case_id':CASE,'candidate':'adaptive_Ca','metal':'La','medium':'vacuum',
            'seed_kind':kind,'charge':-3,'multiplicity':1,'gradient_requested':False,'source':src,'seed_source':row['seed_source'],
            'xyz':record(d/'core.xyz'),'input':record(d/'endpoint.inp'),'output_path':str(d/'endpoint.out'),
            'immutable_seed':record(d/'seed.immutable.xtbw'),'immutable_gbw':record(d/'seed.immutable.gbw'),
            'active_seed_path':str(d/'endpoint.runtime.xtbw'),'active_gbw_path':str(d/'endpoint.runtime.gbw')})
    impl=snapshot(out/'implementation');m={'protocol_id':PROTOCOL,'settings':SETTINGS,'inventory':record(inventory),
        'agreement':record(agreement),'tasks':tasks,'orca':target['orca'],'implementation':impl,
        'execution_resources':{'mpi_ranks':1,'concurrent_tasks':2},'execution_policy':{
        'task_runner':impl['run_orca_task_manifest.py'],'runtime_renderer':impl['render_orca_runtime_input.py']},
        'primary_result_changed':False,'new_MACE_DFT_geometry_calls':0}
    mp=out/'manifest.json';write_new(mp,m);v=validate(mp,True);write_new(out/'PREFLIGHT.json',v);return v


def validate(manifest,fresh=False):
    m=read_json(manifest);inv=read_json(verify(m['inventory']));target=inv['target']['task']
    if m['protocol_id']!=PROTOCOL or m['settings']!=SETTINGS or m['execution_resources']!={'mpi_ranks':1,'concurrent_tasks':2}:
        raise InvalidArtifact('fixed two-attempt protocol differs')
    if [t['task_id']for t in m['tasks']]!=list(SEEDS):raise InvalidArtifact('exact ordered starts required')
    verify(m['agreement']);verify(m['orca'])
    for p in m['implementation'].values():verify(p)
    load_manifest_tasks(Path(manifest))
    for t,r in zip(m['tasks'],inv['seeds']):
        if (t['case_id'],t['candidate'],t['metal'],t['medium'],t['charge'],t['multiplicity'],t['gradient_requested'])!=(CASE,'adaptive_Ca','La','vacuum',-3,1,False):raise InvalidArtifact('state differs')
        if verify(t['xyz']).read_bytes()!=verify(target['xyz']).read_bytes() or verify(t['input']).read_text()!=recipe(-3,'vacuum','cold_seed'):
            raise InvalidArtifact('failed geometry or strict recipe changed')
        expected_source={'electron_count':636,'parameter_export':inv['target']['parameter_export'],'source_energy_hartree':r['actual']['energy_hartree']}
        if t['source']!=expected_source:raise InvalidArtifact('source audit identity differs')
        if t['seed_kind']!=r['seed_kind'] or t['seed_source']!=r['seed_source']:raise InvalidArtifact('seed selection differs')
        for field,k in [('immutable_seed','xtbw'),('immutable_gbw','gbw')]:
            if verify(t[field]).read_bytes()!=verify(r['seed_source'][k]).read_bytes():raise InvalidArtifact('seed bytes changed')
        d=Path(t['output_path']).parent
        if Path(t['active_seed_path'])!=d/'endpoint.runtime.xtbw' or Path(t['active_gbw_path'])!=d/'endpoint.runtime.gbw':raise InvalidArtifact('runtime basename differs')
        if fresh and {p.name for p in d.iterdir()}!={'core.xyz','endpoint.inp','seed.immutable.gbw','seed.immutable.xtbw'}:raise InvalidArtifact('nonfresh attempt directory')
    return {'status':'validated','manifest':record(manifest),'tasks':2,'new_molecular_calls_in_validation':0}


def execute(manifest):
    mp=Path(manifest).resolve();m=read_json(mp)
    if not os.environ.get('SLURM_JOB_ID') or int(os.environ.get('SLURM_CPUS_ON_NODE','0'))!=2:raise InvalidArtifact('exact2CPU allocation required')
    with (mp.parent/'execute.lock').open('a')as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB);validate(mp,True);before=[]
        for t in m['tasks']:
            shutil.copyfile(verify(t['immutable_seed']),t['active_seed_path']);shutil.copyfile(verify(t['immutable_gbw']),t['active_gbw_path'])
            before.append({'task_id':t['task_id'],'active_seed':record(t['active_seed_path']),'active_gbw':record(t['active_gbw_path'])})
        write_new(mp.parent/'SEEDS_BEFORE.json',{'manifest':record(mp),'rows':before});start=time.monotonic();error=None;results=None
        try:results=run_manifest(mp,orca_path=verify(m['orca']),workers=2,nprocs=1)
        except Exception as exc:error=repr(exc)
        finally:
            after=[]
            for t in m['tasks']:
                d=Path(t['output_path']).parent;r={'task_id':t['task_id']}
                for active,k,field in [('active_seed_path','xtbw','preserved_after'),('active_gbw_path','gbw','preserved_gbw_after')]:
                    p=Path(t[active]);dest=d/('seed.after.'+k)
                    if p.exists():shutil.copyfile(p,dest)
                    r[field]=record(dest)if dest.exists()else None
                after.append(r)
            write_new(mp.parent/'SEEDS_AFTER.json',{'manifest':record(mp),'rows':after})
            write_new(mp.parent/'EXECUTION.json',{'manifest':record(mp),'results':results,'error':error,'wall_seconds':time.monotonic()-start,'slurm_job_id':os.environ['SLURM_JOB_ID'],'allocated_cpus':2})
        if error:raise InvalidArtifact(error)
    return {'status':'executed','tasks':2}


def decision(rows):
    complete=len(rows)==2 and all(r['status']=='confirmed_restart'for r in rows)
    diff=(rows[1]['energy_hartree']-rows[0]['energy_hartree'])*HA_TO_KCAL if complete else None
    passed=complete and abs(diff)<=SETTINGS['component_agreement_kcal_mol']
    return {'status':'qualified_recovery_sensitivity'if passed else 'recovery_unavailable',
        'seed_difference_kcal_mol':diff,'agreement_pass':passed,
        'recovery_energy_hartree':rows[0]['energy_hartree']if passed else None,
        'reported_seed':'origin','primary_status':'unavailable_unchanged'}


def collect(manifest,output):
    validate(manifest);mp=Path(manifest).resolve();m=read_json(mp)
    def idx(name):
        p=mp.parent/name;return {r['task_id']:r for r in read_json(p)['rows']}if p.exists()else {}
    before=idx('SEEDS_BEFORE.json');after=idx('SEEDS_AFTER.json');rows=[]
    for t in m['tasks']:
        r=collect_task(mp,t,before.get(t['task_id']),after.get(t['task_id']));r['seed_kind']=t['seed_kind']
        # Seed geometry differs: its energy is not a prior energy at this target.
        r.pop('change_from_archive_kcal_mol',None)
        if r['status']=='confirmed_restart':
            txt=verify(r['actual']['output']).read_text();tol=re.findall(r'Energy Change\s+TolE\s+\.{4}\s+([-+0-9.eE]+)',txt)
            if len(tol)!=1 or float(tol[0])!=1e-10 or read_json(verify(r['actual']['receipt']))['parallelism']['nprocs']!=1:
                r.update(status='audit_failed',energy_hartree=None,reason='strict tolerance/rank absent')
        rows.append(r)
    d={'protocol_id':PROTOCOL,'manifest':record(mp),'rows':rows,'denominator':2,**decision(rows),
       'primary_result_changed':False,'new_molecular_calls_in_collection':0}
    write_new(output,d);return {k:v for k,v in d.items()if k!='rows'}

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);sp=p.add_subparsers(dest='op',required=True)
    for name,fields in {'inventory':('output',),'prepare':('inventory','agreement','output'),'validate':('manifest',),'execute':('manifest',),'collect':('manifest','output')}.items():
        s=sp.add_parser(name)
        for f in fields:s.add_argument('--'+f,required=True,type=Path)
    args=vars(p.parse_args());print(json.dumps(globals()[args.pop('op')](**args),indent=2))
