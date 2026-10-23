#!/usr/bin/env python3
"""Four conditional Cartesian checks; no fit outcomes opened during preparation."""
import argparse,concurrent.futures,copy,hashlib,importlib.util,json,os,resource,shutil,threading,time
from pathlib import Path
H=Path(__file__).resolve().parent;E=H.parents[1];ROOT=E.parents[2];W=ROOT/'workspaces/lanm_selectivity_free_energy_20261021/explicit_model_v1'
CELLS=('aquo8_central_La','aquo8_central_Dy','aquo9_central_La','aquo9_central_Dy')
TEN=tuple(q+'_'+m for q in ('water_minus','water_central','water_plus','aquo8_central','aquo9_central') for m in ('La','Dy'))
ATOMS={g:25 if 'aquo8' in g else 28 for g in CELLS};THREADS=('OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS')
ARGS=['system.xyz','Y','Y','0.00001']

def require(ok,message):
    if not ok:raise ValueError(message)
def pin(path):
    p=Path(path).resolve();return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def verify(p):
    require(pin(p['path'])==p,'Changed pin: '+p['path']);return Path(p['path'])
def read(p):return json.loads(verify(p).read_text())
def module(p,name):
    s=importlib.util.spec_from_file_location(name,verify(p));m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
def save(path,obj):
    p=Path(path);t=p.with_name(p.name+'.tmp');t.write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n');os.replace(t,p)

def default_config():
    fp=pin(W/'paired_water_coordination_joint_v1/run_v1/CONFIG_RUN.json');f=read(fp);ec=read(f['evaluator_config'])
    return {'implementation':pin(__file__),'declaration':pin(H/'DECLARATION.md'),'fit_config':fp,'fit_source':f['implementation'],'fit_evaluator':f['evaluator'],
      'fit_training':f['training_result'],'reviewer_source':pin(H.parent/'selected_review_v1/REVIEW.py'),
      'native_helper':ec['native_helper'],'numerical_helper':pin(E/'paired_water_transfer_v1/RUN.py'),'executables':ec['executables'],
      'testgrad_source':pin(W/'tinker_source/source/testgrad.f'),'previous_transfer':pin(W/'paired_water_six_transfer_v1/run_v1/execution/RESULT.json'),
      'candidate_result':None,'candidate_review':None,'root_release':None,'maximum_native_programs':4,'maximum_internal_FD_energies':636,'maximum_workers':4}

def metadata(cfg):
    base=copy.deepcopy(cfg)
    for k in ('candidate_result','candidate_review','root_release'):base[k]=None
    require(base==default_config(),'Frozen configuration differs')
    for k in ('declaration','fit_source','fit_evaluator','fit_training','reviewer_source','native_helper','numerical_helper','testgrad_source','previous_transfer'):verify(cfg[k])
    for p in cfg['executables'].values():verify(p)
    require(cfg['testgrad_source']['sha256']=='5b3523b67e95441bbc0d45582b0055c263cbda13fd0a8b03dc81d8ce33924333','Derivative source changed')
    require(sum(6*n for n in ATOMS.values())==636,'Internal FD denominator differs')

def candidate_gate(cfg):
    require(cfg['candidate_result'] is not None and cfg['candidate_review'] is not None,'Terminal candidate/review intentionally unbound')
    require(Path(cfg['candidate_result']['path'])==W/'paired_water_coordination_joint_v1/run_v1/execution/RESULT.json','Wrong candidate run')
    r=read(cfg['candidate_result']);v=read(cfg['candidate_review']);selected=r['selected']
    require(r['status']=='bounded_joint_coordination_water_complete' and r['job_id']=='17676' and r['error'] is None and r['candidate_frozen'] and r['development_checks_passed'],'Passing terminal17676 required')
    require(r['implementation']==cfg['fit_source'] and r['config']==cfg['fit_config'],'Fitter source/config differs')
    score=selected['score'];energy=score['energy_controls'];force=score['per_fixture']
    require(sum(len(x) for x in energy.values())==10 and sum(len(x) for x in force.values())==15 and score['development_checks_passed'],'All25 controls required')
    require(all(x['passed'] and abs(x['error_kcal_mol'])<=1 for rows in energy.values() for x in rows.values()) and all(x['passed'] and x['rigid_RMS_kcal_mol_A']<=2 and x['maximum_body_component_error_kcal_mol_A']<=5 for rows in force.values() for x in rows.values()),'A development control failed')
    require(v['status']=='selected_ten_raw_review_passed' and v['result']==cfg['candidate_result'] and v['selected_parameters']==selected['parameters'] and v['implementation']==cfg['reviewer_source'],'Independent selected-ten authority differs')
    require(v['complete']==10 and set(v['selected_cell_ids'])==set(TEN) and v['new_native_calls']==0 and v['error'] is None and all(v[k] is True for k in ('all_selected_raw_replays_exact','exact_selected_score_replayed','all_development_controls_passed')),'Independent review did not fully pass')
    require(selected['status']=='complete' and set(selected['model'])==set(TEN) and len(selected['parameters'])==6,'Exact selected ten required')
    return r

def compare(check,value):
    require(check['fd_passed'] and check['all_atoms_active_and_components_present'],'Cartesian numerical derivative failed')
    de=abs(check['energy']-value['energy']);df=max(abs(x-y) for a,b in zip(check['forces'],value['forces']) for x,y in zip(a,b))
    require(len(check['forces'])==len(value['forces']) and de<=1e-4 and df<=1e-3,'Selected saved E/F replay failed')
    return {'passed':True,'energy_replay_error_kcal_mol':de,'analytic_force_replay_max_error_kcal_mol_A':df,'Cartesian_max_error_kcal_mol_A':check['max_component_error_kcal_mol_A']}

def bind(config_path,out):
    cfg=json.loads(Path(config_path).read_text());metadata(cfg);r=candidate_gate(cfg);data=read(cfg['fit_training']);old=read(cfg['previous_transfer'])
    require(old['status']=='consumed_six_parameter_transfer_complete' and old['job_id']=='17612' and old['error'] is None and old['stage_A_passed'],'Original Cartesian framework unavailable')
    native=module(cfg['native_helper'],'cartesian_bind_native');helper=module(cfg['numerical_helper'],'cartesian_bind_parser');rows={}
    for gid in CELLS:
        row=data['rows'][gid];value=r['selected']['model'][gid];pv=value['provenance'];inp=pv['effective_input'];require(row['geometry']['atoms']==ATOMS[gid],'Consumed atom count differs')
        require(pv['requested_vector']==r['selected']['parameters'] and pv['request_cache_implementation']==cfg['fit_evaluator'],'Selected input context differs')
        verify(pv['native_authority'])
        for k in ('analyze_output','gradient_output'):verify(value['native'][k])
        directory=out/'inputs'/gid;directory.mkdir(parents=True)
        for k,n in (('xyz','system.xyz'),('key','system.key'),('parameters','amoeba09.prm')):shutil.copyfile(verify(inp[k]),directory/n)
        files={k:pin(directory/n) for k,n in (('xyz','system.xyz'),('key','system.key'),('parameters','amoeba09.prm'))}
        require(all(files[k]['sha256']==inp[k]['sha256'] for k in files) and all(files[k]['sha256']==row['input'][k]['sha256'] for k in ('xyz','parameters')),'Selected input copy differs')
        prior=None
        if gid.startswith('aquo8'):
            a=old['rows'][gid]
            if a['status']=='complete' and all(a['input']['files'][k]['sha256']==files[k]['sha256'] for k in files):
                check=helper.gradient(verify(a['cartesian_check']['output']),ATOMS[gid],True,native);require(check==a['cartesian_check'],'Old Cartesian record differs')
                admission=compare(check,value);prior={'authority':cfg['previous_transfer'],'check':check,'admission':admission}
        rows[gid]={'geometry':row['geometry'],'input':dict(row['input'],**files,derivative_test=True),'files':files,'selected_model':value,'native_authority':pv['native_authority'],'reuse':prior,'internal_FD_energies_if_new':6*ATOMS[gid]}
    bound={'status':'cartesian_bound_no_calls','implementation':pin(__file__),'config':pin(config_path),'candidate_result':cfg['candidate_result'],'candidate_review':cfg['candidate_review'],
      'parameters':r['selected']['parameters'],'rows':rows,'maximum_native_programs':4,'maximum_internal_FD_energies':636,
      'planned_new_programs':sum(x['reuse'] is None for x in rows.values()),'planned_new_internal_FD_energies':sum(x['internal_FD_energies_if_new'] for x in rows.values() if x['reuse'] is None),
      'native_calls':0,'new_geometry':False,'raw_selected_ten_review_repeated':False}
    save(out/'BOUND.json',bound);return bound

def bound_inputs(cfg,bound):
    metadata(cfg);require(bound['status']=='cartesian_bound_no_calls' and bound['implementation']==pin(__file__) and set(bound['rows'])==set(CELLS),'Bound scope differs')
    for k in ('candidate_result','candidate_review'):require(cfg[k]==bound[k],'Bound candidate differs');verify(cfg[k])
    require(bound['maximum_native_programs']==4 and bound['maximum_internal_FD_energies']==636,'Bound ceilings differ')
    for row in bound['rows'].values():
        for p in row['files'].values():verify(p)
        verify(row['native_authority'])

def run(bound_path,release_path,out):
    started=time.monotonic();cpu=time.process_time();children=resource.getrusage(resource.RUSAGE_CHILDREN)
    bp=pin(bound_path);b=read(bp);cfg=read(b['config']);bound_inputs(cfg,b);release=json.loads(Path(release_path).read_text())
    require(release['status']=='approved_joint_water_Cartesian_checks' and release['bound']==bp and release['implementation']==pin(__file__) and release['candidate_result']==cfg['candidate_result'],'Root release differs')
    require(release['maximum_native_programs']==4 and release['maximum_internal_FD_energies']==636 and release['maximum_workers']==4,'Release limits differ')
    require(os.environ.get('SLURM_JOB_ID','').isdigit() and int(os.environ.get('SLURM_NTASKS','1'))==1 and int(os.environ.get('SLURM_CPUS_PER_TASK','0'))==4 and all(os.environ.get(k)=='1' for k in THREADS),'Exact4CPU/single-thread runtime required')
    native=module(cfg['native_helper'],'cartesian_runtime_native');helper=module(cfg['numerical_helper'],'cartesian_runtime_parser')
    result={'status':'running','implementation':pin(__file__),'bound':bp,'release':pin(release_path),'candidate_result':cfg['candidate_result'],'candidate_review':cfg['candidate_review'],'job_id':os.environ['SLURM_JOB_ID'],
      'rows':{g:{'status':'not_attempted','receipts':[]} for g in CELLS},'native_program_invocations':0,'internal_FD_energies_requested':0,'internal_FD_energies_completed':0,'error':None,'numerical_checks_passed':False,'physical_accuracy_qualified':False}
    lock=threading.RLock();stop=threading.Event()
    def persist():save(out/'RESULT.json',result)
    def cell(gid):
        if stop.is_set():return
        row=b['rows'][gid];target=result['rows'][gid]
        try:
            if row['reuse'] is not None:
                prior=row['reuse'];verify(prior['authority']);check=helper.gradient(verify(prior['check']['output']),ATOMS[gid],True,native);require(check==prior['check'],'Reused numerical evidence changed');admission=compare(check,row['selected_model'])
                with lock:target.update(status='complete',passed=True,reused=True,authority=prior['authority'],check=check,admission=admission,input=row['files']);persist()
                return
            directory=out/gid;directory.mkdir()
            for k,n in (('xyz','system.xyz'),('key','system.key'),('parameters','amoeba09.prm')):shutil.copyfile(verify(row['files'][k]),directory/n);require(pin(directory/n)['sha256']==row['files'][k]['sha256'],'Runtime copy differs')
            with lock:
                require(not stop.is_set() and result['native_program_invocations']<4 and result['internal_FD_energies_requested']+6*ATOMS[gid]<=636,'Native stop/ceiling')
                result['native_program_invocations']+=1;result['internal_FD_energies_requested']+=6*ATOMS[gid];target.update(status='started',input=row['files'],reused=False);persist()
            receipts=[]
            try:
                native.invoke([cfg['executables']['testgrad']['path']]+ARGS,directory,directory/'testgrad.out',receipts)
                check=helper.gradient(directory/'testgrad.out',ATOMS[gid],True,native)
                with lock:target['check']=check;result['internal_FD_energies_completed']+=6*ATOMS[gid]
                admission=compare(check,row['selected_model'])
            except Exception:
                stop.set();raise
            finally:
                with lock:target['receipts']=receipts;persist()
            with lock:target.update(status='complete',passed=True,check=check,admission=admission);persist()
        except Exception as exc:
            stop.set()
            with lock:target.update(status='failed',passed=False,error=repr(exc));result['error']=repr(exc);persist()
    try:
        persist()
        with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:list(pool.map(cell,CELLS))
        require(not stop.is_set() and all(v['status']=='complete' and v['passed'] for v in result['rows'].values()),'Incomplete or failed Cartesian checks')
        result.update(status='joint_water_Cartesian_checks_complete',numerical_checks_passed=True)
    except Exception as exc:stop.set();result.update(status='native_or_Cartesian_failure',error=repr(exc))
    finally:
        now=resource.getrusage(resource.RUSAGE_CHILDREN);result.update(global_stop=stop.is_set(),elapsed_seconds=time.monotonic()-started,parent_CPU_seconds=time.process_time()-cpu,aggregate_child_user_CPU_seconds=now.ru_utime-children.ru_utime,aggregate_child_system_CPU_seconds=now.ru_stime-children.ru_stime);persist()
    return result

def main():
    p=argparse.ArgumentParser();s=p.add_subparsers(dest='action',required=True);a=s.add_parser('prepare');a=s.add_parser('bind');a.add_argument('--config',required=True);a=s.add_parser('run');a.add_argument('--bound',required=True);a.add_argument('--release',required=True)
    for a in s.choices.values():a.add_argument('--out',required=True)
    a=p.parse_args();out=Path(a.out).resolve();out.mkdir(parents=True,exist_ok=False)
    if a.action=='prepare':
        cfg=default_config();metadata(cfg);save(out/'CONFIG_SOURCE.json',cfg);r={'status':'source_prepared_unbound','implementation':pin(__file__),'config':pin(out/'CONFIG_SOURCE.json'),'cells':ATOMS,'maximum_native_programs':4,'maximum_internal_FD_energies':636,'native_calls':0,'candidate_bound':False};save(out/'MANIFEST.json',r)
    else:r=bind(a.config,out) if a.action=='bind' else run(a.bound,a.release,out)
    print(json.dumps({'status':r['status'],'out':str(out)}))
if __name__=='__main__':main()
