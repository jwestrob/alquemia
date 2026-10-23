#!/usr/bin/env python3
"""Finite conditional model binding/checks/orientation evaluation; no fitting."""
import argparse
import concurrent.futures
import copy
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import resource
import shutil
import threading
import time

HERE=Path(__file__).resolve().parent
EXPLICIT=HERE.parents[1]
ROOT=EXPLICIT.parents[2]
WORK=ROOT/'workspaces/lanm_selectivity_free_energy_20261021/explicit_model_v1'
FAMILIES=('acetate_mono','acetate_bi','acetamide')
CHECKS=tuple(prefix+f+'_central_'+m for prefix in ('hydrated_','matched_dry_') for f in FAMILIES for m in ('La','Dy'))
ORIENTATIONS=('validation_aquo8_orientation',)+tuple('validation_hydrated_'+f+'_orientation' for f in FAMILIES)
RESERVED=tuple(q+'_'+s+'_'+m for q in ORIENTATIONS for s in ('minus','plus') for m in ('La','Dy'))
THREADS=('OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS')

def require(ok,message):
    if not ok:raise ValueError(message)
def pin(path):
    p=Path(path).resolve();return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def verify(p):
    require(pin(p['path'])==p,'Changed pin: '+p['path']);return Path(p['path'])
def read(p):return json.loads(verify(p).read_text())
def module(p,name):
    s=importlib.util.spec_from_file_location(name,verify(p));m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
def save(path,value):
    p=Path(path);tmp=p.with_name(p.name+'.tmp');tmp.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n');os.replace(tmp,p)

def default_config():
    fitpin=pin(WORK/'paired_donor_fit_v1/run_v1/CONFIG_RUN.json');fit=read(fitpin);ecfg=read(fit['evaluator_config'])
    return {'implementation':pin(__file__),'specification':pin(HERE/'SPEC.md'),'root_declaration':pin(HERE.parent/'ROOT_DECLARATION.md'),
        'outline':pin(EXPLICIT/'paired_donor_validation_outline_v1/OUTLINE.md'),'candidate_config':fitpin,'candidate_job_id':'17621',
        'candidate_source':fit['implementation'],'evaluator':fit['evaluator'],'evaluator_config':fit['evaluator_config'],
        'training_result':fit['training_result'],'objective':fit['objective'],'body_helper':fit['body_helper'],
        'reserved_panel':pin(WORK/'matched_dry_v1/reserved_validation_v1/RESERVED_PANEL.json'),
        'reserved_native_panel':pin(WORK/'matched_dry_v1/reserved_validation_v1/NATIVE_PANEL.json'),
        'numerical_helper':pin(EXPLICIT/'paired_water_transfer_v1/RUN.py'),'native_helper':ecfg['native_helper'],'executables':ecfg['executables'],
        'candidate_result':None,'root_release':None,'maximum_workers':8,'maximum_check_programs':12,'maximum_internal_FD_energies':1344,'maximum_reserved_programs':32}

def inputs(cfg):
    require(cfg==default_config(),'Source/input configuration differs')
    require(cfg['candidate_source']['sha256']=='9559a2692042781268567933f805a9e2827ba6e48352e00f12e79bddc4ce4c76' and cfg['evaluator']['sha256']=='c0d3085cce6e4b7d8e74687ee0729922fe57836426e53731e2b526dc4790592a','Donor implementation changed')
    require(cfg['reserved_panel']['sha256']=='eed99b0fb6dbe7b49857d68853fe809481f6b657c928a7e75197d77c61f3432b' and cfg['reserved_native_panel']['sha256']=='75d3fc69703c8a5acca61f547e3fcf5d614cfe6ac653d2f3ecd3f34a93bcaf5a','Reserved inventory changed')
    for key in ('specification','root_declaration','outline','candidate_source','evaluator','evaluator_config','objective','body_helper','numerical_helper','native_helper'):verify(cfg[key])
    for p in cfg['executables'].values():verify(p)
    data=read(cfg['training_result']);panel=read(cfg['reserved_panel']);native=read(cfg['reserved_native_panel'])
    require(panel['status']=='reserved_complete' and panel['declared_cells']==24 and panel['new_molecular_calls']==0 and not panel['outcomes_opened'],'Original reserved preparation unavailable')
    require(native['source_inventory']==cfg['reserved_panel'],'Reserved native/geometry join differs')
    geometries={g['geometry_id']:g for g in panel['geometries'] if '_orientation_' in g['geometry_id']}
    require(set(geometries)==set(RESERVED),'Sixteen reserved orientations required')
    checks={x['id']:x for x in panel['checks']}
    for gid,g in geometries.items():
        require(g['validation_status']=='reserved' and checks[gid]['passed'],'Reserved geometry gate fails');verify(g['xyz'])
        require(any(c['geometry']==g and c['geometry_id']==gid and c['recipe']=='ecp_pbe0' for c in native['cells']),'Native geometry manifest mismatch')
    require(sum(data['rows'][g]['geometry']['atoms'] for g in CHECKS)==224,'Cartesian atom denominator differs')
    return data,geometries

def parent_id(gid):
    metal=gid.rsplit('_',1)[1]
    if gid.startswith('validation_aquo8_'):return 'aquo8_central_'+metal
    family=next(f for f in FAMILIES if gid.startswith('validation_hydrated_'+f+'_orientation_'))
    return 'hydrated_'+family+'_central_'+metal

def typed_xyz(g,row):
    lines=verify(row['input']['xyz']).read_text().splitlines();require(int(lines[0].split()[0])==g['atoms'],'Parent XYZ count')
    require(g['symbols']==row['geometry']['symbols'] and g['charge']==row['geometry']['charge'] and g['physical_multiplicity']==row['geometry']['physical_multiplicity'],'Parent chemical state differs')
    verify(g['parent_geometry']);verify(row['geometry']['xyz'])
    require(g['parent_geometry']['sha256']==row['geometry']['xyz']['sha256'],'Reserved parent geometry bytes differ')
    raw=verify(g['xyz']).read_text().splitlines();require(int(raw[0])==g['atoms'],'Reserved XYZ count')
    require([x.split()[0] for x in raw[2:]]==g['symbols'] and [[float(v) for v in x.split()[1:4]] for x in raw[2:]]==g['positions_A'],'Reserved XYZ coordinate replay differs')
    result=[f"{g['atoms']}  {g['geometry_id']} reserved orientation; exact candidate topology"]
    for i,(line,x) in enumerate(zip(lines[1:],g['positions_A'])):
        fields=line.split();require(int(fields[0])==i+1 and int(fields[5])==row['input']['types'][i] and [int(j)-1 for j in fields[6:]]==row['input']['bonds_zero_based'][i],'Parent Tinker graph/type differs')
        result.append(f"{i+1:5d} {fields[1]:<4s} "+' '.join(f'{v:20.12f}' for v in x)+f" {fields[5]:>5s}"+''.join(f'{int(j):6d}' for j in fields[6:]))
    require(len(result)==g['atoms']+1,'Typed XYZ count differs')
    require([[float(v) for v in x.split()[2:5]] for x in result[1:]]==g['positions_A'],'Typed coordinate serialization changed')
    return '\n'.join(result)+'\n'

def selected_authority(ev,candidate_path,candidate,gid):
    value=candidate['selected']['model'][gid];values=candidate['selected']['parameters'];pv=value['provenance']
    token,identity,key,render=ev.cache_identity(values,gid)
    require(pv['requested_vector']==values and pv['cache_key']==token and pv['request_cache_implementation']==ev.config['implementation'] and pv['effective_parameters']==render['effective_values'],'Selected request provenance differs')
    inp=pv['effective_input']
    for p in inp.values():verify(p)
    require(inp['key']['sha256']==hashlib.sha256(key.encode()).hexdigest() and all(inp[k]['sha256']==ev.data['rows'][gid]['input'][k]['sha256'] for k in ('xyz','parameters')),'Selected physical input differs')
    require(pv['source_input']==ev.config['source_inputs'][gid],'Selected source lineage differs')
    if pv['reused_baseline']:
        require(pv['native_authority']==ev.config['transfer_result'] and pv['native_implementation'] is None and pv['baseline_origin']==ev.baseline_review[gid],'Baseline producer relabeled')
        require(inp==ev.config['baseline_inputs'][gid] and all(value[k]==ev.baselines[gid][k] for k in ('energy','forces','native')),'Selected baseline differs')
    else:
        ar=read(pv['native_authority']);expected=Path(candidate_path).resolve().parent/'cache/cells'/token/'RESULT.json'
        require(pv['native_authority']==pin(expected) and ar['status']=='complete' and ar['cache_identity']==identity and ar['cache_key']==pv['native_cache_key']==token,'Selected native authority differs')
        require(pv['native_implementation']==ev.config['implementation'] and ar['input']==inp and all(ar[k]==value[k] for k in ('energy','forces','native')),'Selected producer/input/value differs')
        require(len(ar['receipts'])==2 and all(x['returncode']==0 for x in ar['receipts']),'Native receipt incomplete')
        for x in ar['receipts']:verify(x['output'])
    n=value['native'];again=ev.native.parse(verify(n['analyze_output']),verify(n['gradient_output']),dict(ev.data['rows'][gid]['input'],derivative_test=n['finite_difference'] is not None))
    require(again==n and value['energy']==n['total_energy_kcal_mol'] and value['forces']==n['force_kcal_mol_A'],'Selected raw native replay differs')
    return key,render

def bind(config_path,candidate_path,out):
    cfg=json.loads(Path(config_path).read_text());data,geometries=inputs(cfg);candidate=json.loads(Path(candidate_path).read_text())
    require(candidate['status']=='bounded_paired_donor_complete' and candidate['error'] is None and candidate['candidate_frozen'] and candidate['development_checks_passed'],'Only a terminal passing donor candidate can bind')
    require(candidate['implementation']==cfg['candidate_source'] and candidate['config']==cfg['candidate_config'] and candidate['job_id']=='17621' and candidate['original_results_unchanged'],'Wrong actual donor producer')
    selected=candidate['selected'];require(selected['status']=='complete' and selected['score']['development_checks_passed'] and set(selected['model'])==set(data['rows']),'All28 selected records and controls required')
    fit=read(cfg['candidate_config']);release=read(fit['root_release']);driver=module(cfg['candidate_source'],'stageA_fit_source')
    require(release['donor_config_sha256']==driver.fingerprint(fit) and release['donor_status']=='approved_bounded_paired_donor','Actual donor release differs')
    evaluator=module(cfg['evaluator'],'stageA_donor_evaluator');ev=evaluator.Evaluator(read(cfg['evaluator_config']),out/'admission_cache')
    ev.checked_values(selected['parameters']);require(candidate['family']==ev.family and candidate['vectors'][evaluator.digest(selected['parameters'])]==selected,'Selected vector/family differs')
    require(candidate['frozen_water_invariants']==ev.water_invariants and ev.water_invariants['passed'],'Frozen water invariants differ')
    objective=module(cfg['objective'],'stageA_objective');body=module(cfg['body_helper'],'stageA_body')
    require(objective.evaluate(data,selected['model'],body,ev.water_baseline)==selected['score'],'Original selected objective replay failed')
    keys={gid:selected_authority(ev,candidate_path,candidate,gid) for gid in data['rows']}
    check_inputs={gid:{'files':copy.deepcopy(selected['model'][gid]['provenance']['effective_input']),'model':copy.deepcopy(selected['model'][gid]),'geometry':data['rows'][gid]['geometry'],'input':data['rows'][gid]['input']} for gid in CHECKS}
    reserved={}
    for gid in RESERVED:
        g=geometries[gid];parent=parent_id(gid);row=data['rows'][parent];key,render=keys[parent];directory=out/'inputs'/gid;directory.mkdir(parents=True)
        (directory/'system.xyz').write_text(typed_xyz(g,row));(directory/'system.key').write_text(key);shutil.copyfile(verify(row['input']['parameters']),directory/'amoeba09.prm')
        files={k:pin(directory/name) for k,name in (('xyz','system.xyz'),('key','system.key'),('parameters','amoeba09.prm'))}
        require(files['key']['sha256']==selected['model'][parent]['provenance']['effective_input']['key']['sha256'] and files['parameters']['sha256']==row['input']['parameters']['sha256'],'Reserved candidate parameters differ')
        meta=copy.deepcopy(row['input']);meta.update(id=gid,derivative_test=False,source_geometry=cfg['reserved_panel'],**files)
        reserved[gid]={'geometry':g,'input':meta,'bodies':copy.deepcopy(row['bodies']),'parent_training_id':parent,'files':files,'render':render}
        require(body.bodies(reserved[gid])==reserved[gid]['bodies'],'Reserved rigid-body graph differs')
    require(ev.snapshot()['native_program_invocations']==0,'Binder made native calls')
    bound={'status':'passing_donor_candidate_bound_no_calls','implementation':pin(__file__),'config':pin(config_path),'candidate_result':pin(candidate_path),'candidate_config':candidate['config'],
        'donor_parameters':selected['parameters'],'water_parameters':candidate['family']['frozen_water_parameters'],'candidate_score':selected['score'],'water_invariants':candidate['frozen_water_invariants'],
        'consumed_checks':check_inputs,'rows':reserved,'anchor_model':{k:copy.deepcopy(selected['model'][k]) for k in ('aquo8_central_La','aquo8_central_Dy')},
        'selected_native_authorities':{gid:v['provenance']['native_authority'] for gid,v in selected['model'].items()},'new_native_calls':0,'reserved_outcomes_opened':False}
    save(out/'BOUND.json',bound);return bound

def admit_checks(receipt,bound_pin,bound,helper,native):
    require(receipt['status']=='consumed_numerical_checks_complete' and receipt['numerical_checks_passed'] and receipt['error'] is None and not receipt['global_stop'],'Consumed numerical qualification unavailable')
    require(receipt['implementation']==pin(__file__) and receipt['stage']=='checks','Consumed check producer differs')
    require(receipt['bound']==bound_pin and receipt['candidate_result']==bound['candidate_result'] and set(receipt['checks'])==set(CHECKS),'Consumed check scope differs')
    require(receipt['native_program_invocations']<=12 and receipt['internal_FD_energies_requested']<=1344 and receipt['internal_FD_energies_requested']==receipt['internal_FD_energies_completed'],'Numerical cost/coverage differs')
    for gid,row in receipt['checks'].items():
        require(row['status']=='complete' and row['passed'],'Consumed check failed')
        require(row['input']['files']==bound['consumed_checks'][gid]['files'],'Consumed checked input differs')
        for p in row['input']['files'].values():verify(p)
        for record in row['receipts']:
            require(record['returncode']==0 and record['argv'][1:]==['system.xyz','Y','Y','0.00001'],'Consumed native receipt differs');verify(record['output'])
        g=helper.gradient(verify(row['cartesian_check']['output']),bound['consumed_checks'][gid]['geometry']['atoms'],True,native)
        require(g==row['cartesian_check'] and g['fd_passed'],'Consumed Cartesian replay failed')
        replay_check(g,bound['consumed_checks'][gid]['model'])

def replay_check(g,old):
    require(g['fd_passed'] and abs(g['energy']-old['energy'])<=1e-4,'Cartesian/saved energy gate failed')
    require(max(abs(x-y) for a,b in zip(g['forces'],old['forces']) for x,y in zip(a,b))<=1e-3,'Saved analytic-force replay failed')

def run(stage,bound_path,release_path,out):
    start=time.monotonic();cpu=time.process_time();child=resource.getrusage(resource.RUSAGE_CHILDREN)
    bp=pin(bound_path);bound=read(bp);cfg=read(bound['config']);data,geometries=inputs(cfg);release=json.loads(Path(release_path).read_text())
    require(bound['implementation']==pin(__file__) and bound['status']=='passing_donor_candidate_bound_no_calls','Unadmitted bound source')
    require(release['status']=='approved_donor_stageA_model_'+stage and release['bound']==bp and release['implementation']==pin(__file__) and release['candidate_result']==bound['candidate_result'],'Exact root stage release required')
    require(release['maximum_workers']==8 and release['maximum_check_programs']==12 and release['maximum_reserved_programs']==32 and release['maximum_internal_FD_energies']==1344,'Release scope differs')
    require(os.environ.get('SLURM_JOB_ID','').isdigit() and int(os.environ.get('SLURM_NTASKS','1'))==1 and int(os.environ.get('SLURM_CPUS_PER_TASK','0'))==8 and all(os.environ.get(k)=='1' for k in THREADS),'Exact8CPU/one-thread runtime required')
    verify(bound['candidate_result']);verify(bound['candidate_config'])
    require(bound['candidate_score']['development_checks_passed'] and bound['water_invariants']['passed'],'Candidate no longer passes')
    native=module(cfg['native_helper'],'stageA_native');helper=module(cfg['numerical_helper'],'stageA_existing_numerical_helpers')
    for p in bound['selected_native_authorities'].values():verify(p)
    if stage=='reserved':admit_checks(read(release['numerical_receipt']),bp,bound,helper,native)
    slots=CHECKS if stage=='checks' else RESERVED;limit=12 if stage=='checks' else 32
    result={'status':'running','stage':stage,'bound':bp,'release':pin(release_path),'implementation':pin(__file__),'candidate_result':bound['candidate_result'],'job_id':os.environ['SLURM_JOB_ID'],
        'checks':{gid:{'status':'not_attempted','passed':False,'receipts':[]} for gid in CHECKS} if stage=='checks' else {},
        'rows':{gid:{'status':'not_attempted','receipts':[]} for gid in RESERVED} if stage=='reserved' else {},
        'native_program_invocations':0,'internal_FD_energies_requested':0,'internal_FD_energies_completed':0,'numerical_checks_passed':stage=='reserved','error':None,'physical_accuracy_qualified':False}
    lock=threading.RLock();stop=threading.Event()
    def persist():save(out/'RESULT.json',result)
    def fail(exc):
        stop.set()
        with lock:result['error']=repr(exc);persist()
    def cell(gid):
        target=result['checks' if stage=='checks' else 'rows'][gid]
        if stop.is_set():return
        try:
            row=bound['consumed_checks'][gid] if stage=='checks' else bound['rows'][gid]
            files=row['files'];directory=out/gid;directory.mkdir()
            for k,name in (('xyz','system.xyz'),('key','system.key'),('parameters','amoeba09.prm')):
                shutil.copyfile(verify(files[k]),directory/name);require(pin(directory/name)['sha256']==files[k]['sha256'],'Runtime input copy differs')
            with lock:target.update(status='started',input={'files':files,'geometry':row['geometry'],'parent_training_id':row.get('parent_training_id')});persist()
            if stage=='checks':
                old=row['model'];helper.replay_native(old,data['rows'][gid],native)
                if old['native']['finite_difference'] is not None:
                    g=helper.gradient(verify(old['native']['gradient_output']),row['geometry']['atoms'],True,native);replay_check(g,old)
                    with lock:target.update(status='complete',passed=True,model=old,cartesian_check=g,reused_existing_check=old['provenance']['native_authority']);persist()
                    return
            programs=[('testgrad',['system.xyz','Y','Y','0.00001'])] if stage=='checks' else [('analyze',['system.xyz','E']),('testgrad',['system.xyz','Y','N'])]
            for program,args in programs:
                with lock:
                    require(not stop.is_set(),'Global native stop active');require(result['native_program_invocations']<limit,'Stage native ceiling')
                    result['native_program_invocations']+=1
                    if stage=='checks':
                        result['internal_FD_energies_requested']+=6*row['geometry']['atoms'];require(result['internal_FD_energies_requested']<=1344,'FD energy ceiling')
                    target['active_program']={'program':program,'argv':args};persist()
                receipts=[];output=directory/(program+'.out')
                try:
                    native.invoke([cfg['executables'][program]['path']]+args,directory,output,receipts);helper.program_admission(output,program,native)
                    if program=='testgrad':
                        if stage=='checks':
                            g=helper.gradient(output,row['geometry']['atoms'],True,native);replay_check(g,old)
                            with lock:target.update(cartesian_check=g,model=old,passed=True);result['internal_FD_energies_completed']+=6*row['geometry']['atoms']
                        else:
                            n=native.parse(directory/'analyze.out',output,dict(row['input'],derivative_test=False))
                            with lock:target['model']={'energy':n['total_energy_kcal_mol'],'forces':n['force_kcal_mol_A'],'native':n,'provenance':{'bound':bp,'input':files,'candidate_result':bound['candidate_result']}}
                except Exception as exc:fail(exc);raise
                finally:
                    with lock:target['receipts'].extend(receipts);persist()
            with lock:target['status']='complete';persist()
        except Exception as exc:
            stop.set()
            with lock:target.update(status='failed',error=repr(exc));result['error']=repr(exc);persist()
    try:
        persist()
        with concurrent.futures.ThreadPoolExecutor(max_workers=8) as pool:list(pool.map(cell,slots))
        records=result['checks' if stage=='checks' else 'rows'];require(not stop.is_set() and all(r['status']=='complete' for r in records.values()),'Incomplete stage; no progression')
        if stage=='checks':require(all(r['passed'] for r in records.values()),'Numerical check failed');result.update(status='consumed_numerical_checks_complete',numerical_checks_passed=True)
        else:result.update(status='reserved_stageA_model_complete',numerical_receipt=release['numerical_receipt'],anchor_model=bound['anchor_model'])
        verify(bound['candidate_result']);result['parent_result_unchanged']=True
    except Exception as exc:stop.set();result.update(status='native_or_admission_failure',error=repr(exc))
    finally:
        end=resource.getrusage(resource.RUSAGE_CHILDREN);result.update(global_stop=stop.is_set(),elapsed_seconds=time.monotonic()-start,parent_CPU_seconds=time.process_time()-cpu,aggregate_child_user_CPU_seconds=end.ru_utime-child.ru_utime,aggregate_child_system_CPU_seconds=end.ru_stime-child.ru_stime,per_program_child_CPU_may_overlap_not_summed=True);persist()
    return result

def main():
    p=argparse.ArgumentParser();s=p.add_subparsers(dest='action',required=True)
    a=s.add_parser('bind');a.add_argument('--config',required=True);a.add_argument('--candidate',required=True)
    for stage in ('checks','reserved'):
        a=s.add_parser(stage);a.add_argument('--bound',required=True);a.add_argument('--release',required=True)
    for a in s.choices.values():a.add_argument('--out',required=True)
    args=p.parse_args();out=Path(args.out).resolve();out.mkdir(parents=True,exist_ok=False)
    try:r=bind(args.config,args.candidate,out) if args.action=='bind' else run(args.action,args.bound,args.release,out)
    except Exception as exc:
        r={'status':'adapter_failed_or_unavailable','error':repr(exc),'new_native_calls':0 if args.action=='bind' else None,'preserve_any_existing_RESULT_and_native_receipts':True,'checks':{g:{'status':'not_attempted'} for g in CHECKS},'rows':{g:{'status':'not_attempted'} for g in RESERVED}};save(out/'FAILURE.json',r)
    print(json.dumps({'status':r['status'],'output':str(out),'error':r.get('error')}))
    return 0 if r['status'] in ('passing_donor_candidate_bound_no_calls','consumed_numerical_checks_complete','reserved_stageA_model_complete') else 1
if __name__=='__main__':raise SystemExit(main())
