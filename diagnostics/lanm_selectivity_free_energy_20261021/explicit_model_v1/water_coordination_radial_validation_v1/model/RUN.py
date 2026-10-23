#!/usr/bin/env python3
"""Conditional four-cell fresh W9 radial model validation; source remains unbound."""
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

HERE=Path(__file__).resolve().parent;EXPLICIT=HERE.parents[1];ROOT=EXPLICIT.parents[2]
WORK=ROOT/'workspaces/lanm_selectivity_free_energy_20261021/explicit_model_v1'
CELLS=tuple('aquo9_'+s+'_'+m for s in ('minus','plus') for m in ('La','Dy'))
TRAINING_CELLS=tuple(q+'_'+m for q in ('water_minus','water_central','water_plus','aquo8_central','aquo9_central') for m in ('La','Dy'))
ANCHORS=('aquo9_central_La','aquo9_central_Dy')
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
def digest(value):return hashlib.sha256(json.dumps(value,sort_keys=True,separators=(',',':'),allow_nan=False).encode()).hexdigest()

def default_config():
    # Static submitted configuration only. Never open the running fit result.
    fit_config=pin(WORK/'paired_water_coordination_joint_v1/run_v1/CONFIG_RUN.json');fcfg=read(fit_config);ecfg=read(fcfg['evaluator_config'])
    return {'implementation':pin(__file__),'specification':pin(HERE/'SPEC.md'),'root_declaration':pin(HERE.parent/'ROOT_DECLARATION.md'),
        'geometry_manifest':pin(WORK/'water_coordination_radial_validation_v1/preparation_v1/MANIFEST.json'),
        'fit_config':fit_config,'fit_implementation':fcfg['implementation'],'fit_evaluator':fcfg['evaluator'],'fit_evaluator_config':fcfg['evaluator_config'],
        'fit_training':fcfg['training_result'],'fit_objective':fcfg['objective'],'body_helper':fcfg['body_helper'],
        'numerical_helper':pin(EXPLICIT/'paired_water_transfer_v1/RUN.py'),'native_helper':ecfg['native_helper'],
        'resolver':ecfg['resolver'],'executables':ecfg['executables'],'candidate_result':None,'candidate_review':None,
        'maximum_native_programs':8,'maximum_workers':4,'new_FD_energies':0,'root_release':None}

def metadata(cfg):
    base=copy.deepcopy(cfg)
    for k in ('candidate_result','candidate_review','root_release'):base[k]=None
    require(base==default_config(),'Frozen source/input configuration differs')
    for k in ('specification','root_declaration','fit_implementation','fit_evaluator','fit_evaluator_config','fit_training','fit_objective','body_helper','native_helper','numerical_helper','resolver'):verify(cfg[k])
    for p in cfg['executables'].values():verify(p)
    manifest=read(cfg['geometry_manifest']);require(set(manifest['rows'])==set(CELLS),'Four exact radial geometries required')
    for gid,row in manifest['rows'].items():
        g=row['geometry'];metal=gid.rsplit('_',1)[1];step=-.1 if '_minus_' in gid else .1
        require(g['atoms']==28 and g['water_count']==9 and g['metal']==metal and g['charge']==3 and g['physical_multiplicity']==(1 if metal=='La' else 6) and g['radial_step_A']==step,'Fresh geometry/state differs')
        verify(g['xyz']);verify(row['input']['xyz'])
        require(row['input']['types'][0]==(500 if metal=='La' else 501) and set(row['input']['types'][1:])=={36,37},'Fresh water typing differs')
        t=verify(row['input']['xyz']).read_text().splitlines();require(len(t)==29 and int(t[0].split()[0])==28,'Typed geometry count differs')
        for i,line in enumerate(t[1:]):
            a=line.split();require(int(a[0])==i+1 and a[1]==g['symbols'][i] and [float(v) for v in a[2:5]]==g['positions_A'][i] and int(a[5])==row['input']['types'][i] and [int(v)-1 for v in a[6:]]==row['input']['bonds_zero_based'][i],'Typed geometry/graph differs')
    return manifest

def candidate_gate(cfg):
    require(cfg['candidate_result'] is not None and cfg['candidate_review'] is not None,'Passing terminal candidate and independent review are intentionally unbound')
    require(Path(cfg['candidate_result']['path'])==WORK/'paired_water_coordination_joint_v1/run_v1/execution/RESULT.json','Wrong candidate run')
    candidate=read(cfg['candidate_result']);review=read(cfg['candidate_review'])
    require(candidate['status']=='bounded_joint_coordination_water_complete' and candidate['job_id']=='17676' and candidate['error'] is None and candidate['candidate_frozen'] and candidate['development_checks_passed'],'Passing terminal17676 required')
    require(candidate['config']==cfg['fit_config'] and candidate['implementation']==cfg['fit_implementation'],'Actual fit source/config differs')
    require(review['status']=='selected_ten_raw_review_passed' and review['result']==cfg['candidate_result'] and review['selected_parameters']==candidate['selected']['parameters'],'Independent selected-ten review differs')
    require(set(review['selected_cell_ids'])==set(TRAINING_CELLS) and review['complete']==10 and review['new_native_calls']==0 and review['error'] is None,'Independent review scope differs')
    require(all(review[k] is True for k in ('all_selected_raw_replays_exact','exact_selected_score_replayed','all_development_controls_passed')),'Independent review did not pass');verify(review['implementation'])
    require(review['implementation']!=pin(__file__),'Independent reviewer source required')
    selected=candidate['selected'];require(set(selected['model'])==set(TRAINING_CELLS) and selected['score']['development_checks_passed'],'Selected ten cells incomplete')
    require(candidate['vectors'][digest(selected['parameters'])]==selected,'Selection absent from frozen candidate inventory')
    data=read(cfg['fit_training']);body=module(cfg['body_helper'],'radial_bind_body');obj=module(cfg['fit_objective'],'radial_bind_objective')
    require(obj.evaluate(data,selected['model'],body)==selected['score'],'Selected full objective does not replay')
    evaluator=module(cfg['fit_evaluator'],'radial_bind_evaluator');ecfg=read(cfg['fit_evaluator_config'])
    # The renderer is a pure method. Do not construct an evaluator or execute
    # the fit; bind its exact data/family and preserve the original intervals.
    renderer=evaluator.Evaluator.__new__(evaluator.Evaluator);renderer.config=ecfg;renderer.data=data;renderer.family=candidate['family']
    require([p['name'] for p in renderer.family['parameters']]==evaluator.NAMES,'Six-parameter order changed')
    for j,p in enumerate(renderer.family['parameters']):
        origin=evaluator.ORIGINAL_INITIAL[j] if j<4 else evaluator.EPSILON['La' if j==4 else 'Dy']
        bounds=(.5*origin,1.5*origin) if j>=4 else ((.05,.5) if j%2 else (.85*origin,1.15*origin))
        require((p['lower'],p['upper'])==bounds,'Original interval changed')
    values=renderer.checked_values(selected['parameters']);native=module(cfg['native_helper'],'radial_bind_native');helper=module(cfg['numerical_helper'],'radial_bind_numeric')
    old=read(ecfg['prior_result']);coord=read(ecfg['coordination_result'])
    inherited={}
    for vector in old['vectors'].values():
        for gid,value in vector['model'].items():
            p=value['provenance'];inherited[(gid,p['native_authority']['path'],p['native_authority']['sha256'])]=value
    authorities={}
    for gid,value in selected['model'].items():
        pv=value['provenance'];ar=read(pv['native_authority']);files=pv['effective_input'];rendered,r=renderer.render(values,gid)
        require(pv['requested_vector']==values and pv['effective_parameters']==r['effective_values'],'Selected requested/effective vector differs')
        for p in files.values():verify(p)
        require(files['key']['sha256']==hashlib.sha256(rendered.encode()).hexdigest() and all(files[k]['sha256']==data['rows'][gid]['input'][k]['sha256'] for k in ('xyz','parameters')),'Selected physical input differs')
        if not pv['inherited_native']:
            require(ar['status']=='complete' and ar['cell_id']==gid and ar['cache_identity']['implementation']==cfg['fit_evaluator'] and ar['cache_key']==pv['cache_key']==digest(ar['cache_identity']),'Fresh candidate native authority differs')
            require(pv['native_authority']==pin(Path(cfg['candidate_result']['path']).parent/'cache/cells'/pv['cache_key']/'RESULT.json'),'New native authority path differs')
            require(ar['input']==files and all(ar[k]==value[k] for k in ('energy','forces','native')) and len(ar['receipts'])==2 and all(x['returncode']==0 for x in ar['receipts']),'New native values/receipts differ')
            for receipt in ar['receipts']:verify(receipt['output'])
        elif pv['native_source_result']==ecfg['prior_result']:
            member=inherited.get((gid,pv['native_authority']['path'],pv['native_authority']['sha256']))
            require(member is not None and all(member[k]==value[k] for k in ('energy','forces','native')) and member['provenance']['effective_input']==files,'Inherited17605 authority absent or relabeled')
        else:
            require(gid in ANCHORS and pv['native_source_result']==ecfg['coordination_result'] and pv['native_authority']==ecfg['coordination_result'],'Unknown inherited authority')
            require(all(coord['rows'][gid]['model'][k]==value[k] for k in ('energy','forces','native')),'Inherited17670 values differ')
        helper.replay_native(value,data['rows'][gid],native);authorities[gid]=copy.deepcopy(pv)
    return candidate,data,values,authorities

def bind(config_path,out):
    cfg=json.loads(Path(config_path).read_text());manifest=metadata(cfg);candidate,data,values,authorities=candidate_gate(cfg)
    body=module(cfg['body_helper'],'radial_binding_body');rows={}
    for gid in CELLS:
        new=manifest['rows'][gid];g=new['geometry'];metal=g['metal'];parent='aquo9_central_'+metal;source=data['rows'][parent];old=candidate['selected']['model'][parent]['provenance']['effective_input']
        require(new['input']['types']==source['input']['types'] and new['input']['bonds_zero_based']==source['input']['bonds_zero_based'],'Validation source graph differs')
        directory=out/'inputs'/gid;directory.mkdir(parents=True)
        for key,name in (('xyz','system.xyz'),('key','system.key'),('parameters','amoeba09.prm')):
            sourcepin=new['input']['xyz'] if key=='xyz' else old[key];shutil.copyfile(verify(sourcepin),directory/name)
        files={k:pin(directory/n) for k,n in (('xyz','system.xyz'),('key','system.key'),('parameters','amoeba09.prm'))}
        require(files['xyz']['sha256']==new['input']['xyz']['sha256'] and all(files[k]['sha256']==old[k]['sha256'] for k in ('key','parameters')),'Validation input copy differs')
        inp=copy.deepcopy(source['input']);inp.update(id=gid,source_geometry=g,derivative_test=False,**files)
        row={'geometry':g,'input':inp,'files':files,'parent_training_id':parent};row['bodies']=body.bodies(row);require(row['bodies']==source['bodies'],'Body graph differs');rows[gid]=row
    bound={'status':'fresh_W9_radial_bound_no_calls','implementation':pin(__file__),'config':pin(config_path),'candidate_result':cfg['candidate_result'],'candidate_review':cfg['candidate_review'],
        'geometry_manifest':cfg['geometry_manifest'],'parameters':values,'rows':rows,'anchor_model':{m:copy.deepcopy(candidate['selected']['model']['aquo9_central_'+m]) for m in ('La','Dy')},
        'anchor_reference':{m:copy.deepcopy(data['rows']['aquo9_central_'+m]['reference']) for m in ('La','Dy')},'selected_ten_authorities':authorities,
        'all_ten_selected_raw_replayed':True,'new_native_calls':0,'validation_outcomes_opened':False,'physical_accuracy_qualified':False}
    save(out/'BOUND.json',bound);return bound

def bound_inputs(cfg,bound):
    metadata(cfg);require(bound['implementation']==pin(__file__) and bound['status']=='fresh_W9_radial_bound_no_calls','Bound source differs')
    require(set(bound['rows'])==set(CELLS) and bound['all_ten_selected_raw_replayed'] and set(bound['selected_ten_authorities'])==set(TRAINING_CELLS),'Bound scope differs')
    for k in ('candidate_result','candidate_review','geometry_manifest'):require(bound[k]==cfg[k],'Bound authority differs');verify(cfg[k])
    for row in bound['rows'].values():
        for p in row['files'].values():verify(p)

def admit_result(config,bound,result):
    """Pure saved admission. Dict inputs; four availability rows, no engines."""
    bound_inputs(config,bound);require(result['implementation']==pin(__file__) and read(result['bound'])==bound,'Actual result/bound mismatch')
    require(set(result['rows'])==set(CELLS) and 0<=result['native_program_invocations']<=8 and result['new_FD_energies']==0,'Actual finite scope mismatch')
    release=read(result['release']);require(release['implementation']==pin(__file__) and release['bound']==result['bound'] and release['status']=='approved_fresh_W9_radial_model','Actual release differs')
    native=module(config['native_helper'],'orientation_saved_native');helper=module(config['numerical_helper'],'orientation_saved_numeric');admitted={}
    for gid in CELLS:
        row=result['rows'][gid]
        try:
            require(row['status']=='complete' and row['input']==bound['rows'][gid]['files'],'Native cell unavailable/input differs')
            require(len(row['receipts'])==2 and all(x['returncode']==0 for x in row['receipts']),'Native receipt incomplete')
            for record,program,args in zip(row['receipts'],('analyze','testgrad'),(['system.xyz','E'],['system.xyz','Y','N'])):
                require(record['argv']==[config['executables'][program]['path']]+args,'Native command differs');verify(record['output'])
            value=row['model'];helper.replay_native(value,bound['rows'][gid],native)
            require(value['native']['finite_difference'] is None,'Unexpected reserved finite differences')
            admitted[gid]={'status':'complete','model':value,'input':row['input'],'receipts':row['receipts']}
        except Exception as exc:admitted[gid]={'status':'unavailable','model':None,'reason':repr(exc)}
    complete=sum(x['status']=='complete' for x in admitted.values())
    require(complete<4 or result['status']=='water_coordination_radial_model_complete' and result['error'] is None and not result['global_stop'],'Contradictory terminal labels')
    return {'status':'saved_model_admission_complete' if complete==4 else 'saved_model_admission_partial','rows':admitted,'complete':complete,
        'anchor_model':bound['anchor_model'],'new_native_calls':0,'full12_validation_complete':False,'physical_accuracy_qualified':False}

def run(bound_path,release_path,out):
    started=time.monotonic();cpu=time.process_time();children=resource.getrusage(resource.RUSAGE_CHILDREN)
    bp=pin(bound_path);bound=read(bp);cfg=read(bound['config']);bound_inputs(cfg,bound);release=json.loads(Path(release_path).read_text())
    require(release['status']=='approved_fresh_W9_radial_model' and release['bound']==bp and release['implementation']==pin(__file__) and release['candidate_result']==cfg['candidate_result'],'Exact root release required')
    require(release['maximum_native_programs']==8 and release['maximum_workers']==4 and release['new_FD_energies']==0,'Release ceiling mismatch')
    require(os.environ.get('SLURM_JOB_ID','').isdigit() and int(os.environ.get('SLURM_NTASKS','1'))==1 and int(os.environ.get('SLURM_CPUS_PER_TASK','0'))==4 and all(os.environ.get(k)=='1' for k in THREADS),'Exact4CPU/single-thread runtime required')
    native=module(cfg['native_helper'],'orientation_runtime_native');helper=module(cfg['numerical_helper'],'orientation_runtime_numeric')
    result={'status':'running','implementation':pin(__file__),'bound':bp,'release':pin(release_path),'candidate_result':cfg['candidate_result'],'candidate_review':cfg['candidate_review'],
        'job_id':os.environ['SLURM_JOB_ID'],'rows':{g:{'status':'not_attempted','receipts':[]} for g in CELLS},'native_program_invocations':0,'new_FD_energies':0,
        'error':None,'anchor_model':bound['anchor_model'],'full12_validation_complete':False,'physical_accuracy_qualified':False}
    lock=threading.RLock();stop=threading.Event()
    def persist():save(out/'RESULT.json',result)
    def fail(exc):
        stop.set()
        with lock:result['error']=repr(exc);persist()
    def cell(gid):
        if stop.is_set():return
        target=result['rows'][gid];row=bound['rows'][gid]
        try:
            directory=out/gid;directory.mkdir()
            for key,name in (('xyz','system.xyz'),('key','system.key'),('parameters','amoeba09.prm')):
                shutil.copyfile(verify(row['files'][key]),directory/name);require(pin(directory/name)['sha256']==row['files'][key]['sha256'],'Runtime copy differs')
            with lock:target.update(status='started',input=row['files']);persist()
            for program,args in (('analyze',['system.xyz','E']),('testgrad',['system.xyz','Y','N'])):
                with lock:require(not stop.is_set(),'Global stop active');require(result['native_program_invocations']<8,'Native ceiling');result['native_program_invocations']+=1;target['active_program']=program;persist()
                receipts=[];output=directory/(program+'.out')
                try:
                    native.invoke([cfg['executables'][program]['path']]+args,directory,output,receipts);helper.program_admission(output,program,native)
                    if program=='testgrad':
                        n=native.parse(directory/'analyze.out',output,dict(row['input'],derivative_test=False))
                        with lock:target['model']={'energy':n['total_energy_kcal_mol'],'forces':n['force_kcal_mol_A'],'native':n,'provenance':{'bound':bp,'input':row['files'],'candidate_result':cfg['candidate_result']}}
                except Exception as exc:fail(exc);raise
                finally:
                    with lock:target['receipts'].extend(receipts);persist()
            with lock:target['status']='complete';persist()
        except Exception as exc:
            stop.set()
            with lock:target.update(status='failed',error=repr(exc));result['error']=repr(exc);persist()
    try:
        persist()
        with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:list(pool.map(cell,CELLS))
        require(not stop.is_set() and all(r['status']=='complete' for r in result['rows'].values()),'Incomplete water orientation cells')
        result['status']='water_coordination_radial_model_complete'
        verify(cfg['candidate_result']);verify(cfg['candidate_review']);result['original_results_unchanged']=True
    except Exception as exc:stop.set();result.update(status='native_or_admission_failure',error=repr(exc))
    finally:
        ch=resource.getrusage(resource.RUSAGE_CHILDREN);result.update(global_stop=stop.is_set(),elapsed_seconds=time.monotonic()-started,parent_CPU_seconds=time.process_time()-cpu,aggregate_child_user_CPU_seconds=ch.ru_utime-children.ru_utime,aggregate_child_system_CPU_seconds=ch.ru_stime-children.ru_stime,per_program_child_CPU_may_overlap_not_summed=True);persist()
    return result

def main():
    p=argparse.ArgumentParser();s=p.add_subparsers(dest='action',required=True)
    a=s.add_parser('bind');a.add_argument('--config',required=True)
    a=s.add_parser('run');a.add_argument('--bound',required=True);a.add_argument('--release',required=True)
    for a in s.choices.values():a.add_argument('--out',required=True)
    args=p.parse_args();out=Path(args.out).resolve();out.mkdir(parents=True,exist_ok=False)
    try:r=bind(args.config,out) if args.action=='bind' else run(args.bound,args.release,out)
    except Exception as exc:r={'status':'adapter_failed_or_unavailable','error':repr(exc),'new_native_calls':0 if args.action=='bind' else None,'rows':{g:{'status':'not_attempted'} for g in CELLS}};save(out/'FAILURE.json',r)
    print(json.dumps({'status':r['status'],'out':str(out),'error':r.get('error')}));return 0 if r['status'] in ('fresh_W9_radial_bound_no_calls','water_coordination_radial_model_complete') else 1
if __name__=='__main__':raise SystemExit(main())
