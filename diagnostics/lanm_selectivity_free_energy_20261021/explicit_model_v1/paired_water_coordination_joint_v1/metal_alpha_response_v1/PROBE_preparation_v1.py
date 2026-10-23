#!/usr/bin/env python3
"""Eight fixed alpha probes; original science failure is evidence, not a gate."""
import argparse
import concurrent.futures
import copy
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import re
import resource
import shutil
import threading
import time
import numpy as np

HERE=Path(__file__).resolve().parent
THREADS=('OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS')
GEOMETRIES=('water_minus','water_central','water_plus','aquo8_central','aquo9_central')
CELLS=tuple(q+'_'+m for q in GEOMETRIES for m in ('La','Dy'))
ALPHA={'La':1.6597,'Dy':.7280}
POINT=[3.919394849221829,.17662704152555941,3.640166322240097,.15457541153280463,.15863887317162664,.202479823314524]
FRACTIONS=((-0.05,'minus5'),(-0.01,'minus1'),(0.01,'plus1'),(0.05,'plus5'))

def require(ok,message):
    if not ok:raise ValueError(message)
def pin(path):
    p=Path(path).resolve();return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def verify(p):require(pin(p['path'])==p,'Changed pin: '+p['path']);return Path(p['path'])
def read(p):return json.loads(verify(p).read_text())
def token(x):return hashlib.sha256(json.dumps(x,sort_keys=True,separators=(',',':'),allow_nan=False).encode()).hexdigest()
def fingerprint(c):return token({k:v for k,v in c.items() if k not in ('root_release','preparation_result')})
def module(p,name):
    s=importlib.util.spec_from_file_location(name,verify(p));m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
def save(p,x):
    p=Path(p);tmp=p.with_name(p.name+'.tmp');tmp.write_text(json.dumps(x,indent=2,allow_nan=False)+'\n');os.replace(tmp,p)

def alpha_key(base,prm,metal,alpha,resolver):
    require(metal in ALPHA and isinstance(alpha,(int,float)) and not isinstance(alpha,bool) and math.isfinite(alpha) and alpha>0,'Invalid alpha')
    mt=500 if metal=='La' else 501;before=resolver.parse_parameters(prm,base)
    require(before['polarize'][mt]['alpha_A3']==ALPHA[metal],'Original alpha differs')
    lines=base.splitlines(keepends=True);changed=[];old_token=None
    for j,line in enumerate(lines):
        fields=list(re.finditer(r'\S+',line))
        if fields and fields[0].group().lower()=='polarize' and int(fields[1].group())==mt:
            require(len(fields)>=4,'Incomplete polarize line');a=fields[2];old_token=a.group()
            require(float(old_token)==ALPHA[metal],'Unexpected alpha token')
            lines[j]=line[:a.start()]+format(alpha,'.17g')+line[a.end():];changed.append(j)
    require(len(changed)==1,'Exactly one metal polarize record required')
    key=''.join(lines);after=resolver.parse_parameters(prm,key);expected=copy.deepcopy(before);expected['polarize'][mt]['alpha_A3']=alpha
    require(after==expected,'Non-alpha resolved parameter changed')
    for t in (36,37):require(resolver.resolve_pair(before,mt,t)==resolver.resolve_pair(after,mt,t),'Numerical pair coefficients changed')
    return key,{'metal':metal,'metal_type':mt,'alpha_origin_A3':ALPHA[metal],'alpha_A3':alpha,'changed_line_index':changed[0],
        'original_alpha_token':old_token,'new_alpha_token':format(alpha,'.17g'),'one_token_only':True,
        'pdamp_origin':ALPHA[metal]**(1/6),'pdamp_probe':alpha**(1/6),'damping_exponent_ratio_at_fixed_r':math.sqrt(ALPHA[metal]/alpha),
        'polarizability_and_damping_length_both_change':True}

def admitted_terminal(result,review,result_pin,reviewer_pin):
    require(result_pin['sha256']=='5001cb5798571245d5dcff7c6afbcc8c482e1525f0e8e09e2c4c1b124eb03832','Wrong terminal result')
    require(result['status']=='bounded_joint_coordination_water_complete' and result['job_id']=='17676' and result['error'] is None and result['stopping_reason']=='no_actual_improvement_at_complete_trial_batch','Wrong terminal state')
    require(result['selected']['parameters']==POINT and result['selected']['status']=='complete' and set(result['selected']['model'])==set(CELLS),'Selected point differs')
    require(result['development_checks_passed'] is False and result['selected']['score']['development_checks_passed'] is False,'Expected preserved scientific failure')
    require(review['status']=='selected_ten_raw_review_passed' and review['result']==result_pin and review['implementation']==reviewer_pin and review['error'] is None,'Raw selected review unavailable')
    require(review['complete']==10 and set(review['selected_cell_ids'])==set(CELLS) and review['selected_parameters']==POINT and review['new_native_calls']==0,'Review scope differs')
    require(review['all_selected_raw_replays_exact'] is True and review['exact_selected_score_replayed'] is True and review['all_development_controls_passed'] is False,'Raw/scientific review distinction differs')

def default_config():
    root=read(pin(HERE/'ROOT_INPUTS.json'));result=read(root['fit_result']);c=read(result['config']);ecfg=read(c['evaluator_config'])
    src=Path(ecfg['executables']['analyze']['path']).parents[1]/'source'
    # The executable lives under the isolated Tinker source tree's bin directory.
    if not (src/'kpolar.f').exists():
        src=HERE.parents[4]/'workspaces/lanm_selectivity_free_energy_20261021/explicit_model_v1/tinker_source/source'
    return {'protocol':'fixed_metal_alpha_response_v1','implementation':pin(__file__),'specification':pin(HERE/'SPEC.md'),
        'root_inputs':pin(HERE/'ROOT_INPUTS.json'),'root_declaration':root['declaration'],'candidate_result':root['fit_result'],'selected_review':root['selected_review'],
        'reviewer':read(root['selected_review'])['implementation'],'candidate_implementation':c['implementation'],'candidate_evaluator':c['evaluator'],
        'objective':c['objective'],'body_helper':c['body_helper'],'training':c['training_result'],'native_helper':ecfg['native_helper'],'resolver':ecfg['resolver'],
        'executables':ecfg['executables'],'native_physics_sources':{n:pin(src/n) for n in ('kpolar.f','damping.f','induce.f')},
        'maximum_new_physical_cells':40,'maximum_new_native_programs':80,'workers':8,'root_release':None,'preparation_result':None}

def checked_config(cfg):
    require(cfg['implementation']==pin(__file__) and all(os.environ.get(k)=='1' for k in THREADS),'Source/thread controls differ')
    require([cfg[k] for k in ('maximum_new_physical_cells','maximum_new_native_programs','workers')]==[40,80,8],'Finite scope differs')
    for p in cfg.values():
        if isinstance(p,dict) and set(p)=={'path','sha256'}:verify(p)
    for p in (*cfg['executables'].values(),*cfg['native_physics_sources'].values()):verify(p)

def prepare(cfg,out):
    checked_config(cfg);r=read(cfg['candidate_result']);review=read(cfg['selected_review']);admitted_terminal(r,review,cfg['candidate_result'],cfg['reviewer'])
    root=read(cfg['root_inputs']);require(root['declaration']==cfg['root_declaration'] and root['fit_result']==cfg['candidate_result'] and root['selected_review']==cfg['selected_review'],'Root input join differs')
    c=read(r['config']);ecfg=read(c['evaluator_config'])
    require(c['implementation']==cfg['candidate_implementation'] and c['evaluator']==cfg['candidate_evaluator'] and c['objective']==cfg['objective'] and c['body_helper']==cfg['body_helper'] and c['training_result']==cfg['training'],'Consumed source join differs')
    require(all(ecfg[k]==cfg[k] for k in ('native_helper','resolver','executables')),'Native context differs')
    data=read(cfg['training']);require(set(data['rows'])==set(CELLS) and not data['orientation_in_training'],'Training scope differs')
    body=module(cfg['body_helper'],'alpha_body');obj=module(cfg['objective'],'alpha_objective');native=module(cfg['native_helper'],'alpha_native');resolver=module(cfg['resolver'],'alpha_resolver')
    auth=module(cfg['reviewer'],'alpha_authority');ev=module(cfg['candidate_evaluator'],'alpha_saved_renderer')
    renderer=object.__new__(ev.Evaluator);renderer.data=data;renderer.family=r['family'];renderer.resolver=resolver
    admitted={}
    def admit(vector):
        require(vector['status']=='complete' and set(vector['model'])==set(CELLS),'Incomplete saved vector')
        require(obj.evaluate(data,vector['model'],body)==vector['score'],'Saved score differs')
        for gid,value in vector['model'].items():
            row=data['rows'][gid];pv=value['provenance'];inp=pv['effective_input'];key,_=renderer.render(vector['parameters'],gid)
            require(pv['requested_vector']==vector['parameters'] and verify(inp['key']).read_text()==key,'Saved effective parameter rendering differs')
            require(all(inp[k]['sha256']==row['input'][k]['sha256'] for k in ('xyz','parameters')),'Saved physical input differs')
            identity=(gid,pv['native_authority']['path'],pv['native_authority']['sha256'])
            if identity not in admitted:
                producer=auth.native_authority(gid,value)
                parsed=native.parse(verify(value['native']['analyze_output']),verify(value['native']['gradient_output']),dict(row['input'],derivative_test=False))
                require(parsed==value['native'] and value['energy']==parsed['total_energy_kcal_mol'] and value['forces']==parsed['force_kcal_mol_A'],'Saved raw E/F differs')
                admitted[identity]={'cell':gid,'producer':producer,'input':inp}
    baseline=r['selected'];admit(baseline);last=r['iterations'][-1]
    require(last['status']=='complete' and last['origin_key']==last['selected_key']==r['current_key']==token(POINT),'Last Jacobian is not at selected point')
    require(len(last['derivative_keys'])==12,'Twelve saved derivative points required');derivatives=[r['vectors'][k] for k in last['derivative_keys']]
    lo=np.asarray([p['lower'] for p in r['family']['parameters']]);hi=np.asarray([p['upper'] for p in r['family']['parameters']]);width=hi-lo;z=(np.asarray(POINT)-lo)/width;columns=[];spans=[]
    for j in range(6):
        for offset,sign in enumerate((-1,1)):
            v=derivatives[2*j+offset];target=z.copy();target[j]=np.clip(target[j]+sign*.005,0.,1.);expected=np.asarray(POINT).copy();expected[j]=np.clip(expected[j]+(target[j]-z[j])*width[j],lo[j],hi[j])
            require(v['parameters']==expected.tolist(),'Saved parameter FD step differs');admit(v)
        a,b=derivatives[2*j:2*j+2];span=b['normalized_parameters'][j]-a['normalized_parameters'][j];spans.append(span)
        columns.append((np.asarray(b['score']['weighted_residuals'])-np.asarray(a['score']['weighted_residuals']))/span)
    J=np.column_stack(columns);require(J.tolist()==last['Jacobian_normalized'] and spans==last['derivative_actual_widths'],'Saved J6 differs')
    compact={'rows':data['rows'],'baseline_model':baseline['model'],'baseline_score':baseline['score'],'baseline_parameters':POINT,
        'selected_normalized_parameters':baseline['normalized_parameters'],'family':r['family'],'six_column_Jacobian':J.tolist(),'six_column_actual_spans':spans,
        'derivative_keys':last['derivative_keys'],'admitted_native_authorities':list(admitted.values()),'provenance':{k:cfg[k] for k in ('candidate_result','selected_review','reviewer','training','candidate_implementation','candidate_evaluator')},
        'historical17670_review_loss_preserved':True,'new_native_calls':0}
    save(out/'DATA.json',compact);probes=[]
    for metal in ('La','Dy'):
        for fraction,label in FRACTIONS:
            alpha=ALPHA[metal]*(1+fraction);probe={'id':metal+'_'+label,'metal':metal,'fraction':fraction,'eta':(alpha-ALPHA[metal])/ALPHA[metal],'alpha':alpha,'inputs':{},'reused_cells':[]}
            for gid in CELLS:
                if not gid.endswith('_'+metal):probe['reused_cells'].append(gid);continue
                inp=baseline['model'][gid]['provenance']['effective_input'];key,render=alpha_key(verify(inp['key']).read_text(),verify(inp['parameters']).read_text(),metal,alpha,resolver)
                directory=out/'inputs'/probe['id']/gid;directory.mkdir(parents=True)
                shutil.copyfile(verify(inp['xyz']),directory/'system.xyz');shutil.copyfile(verify(inp['parameters']),directory/'amoeba09.prm');(directory/'system.key').write_text(key)
                probe['inputs'][gid]={'files':{k:pin(directory/n) for k,n in (('xyz','system.xyz'),('key','system.key'),('parameters','amoeba09.prm'))},'baseline_input':inp,'render':render}
            require(len(probe['inputs'])==len(probe['reused_cells'])==5,'Active metal inventory differs');probes.append(probe)
    protected=[cfg['candidate_result'],cfg['selected_review'],r['config'],r['evaluator']['ledger'],cfg['training']]
    for p in protected:verify(p)
    manifest={'status':'prepared','implementation':pin(__file__),'scientific_config_sha256':fingerprint(cfg),'data':pin(out/'DATA.json'),'probes':probes,
        'maximum_new_native_programs':80,'maximum_new_physical_cells':40,'workers':8,'new_native_calls':0,'protected_parents':protected,'candidate_result':cfg['candidate_result'],'selected_review':cfg['selected_review']}
    save(out/'MANIFEST.json',manifest);return manifest

class Counter:
    def __init__(self):self.lock=threading.RLock();self.stop=threading.Event();self.count=0
    def reserve(self):
        with self.lock:
            require(not self.stop.is_set() and self.count<80,'Global stop or eighty-program ceiling');self.count+=1;return self.count

def run(cfg,out):
    start=time.monotonic();cpu=time.process_time();ch0=resource.getrusage(resource.RUSAGE_CHILDREN);checked_config(cfg)
    require(os.environ.get('SLURM_JOB_ID','').isdigit() and int(os.environ.get('SLURM_NTASKS','1'))==1 and int(os.environ.get('SLURM_CPUS_PER_TASK','0'))==8,'Exactly one eight-CPU Slurm task required')
    prep=read(cfg['preparation_result']);manifest=read(prep['manifest']);data=read(manifest['data']);release=read(cfg['root_release'])
    require(prep['status']=='prepared' and prep['new_native_calls']==0 and manifest['scientific_config_sha256']==fingerprint(cfg),'Prepared binding differs')
    require(release['status']=='approved_fixed_metal_alpha_diagnostic' and release['implementation']==cfg['implementation'] and release['preparation_result']==cfg['preparation_result'] and release['scientific_config_sha256']==fingerprint(cfg),'Root release binding differs')
    require(release['maximum_new_native_programs']==80 and release['maximum_new_physical_cells']==40 and release['workers']==8 and release['candidate_selection_authorized'] is False,'Release scope differs')
    require([p['id'] for p in manifest['probes']]==[m+'_'+l for m in ('La','Dy') for _,l in FRACTIONS] and len(data['baseline_model'])==10,'Exact eight-probe scope differs')
    for p in manifest['protected_parents']:verify(p)
    native=module(cfg['native_helper'],'alpha_run_native');ev=module(cfg['candidate_evaluator'],'alpha_native_admission');admission=object.__new__(ev.Evaluator);admission.native=native
    resolver=module(cfg['resolver'],'alpha_run_resolver');body=module(cfg['body_helper'],'alpha_run_body');obj=module(cfg['objective'],'alpha_run_objective')
    count=Counter();result={'status':'running','implementation':pin(__file__),'config':cfg,'manifest':prep['manifest'],'preparation_result':cfg['preparation_result'],
        'job_id':os.environ['SLURM_JOB_ID'],'native_program_invocations':0,'global_stop':False,'error':None,'new_candidate_selected':False,'reserved_outcomes_opened':False,
        'records':{p['id']:{k:p[k] for k in ('metal','alpha','eta','fraction')}|{'status':'not_attempted','cells':{g:{'status':'not_attempted'} for g in CELLS}} for p in manifest['probes']}}
    def persist():
        with count.lock:
            result.update(native_program_invocations=count.count,global_stop=count.stop.is_set());save(out/'RESULT.json',result)
    def worker(probe):
        record=result['records'][probe['id']]
        try:
            with count.lock:
                record['status']='started';record['model']={g:copy.deepcopy(data['baseline_model'][g]) for g in probe['reused_cells']}
                for g,v in record['model'].items():record['cells'][g]={'status':'complete','reused':True,'authority':v['provenance']['native_authority']}
            for gid,inp in probe['inputs'].items():
                require(not count.stop.is_set(),'Global stop active');row=data['rows'][gid];directory=out/probe['id']/gid;directory.mkdir(parents=True)
                cell={'status':'started','cell_id':gid,'input':inp,'receipts':[]}
                with count.lock:record['cells'][gid]=cell
                try:
                    key,render=alpha_key(verify(inp['baseline_input']['key']).read_text(),verify(inp['baseline_input']['parameters']).read_text(),probe['metal'],probe['alpha'],resolver)
                    require(key==verify(inp['files']['key']).read_text() and render==inp['render'],'Prepared alpha key differs')
                    for k,n in (('xyz','system.xyz'),('key','system.key'),('parameters','amoeba09.prm')):shutil.copyfile(verify(inp['files'][k]),directory/n)
                    cell['runtime_input']={k:pin(directory/n) for k,n in (('xyz','system.xyz'),('key','system.key'),('parameters','amoeba09.prm'))};save(directory/'RESULT.json',cell)
                    for program,args in (('analyze',['system.xyz','E']),('testgrad',['system.xyz','Y','N'])):
                        with count.lock:
                            count.reserve();save(out/'CALL_COUNT.json',{'native_program_invocations':count.count,'maximum':80})
                        receipts=[]
                        try:
                            native.invoke([cfg['executables'][program]['path']]+args,directory,directory/(program+'.out'),receipts);admission._program_admission(directory/(program+'.out'),program)
                        except Exception:count.stop.set();raise
                        finally:
                            with count.lock:cell['receipts'].extend(receipts);save(directory/'RESULT.json',cell)
                    parsed=native.parse(directory/'analyze.out',directory/'testgrad.out',dict(row['input'],derivative_test=False))
                    value={'energy':parsed['total_energy_kcal_mol'],'forces':parsed['force_kcal_mol_A'],'native':parsed,
                        'provenance':{'probe_id':probe['id'],'alpha':probe['alpha'],'baseline':data['provenance']['candidate_result'],'effective_input':cell['runtime_input'],'native_implementation':pin(__file__)}}
                    with count.lock:cell.update(status='complete',model=value);save(directory/'RESULT.json',cell);record['model'][gid]=value
                except Exception as exc:
                    count.stop.set()
                    with count.lock:cell.update(status='failed',error=repr(exc));save(directory/'RESULT.json',cell)
                    raise
            score=obj.evaluate({'rows':data['rows']},record['model'],body)
            with count.lock:record.update(status='complete',score=score)
        except Exception as exc:
            count.stop.set()
            with count.lock:record.update(status='failed_or_unavailable',error=repr(exc));result['error']=repr(exc);persist()
    try:
        persist()
        with concurrent.futures.ThreadPoolExecutor(max_workers=8) as pool:list(pool.map(worker,manifest['probes']))
        require(not count.stop.is_set() and all(v['status']=='complete' for v in result['records'].values()) and count.count==80,'Native probe denominator incomplete')
        result['status']='metal_alpha_diagnostic_complete'
        for p in manifest['protected_parents']:verify(p)
        result['protected_parent_pins_unchanged']=True
    except Exception as exc:count.stop.set();result.update(status='failed_or_unavailable',error=repr(exc))
    finally:
        ch=resource.getrusage(resource.RUSAGE_CHILDREN);result.update(elapsed_seconds=time.monotonic()-start,parent_cpu_seconds=time.process_time()-cpu,
            aggregate_child_user_cpu_seconds=ch.ru_utime-ch0.ru_utime,aggregate_child_system_cpu_seconds=ch.ru_stime-ch0.ru_stime);persist()
    return result

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=('prepare','run','config'));p.add_argument('--config');p.add_argument('--output',required=True);a=p.parse_args()
    if a.stage=='config':save(a.output,default_config());raise SystemExit(0)
    cfg=json.loads(Path(a.config).read_text());out=Path(a.output).resolve();out.mkdir(parents=True,exist_ok=False);start=time.monotonic();cpu=time.process_time()
    try:
        if a.stage=='prepare':
            require(cfg['root_release'] is None and cfg['preparation_result'] is None,'Preparation must be unreleased');m=prepare(cfg,out)
            r={'status':'prepared','implementation':pin(__file__),'config':pin(a.config),'manifest':pin(out/'MANIFEST.json'),'new_native_calls':0,'wall_seconds':time.monotonic()-start,'CPU_seconds':time.process_time()-cpu};save(out/'RESULT.json',r)
        else:r=run(cfg,out)
    except Exception as exc:r={'status':'precall_or_preparation_failed','error':repr(exc),'implementation':pin(__file__)};save(out/'FAILURE.json',r)
    print(json.dumps({'status':r['status'],'output':str(out)}));raise SystemExit(0 if r['status'] in ('prepared','metal_alpha_diagnostic_complete') else 1)
