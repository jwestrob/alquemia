#!/usr/bin/env python3
"""One frozen consumed-transfer panel; zero-call binding or root-released Slurm."""
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

HERE=Path(__file__).resolve().parent
sp=importlib.util.spec_from_file_location('transfer_data_adapter',HERE/'DATA.py');d=importlib.util.module_from_spec(sp);sp.loader.exec_module(d)
pin,verify,read,require,module=d.pin,d.verify,d.read,d.require,d.module
save=d.e.save

def loaded(path):
    data=json.loads(Path(path).read_text());require(data['status']=='consumed32_prepared','Prepared data unavailable')
    cfg=data['config'];require(cfg['data_source']==pin(HERE/'DATA.py'),'Prepared data source changed')
    for key in ('specification','evaluator_helper','native_helper','body_helper','resolver','reference_helper','testgrad_source'):verify(cfg[key])
    for p in cfg['executables'].values():verify(p)
    require(set(data['rows'])==set(d.CELLS),'Prepared cell denominator')
    verify(cfg['training']);verify(cfg['original_panel'])
    for p in cfg['radial_references'].values():verify(p)
    for row in data['rows'].values():
        for p in row['input']['xyz'],row['input']['key'],row['input']['parameters'],row['geometry']['xyz'],row['baseline_result'],row['reference_result']:verify(p)
        for k in ('output','engrad','input','panel'):verify(row['reference'][k])
        for k in ('analyze_output','gradient_output'):verify(row['baseline']['native'][k])
    return data

def gradient(path,natoms,fd,native):
    text=Path(path).read_text();program_admission(path,'testgrad',native)
    require('Cartesian Gradient Breakdown over Individual Atoms :'in text,'Missing atomic gradients')
    tail=text.split('Cartesian Gradient Breakdown over Individual Atoms :',1)[1];arrays={'Anlyt':[],'Numer':[]}
    for line in tail.splitlines():
        m=re.fullmatch(r'\s*(Anlyt|Numer)\s+(\d+)\s+('+native.NUMBER+r')\s+('+native.NUMBER+r')\s+('+native.NUMBER+r')\s+('+native.NUMBER+r')\s*',line)
        if m:
            require(int(m[2])==len(arrays[m[1]])+1,'Atomic gradient ordering');arrays[m[1]].append([float(m[k])for k in (3,4,5)])
    require(len(arrays['Anlyt'])==natoms and len(arrays['Numer'])==(natoms if fd else 0),'All atom/component gradients required')
    require(all(math.isfinite(v)for ar in arrays.values()for xyz in ar for v in xyz),'Nonfinite gradients')
    energy=float(re.findall(r'Total Potential Energy\s*:\s*('+native.NUMBER+')',text)[0]);error=max(abs(x-y)for a,b in zip(arrays['Anlyt'],arrays['Numer'])for x,y in zip(a,b))if fd else None
    return {'energy':energy,'forces':[[-v for v in xyz]for xyz in arrays['Anlyt']],
      'gradient':arrays['Anlyt'],'numerical_gradient':arrays['Numer'],'fd_total_separation_A':1e-5 if fd else None,
      'fd_each_side_displacement_A':5e-6 if fd else None,'max_component_error_kcal_mol_A':error,
      'fd_passed':error<=.001 if fd else None,'all_atoms_active_and_components_present':True,'output':pin(path)}

def program_admission(path,program,native):
    text=Path(path).read_text();require('Tinker'in text and 'Version 26.2'in text,'Native identity')
    require(not re.search(r'(?i)not converged|undefined|fatal|error termination|\bNaN\b|\bInfinity\b',text),'Native failure/nonfinite marker')
    es=re.findall(r'Total Potential Energy\s*:\s*('+native.NUMBER+')',text)
    require(len(es)==1 and math.isfinite(float(es[0])),'Finite native energy missing')
    if program=='analyze':require('Energy Component Breakdown :'in text,'Native energy components missing')

def replay_native(value,row,native):
    n=value['native'];meta=dict(row['input'],derivative_test=n['finite_difference'] is not None)
    replay=native.parse(verify(n['analyze_output']),verify(n['gradient_output']),meta)
    require(replay==n and value['energy']==n['total_energy_kcal_mol'] and value['forces']==n['force_kcal_mol_A'],'Saved native E/F differs')
    return n

def bind(data_path,candidate_path,out):
    data=loaded(data_path);cfg=data['config'];candidate=json.loads(Path(candidate_path).read_text())
    require(candidate['status']=='bounded_joint_refinement_complete' and candidate['error'] is None,'Incomplete/failed candidate cannot execute')
    require(candidate['implementation']==pin(HERE.parent/'paired_water_repair_v1/joint_refinement_v1/REFINE.py'),'Unexpected candidate source')
    ccfg=read(candidate['config']);require(ccfg['implementation']==candidate['implementation'] and ccfg['evaluator']==cfg['evaluator_helper'],'Candidate source/evaluator mismatch')
    for key in ('specification','objective','body_helper','evaluator','evaluator_config'):verify(ccfg[key])
    require(ccfg['body_helper']==cfg['body_helper'],'Candidate body weights differ')
    ecfg=read(ccfg['evaluator_config'])
    for key in ('native_helper','resolver','executables'):require(ecfg[key]==cfg[key],'Candidate native authority differs')
    selected=candidate['selected'];values=selected['parameters'];require(set(selected['model'])==set(d.CENTERS),'Four complete selected centers required')
    native=module(cfg['native_helper'],'transfer_bind_native');body=module(cfg['body_helper'],'transfer_bind_body');objective=module(ccfg['objective'],'transfer_bind_objective')
    score=objective.evaluate({'rows':{k:data['rows'][k]for k in d.CENTERS}},selected['model'],body)
    require(score==selected['score'],'Selected objective replay differs')
    inputs={};reusable={};authority=[]
    for gid in d.CELLS:
        row=data['rows'][gid];key,render=d.render(cfg,row,values);directory=out/'inputs'/gid;directory.mkdir(parents=True)
        shutil.copyfile(verify(row['input']['xyz']),directory/'system.xyz');shutil.copyfile(verify(row['input']['parameters']),directory/'amoeba09.prm');(directory/'system.key').write_text(key)
        inp={k:pin(directory/file)for k,file in (('xyz','system.xyz'),('key','system.key'),('parameters','amoeba09.prm'))}
        require(inp['xyz']['sha256']==row['input']['xyz']['sha256'] and inp['parameters']['sha256']==row['input']['parameters']['sha256'],'Input copy differs')
        inputs[gid]={'files':inp,'render':render,'input_identity':d.e.digest({'cell':gid,'xyz_sha256':inp['xyz']['sha256'],'key_sha256':inp['key']['sha256'],'prm_sha256':inp['parameters']['sha256'],'native':cfg['native_helper'],'executables':cfg['executables']})}
        if gid in d.CENTERS:
            value=copy.deepcopy(selected['model'][gid]);pv=value['provenance'];ar=read(pv['native_authority']);authority.append(pv['native_authority'])
            require(pv['source_input']=={k:row['input'][k]for k in ('xyz','key','parameters')},'Selected source input differs')
            if pv['reused_baseline']:
                require(render['original_key_reused'] and value['native']==row['baseline']['native'],'Selected baseline reuse differs')
            else:
                require(ar['status']=='complete' and ar['cell_id']==gid and ar['cache_key']==pv['cache_key'],'Selected cache authority differs')
                for k,p in ar['input'].items():verify(p);require(p['sha256']==inp[k]['sha256'],'Selected effective input differs')
                require(ar['cache_identity']['implementation']==cfg['evaluator_helper'] and ar['cache_identity']['native_helper']==cfg['native_helper'] and ar['cache_identity']['executables']==cfg['executables'],'Selected native engine differs')
                require(ar['energy']==value['energy'] and ar['forces']==value['forces'] and ar['native']==value['native'],'Selected native values differ')
                require(len(ar['receipts'])==2 and all(x['returncode']==0 for x in ar['receipts']),'Selected native receipts incomplete')
            replay_native(value,row,native);reusable[gid]={'source':'selected_candidate','authority':pv['native_authority'],'value':value}
        elif render['original_key_reused']:
            value=copy.deepcopy(row['baseline']);replay_native(value,row,native);reusable[gid]={'source':'unchanged_original_input','authority':row['baseline_result'],'value':value}
    require(all(k in reusable and not inputs[k]['render']['water_recipient']for k in d.DRY),'Dry control reuse not established')
    flags={'development_checks_passed':candidate['development_checks_passed'],'selected_score_development_passed':selected['score']['development_checks_passed'],'stopping_reason':candidate['stopping_reason'],'completed_outer_steps':candidate['completed_outer_steps'],'rank_records_at_iteration_origins':[{'origin_key':x['origin_key'],'linear_directions':x['linear_directions']}for x in candidate['iterations']],'final_point_rank_not_inferred':True,'at_parameter_bound':[v==(.05 if j%2 else .85*d.e.INITIAL[j])or v==(.5 if j%2 else 1.15*d.e.INITIAL[j])for j,v in enumerate(values)],'selected_method_sensitivity':candidate['selected_method_sensitivity']}
    bound={'status':'candidate_bound_no_native_calls','implementation':pin(__file__),'data_source':pin(HERE/'DATA.py'),'score_source':pin(HERE/'SCORE.py'),'specification':cfg['specification'],'data':pin(data_path),'candidate_result':pin(candidate_path),'candidate_config':candidate['config'],'parameters':values,'selected_center_model':selected['model'],'candidate_flags':flags,'inputs':inputs,'reuse':reusable,'selected_native_authorities':authority,'maximum_native_programs':48,'maximum_internal_FD_energies':1380,'maximum_workers':8,'new_native_calls':0,'reserved_validation_opened':False}
    save(out/'BOUND.json',bound);return bound

def run(bound_path,release_path,out):
    start=time.monotonic();parent=time.process_time();child0=resource.getrusage(resource.RUSAGE_CHILDREN)
    bound=json.loads(Path(bound_path).read_text());data=loaded(verify(bound['data']));cfg=data['config'];release=json.loads(Path(release_path).read_text())
    require(bound['implementation']==pin(__file__) and bound['data_source']==pin(HERE/'DATA.py') and bound['score_source']==pin(HERE/'SCORE.py'),'Bound source changed')
    require(release['status']=='approved_consumed_water_transfer' and release['bound']==pin(bound_path) and release['implementation']==pin(__file__),'Exact root release required')
    require(release['maximum_native_programs']==48 and release['maximum_workers']==8 and release['reserved_validation_authorized'] is False,'Release scope differs')
    if not all(bound['candidate_flags'][k] for k in ('development_checks_passed','selected_score_development_passed')):
        require(release.get('explicit_diagnostic_candidate_decision') is True,'Development-failing candidate needs explicit diagnostic decision')
    require(os.environ.get('SLURM_JOB_ID','').isdigit() and int(os.environ.get('SLURM_NTASKS','1'))==1 and int(os.environ.get('SLURM_CPUS_PER_TASK','0'))==8,'One root-owned eight-CPU Slurm task required')
    require(all(os.environ.get(k)=='1' for k in d.e.THREADS),'One numerical thread per native process required')
    verify(bound['candidate_result']);verify(bound['candidate_config'])
    for p in bound['selected_native_authorities']:verify(p)
    native=module(cfg['native_helper'],'transfer_run_native');body=module(cfg['body_helper'],'transfer_run_body');score=module(bound['score_source'],'transfer_score')
    for gid,i in bound['inputs'].items():
        for p in i['files'].values():verify(p)
        key,render=d.render(cfg,data['rows'][gid],bound['parameters']);require(key==Path(i['files']['key']['path']).read_text() and render==i['render'],'Rendered input changed')
    for gid,reuse in bound['reuse'].items():verify(reuse['authority']);replay_native(reuse['value'],data['rows'][gid],native)
    result={'status':'running','bound':pin(bound_path),'release':pin(release_path),'implementation':pin(__file__),'job_id':os.environ['SLURM_JOB_ID'],'candidate_flags':bound['candidate_flags'],'rows':{k:{'status':'not_attempted','receipts':[]}for k in d.CELLS},'native_program_invocations':0,'internal_FD_energies_requested':0,'internal_FD_energies_completed':0,'stage_A_passed':False,'error':None,'reserved_validation_opened':False,'physical_accuracy_qualified':False}
    lock=threading.RLock();stop=threading.Event()
    def persist():save(out/'RESULT.json',result)
    def fail(exc):
        stop.set()
        with lock:result['error']=repr(exc);persist()
    def cell(gid):
        row=data['rows'][gid];target=result['rows'][gid];fd=gid in d.FD;center=gid in d.CENTERS
        if stop.is_set():return
        try:
            directory=out/gid;directory.mkdir();meta=dict(row['input'],derivative_test=fd)
            for k,name in (('xyz','system.xyz'),('key','system.key'),('parameters','amoeba09.prm')):shutil.copyfile(verify(bound['inputs'][gid]['files'][k]),directory/name)
            with lock:target.update(status='started',input=bound['inputs'][gid]);persist()
            reuse=bound['reuse'].get(gid);old=reuse['value']if reuse else None
            if old:
                ap=verify(old['native']['analyze_output']);gp=verify(old['native']['gradient_output'])
                if not fd or old['native']['finite_difference'] is not None:
                    extra={}
                    if fd:
                        g=gradient(gp,row['geometry']['atoms'],True,native);require(g['fd_passed'],'Reused Cartesian force check fails');extra['cartesian_check']=g
                    with lock:target.update(status='complete',model=old,reused_from=reuse['authority'],**extra);persist()
                    return
            else:ap=directory/'analyze.out'
            programs=[]if old else[('analyze',['system.xyz','E'])]
            programs.append(('testgrad',['system.xyz','Y','Y','0.00001']if fd else['system.xyz','Y','N']))
            for program,args in programs:
                with lock:
                    require(not stop.is_set(),'Global native stop active');require(result['native_program_invocations']<48,'48-program ceiling')
                    result['native_program_invocations']+=1
                    if program=='testgrad' and fd:result['internal_FD_energies_requested']+=6*row['geometry']['atoms'];require(result['internal_FD_energies_requested']<=1380,'FD energy ceiling')
                    target['active_program']={'program':program,'args':args,'status':'started'}
                    persist()
                output=directory/(program+'.out');local_receipts=[];updates={}
                try:
                    native.invoke([cfg['executables'][program]['path']]+args,directory,output,local_receipts);program_admission(output,program,native)
                    if program=='testgrad':
                        g=gradient(output,row['geometry']['atoms'],fd,native)
                        if fd:require(g['fd_passed'],'Cartesian force consistency failed')
                        if center:
                            require(abs(g['energy']-old['energy'])<=.0001,'Frozen-center energy replay failed')
                            require(max(abs(x-y)for a,b in zip(g['forces'],old['forces'])for x,y in zip(a,b))<=.001,'Frozen-center force replay failed')
                            updates.update(cartesian_check=g,model=old)
                        else:
                            parsed=native.parse(ap,output,meta);updates['model']={'energy':parsed['total_energy_kcal_mol'],'forces':parsed['force_kcal_mol_A'],'native':parsed}
                            if fd:updates['cartesian_check']=g
                        if fd:
                            with lock:result['internal_FD_energies_completed']+=6*row['geometry']['atoms']
                except Exception as exc:fail(exc);raise
                finally:
                    with lock:
                        target['receipts'].extend(local_receipts);target['active_program']['status']='returned';target.update(updates);persist()
            with lock:target.update(status='complete');persist()
        except Exception as exc:
            with lock:target.update(status='failed',error=repr(exc));fail(exc)
    try:
        for gid in d.DRY:result['rows'][gid].update(status='complete',model=bound['reuse'][gid]['value'],reused_from=bound['reuse'][gid]['authority'])
        persist()
        with concurrent.futures.ThreadPoolExecutor(max_workers=4)as pool:list(pool.map(cell,d.CENTERS))
        require(not stop.is_set() and all(result['rows'][k]['status']=='complete' for k in d.CENTERS),'Fitted-center gate failed')
        result['stage_A_passed']=True;persist()
        with concurrent.futures.ThreadPoolExecutor(max_workers=8)as pool:list(pool.map(cell,d.RADIAL+d.HYDRATED))
        require(not stop.is_set() and all(r['status']=='complete'for r in result['rows'].values()),'Transfer native/Cartesian gate failed')
        result['comparison']=score.evaluate(data,{k:r['model']for k,r in result['rows'].items()},body);result['status']='consumed_transfer_complete'
    except Exception as exc:result.update(status='native_or_cartesian_failure',error=repr(exc));stop.set()
    finally:
        ch=resource.getrusage(resource.RUSAGE_CHILDREN);elapsed=time.monotonic()-start
        result.update(elapsed_seconds=elapsed,allocated_cpu_seconds=elapsed*8,parent_cpu_seconds=time.process_time()-parent,aggregate_child_user_cpu_seconds=ch.ru_utime-child0.ru_utime,aggregate_child_system_cpu_seconds=ch.ru_stime-child0.ru_stime,per_invocation_child_cpu_may_overlap_not_summed=True,global_stop=stop.is_set());persist()
    return result

def main():
    ap=argparse.ArgumentParser();sub=ap.add_subparsers(dest='action',required=True)
    a=sub.add_parser('bind');a.add_argument('--data',required=True);a.add_argument('--candidate',required=True);a.add_argument('--out',required=True)
    a=sub.add_parser('run');a.add_argument('--bound',required=True);a.add_argument('--release',required=True);a.add_argument('--out',required=True)
    args=ap.parse_args();out=Path(args.out).resolve();out.mkdir(parents=True,exist_ok=False)
    try:
        if args.action=='bind':r=bind(args.data,args.candidate,out)
        else:r=run(args.bound,args.release,out)
    except Exception as exc:
        r={'status':'precall_admission_failed','error':repr(exc),'implementation':pin(__file__),'rows':{k:{'status':'not_attempted'}for k in d.CELLS},'new_native_calls':0};save(out/'FAILURE.json',r)
        print(json.dumps({'status':r['status'],'result':pin(out/'FAILURE.json'),'error':r['error']}));return 1
    print(json.dumps({'status':r['status'],'result':pin(out/('BOUND.json'if args.action=='bind'else'RESULT.json'))}));return 0 if r['status']in ('candidate_bound_no_native_calls','consumed_transfer_complete')else 1

if __name__=='__main__':raise SystemExit(main())
