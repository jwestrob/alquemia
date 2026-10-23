"""Nine saved cutoff12 baseline gradients; inherited six derivative checks."""
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import sys
import time

HERE=Path(__file__).resolve().parent
IDS=('central',)+tuple(d+'__'+h+'__'+s for d in ('whole_box','strain_1') for h in ('1e-04','5e-05') for s in ('minus','plus'))
THREADS=('OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','BLIS_NUM_THREADS')
COUNTS={'native_programs':1,'gradient_calls':1,'energy_only_calls':0,'extra_operator_applications':0}
def require(ok,message):
    if not ok:raise ValueError(message)
def pin(path):
    p=Path(path).resolve()
    with p.open('rb') as f:sha=hashlib.file_digest(f,'sha256').hexdigest()
    return {'path':str(p),'sha256':sha}
def verify(p):require(isinstance(p,dict) and set(p)=={'path','sha256'} and pin(p['path'])==p,'Changed pin');return Path(p['path'])
def read(p):return json.loads(verify(p).read_text())
def save(path,obj):
    with Path(path).open('x') as f:json.dump(obj,f,indent=2,allow_nan=False);f.write('\n')
def module(p,name):
    spec=importlib.util.spec_from_file_location(name,verify(p));m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m

def receipt_identity(rec,sid,row,cfg,m,sub,build):
    require(rec['status']=='stage2_gradient_complete' and rec['error'] is None and rec['state_id']==sid and rec['role']=='BYPASS_OLD_ABSENT','Exact terminal inactive state')
    require(rec['implementation']==m['execution_sources']['RUN.py'] and rec['manifest']==cfg['manifest'] and rec['build']==cfg['build'] and rec['source_manifest']==m['native_sources'],'Exact source/build/manifest')
    require(rec['configuration']==sub['configuration'] and rec['release']==sub['release'] and rec['inputs']==row['inputs'] and rec['source_inputs']==row['source_inputs'],'Exact submitted inputs')
    require(rec['certificate'] is None and rec['native_attempts']==1 and rec['counts']==COUNTS,'One gradient without initializer/factor/extra action')
    require(str(rec['allocation']['job_id'])==str(sub['job_id']) and rec['allocation']['allocated_CPUs']==rec['allocation']['tasks']==rec['allocation']['nodes']==1,'Actual submitted one-CPU allocation')
    cmd=rec['commands']['native'];require(cmd['status']=='complete' and cmd['returncode']==0 and cmd['argv']==[build['executables']['gradient_old']['path'],'system.xyz'],'Actual admitted executable')
    require(cmd['wall_seconds_limit']==cmd['CPU_seconds_limit']==3600 and cmd['AS_bytes']==3*1024**3,'Unchanged child limits')

def scheduler(text,submissions):
    expected={str(x['job_id']) for x in submissions};require(len(expected)==9 and len(submissions)==9,'Nine distinct submitted allocations including center')
    records=[];parents={}
    for line in text.splitlines():
        if not line.strip():continue
        f=line.split('|');require(len(f)==7,'Seven-column sacct schema')
        job,name,state,code,seconds,cpus,totalcpu=f;root=job.split('.')[0]
        require(root in expected,'No unrelated scheduler records')
        rec=dict(job_id=job,job_name=name,state=state,exit_code=code,elapsed_seconds=int(seconds),allocated_CPUs=int(cpus),allocated_CPU_seconds=int(seconds)*int(cpus),measured_TotalCPU=totalcpu)
        require(rec['elapsed_seconds']>=0 and rec['allocated_CPUs']>=1,'Finite allocation');records.append(rec)
        if job==root:require(job not in parents,'Duplicate parent');parents[job]=rec
    require(set(parents)==expected,'All nine parent records')
    require(all(x['state']=='COMPLETED' and x['exit_code']=='0:0' and x['allocated_CPUs']==1 for x in parents.values()),'Nine successful one-CPU allocations')
    return {'records':records,'parent_jobs':sorted(parents),'allocation_count':9,'allocated_CPU_seconds':sum(x['allocated_CPU_seconds'] for x in parents.values()),'batch_or_extern_not_added':True,'scope':'Nine new cutoff12 gradient programs including new central; build and historical cutoff9 costs separate.'}

def run(config,release,out):
    out=Path(out);out.mkdir(parents=True,exist_ok=False);tick,cpu=time.monotonic(),time.process_time()
    result={'status':'cutoff12_baseline_collection_incomplete','implementation':pin(__file__),'configuration':pin(config),'release':pin(release),'complete':False,'rows':{},'directions':{},'failures':[],'error':None,'counts':{'expected_new_gradients':9,'admitted_new_gradients':0,'reused_center_gradients':0,'checks_completed':0},'new_native_calls':0,'new_FFTs':0,'new_factorizations':0,'new_scheduler_queries':0,'physical_accuracy_qualified':False,'active_response_qualified':False,'trajectory_stability_qualified':False,'e_app':None}
    try:
        cfg=json.loads(Path(config).read_text());rel=json.loads(Path(release).read_text());inp=read(cfg['inputs'])
        require(cfg['status']=='bound_cutoff12_baseline_collection' and cfg['implementation']==pin(__file__) and cfg['inputs']==pin(HERE/'INPUTS.json'),'Exact collector configuration')
        require(rel['status']=='released_cutoff12_baseline_collection' and rel['implementation']==pin(__file__) and rel['configuration']==pin(config) and Path(rel['output']).resolve()==out.resolve(),'Exact collection release')
        review=read(rel['source_review']);require(review['status']=='passed_cutoff12_baseline_collector_source_review' and review['implementation']==pin(__file__) and review['inputs']==cfg['inputs'],'Exact collector source review')
        require(all(os.environ.get(k)=='1' for k in THREADS),'One numerical thread')
        m=read(cfg['manifest']);require(cfg['manifest']==inp['manifest'] and tuple(m['rows'])==IDS and m['maximum_native_programs']==9,'Exact prepared nine states')
        decl=read(inp['root_declaration']);require(m['root_declaration']==inp['root_declaration'] and tuple(decl['states'])==IDS,'Prospective panel identity')
        require(inp['parent_collector']==m['collector'] and inp['parent_collector']['sha256']=='66c6259d428511b6b399d30aae4dfb9c65ef470135a1324812dc5c0dbc5218a8' and inp['component_helper']['sha256']=='19f9f3eb9483b3d28c3b850897716365039b09e9d5b094512c3f82e2ac2140e0','Frozen numerical helper identities')
        mod=module(inp['parent_collector'],'cutoff12_inherited_numerics');base=module(inp['component_helper'],'cutoff12_inherited_components')
        require(float(mod.GATE)==decl['numerical_criteria']['gate_kcal_mol_A']==.001 and list(map(float,mod.STEPS))==decl['numerical_criteria']['half_steps'],'Unchanged numerical gates')
        old=read(inp['cutoff9_collection']);oldreview=read(inp['cutoff9_actual_review'])
        require(inp['cutoff9_collection']==m['baseline_failed_collection'] and inp['cutoff9_actual_review']==decl['evidence'] and oldreview['status']=='reviewed_complete_baseline_reproduction_of_periodic_failure' and oldreview['result']==inp['cutoff9_collection'],'Admitted cutoff9 authority')
        require(old['complete'] and old['status']=='baseline_broad_controls_numerical_failed' and old['counts']['baseline_checks_completed']==6,'Preserve complete cutoff9 numerical failure')
        build=read(cfg['build']);br=read(cfg['build_root_review'])
        require(build['status']=='stage2_build_complete' and build['error'] is None and build['manifest']==cfg['manifest'] and build['source_manifest']==m['native_sources'] and build['compile_attempts']==build['link_attempts']==1 and build['native_programs']==0,'Exact successful isolated build')
        require(br['status']=='admitted_cutoff12_baseline_build' and br['result']==cfg['build'] and br['manifest']==cfg['manifest'],'Root build admission')
        verify(build['executables']['gradient_old']);verify(build['library'])
        subs=read(cfg['submissions']);require(len(subs)==9 and tuple(x['state_id'] for x in subs)==IDS and set(cfg['results'])==set(IDS),'Exact nine submitted/result slots')
        submap={x['state_id']:x for x in subs};parsed={}
        result.update(inputs=cfg['inputs'],manifest=cfg['manifest'],build=cfg['build'],build_root_review=cfg['build_root_review'],cutoff9_collection=inp['cutoff9_collection'],cutoff9_actual_review=inp['cutoff9_actual_review'],source_review=rel['source_review'])
        for sid in IDS:
            try:
                rec=read(cfg['results'][sid]);row=m['rows'][sid];sub=submap[sid]
                receipt_identity(rec,sid,row,cfg,m,sub,build)
                require(Path(cfg['results'][sid]['path']).resolve()==Path(row['output']).resolve()/'RESULT.json','Fixed output path')
                for p in list(row['inputs'].values())+list(row['source_inputs'].values()):verify(p)
                native_cfg=read(rec['configuration']);native_rel=read(rec['release'])
                require(native_cfg['state_id']==sid and native_cfg['operation']=='gradient' and native_cfg['manifest']==cfg['manifest'] and native_cfg['build']==cfg['build'] and native_cfg['build_root_review']==cfg['build_root_review'],'Native configuration identity')
                require(native_rel['status']=='released_cutoff12_native' and native_rel['configuration']==rec['configuration'] and Path(native_rel['output']).resolve()==Path(row['output']).resolve(),'Native root release')
                ns=read(native_rel['source_review']);require(ns['status']=='passed_cutoff12_baseline_source_review' and ns['manifest']==cfg['manifest'] and ns['binding']==m['execution_sources']['BINDING.py'] and ns['implementations']['native']==rec['implementation'],'Native source review')
                for key in ('stdout','stderr'):verify(rec['commands']['native'][key])
                data=mod.parse_gradient(verify(rec['gradient_output']));require(data['role']=='BYPASS_OLD_ABSENT' and data['state_id']==sid and data['input_binding']==rec['configuration']['sha256'],'Raw state/config binding')
                parsed[sid]=data;result['rows'][sid]={'status':'admitted_saved_cutoff12_gradient','result':cfg['results'][sid],'gradient_output':rec['gradient_output'],'energy_token':data['energy_token'],'polarization_token':data['polarization_token'],'job_id':sub['job_id'],'commands':rec['commands'],'cost':rec['cost'],'internal_counters':'unavailable through unchanged old archive'}
            except Exception as exc:
                result['failures'].append({'state_id':sid,'error':repr(exc)});result['rows'][sid]={'status':'missing_or_failed_saved_gradient','result':cfg['results'].get(sid),'error':repr(exc)}
        result['counts']['admitted_new_gradients']=len(parsed)
        directions=read(m['directions']);require(m['directions']==inp['directions'],'Frozen direction/tangent authority')
        chosen=[d for d in directions['directions'] if d['id'] in ('whole_box','strain_1')];require([d['id'] for d in chosen]==['whole_box','strain_1'],'Two fixed directions')
        for d in chosen:
            name=d['id'];needed=['central']+[s for s in IDS if s.startswith(name+'__')]
            if any(s not in parsed for s in needed):result['directions'][name]={'status':'incomplete'};continue
            tangent=mod.load_tangent(d['coefficients_dr_dq']);analytic,printing=mod.projection(parsed['central'],d,tangent)
            energies={h:tuple(parsed[name+'__'+label+'__'+s] for s in ('minus','plus')) for h,label in zip(mod.STEPS,('1e-04','5e-05'))}
            slopes=mod.slope_checks(name,analytic,printing,{h:tuple(x['energy_token'] for x in pair) for h,pair in energies.items()})
            prior=old['directions'][name]
            comparison={'status':'descriptive_no_new_gate','central_analytic_12_minus_9':float(analytic)-prior['central_analytic_baseline'],'steps':[{'half_step_A':a['half_step_A'],'FD_12_minus_9':a['FD_kcal_mol_A']-b['FD_kcal_mol_A'],'absolute_error_12':a['check']['error'],'absolute_error_9':b['check']['error']} for a,b in zip(slopes['steps'],prior['baseline_original_checks']['steps'])]}
            components=base.component_slopes(mod,energies)
            result['directions'][name]={'status':'complete','cutoff12_original_checks':slopes,'cutoff9_original_checks':prior['baseline_original_checks'],'central_analytic_cutoff12':float(analytic),'central_analytic_cutoff9':prior['central_analytic_baseline'],'cutoff12_components':components,'cutoff9_components':prior['baseline_components'],'comparison':comparison}
            result['counts']['checks_completed']+=3
        result['scheduler_accounting']=cfg['scheduler_accounting'];result['new_gradient_scheduler_cost']=scheduler(verify(cfg['scheduler_accounting']).read_text(),subs)
        result['build_cost_separate']={'build':cfg['build'],'cost':build['cost'],'scheduler_accounting':cfg['build_scheduler_accounting']}
        verify(cfg['build_scheduler_accounting'])
        result['historical_cutoff9_cost_separate']=old['new_baseline_scheduler_cost']
        result['complete']=len(parsed)==9 and not result['failures'] and result['counts']['checks_completed']==6
        result['numerical_checks_passed']=result['complete'] and all(x['cutoff12_original_checks']['passed'] for x in result['directions'].values())
        result['status']='cutoff12_baseline_checks_passed' if result['numerical_checks_passed'] else ('cutoff12_baseline_numerical_failed' if result['complete'] else 'cutoff12_baseline_collection_incomplete')
    except Exception as exc:result['error']=repr(exc)
    result['collector_cost']={'wall_seconds':time.monotonic()-tick,'CPU_seconds':time.process_time()-cpu,'threads':1}
    save(out/'RESULT.json',result);return 0 if result['complete'] else 1

if __name__=='__main__':raise SystemExit(run(*sys.argv[1:]))
