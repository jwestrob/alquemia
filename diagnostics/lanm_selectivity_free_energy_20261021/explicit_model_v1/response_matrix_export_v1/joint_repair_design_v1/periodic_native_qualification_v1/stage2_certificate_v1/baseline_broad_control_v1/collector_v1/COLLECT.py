"""Finite saved eight-baseline comparison; no native/scheduler/FFT/factor invocation."""
import importlib.util
import json
import math
import os
import sys
import time
from decimal import Decimal,localcontext
from pathlib import Path
import hashlib

HERE=Path(__file__).resolve().parent
IDS=tuple(d+'__'+h+'__'+s for d in ('whole_box','strain_1') for h in ('1e-04','5e-05') for s in ('minus','plus'))
THREADS=('OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','BLIS_NUM_THREADS')
def require(ok,msg):
    if not ok:raise ValueError(msg)
def pin(path):
    p=Path(path).resolve()
    with p.open('rb') as f:h=hashlib.file_digest(f,'sha256').hexdigest()
    return {'path':str(p),'sha256':h}
def verify(item):
    require(isinstance(item,dict) and set(item)=={'path','sha256'} and pin(item['path'])==item,'Missing/changed saved pin')
    return Path(item['path'])
def read(item):return json.loads(verify(item).read_text())
def save(path,obj):
    with Path(path).open('x') as f:json.dump(obj,f,indent=2,allow_nan=False);f.write('\n')
def parent(item):
    require(item['sha256']=='66c6259d428511b6b399d30aae4dfb9c65ef470135a1324812dc5c0dbc5218a8','Exact original numerical source')
    spec=importlib.util.spec_from_file_location('unchanged_stage2_saved_parser',verify(item))
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod);return mod

def scheduler(text,submissions):
    """Exact existing seven-column sacct file; retain steps, count parents once."""
    expected={str(x['job_id']) for x in submissions};require(len(expected)==8,'Eight unique submitted jobs')
    records=[];parents={}
    for line in text.splitlines():
        if not line.strip():continue
        f=line.split('|');require(len(f)==7,'Exact sacct schema')
        job,name,state,code,seconds,cpus,totalcpu=f
        root=job.split('.')[0];require(root in expected,'No unrelated jobs in accounting')
        rec=dict(job_id=job,job_name=name,state=state,exit_code=code,elapsed_seconds=int(seconds),allocated_CPUs=int(cpus),
                 allocated_CPU_seconds=int(seconds)*int(cpus),measured_TotalCPU=totalcpu)
        require(rec['elapsed_seconds']>=0 and rec['allocated_CPUs']>=1,'Allocation arithmetic')
        records.append(rec)
        if job==root:
            require(job not in parents,'Duplicate parent allocation');parents[job]=rec
    require(set(parents)==expected,'All eight parent accounting rows')
    require(all(x['state']=='COMPLETED' and x['exit_code']=='0:0' and x['allocated_CPUs']==1 for x in parents.values()),
            'Actual eight successful one-CPU allocations')
    return {'records':records,'parent_jobs':sorted(parents),'allocated_CPU_seconds':sum(x['allocated_CPU_seconds'] for x in parents.values()),
            'allocation_count':8,'batch_or_extern_not_added':True,'scope':'Only eight new baseline endpoints; central and active panel are reused, separately pinned historical costs.'}

def component_slopes(mod,energies):
    out=[]
    with localcontext() as ctx:
        ctx.prec=70
        for h in mod.STEPS:
            minus,plus=energies[h]
            values={}
            for kind in ('total','polarization','nonpolar'):
                def value(x):
                    e=mod.number(x['energy_token']);ep=mod.number(x['polarization_token'])
                    return e if kind=='total' else ep if kind=='polarization' else e-ep
                def rb(x):
                    e=mod.rounding(x['energy_token']);ep=mod.rounding(x['polarization_token'])
                    return e if kind=='total' else ep if kind=='polarization' else e+ep
                values[kind]={'FD':float((value(plus)-value(minus))/(2*h)),
                              'printing_only_bound':float((rb(plus)+rb(minus))/(2*h))}
            out.append({'half_step_A':float(h),'components':values})
    return {'steps':out,'half_minus_full_step_slope_changes':{k:out[1]['components'][k]['FD']-out[0]['components'][k]['FD'] for k in out[0]['components']},
            'component_analytic_derivatives_unavailable':True,'scope':'Raw Etotal, Epolar and Etotal-Epolar slope arithmetic only; no extra component-force attribution.'}

def increment(mod,active,baseline,analytic_active,analytic_baseline,print_active,print_baseline):
    rows=[]
    with localcontext() as ctx:
        ctx.prec=70;analytic=analytic_active-analytic_baseline;fds=[]
        for h in mod.STEPS:
            am,ap=active[h];bm,bp=baseline[h]
            dm=mod.number(am['energy_token'])-mod.number(bm['energy_token'])
            dp=mod.number(ap['energy_token'])-mod.number(bp['energy_token'])
            fd=(dp-dm)/(2*h);fds.append(fd)
            bound=sum((mod.rounding(x['energy_token']) for x in (am,ap,bm,bp)),Decimal(0))/(2*h)+print_active+print_baseline
            rows.append({'half_step_A':float(h),'increment_minus_kcal_mol':float(dm),'increment_plus_kcal_mol':float(dp),
                         'FD_increment':float(fd),'analytic_increment':float(analytic),'signed_residual':float(fd-analytic),
                         'absolute_residual':float(abs(fd-analytic)),'printing_only_bound':float(bound)})
        return {'status':'descriptive_no_pass_criterion','steps':rows,'half_minus_full_step_FD_increment_change':float(fds[1]-fds[0]),
                'interpretation':'Difference of the two separately computed Hamiltonians at identical saved coordinates; does not repair either original gate.',
                'printing_bound_scope':'Printed token uncertainty only, not finite-step, induction or floating-application error.'}

def run(config,release,out):
    out=Path(out);out.mkdir(parents=True,exist_ok=False);tick=time.monotonic();cpu=time.process_time()
    result={'status':'baseline_saved_collection_failed','implementation':pin(__file__),'configuration':pin(config),'release':pin(release),
            'complete':False,'failures':[],'rows':{},'directions':{},'counts':{'expected_new_baseline_gradients':8,'admitted_new_baseline_gradients':0,'reused_central_gradients':1,'baseline_checks_completed':0},
            'new_native_calls':0,'new_FFTs':0,'new_factorizations':0,'new_scheduler_queries':0,'e_app':None,
            'physical_accuracy_qualified':False,'trajectory_stability_qualified':False,'error':None}
    try:
        cfg=json.loads(Path(config).read_text());rel=json.loads(Path(release).read_text());inp=read(cfg['inputs'])
        require(cfg['status']=='bound_baseline_saved_collection' and cfg['implementation']==pin(__file__) and cfg['inputs']==pin(HERE/'INPUTS.json'),'Bound exact collection')
        require(cfg['scheduler_accounting']==inp['scheduler_accounting'],'Exact saved eight-job accounting')
        require(rel['status']=='released_baseline_saved_collection' and rel['implementation']==pin(__file__) and rel['configuration']==pin(config) and Path(rel['output']).resolve()==out.resolve(),'Exact root release')
        review=read(rel['source_review']);require(review['status']=='passed_baseline_saved_collection_source_review' and review['implementation']==pin(__file__) and review['inputs']==cfg['inputs'],'Source review')
        require(all(os.environ.get(k)=='1' for k in THREADS),'One numerical thread')
        m=read(inp['baseline_manifest']);mod=parent(inp['parent_collector']);decl=read(inp['root_declaration']);subs=read(inp['submissions'])
        require(tuple(m['rows'])==IDS and set(cfg['baseline_results'])==set(IDS) and tuple(x['state_id'] for x in subs)==IDS,'Exact eight slots and submissions')
        require(decl['directions']==['whole_box','strain_1'] and decl['original_slope_gate_kcal_mol_A']==float(mod.GATE)==0.001,'Original six gates')
        require(read(inp['baseline_source_review'])['status']=='passed_stage2_baseline_broad_control_source_review','Baseline execution review')
        active=read(inp['active_collection']);admit=read(inp['active_root_admission'])
        require(admit['result']==inp['active_collection'] and active['status']=='stage2_finite_numerical_qualification_failed' and active['complete'],'Original complete failure retained')
        result.update(inputs=cfg['inputs'],active_collection=inp['active_collection'],active_original_status=active['status'],baseline_manifest=inp['baseline_manifest'],
                      reused_active_historical_cost=active['saved_scheduler_accounting'],source_review=rel['source_review'])
        central=read(inp['central_baseline']);require(central['status']=='stage2_gradient_complete' and central['state_id']=='BYPASS_OLD_ABSENT' and central['build']==m['build'],'Exact old baseline center')
        baseline_center=mod.parse_gradient(verify(central['gradient_output']));active_center=mod.parse_gradient(verify(active['rows']['central']['gradient_output']))
        require(baseline_center['role']=='BYPASS_OLD_ABSENT' and active_center['role']=='ACTIVE','Center role identity')
        result['central_baseline']=inp['central_baseline'];result['central_active']=active['rows']['central']['result']
        parsed={};active_tokens={};submap={x['state_id']:x for x in subs}
        for sid in IDS:
            try:
                receipt=read(cfg['baseline_results'][sid]);row=m['rows'][sid];submission=submap[sid]
                require(receipt['status']=='stage2_gradient_complete' and receipt['error'] is None and receipt['state_id']==sid
                        and receipt['role']=='BYPASS_OLD_ABSENT' and receipt['implementation']==inp['baseline_implementation']
                        and receipt['manifest']==inp['baseline_manifest'] and receipt['build']==m['build'],'Exact complete baseline receipt')
                require(receipt['configuration']==submission['configuration'] and receipt['release']==submission['release']
                        and receipt['inputs']==row['inputs'] and receipt['source_inputs']==row['source_inputs'],'Actual submission/input joins')
                require(receipt['certificate'] is None and receipt['native_attempts']==1 and receipt['counts']=={'native_programs':1,'gradient_calls':1,'energy_only_calls':0,'extra_operator_applications':0},'One old baseline gradient only')
                command=receipt['commands']['native'];require(command['returncode']==0 and command['status']=='complete','Native completion')
                require(command['argv']==[read(m['build'])['executables']['gradient_old']['path'],'system.xyz'],'Unchanged admitted old executable')
                data=mod.parse_gradient(verify(receipt['gradient_output']));require(data['role']=='BYPASS_OLD_ABSENT' and data['state_id']==sid and data['input_binding']==receipt['configuration']['sha256'],'Raw baseline identity')
                parsed[sid]=data;arow=active['rows'][sid]
                with verify(arow['gradient_output']).open() as stream:head=[stream.readline().split() for _ in range(4)]
                require(head[2]==['TOTAL_ENERGY',arow['energy_token']] and head[3]==['POLARIZATION_ENERGY',arow['polarization_token']],'Admitted active raw energies')
                active_tokens[sid]={'energy_token':arow['energy_token'],'polarization_token':arow['polarization_token']}
                result['rows'][sid]={'status':'admitted_saved_baseline_gradient','result':cfg['baseline_results'][sid],'gradient_output':receipt['gradient_output'],
                                    'energy_token':data['energy_token'],'polarization_token':data['polarization_token'],'job_id':submission['job_id'],
                                    'commands':receipt['commands'],'cost':receipt['cost'],'internal_counters':'unavailable through unchanged old archive'}
            except Exception as exc:
                result['failures'].append({'state_id':sid,'error':repr(exc)});result['rows'][sid]={'status':'missing_or_failed_saved_baseline','result':cfg['baseline_results'].get(sid),'error':repr(exc)}
        result['counts']['admitted_new_baseline_gradients']=len(parsed)
        directions=read(inp['directions'])
        for direction in directions['directions']:
            name=direction['id']
            if name not in inp['directions_in_order']:continue
            required=[sid for sid in IDS if sid.startswith(name+'__')]
            if any(sid not in parsed for sid in required):result['directions'][name]={'status':'incomplete'};continue
            tangent=mod.load_tangent(direction['coefficients_dr_dq']);ba,bp=mod.projection(baseline_center,direction,tangent);aa,ap=mod.projection(active_center,direction,tangent)
            energies={};aenergies={}
            for h,label in zip(mod.STEPS,('1e-04','5e-05')):
                keys=[name+'__'+label+'__'+s for s in ('minus','plus')];energies[h]=tuple(parsed[k] for k in keys);aenergies[h]=tuple(active_tokens[k] for k in keys)
            slope=mod.slope_checks(name,ba,bp,{h:tuple(x['energy_token'] for x in pair) for h,pair in energies.items()})
            result['directions'][name]={'status':'complete','baseline_original_checks':slope,'active_original_checks':active['numerical']['directions'][name],
                 'baseline_components':component_slopes(mod,energies),'active_components':component_slopes(mod,aenergies),
                 'central_analytic_baseline':float(ba),'central_analytic_active':float(aa),
                 'response_increment':increment(mod,aenergies,energies,aa,ba,ap,bp)}
            result['counts']['baseline_checks_completed']+=3
        result['scheduler_accounting']=cfg['scheduler_accounting'];result['new_baseline_scheduler_cost']=scheduler(verify(cfg['scheduler_accounting']).read_text(),subs)
        result['complete']=len(parsed)==8 and not result['failures'] and result['counts']['baseline_checks_completed']==6
        result['baseline_checks_passed']=result['complete'] and all(x['baseline_original_checks']['passed'] for x in result['directions'].values())
        result['status']='baseline_broad_controls_passed' if result['baseline_checks_passed'] else ('baseline_broad_controls_numerical_failed' if result['complete'] else 'baseline_broad_controls_incomplete')
    except Exception as exc:result['error']=repr(exc)
    result['collector_cost']={'wall_seconds':time.monotonic()-tick,'CPU_seconds':time.process_time()-cpu,'threads':1}
    save(out/'RESULT.json',result)
    return 0 if result['complete'] else 1

if __name__=='__main__':raise SystemExit(run(*sys.argv[1:]))
