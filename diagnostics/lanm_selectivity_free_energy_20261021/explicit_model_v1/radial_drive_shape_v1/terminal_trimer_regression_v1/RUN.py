#!/usr/bin/env python3
"""Six fresh serial radial14 analyzes; admitted fragments and unchanged cycle equations."""
import argparse,ast,copy,hashlib,importlib.util,json,math,os,re,resource,shutil,socket,sys,threading,time,types
from pathlib import Path

THREADS=('OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','BLIS_NUM_THREADS')
IDS=tuple(f'aquo{n}_central_{metal}__{pair}__metal_two_water' for n in (8,9) for metal in ('La','Dy') for pair in (('near','far') if n==8 else ('near',)))
JOINS=('implementation','scientific_ready','inventory','method_inputs','candidate','source_review','hostname','output_directory')

def require(ok,message):
    if not ok:raise ValueError(message)
def pin(path):
    p=Path(path).resolve();return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def verify(q):
    require(isinstance(q,dict) and pin(q['path'])==q,'Changed or unbound pin: '+str(q));return Path(q['path'])
def read(q):return json.loads(verify(q).read_text())
def save(path,x):
    path=Path(path);tmp=path.with_name(path.name+'.tmp');tmp.write_text(json.dumps(x,indent=2,allow_nan=False)+'\n');os.replace(tmp,path)
def module(q,name):
    spec=importlib.util.spec_from_file_location(name,verify(q));m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m

def extract(q,names,namespace,container=None):
    """Compile exact named AST bodies; never execute their producer's setup/run path."""
    tree=ast.parse(verify(q).read_text());body=tree.body
    if container:
        for name in container:
            body=next(n for n in body if isinstance(n,(ast.FunctionDef,ast.ClassDef)) and n.name==name).body
    nodes=[next(n for n in body if isinstance(n,ast.FunctionDef) and n.name==name) for name in names]
    exec(compile(ast.Module(body=nodes,type_ignores=[]),q['path'],'exec'),namespace)
    return {name:namespace[name] for name in names}

def methods(s,extra):
    import numpy as np
    base=module(s['parent_evaluator'],'trimer_original_evaluator')
    ns={'require':require,'copy':copy,'math':math,'re':re,'np':np,'Path':Path}
    checks=types.SimpleNamespace(**extract(extra['radial_checks'],('gate','delta','readback'),ns))
    radial_ns={'require':require,'copy':copy,'math':math,'base':base,'checks':checks}
    extract(s['radial_adapter'],('radial',),radial_ns)
    funcs=extract(s['radial_adapter'],('checked_values','render','_parse'),radial_ns,('evaluator_class','RadialEvaluator'))
    cls=type('FiniteRadialMethods',(),funcs)
    scalar=extract(s['parent_adapter'],('cycles','comparisons'),{'require':require,'math':math})
    return base,cls,scalar

def validate_partition(inv):
    require(tuple(inv['new_analyze_order'])==IDS and set(inv['rows'])==set(IDS),'Exact six-call order')
    require(len(inv['cycles'])==8 and len(inv['contrasts'])==4,'Eight consumed cycles/four pairs')
    a,b,c=map(set,(inv['reused_selected'],inv['static_terms'],inv['rows']))
    require(len(a)==14 and len(b)==22 and len(c)==6 and not(a&b or a&c or b&c),'Exact disjoint 14+22+6 partition')
    require(a|b|c==set().union(*(set(v['model_terms']) for v in inv['cycles'].values())),'Complete literal cycle dependency set')
    require(inv['historical_six_row_import_authorized'] is False,'No historical fresh-row import')
    require(all(not r['needs_export'] and not r['needs_gradient'] for r in inv['rows'].values()),'Analyze-only fresh rows; no matrix/gradient branch')

def context(cfg):
    require(cfg['implementation']==pin(__file__),'Exact implementation')
    ready=read(cfg['scientific_ready']);inv=read(cfg['inventory']);extra=read(cfg['method_inputs']);s=inv['sources']
    require(ready['status']=='frozen_radial14_terminal_trimer_scientific_inventory' and ready['inventory']==cfg['inventory'] and ready['candidate']==cfg['candidate'],'Frozen scientific inventory joins')
    require(read(ready['source_inputs'])==s and inv['source_inputs']==ready['source_inputs'],'Exact source inventory')
    verify(ready['scientific_specification']);verify(ready['metadata_result']);decl=read(s['root_declaration'])
    require(ready['root_declaration']==s['root_declaration']==inv['root_declaration'],'Declared scope')
    require(cfg['candidate']==s['terminal_candidate']==inv['candidate']['artifact']==decl['parents']['terminal_candidate'],'Fixed terminal artifact')
    require(inv['candidate']['key']==decl['candidate_key'],'Fixed terminal key')
    admission=read(s['terminal_root_admission'])
    require(admission['status']=='admitted_bounded_radial14_development_fit' and admission['selected_artifact']==cfg['candidate'] and admission['selected_key']==inv['candidate']['key'] and admission['producer']==s['terminal_result'] and admission['comparison']==s['terminal_comparison'],'Admitted exact terminal; no new selected-cell admission')
    selected=read(cfg['candidate']);require(selected['status']=='complete' and selected['parameters']==inv['candidate']['parameters'] and len(selected['parameters'])==14 and len(selected['cells'])==46,'Complete radial14 selected46')
    validate_partition(inv)
    build=read(s['build']);qual=read(s['qualification_root_review']);review=read(s['qualification_review'])
    require(build['status']=='complete' and build['executables']['analyze']==s['analyze'],'Qualified radial analyze binary')
    require(qual['status']=='admitted_radial_drive_numerical_qualification' and qual['independent_actual_review']==s['qualification_review'] and review['status']=='passed_radial_drive_qualification_actual_review' and qual['result']==review['result'] and qual['numerical_checks_passed'] and qual['zero_parity_passed'],'Actual radial qualification')
    verify(s['native_source_manifest']);verify(s['analyze'])
    old=read(s['feca_root_admission']);require(old['status']=='admitted_consumed_physical_gate_terminal_trimer_regression' and old['result']==s['feca_result'],'Explicit feca12 predecessor')
    previous=read(s['feca_result']);require(all(x['status']=='complete' for x in previous['cycles'].values()),'Admitted prior eight cycles')
    env=read(s['environment']);import numpy as np
    require(sys.version==env['python_version'] and np.__version__==env['versions']['numpy'],'Pinned Python/NumPy runtime')
    require(pin(sys.executable)==env['files'][0],'Pinned interpreter binary')
    for q in env['files']:verify(q)
    base,cls,scalar=methods(s,extra);ec=read(s['evaluator_configuration']);em=read(s['evaluator_manifest'])
    require(ec['manifest']==s['evaluator_manifest'] and ec['native_helper']==s['native_helper'],'Inherited evaluator/helper identity')
    proxy=cls();proxy.domain=read(s['original_domain']);proxy.m=copy.deepcopy(em);proxy.m['rows'].update(inv['rows']);proxy.m['units']=inv['units']
    proxy.h={k:base.module(ec[k],'trimer_'+k) for k in base.HELPERS};proxy.lock=threading.Lock();proxy.rendered={};proxy.batch={'adapter':s['radial_adapter']}
    values=proxy.checked_values(inv['candidate']['parameters']);require(values==selected['parameters'],'No parameter coercion/change')
    joins={}
    for cid,authority in inv['reused_selected'].items():
        record=selected['cells'][cid];row=authority['definition']
        require(authority['artifact']==cfg['candidate'] and authority['json_path']==['cells',cid] and authority['native_admission_reused']==s['terminal_root_admission'],'Selected raw authority')
        require(record['status']=='complete' and record['program_inputs']==authority['program_inputs'] and record['receipts']==authority['receipts'] and record['cache_identity']==authority['cache_identity'],'Literal selected raw evidence')
        key,meta=proxy.render(values,cid)
        inp=record['program_inputs']['analyze'];require(verify(inp['key']).read_text()==key,'Exact radial14 selected rendered key')
        for k in ('xyz','parameters'):
            require(inp[k]['sha256']==row['source_input'][k]['sha256'],'Selected literal geometry/parameters');verify(inp[k])
        ident=authority['cache_identity'];require(ident['source_manifest']==s['native_source_manifest'] and ident['executables']['analyze']==s['analyze'] and ident['radial_drive']==meta['radial_drive'],'Exact radial source/binary/e identity')
        for receipt in record['receipts']:require(receipt['returncode']==0,'Admitted child success');verify(receipt['output'])
        joins[cid]={'artifact':cfg['candidate'],'json_path':['cells',cid],'raw_inputs':record['program_inputs'],'receipts':record['receipts'],'native_admission':s['terminal_root_admission'],'raw_numerical_replay_performed':False}
    for cid,row in inv['rows'].items():
        for q in row['source_input'].values():verify(q)
        key,meta=proxy.render(values,cid,exports=False)
        require(meta['radial_drive']==values[12 if row['metal']=='La' else 13] and not meta['export_enabled'],'Fresh fourteen-coordinate analyze render')
    vals={cid:{'status':'complete','energy':r['energy_kcal_mol'],'mu_Debye':r['mu_Debye']} for cid,r in inv['static_terms'].items()}
    vals.update({cid:copy.deepcopy(selected['cells'][cid]) for cid in inv['reused_selected']})
    return types.SimpleNamespace(inv=inv,s=s,proxy=proxy,values=values,vals=vals,scalar=scalar,previous=previous,reuse_joins=joins)

def execution_gate(cfg,out):
    require(not os.environ.get('SLURM_JOB_ID') and not os.environ.get('SLURM_STEP_ID'),'Truthful login execution')
    require(socket.gethostname()==cfg['hostname']=='biotite' and str(out)==cfg['output_directory'],'Exact login host/output')
    require(all(os.environ.get(k)=='1' for k in THREADS),'One numerical thread per process')
    rel=read(cfg['root_release']);review=read(cfg['source_review'])
    require(rel['status']=='released_radial14_terminal_trimer_regression' and all(rel[k]==cfg[k] for k in JOINS),'Exact root release')
    require(review['status']=='passed_radial14_terminal_trimer_source_review' and all(review[k]==cfg[k] for k in ('implementation','inventory','scientific_ready','method_inputs')),'Exact independent source review')
    require(rel['maximum_native_calls']==6 and rel['native_workers']==1 and rel['maximum_aggregate_numerical_threads']==2 and rel['aggregate_login_limit']==16 and rel['execution_mode']=='login' and rel['active_process_check'] is not None,'Finite serial release/resource check')
    verify(rel['active_process_check']);return rel

def relabel(comparison):
    mapping={'old_component_terminal':'feca12_physical_gate_terminal','new_terminal':'frozen_radial14_terminal5d992'}
    for group in ('cycles','paired_Dy_minus_La'):
        for row in comparison[group].values():row['models']={mapping.get(k,k):v for k,v in row['models'].items()}
    return comparison

def execute_six(ctx,result,out,persist):
    """One fixed order, persisted intent before each invoke, stop on first failure."""
    for cid in IDS:
        row=ctx.inv['rows'][cid];folder=out/'cells'/cid;folder.mkdir(parents=True,exist_ok=False)
        key,meta=ctx.proxy.render(ctx.values,cid,exports=False)
        for k,name in (('xyz','system.xyz'),('parameters','amoeba09.prm')):shutil.copyfile(verify(row['source_input'][k]),folder/name)
        (folder/'system.key').write_text(key);inp={k:pin(folder/name) for k,name in (('xyz','system.xyz'),('key','system.key'),('parameters','amoeba09.prm'))}
        rec={'status':'started','cell_id':cid,'inherited_native':False,'render':meta,'program_inputs':{'analyze':inp},'receipts':[],'executor':result['executor']};result['rows'][cid]=rec
        require(result['native_program_invocations']<6,'Six-call hard ceiling');result['native_program_invocations']+=1;persist()
        try:
            ctx.proxy.h['native_helper'].invoke([ctx.s['analyze']['path'],'system.xyz','EM'],folder,folder/'analyze.out',rec['receipts'])
            rec.update(ctx.proxy._parse(rec,row));rec['status']='complete'
        except Exception as exc:rec.update(status='failed',error=repr(exc));raise
        finally:
            for receipt in rec['receipts']:receipt['executor']=result['executor']
            save(folder/'RESULT.json',rec);persist()

def run(cfg,out,preflight=False):
    tick,cpu=time.monotonic(),time.process_time();child=resource.getrusage(resource.RUSAGE_CHILDREN)
    out=Path(out).resolve();out.mkdir(parents=True,exist_ok=False)
    executor={'mode':'login','job_id':None,'hostname':socket.gethostname(),'controller_pid':os.getpid(),'native_workers':1,'threads_per_process':1,'maximum_owned_numerical_threads':2,'thread_environment':{k:os.environ.get(k) for k in THREADS}}
    result={'status':'started','implementation':pin(__file__),'configuration':cfg,'inventory':cfg.get('inventory'),'candidate':cfg.get('candidate'),'root_release':cfg.get('root_release'),'executor':executor,'job_id':None,'native_program_invocations':0,'rows':{cid:{'status':'not_attempted'} for cid in IDS},'cycles':{},'comparison':{},'error':None,'preflight_only':preflight,'fit':False,'physical_accuracy_qualified':False,'reserved_outcomes_opened':False}
    save(out/'RESULT.json',result);ctx=None
    def persist():
        result['actual_counts']={'new_analyze_attempts':result['native_program_invocations'],'complete_new_analyze':sum(r['status']=='complete' for r in result['rows'].values()),'admitted_native_induction_blocks':sum(v['native_solve_blocks'] for r in result['rows'].values() for v in r.get('radial_response_readback',{}).values()),'new_optimizer_or_conic_solves':0,'new_testgrad':0,'new_reference_calls':0,'inherited_costs_or_calls_added':False}
        if ctx is not None and not preflight:
            result['cycles']=ctx.scalar['cycles'](ctx.inv['cycles'],dict(ctx.vals,**result['rows']))
            result['comparison']=relabel(ctx.scalar['comparisons'](result['cycles'],ctx.previous['cycles']))
        save(out/'RESULT.json',result)
    try:
        require(all(os.environ.get(k)=='1' for k in THREADS),'One thread before scientific imports')
        if not preflight:execution_gate(cfg,out)
        ctx=context(cfg)
        result.update(parameters=ctx.values,selected_key=ctx.inv['candidate']['key'],terminal_metrics=ctx.inv['terminal_metrics'],reused_selected=ctx.reuse_joins,static_terms=ctx.inv['static_terms'],reference=ctx.inv['reference'],historical_comparators=ctx.inv['historical_comparators'],source_authorities=ctx.s,counts=ctx.inv['counts'])
        if preflight:result['status']='passed_radial14_terminal_trimer_metadata_preflight'
        else:
            persist();execute_six(ctx,result,out,persist)
            require(result['native_program_invocations']==6 and all(r['status']=='complete' for r in result['rows'].values()) and all(r['status']=='complete' for r in result['cycles'].values()),'Complete six/eight denominator')
            require(result['actual_counts']['admitted_native_induction_blocks']==12,'Two admitted induction blocks per new analyze');result['status']='radial14_consumed_trimer_regression_complete'
    except Exception as exc:result.update(status='radial14_consumed_trimer_regression_failed',error=repr(exc))
    finally:
        end=resource.getrusage(resource.RUSAGE_CHILDREN)
        result['cost']={'wall_seconds':time.monotonic()-tick,'parent_CPU_seconds':time.process_time()-cpu,'child_user_CPU_seconds':end.ru_utime-child.ru_utime,'child_system_CPU_seconds':end.ru_stime-child.ru_stime,'maximum_owned_numerical_threads':2,'native_workers':1,'threads_per_process':1,'scheduler_allocation':None,'scheduler_allocated_CPU_seconds':0,'inherited_costs_not_added':True}
        try:persist()
        except Exception as exc:result['persistence_error']=repr(exc);result['status']='radial14_consumed_trimer_regression_failed';save(out/'RESULT.json',result)
    return result

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--config',required=True);p.add_argument('--out',required=True);p.add_argument('--preflight',action='store_true');a=p.parse_args()
    config_path=Path(a.config);cfg=json.loads(config_path.read_text());cfg['_configuration_file']=pin(config_path)
    r=run(cfg,Path(a.out),a.preflight);print(json.dumps({'status':r['status'],'error':r['error'],'native_program_invocations':r['native_program_invocations'],'cost':r['cost']}));sys.exit(0 if r['status'] in ('passed_radial14_terminal_trimer_metadata_preflight','radial14_consumed_trimer_regression_complete') else 1)
