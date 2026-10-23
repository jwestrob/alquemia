#!/usr/bin/env python3
"""Six saved native parses and independent eight-cycle arithmetic; no new native."""
import copy,datetime,hashlib,importlib.util,json,math,os,re,sys,time,threading,subprocess
from decimal import Decimal,getcontext
from pathlib import Path
getcontext().prec=50
HERE=Path(__file__).resolve().parent
PINS={};CHECKS=[]
def pin(p):
    p=Path(p).resolve();return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def verify(q):
    assert pin(q['path'])==q,q;PINS[q['path']]=q;return Path(q['path'])
def read(q):return json.loads(verify(q).read_text())
def check(ok,name):
    assert ok,name;CHECKS.append(name)
def D(v):return Decimal(str(v).replace('D','E'))
def arithmetic(*v):return 128*24*sys.float_info.epsilon*max([1.]+[abs(float(x)) for x in v])
def norm(v):return math.sqrt(math.fsum(float(x)**2 for x in v))
def close(a,b,name):
    aa=a if isinstance(a,(list,tuple)) else [a];bb=b if isinstance(b,(list,tuple)) else [b]
    check(len(aa)==len(bb) and max(abs(float(x)-float(y)) for x,y in zip(aa,bb))<=arithmetic(*aa,*bb),name)
def module(q,name):
    spec=importlib.util.spec_from_file_location(name,verify(q));m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m
def forbidden(*a,**kw):raise AssertionError('Saved-only review forbids native/linear-solve execution')

def main(binding,out):
    tick,cpu=time.monotonic(),time.process_time()
    check(all(os.environ.get(k)=='1' for k in ('OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','BLIS_NUM_THREADS')),'One-thread saved review')
    check(binding['implementation']==pin(__file__) and binding['output_directory']==str(out),'Exact reviewer binding')
    verify(binding['declaration']);verify(binding['arithmetic_parent'])
    r=read(binding['result']);cfg=read(binding['configuration']);ready=read(binding['ready'])
    check(r['status']=='radial14_consumed_trimer_regression_complete' and r['error'] is None and not r.get('persistence_error'),'Complete terminal result without hidden persistence failure')
    check(r['configuration']==dict(cfg,_configuration_file=binding['configuration']),'Exact actual configuration')
    check(r['implementation']==cfg['implementation']==ready['implementation'],'Exact executed reviewed source')
    source=read(cfg['source_review']);release=read(cfg['root_release']);root_source=read(release['root_source_review'])
    check(source['status']=='passed_radial14_terminal_trimer_source_review' and root_source['status']=='admitted_radial14_terminal_trimer_execution_source','Independent/root source admissions')
    check(root_source['independent_review']==cfg['source_review'] and root_source['implementation']==cfg['implementation'] and root_source['ready']==binding['ready'],'Root source pin joins')
    check(release['status']=='released_radial14_terminal_trimer_regression','Exact finite actual release')
    for k in ('implementation','scientific_ready','inventory','method_inputs','candidate','source_review','hostname','output_directory'):
        check(release[k]==cfg[k],'Release join '+k)
    check(all(source[k]==cfg[k]==ready[k] for k in ('implementation','scientific_ready','inventory','method_inputs')),'Source receipt scientific and method joins')
    resource_check=read(release['active_process_check']);intent=read(binding['execution_intent'])
    check(intent['configuration']==binding['configuration'] and intent['root_release']==cfg['root_release'] and intent['no_retry'] and intent['maximum_native_calls']==6,'Pre-execution exact intent')
    inv=read(cfg['inventory']);s=inv['sources'];selected=read(cfg['candidate']);adm=read(s['terminal_root_admission'])
    check(adm['status']=='admitted_bounded_radial14_development_fit' and adm['selected_artifact']==cfg['candidate']==inv['candidate']['artifact'],'Existing selected46 admission reused')
    check(r['parameters']==selected['parameters']==inv['candidate']['parameters'] and r['selected_key']==inv['candidate']['key']==adm['selected_key'],'Exact fixed terminal14 parameters')
    check(r['inventory']==cfg['inventory'] and r['reference']==inv['reference'] and r['static_terms']==inv['static_terms'] and r['terminal_metrics']==inv['terminal_metrics'],'Literal inventory/reference/static/terminal metrics')
    check(r['historical_comparators']==inv['historical_comparators'] and r['source_authorities']==s,'Retained comparison/source authorities')
    check(r['fit'] is False and r['physical_accuracy_qualified'] is False and r['reserved_outcomes_opened'] is False,'No fit, heldout or physical accuracy claim')
    check(r['counts']==inv['counts'] and r['native_program_invocations']==6 and len(r['rows'])==6 and len(r['cycles'])==8,'Six/eight denominator')
    check(r['actual_counts']=={'new_analyze_attempts':6,'complete_new_analyze':6,'admitted_native_induction_blocks':12,'new_optimizer_or_conic_solves':0,'new_testgrad':0,'new_reference_calls':0,'inherited_costs_or_calls_added':False},'Exact actual new-call accounting')
    check(release['execution_mode']=='login' and release['native_workers']==1 and release['maximum_aggregate_numerical_threads']==2 and release['aggregate_login_limit']==16,'Released finite login envelope')
    ex=r['executor'];check(ex['mode']=='login' and ex['job_id'] is None and ex['hostname']==cfg['hostname']=='biotite' and ex['controller_pid']>0 and ex['native_workers']==ex['threads_per_process']==1 and ex['maximum_owned_numerical_threads']==2 and all(v=='1' for v in ex['thread_environment'].values()),'Truthful actual login/thread metadata')
    cost=r['cost'];check(cost['scheduler_allocation'] is None and cost['scheduler_allocated_CPU_seconds']==0 and cost['inherited_costs_not_added'] is True,'Local costs exclude old costs and scheduler allocation')

    # Exact qualified methods only; no producer context or constructor, no old-cell replay.
    producer=module(cfg['implementation'],'trimer_saved_producer');extra=read(cfg['method_inputs'])
    base,cls,_=producer.methods(s,extra);ec=read(s['evaluator_configuration']);em=read(s['evaluator_manifest'])
    import numpy as np
    np.linalg.solve=forbidden;subprocess.run=forbidden;subprocess.Popen=forbidden
    proxy=cls();proxy.domain=read(s['original_domain']);proxy.m=copy.deepcopy(em);proxy.m['rows'].update(inv['rows']);proxy.m['units']=inv['units']
    proxy.h={k:module(ec[k],'actual_trimer_'+k) for k in base.HELPERS};proxy.h['native_helper'].invoke=forbidden
    proxy.lock=threading.Lock();proxy.rendered={};proxy.batch={'adapter':s['radial_adapter']}
    for q in (s['analyze'],s['native_source_manifest'],s['radial_adapter'],extra['radial_checks']):verify(q)
    values={};terms={};raw_rows={};last_end=None;receipts=[];maxterm=0.
    for cid in inv['new_analyze_order']:
        cell=r['rows'][cid];row=inv['rows'][cid]
        check(cell['status']=='complete' and cell['cell_id']==cid and not cell['inherited_native'] and len(cell['receipts'])==1,'One successful fresh cell '+cid)
        rec=cell['receipts'][0];receipts.append(rec);folder=Path(cfg['output_directory'])/'cells'/cid
        check(read(pin(folder/'RESULT.json'))==cell,'Immutable per-cell result '+cid)
        check(rec['argv']==[s['analyze']['path'],'system.xyz','EM'] and rec['cwd']==str(folder) and rec['returncode']==0 and rec['timeout_seconds']==300,'Exact analyze receipt '+cid)
        check(rec['executor']==cell['executor']==ex and all(v=='1' for v in rec['threads'].values()),'Same actual executor '+cid)
        start=datetime.datetime.fromisoformat(rec['started_utc']);end=datetime.datetime.fromisoformat(rec['ended_utc'])
        check(end>=start and (last_end is None or start>=last_end),'Strict serial program intervals '+cid);last_end=end
        inp=cell['program_inputs']['analyze'];key,meta=proxy.render(r['parameters'],cid,exports=False)
        check(set(cell['program_inputs'])=={'analyze'} and verify(inp['key']).read_text()==key and cell['render']==meta,'Exact fourteen-coordinate rendered key '+cid)
        for k in ('xyz','parameters'):check(verify(inp[k]).read_bytes()==verify(row['source_input'][k]).read_bytes(),'Unchanged literal '+k+' '+cid)
        check(not list(folder.glob('response_matrix_*.dat')) and not (folder/'testgrad.out').exists(),'No export/testgrad side branch '+cid)
        parsed=proxy._parse(cell,row);check(all(cell[k]==v for k,v in parsed.items()),'Qualified saved native parser exact '+cid)
        check(cell['radial_response_readback']['analyze']['native_solve_blocks']==2,'Two internal induction blocks '+cid)
        text=verify(rec['output']).read_text();check(rec['output']==cell['moments']['analyze_output'],'Native output authority '+cid);verify(cell['moments']['induced_output'])
        tokens=re.findall(r'^\s*Total Potential Energy\s*:\s*([-+\d.EeDd]+)\s+Kcal/mole',text,re.M)
        check(len(tokens)==1 and float(D(tokens[0]))==cell['energy']==cell['moments']['energy_kcal_mol'],'Raw printed energy '+cid)
        n=row['atoms'];dual=[line.split()[2:] for line in text.splitlines() if line.startswith('LOCAL_RESPONSE_INDUCED')][-n:]
        check(len(dual)==n and all(len(x)==6 for x in dual),'Raw last dual vector denominator '+cid)
        mu=[D(row['permanent_lab_Debye'][j])+D(inv['units']['debye_per_eA'])*sum((D(x[j])+D(x[j+3]))/2 for x in dual) for j in range(3)]
        close(mu,cell['mu_Debye'],'Independent raw dual dipole '+cid)
        delta=max(abs(float(a)-b) for a,b in zip(mu,cell['mu_Debye']));maxterm=max(maxterm,delta)
        values[cid]=(D(tokens[0]),mu);terms[cid]={'authority':'six_new_raw_tokens','energy_decimal':tokens[0],'mu_Debye_decimal':list(map(str,mu)),'raw_output':rec['output']}
        raw_rows[cid]={'status':'passed_saved_radial_native_parse','receipt':rec,'input':inp,'maximum_dual_mu_arithmetic_difference_Debye':delta}
    for cid,authority in inv['reused_selected'].items():
        cell=selected['cells'][cid];got=r['reused_selected'][cid]
        check(cell['status']=='complete' and authority['artifact']==cfg['candidate'] and got['artifact']==cfg['candidate'] and got['json_path']==['cells',cid],'Admitted selected scalar authority '+cid)
        check(got['raw_inputs']==cell['program_inputs']==authority['program_inputs'] and got['receipts']==cell['receipts']==authority['receipts'] and got['native_admission']==s['terminal_root_admission'] and got['raw_numerical_replay_performed'] is False,'Selected provenance reused without replay '+cid)
        values[cid]=(D(cell['energy']),list(map(D,cell['mu_Debye'])));terms[cid]={'authority':'admitted_selected_scalar','artifact':cfg['candidate'],'json_path':['cells',cid]}
    for cid,cell in inv['static_terms'].items():
        values[cid]=(D(cell['energy_kcal_mol']),list(map(D,cell['mu_Debye'])));terms[cid]={'authority':'admitted_static_scalar','inventory':cfg['inventory'],'json_path':['static_terms',cid],'historical_source':cell['source']}
    check(len(values)==42 and len(raw_rows)==6 and len(inv['reused_selected'])==14 and len(inv['static_terms'])==22,'Six raw plus36 admitted terms')
    for k,receipt_key in (('child_user_CPU_seconds','child_user_cpu_seconds'),('child_system_CPU_seconds','child_system_cpu_seconds')):close(math.fsum(x[receipt_key] for x in receipts),cost[k],'Exact summed six-child cost '+k)
    check(math.fsum(x['elapsed_seconds'] for x in receipts)<=cost['wall_seconds'] and cost['parent_CPU_seconds']>=0,'Actual parent/wall cost scope')
    feca=read(s['feca_result']);feca_adm=read(s['feca_root_admission'])
    check(feca_adm['result']==s['feca_result'] and feca_adm['actual_review']==s['feca_actual_review'] and feca['status']=='consumed_trimer_regression_complete','Existing feca trimer admission reused');verify(s['feca_actual_review'])
    cycles={};maxE=maxM=0.
    labels=('original_pre_response_baseline','feca12_physical_gate_terminal','frozen_radial14_terminal5d992')
    for cid,c in inv['cycles'].items():
        got=r['cycles'][cid]
        check(got['status']=='complete' and got['missing']==[] and got['model_terms']==c['model_terms'] and len(c['model_terms'])==7 and sorted(c['model_terms'].values())==[-1,-1,-1,1,1,1,1],'Exact seven-term cycle '+cid)
        E=sum(D(co)*values[k][0] for k,co in c['model_terms'].items());mu=[sum(D(co)*values[k][1][j] for k,co in c['model_terms'].items()) for j in range(3)]
        close(E,got['candidate_E3_kcal_mol'],'Independent cycle energy '+cid);close(mu,got['candidate_mu3_Debye'],'Independent cycle dipole '+cid)
        maxE=max(maxE,abs(float(E)-got['candidate_E3_kcal_mol']));maxM=max(maxM,max(abs(float(a)-b) for a,b in zip(mu,got['candidate_mu3_Debye'])))
        check(got['reference_electronic_E3_kcal_mol']==c['reference_electronic_E3_kcal_mol'] and got['reference_mu3_Debye']==c['reference_mu3_Debye'],'Unchanged reference '+cid)
        check(got['original_reference_print_bound_Debye_per_component']==c['reference_print_bound_Debye_per_component'] and got['original_model_print_bound_Debye_per_component']==c['model_print_bound_Debye_per_component'],'Original printing envelopes '+cid)
        old=feca['cycles'][cid];check(old['model_terms']==c['model_terms'] and old['reference_mu3_Debye']==c['reference_mu3_Debye'] and old['reference_electronic_E3_kcal_mol']==c['reference_electronic_E3_kcal_mol'],'Matched feca comparison '+cid)
        models={labels[0]:(D(c['model_E3_kcal_mol']),list(map(D,c['model_mu3_Debye']))),labels[1]:(D(old['candidate_E3_kcal_mol']),list(map(D,old['candidate_mu3_Debye']))),labels[2]:(E,mu)}
        row={'models':{},'reference_electronic_E3_kcal_mol':c['reference_electronic_E3_kcal_mol'],'reference_mu3_Debye':c['reference_mu3_Debye'],'model_terms':c['model_terms']}
        for label,(en,mv) in models.items():
            err=en-D(row['reference_electronic_E3_kcal_mol']);me=[x-D(y) for x,y in zip(mv,row['reference_mu3_Debye'])]
            row['models'][label]={'E3_kcal_mol':float(en),'mu3_Debye':list(map(float,mv)),'energy_error_kcal_mol':float(err),'dipole_error_Debye':list(map(float,me)),'dipole_error_norm_Debye':norm(me)}
        close(row['models'][labels[2]]['energy_error_kcal_mol'],got['energy_error_kcal_mol'],'Cycle signed energy error '+cid);close(row['models'][labels[2]]['dipole_error_Debye'],got['dipole_error_Debye'],'Cycle signed vector error '+cid)
        cycles[cid]=row
    paired={}
    for cid,c in inv['contrasts'].items():
        la=cycles[c['La_cycle']];dy=cycles[c['Dy_cycle']];refE=D(dy['reference_electronic_E3_kcal_mol'])-D(la['reference_electronic_E3_kcal_mol']);refmu=[D(d)-D(l) for l,d in zip(la['reference_mu3_Debye'],dy['reference_mu3_Debye'])]
        row={'models':{},'reference_electronic_Dy_minus_La_E3_kcal_mol':float(refE),'reference_Dy_minus_La_mu3_Debye':list(map(float,refmu))}
        for label in labels:
            l=la['models'][label];d=dy['models'][label];en=D(d['E3_kcal_mol'])-D(l['E3_kcal_mol']);mv=[D(y)-D(x) for x,y in zip(l['mu3_Debye'],d['mu3_Debye'])];me=[x-y for x,y in zip(mv,refmu)]
            row['models'][label]={'Dy_minus_La_E3_kcal_mol':float(en),'Dy_minus_La_mu3_Debye':list(map(float,mv)),'energy_error_kcal_mol':float(en-refE),'dipole_error_Debye':list(map(float,me)),'dipole_error_norm_Debye':norm(me)}
        paired[cid]=row
    comp=r['comparison'];check(comp['cycle_count']==8 and comp['paired_count']==4 and comp['scientific_threshold_added'] is False,'Complete descriptive comparison denominator')
    for group,rows in (('cycles',cycles),('paired_Dy_minus_La',paired)):
        check(set(comp[group])==set(rows),'Exact comparison '+group)
        for cid,row in rows.items():
            got=comp[group][cid];check(got['status']=='complete' and set(got['models'])==set(labels),'Exact model labels '+cid)
            for label,model in row['models'].items():
                for k,v in model.items():close(v,got['models'][label][k],'Independent comparison '+cid+'/'+label+'/'+k)
            for k,v in row.items():
                if k not in ('models','model_terms'):close(v,got[k],'Independent comparison reference '+cid+'/'+k)
    review={'status':'passed_radial14_terminal_trimer_actual_review','implementation':pin(__file__),'binding':pin(HERE/'CONFIG.json'),'declaration':binding['declaration'],'result':binding['result'],'configuration':binding['configuration'],'producer':r['implementation'],'source_review':cfg['source_review'],'root_source_review':release['root_source_review'],'root_release':cfg['root_release'],'candidate':cfg['candidate'],'inventory':cfg['inventory'],'new_native_cell_admission_passed':True,'cycle_arithmetic_passed':True,'new_native_programs_reviewed':6,'native_induction_blocks_reviewed':12,'reused_selected_scalar_terms':14,'reused_static_scalar_terms':22,'raw_new_rows':raw_rows,'terms':terms,'cycles':cycles,'paired_Dy_minus_La':paired,'maximum_new_term_mu_reconstruction_difference_Debye':maxterm,'maximum_cycle_energy_difference_kcal_mol':maxE,'maximum_cycle_mu_difference_Debye':maxM,'terminal_metrics':r['terminal_metrics'],'checks':CHECKS,'checks_passed':len(CHECKS),'verified_pins':list(PINS.values()),'verified_pin_count':len(PINS),'new_native_calls':0,'new_optimizer_calls':0,'new_reference_calls':0,'scientific_threshold_added':False,'physical_accuracy_qualified':False,'model_promotion':False,'limits':['Only six new rows repeat the qualified saved parser; the14selected+22static scalar admissions are reused without fullfit/raw-old-cell replay.','Independent50-digit cycle arithmetic combines six raw-token reconstructions with36admitted stored scalars; no claim of42fresh raw re-admissions.','Consumed eight-cycle development comparison; no heldout, new physical threshold, unique mechanism or free-energy validation.','Original printing envelopes and PBE0/frozen-state/reference limitations remain.'],'cost':{'wall_seconds':time.monotonic()-tick,'CPU_seconds':time.process_time()-cpu,'threads':1,'producer_application':cost,'inherited_costs_not_added':True}}
    lines=['Terminal radial14 consumed trimer actual review: PASSED','Six new analyze rows and12internal induction blocks passed unchanged saved native admission.','Eight seven-term E3/mu3 cycles and four signed Dy-minus-La comparisons match independent arithmetic.','Fourteen selected and22static scalar terms reused their existing admissions; no selected46 re-audit.','Signed energy error is model minus electronic reference in kcal/mol; vector-error norm is Debye.','Model order: original baseline / feca12 / terminal radial14.','','Cycle | signed energy errors | dipole-vector error norms']
    for cid,row in cycles.items():lines.append(cid+' | '+' / '.join(f"{row['models'][l]['energy_error_kcal_mol']:.9f}" for l in labels)+' | '+' / '.join(f"{row['models'][l]['dipole_error_norm_Debye']:.9f}" for l in labels))
    lines+=['','Dy-minus-La contrasts | signed energy errors | dipole-vector error norms']
    for cid,row in paired.items():lines.append(cid+' | '+' / '.join(f"{row['models'][l]['energy_error_kcal_mol']:.9f}" for l in labels)+' | '+' / '.join(f"{row['models'][l]['dipole_error_norm_Debye']:.9f}" for l in labels))
    lines+=['','No new scientific acceptance threshold, accuracy admission or model promotion.','Producer actual cost: '+json.dumps(cost,sort_keys=True),'Saved review cost: '+json.dumps(review['cost'],sort_keys=True),'Maximum cycle arithmetic differences (kcal/mol, Debye): '+str([maxE,maxM])]
    (out/'REPORT.txt').write_text('\n'.join(lines)+'\n');review['report']=pin(out/'REPORT.txt')
    (out/'REVIEW.json').write_text(json.dumps(review,indent=2,allow_nan=False)+'\n')
    print(json.dumps({'review':pin(out/'REVIEW.json'),'report':review['report'],'checks_passed':len(CHECKS),'maximum_cycle_differences':[maxE,maxM],'cost':review['cost']}))

if __name__=='__main__':
    tick,cpu=time.monotonic(),time.process_time();binding=read(pin(HERE/'CONFIG.json'));out=Path(binding['output_directory']);out.mkdir(parents=True,exist_ok=False)
    try:main(binding,out)
    except BaseException as exc:
        (out/'FAILED_REVIEW.json').write_text(json.dumps({'status':'failed_radial14_terminal_trimer_actual_review','implementation':pin(__file__),'binding':pin(HERE/'CONFIG.json'),'error':{'type':type(exc).__name__,'message':str(exc)},'checks_completed':CHECKS,'new_native_calls':0,'cost':{'wall_seconds':time.monotonic()-tick,'CPU_seconds':time.process_time()-cpu,'threads':1}},indent=2)+'\n')
        raise
