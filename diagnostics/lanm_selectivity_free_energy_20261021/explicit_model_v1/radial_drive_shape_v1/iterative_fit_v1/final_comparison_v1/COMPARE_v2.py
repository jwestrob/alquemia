#!/usr/bin/env python3
"""Terminal-only scalar comparison and durable accounting. No scientific executor imports."""
import argparse
from collections import Counter
import datetime
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import time

THREADS=('OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','BLIS_NUM_THREADS')
FRACTIONS=(1.,.5,.25,.125)
COMPONENTS=('P8','B8','P9','B9')
ELECTRICAL=('pair_mu','C_mu','chi')
CHECKS=0
VERIFIED={}


def require(ok,message):
    global CHECKS
    CHECKS+=1
    if not ok: raise ValueError(message)


def pin(path):
    path=Path(path).resolve(); h=hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda:f.read(1024*1024),b''): h.update(chunk)
    return dict(path=str(path),sha256=h.hexdigest())


def verify(p):
    require(isinstance(p,dict) and set(p)=={'path','sha256'},'Exact file pin')
    key=(p['path'],p['sha256'])
    if key not in VERIFIED: require(pin(p['path'])==p,'Changed file '+p['path']); VERIFIED[key]=True
    return Path(p['path'])


def read(p): return json.loads(verify(p).read_text())
def token(v): return hashlib.sha256(json.dumps(v,sort_keys=True,separators=(',',':'),allow_nan=False).encode()).hexdigest()
def save(path,value): Path(path).write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')
def close(a,b,message):
    require(math.isfinite(a) and math.isfinite(b) and abs(a-b)<=1e-10*max(1.,abs(a),abs(b)),message+' arithmetic only')


def decisions(metrics,fixed):
    """Literal signs decide gates; closeness is never an acceptance tolerance."""
    e=metrics['component_errors_kcal_mol']
    jp=math.fsum(e[k]**2 for k in ('P8','P9'))/2; jb=math.fsum(e[k]**2 for k in ('B8','B9'))/2
    q=max(jp/fixed['q_denominators']['J_P'],jb/fixed['q_denominators']['J_B'])
    for a,b,n in ((jp,metrics['J_P'],'JP'),(jb,metrics['J_B'],'JB'),(q,metrics['q'],'q')):close(a,b,n)
    require(len(metrics['original25'])==25,'Exactly25 physical gates')
    physical={}
    for name,row in metrics['original25'].items():
        if name.startswith('energy/'):
            margin=1-abs(row['error_kcal_mol']); flag=margin>=0
            physical[name]=dict(passed=flag,margin_kcal_mol=margin)
        else:
            require(name.startswith('force/'),'Physical gate family')
            a=2-row['rigid_RMS_kcal_mol_A'];b=5-row['maximum_body_component_error_kcal_mol_A'];flag=a>=0 and b>=0
            physical[name]=dict(passed=flag,RMS_margin_kcal_mol_A=a,maximum_margin_kcal_mol_A=b)
        require(row['passed']==flag,'Literal original gate flag '+name)
    ph=all(x['passed'] for x in physical.values())
    require(ph==metrics['all25_original_physical_gates']==metrics['physical_constraints_pass'],'Physical25 meaning stays separate')
    electrical={}
    for name in ELECTRICAL:
        score=metrics['original_block_scores'][name];ceiling=fixed['electrical_actual_scores'][name]
        require(math.isfinite(score) and score>=0,'Finite nonnegative score')
        row=metrics['fixed_electrical_safeguards'][name];flag=score<=ceiling
        require(row==dict(score=score,ceiling=ceiling,passed=flag),'Fixed electrical scalar authority')
        electrical[name]=dict(score=score,ceiling=ceiling,passed=flag,score_margin=ceiling-score,RMS=math.sqrt(score),units='bohr^3' if name=='chi' else 'Debye')
    s=metrics['signed_W9_z_safeguard'];ceiling=fixed['signed_W9_absolute_ceiling_Debye'];r=s['model_minus_reference_Debye']
    require(s['reference_Dy_minus_La_z_Debye']==fixed['signed_W9_reference_Debye'] and s['absolute_ceiling_Debye']==ceiling,'Fixed signed W9 authority')
    signed=abs(r)<=ceiling
    require(s['passed']==signed,'Literal signed W9 gate')
    full=ph and all(x['passed'] for x in electrical.values()) and signed
    require(full==metrics['all_declared_safeguards_pass'],'Combined safeguards distinct from physical25')
    return dict(q=metrics['q'],J_P=metrics['J_P'],J_B=metrics['J_B'],
                independent_reconstruction=dict(q=q,J_P=jp,J_B=jb,decision_authority='producer literals after bounded consistency check'),
                component_RMS_kcal_mol=dict(P=math.sqrt(metrics['J_P']),B=math.sqrt(metrics['J_B'])),component_errors_kcal_mol=e,
                physical25=physical,physical25_pass=ph,physical25_pass_count=sum(x['passed'] for x in physical.values()),electrical=electrical,
                signed_W9=dict(s,absolute_margin_Debye=ceiling-abs(r)),all_declared_safeguards_pass=full,
                old_EF_score_descriptive=metrics['original_block_scores']['old_EF'])


def recorded_choice(current_q,slots):
    """Audit the existing choice; this never writes/changes a candidate."""
    best=current_q;chosen=None;epsilon=1e-10*max(1.,abs(current_q))
    require([x['fraction'] for x in slots]==list(FRACTIONS),'Four declared ordered fractions')
    require(all(not x['eligible'] or x.get('execution_status')=='complete' for x in slots),'No choice from incomplete eligible batch')
    for i,x in enumerate(slots):
        if x['eligible'] and x['actual']['all_declared_safeguards_pass'] and x['actual']['q']<best-epsilon: chosen=i;best=x['actual']['q']
    return chosen,best,epsilon


def parameter_report(values,start,domain):
    require(len(values)==len(start)==len(domain)==14,'Fourteen parameter coordinates')
    rows=[]
    for v,b,p in zip(values,start,domain):
        low,high=p['lower'],p['upper'];width=high-low
        require(math.isfinite(v) and low<=v<=high,'Literal original parameter box')
        rows.append(dict(name=p['name'],value=v,start=b,change=v-b,normalized_change=(v-b)/width,normalized_lower_distance=(v-low)/width,
                         normalized_upper_distance=(high-v)/width,at_literal_boundary=v in (low,high),units=p.get('units','dimensionless')))
    S=math.log(4.)
    native={metal:dict(c=S*(values[i+2]+values[i+3])/4,d=S*(values[i+2]-values[i+3])/4) for metal,i in (('La',0),('Dy',4))}
    return dict(box_coordinates=rows,mapped_native_collective_coefficients=native,box_is_development_domain_not_accuracy_threshold=True)


def vector_metrics(v,parameters,metrics,fixed,manifest,decomposition):
    require(v['status']=='complete' and (v['parameters']==parameters or v['parameters']+[0.,0.]==parameters),'Literal complete vector or admitted zero view')
    require(metrics['key']==token(parameters) and metrics['parameters']==parameters,'Metric/vector identity')
    require(len(v['cells'])==46 and set(v['cells'])==set(manifest['rows']) and all(x['status']=='complete' for x in v['cells'].values()),'Complete46 stored native-admission records')
    energies={k:x['energy_kcal_mol'] for k,x in manifest['static'].items()};energies.update({k:x['energy'] for k,x in v['cells'].items()})
    parents={k:{c:math.fsum(a*energies[n] for n,a in terms.items()) for c,terms in p['model_terms'].items()} for k,p in decomposition['parents'].items()}
    for n in (8,9):
        for c in ('P','B'):
            error=parents[f'aquo{n}_central_Dy'][c]-parents[f'aquo{n}_central_La'][c]-decomposition['paired_Dy_minus_La']['W'+str(n)]['reference_total_kcal_mol'][c]
            close(error,metrics['component_errors_kcal_mol'][c+str(n)],'Saved literal energy cycle')
    score=v['score'];old=score['old_EF'];gates={}
    for group,rows in old['energy_controls'].items():
        for name,r in rows.items():gates['energy/'+group+'/'+name]=dict(error_kcal_mol=r['error_kcal_mol'],passed=r['passed'])
    for kind,rows in old['per_fixture'].items():
        for name,r in rows.items():
            raw=r['raw_body_components_kcal_mol_A'];nw=8 if name=='aquo8_central' else 9 if name=='aquo9_central' else 1
            require(len(raw)==3+6*nw,'Original body-component denominator')
            weights=[1/6]*3+[1/(12*nw)]*(6*nw)
            close(math.sqrt(math.fsum(w*x*x for w,x in zip(weights,raw))),r['rigid_RMS_kcal_mol_A'],'Saved body RMS arithmetic')
            close(max(abs(x) for x in raw),r['maximum_body_component_error_kcal_mol_A'],'Saved body maximum')
            gates['force/'+kind+'/'+name]={k:r[k] for k in ('rigid_RMS_kcal_mol_A','maximum_body_component_error_kcal_mol_A','passed')}
    require(gates==metrics['original25'],'Exact stored25 quantities and metric flags')
    for a,b in (('old_EF','old_EF'),('pair_mu','pair'),('C_mu','C'),('chi','chi')):require(score['block_objectives'][b]==metrics['original_block_scores'][a],'Absolute original score')
    dip=v['response']['C_dipoles_Debye'];s=metrics['signed_W9_z_safeguard']
    require(s['model_Dy_minus_La_z_Debye']==math.fsum([dip['aquo9_central_Dy'][2],-dip['aquo9_central_La'][2]]),'Signed native stored W9 contrast')
    require(s['model_minus_reference_Debye']==math.fsum([dip['aquo9_central_Dy'][2],-dip['aquo9_central_La'][2],-fixed['signed_W9_reference_Debye']]),'Signed residual orientation')
    answer=decisions(metrics,fixed);answer['response']=v['response'];answer['component_parents_kcal_mol']=parents
    return answer


def accounting(text,job):
    rows=[x.split('|') for x in text.strip().splitlines() if x.strip()];header=rows.pop(0)
    if header[-1]=='':header=header[:-1]
    require(all(k in header for k in ('JobIDRaw','State','ExitCode','ElapsedRaw','AllocCPUS','TotalCPU','MaxRSS')),'Declared sacct columns')
    docs=[]
    for row in rows:
        if len(row)==len(header)+1 and row[-1]=='':row=row[:-1]
        require(len(row)==len(header),'Exact sacct row width');docs.append(dict(zip(header,row)))
    parents=[r for r in docs if r['JobIDRaw']==str(job)];require(len(parents)==1,'Exactly one parent allocation row')
    parent=parents[0];state=parent['State'].split()[0].rstrip('+')
    require(state in ('COMPLETED','FAILED','CANCELLED','TIMEOUT','OUT_OF_MEMORY','NODE_FAIL','PREEMPTED','BOOT_FAIL','DEADLINE','REVOKED'),'Terminal scheduler state')
    require(int(parent['AllocCPUS'])==13,'Actual13-CPU parent allocation')
    elapsed=int(parent['ElapsedRaw']);require(elapsed>=0,'Nonnegative parent elapsed')
    return dict(parent=parent,children=[r for r in docs if r is not parent],parent_allocated_CPU_seconds=elapsed*13,
                parent_cost_only_no_child_double_counting=True,total_CPU_literal=parent['TotalCPU'],MaxRSS_scope='Retain parent/step rows separately; blank is unavailable, not zero.')


def stream_journal(path,J,cfg,result,successful):
    state=J.State();previous='0'*64;seq=0;whole=hashlib.sha256();partial=b'';intents={};terminal=set();events=[];counts=Counter();request_events=Counter();cells={};initial=None
    with Path(path).open('rb') as f:
        for line in f:
            whole.update(line)
            if not line.endswith(b'\n'):partial=line;break
            row=json.loads(line);body={k:v for k,v in row.items() if k!='sha256'}
            require(set(row)=={'sequence','previous_sha256','event','sha256'} and row['sequence']==seq+1 and row['previous_sha256']==previous and J.sha(body)==row['sha256'],'Durable ordered journal hash')
            event=row['event'];kind=event['kind'];d=event['data'];state.apply(event);counts[kind]+=1
            if kind=='initialize':initial=d
            elif kind=='attempt_intent':intents[d['attempt']['number']]=d
            elif kind=='attempt_return':
                a=d['attempt'];n=a['number'];require(n not in terminal and d['program_inputs']==intents[n]['program_inputs'],'Unique return/intent inputs');terminal.add(n)
                for p in d['program_inputs'].values():verify(p)
                for receipt in d['receipt_delta']:
                    if 'output' in receipt:verify(receipt['output'])
                    if a['status']=='returned_success':
                        require(receipt['returncode']==0 and receipt['timeout_seconds']==300,'Inherited successful invocation gates')
                        require(receipt['executor']['job_id']==str(result['job_id']),'Actual receipt Slurm job')
                        require(receipt['argv']==[cfg['executables'][a['program']]['path']]+(['system.xyz','EM'] if a['program']=='analyze' else ['system.xyz','Y','N']) and receipt['cwd']==a['folder'],'Exact native program invocation')
                        require(all(v=='1' for v in receipt['threads'].values()),'One numerical thread per native program')
                    if 'started_utc' in receipt and 'ended_utc' in receipt:
                        start=datetime.datetime.fromisoformat(receipt['started_utc']);end=datetime.datetime.fromisoformat(receipt['ended_utc']);require(end>=start,'Receipt chronology');events.extend([(start,1),(end,-1)])
            elif kind=='request_cell':
                request_events[d['index']]+=1;ev=d['evidence']
                if 'result' in ev:cells[(ev['result']['path'],ev['result']['sha256'])]=ev['result']
            seq+=1;previous=row['sha256']
    require(seq>0 and initial['policy']==J.POLICY,'Initialized inherited journal')
    require(initial['authority']==dict(implementation=cfg['adapter'],configuration_fingerprint=result['configuration_fingerprint'],root_release=cfg['root_release']),'Exact journal initialization authority')
    ledger=state.ledger;attempts=ledger.get('native_attempts',[])
    require(ledger['limits']==cfg['limits'] and ledger['limits']['new_parameter_vectors']==228 and ledger['limits']['new_native_programs']==7280,'Current228/7280 journal limits')
    require(len(attempts)==len(intents)==ledger['native_program_invocations']==result['native_programs']<=7280,'Unique charged native attempts')
    require(len({(x['cache_key'],x['program']) for x in attempts})==len(attempts),'No repeated native program identity')
    require(len(ledger['vector_keys'])==len(set(ledger['vector_keys']))==len(ledger['requests'])<=228,'Unique new-vector request count')
    require(ledger['native_returncodes_observed']==sum(bool(a.get('returncode_observed')) for a in attempts),'Observed returns separate from attempts')
    for i,r in enumerate(ledger['requests']):
        require(r['vector_key']==ledger['vector_keys'][i] and token(r['values'])==r['vector_key'] and len(r['cells'])==46,'Exact unique46-cell request')
        charged=sum(x['request_index']==i for x in intents.values())
        require(charged<=r['required_new_programs'],'Whole-vector planned ceiling')
        if r['status']=='complete':require(charged==r['required_new_programs'] and request_events[i]==46,'Completed vector accounting')
    for p in cells.values():verify(p)
    overlap=maximum=0
    for _,d in sorted(events,key=lambda x:(x[0],x[1])):overlap+=d;maximum=max(maximum,overlap)
    require(overlap==0 and maximum<=12,'Receipt native overlap at most12')
    unresolved=[n for n,a in state.attempts.items() if a['status']=='invoke_boundary_entered']
    if successful:
        require(not partial and not unresolved and set(intents)==terminal and not ledger['stopped'],'Successful complete durable ledger')
        require(all(a['status']=='returned_success' for a in attempts) and all(r['status']=='complete' for r in ledger['requests']),'Successful invocation/request statuses')
    summary=dict(journal=dict(path=str(Path(path).resolve()),sha256=whole.hexdigest()),events=seq,event_counts=dict(counts),last_event_sha256=previous,partial_suffix_bytes=len(partial),
                 unresolved_attempts=unresolved,charged_program_attempts=len(attempts),observed_returns=ledger['native_returncodes_observed'],
                 returned_success=sum(a['status']=='returned_success' for a in attempts),failed_attempts=[a for a in attempts if a['status']=='failed_native_attempt'],
                 unique_new_vectors=len(ledger['vector_keys']),request_status_counts=dict(Counter(r['status'] for r in ledger['requests'])),
                 cell_status_counts=dict(Counter(c['status'] for r in ledger['requests'] for c in r['cells'].values())),
                 immutable_import_cell_occurrences=sum(c.get('immutable_import',False) for r in ledger['requests'] for c in r['cells'].values()),
                 unique_new_cell_artifacts=len(cells),unique_requested_cell_keys=len({x['cache_key'] for r in ledger['requests'] for x in r['cells'].values() if 'cache_key' in x}),
                 unique_immutable_import_cell_keys=len({x['cache_key'] for r in ledger['requests'] for x in r['cells'].values() if x.get('immutable_import')}),maximum_native_overlap=maximum,automatic_retry_permitted=False)
    return ledger,summary


def run(c,out):
    require(c['implementation']==pin(__file__),'Frozen comparator')
    declaration=read(c['declaration']);cfg=read(c['campaign_configuration']);root=read(c['root_release']);review=read(c['source_review'])
    for k in ('implementation','declaration','dependencies','campaign_configuration','result','journal','ledger','ledger_sidecar','scheduler_accounting','output_directory','source_review'):
        require(root[k]==c[k],'Exact terminal release '+k)
    require(root['status']=='released_one_saved_radial14_terminal_comparison' and review['status']=='passed_radial14_terminal_comparison_source_review' and review['implementation']==c['implementation'],'Source review and one terminal release')
    require(str(out)==c['output_directory'],'Fresh declared comparison output')
    require(c['campaign_configuration']==declaration['inputs']['CONFIG_RUN.json'],'Declared campaign config')
    submission=read(declaration['inputs']['SUBMISSION.json']);require(submission['configuration']==c['campaign_configuration'] and submission['job_id']=='20605','Single declared job20605')
    scheduler=accounting(verify(c['scheduler_accounting']).read_text(),submission['job_id'])
    r=read(c['result']);require(c['result']['path']==declaration['terminal_result_path'] and r['status']!='started','Exact terminal producer result')
    require(r['configuration']==c['campaign_configuration'] and r['implementation']==cfg['implementation'] and r['limits']==cfg['limits'],'Producer/source/limits joins')
    require(r['configuration_fingerprint']==token({k:v for k,v in cfg.items() if k not in ('root_release','source_review','preparation','fingerprint')}),'Exact campaign fingerprint')
    require(cfg['root_release']==declaration['inputs']['ROOT_RELEASE.json'],'Declared producer release')
    for k in ('implementation','adapter','core','settings','root_release','source_review','reuse_inventory','root_cache_bridge'):verify(cfg[k])
    for executable in cfg['executables'].values():verify(executable)
    require(r['job_id']==submission['job_id'] and r['cost']['execution_mode']=='slurm','Actual scheduler producer identity')
    successful=r['status']=='bounded_radial14_iterative_fit_complete' and r['error'] is None
    if successful:require(scheduler['parent']['State']=='COMPLETED' and scheduler['parent']['ExitCode']=='0:0','Successful result and scheduler agree')
    inv=read(cfg['reuse_inventory']);fixed=inv['fixed_authorities'];registry={x['key']:x for x in inv['entries']}
    dependencies=read(c['dependencies']);authorities=read(cfg['inputs'])['sources'];parent=read(cfg['parent_configuration']);manifest=read(parent['manifest']);dec=read(dependencies['decomposition'])
    require(dependencies['journal']==authorities['journal'] and dependencies['decomposition']==authorities['decomposition'],'Exact inherited scalar/journal authorities')
    domain=read(parent['domain'])['parameters']+[dict(name=n,lower=-1.,upper=1.,units='dimensionless') for n in ('La_radial_drive','Dy_radial_drive')]
    spec=importlib.util.spec_from_file_location('saved_only_journal',verify(dependencies['journal']));J=importlib.util.module_from_spec(spec);spec.loader.exec_module(J)
    ledger,journal=stream_journal(verify(c['journal']),J,cfg,r,successful)
    disk=read(c['ledger']);side=read(c['ledger_sidecar']);reconciled=disk==ledger
    require(side['policy']==J.POLICY and side['journal_path']==c['journal']['path'],'Exact sidecar policy/path')
    require(journal['journal']==c['journal'],'Journal bytes unchanged during streaming review')
    if successful or r.get('evaluator_snapshot_available') is not False:
        require(reconciled and side['ledger']==c['ledger'] and side['sequence']==journal['events'] and side['last_event_sha256']==journal['last_event_sha256'],'Exact terminal journal/ledger/sidecar equality')
        require(r['evaluator']['ledger']==c['ledger'] and r['evaluator']['active_limits']==cfg['limits'],'Active current snapshot limits')
    journal['terminal_compatibility_ledger_reconciled']=reconciled
    journal['legacy_snapshot_labels_not_current_budget']={k:r.get('evaluator',{}).get(k) for k in ('max_parameter_vectors','max_native_programs')}
    require(set(r['vectors'])-set(registry)>=set(ledger['vector_keys']),'Charged vector identities retained in producer')
    require(set(k for k,v in r['vectors'].items() if not v.get('inherited',False) and v['status']=='complete')<=set(ledger['vector_keys']),'No uncharged new complete vectors')
    for request in ledger['requests']:
        row=r['vectors'][request['vector_key']];require(row['parameters']==request['values'] and request['role']['parameters']==request['values'],'Request/vector parameter identity')
        role={k:request['role'][k] for k in ('outer_step','role')};require(role in row['requested_roles'] and role['role'] in ('probe','trial') and 1<=role['outer_step']<=8,'Explicit bounded request role')
        if request['status']=='complete':require(row['status']=='complete','Committed request has completed artifact')
    for key,row in r['vectors'].items():
        require(token(row['parameters'])==key,'All reported vector identities')
        if 'artifact' in row:verify(row['artifact'])
        if key in registry:require(row['artifact']==registry[key]['artifact'] and row['parameters']==registry[key]['parameters'] and row['inherited'],'Original seed provenance preserved')
    start=registry[inv['origin_key']]['parameters'];cache={}
    def point(key,metrics):
        row=r['vectors'][key]
        if key not in cache:
            v=read(row['artifact']);m=vector_metrics(v,row['parameters'],metrics,fixed,manifest,dec)
            m.update(key=key,artifact=row['artifact'],inherited=row.get('inherited',False),parameters=parameter_report(row['parameters'],start,domain));cache[key]=m
        else:require(cache[key]['q']==metrics['q'],'Consistent repeated point metrics')
        return cache[key]
    current=inv['origin_key'];initial=point(current,r['initial_metrics']);trajectory=[dict(step=0,point=initial)];trials=[];proposals=[]
    durable=[read(p) for p in r.get('proposals',[])]
    require(sum(p['solver_calls'] for p in durable)==r['solver_calls']<=7 and len(durable)==r['logical_proposals']<=8 and sum('imported_saved_proposal' in p for p in durable)==r['imported_proposals']<=1,'Imported/new proposal count')
    for i,it in enumerate(r['iterations'],1):
        require(it['outer_step']==i and i<=8 and it['origin_key']==current,'Sequential accepted-origin chain')
        if 'proposal' not in it:
            require(not successful or it['status']=='hard_budget_exhausted','Unfinished proposal is failure/budget only');break
        p=read(it['proposal']);require(p['origin_key']==current and p['derivative_keys']==it['derivative_keys'],'Proposal current/stencil identities')
        require(p['solver_calls']==(0 if i==1 else 1),'Exactly first proposal reuse')
        if i==1:require(p['imported_saved_proposal']==inv['saved_proposal']['result'] and it['derivative_keys']==inv['initial_derivative_keys'],'Exact first28/saved solution authority')
        for k in ('linearization','conic_data','row_layout'):verify(p[k])
        if 'attempt' in p:verify(p['attempt'])
        require(p['solver_success']==(p['solver_status']=='Solved'),'Raw solver status not upgraded')
        require(len(it['trial_slots'])==len(p['trial_slots']) and (not p['trial_slots'] or [x['fraction'] for x in p['trial_slots']]==list(FRACTIONS)),'Complete declared fraction denominator')
        proposals.append(dict(step=i,status=p['status'],solver_status=p['solver_status'],solver_success=p['solver_success'],solver_calls=p['solver_calls'],
                              imported_saved_proposal=p.get('imported_saved_proposal'),independent_residuals=p.get('independent_residuals'),info=p.get('info'),new_solver_cost=p.get('new_solver_cost')))
        for slot,declared in zip(it['trial_slots'],p['trial_slots']):
            require(all(slot[k]==declared[k] for k in ('key','parameters','fraction','eligible')),'Literal prescribed trial slot')
            record=r['vectors'].get(slot['key'],{})
            report=dict(step=i,key=slot['key'],fraction=slot['fraction'],eligible=slot['eligible'],execution_status=slot['execution_status'],
                        vector_record_status=record.get('status','not_requested'),vector_record_error=record.get('error'),prediction=declared['prediction'],producer_selected=False)
            if slot['execution_status']=='complete':
                require(slot['eligible'] and slot['native_admitted'] and slot['artifact']==r['vectors'][slot['key']]['artifact'],'Complete eligible trial lineage')
                actual=point(slot['key'],slot['actual']);pred=declared['prediction'];delta=dict(q=actual['q']-pred['actual_component_q'],component_errors_kcal_mol={k:actual['component_errors_kcal_mol'][k]-pred['component_errors_kcal_mol'][k] for k in COMPONENTS},signed_W9_residual_Debye=actual['signed_W9']['model_minus_reference_Debye']-pred['signed_W9_z_residual_Debye'])
                require(delta==slot['actual_minus_predicted'],'Literal actual-minus-predicted report');report.update(actual=actual,actual_minus_predicted=delta)
            trials.append(report)
        if 'selection' in it:
            chosen,q,eps=recorded_choice(cache[current]['q'],it['trial_slots']);s=it['selection']
            require(s['current_metrics']['key']==current and s['current_metrics']['q']==cache[current]['q'],'Selection uses current incumbent')
            require(chosen==s['selected_trial_index'] and eps==s['strict_q_improvement_tolerance'],'Recorded complete-batch choice/epsilon')
            if chosen is not None:
                slot=it['trial_slots'][chosen];require(it['status']=='complete' and it['selected_key']==s['selected_key']==slot['key'],'Producer accepted point')
                current=slot['key'];require(cache[current]['q']==q and cache[current]['all_declared_safeguards_pass'],'Accepted safeguards/q')
                trajectory.append(dict(step=i,point=cache[current]));next(x for x in reversed(trials) if x['step']==i and x['key']==current)['producer_selected']=True
            else:require(it['status']=='no_actual_improvement','No-improvement terminal status')
        elif it['status']!='no_eligible_fraction':require(not successful or it['status']=='hard_budget_exhausted','No selection from incomplete batch')
    endpoint_key=r.get('selected_key') if successful else r.get('last_complete_incumbent_key',r.get('current_key'))
    require(endpoint_key==current,'Final/last complete incumbent preserves producer chain')
    endpoint_metrics=r['selected_metrics'] if successful else r.get('last_complete_incumbent_metrics')
    if endpoint_metrics is not None:point(current,endpoint_metrics)
    require(len(trajectory)-1==sum(it['status']=='complete' for it in r['iterations']),'Accepted-step count')
    stop=r['stopping_reason']
    if successful:
        require(r['selected']==r['vectors'][current]['artifact'] and r['selected_parameters']==r['vectors'][current]['parameters'],'Final artifact/parameter pin')
        require(r['completed_outer_steps']==len(trajectory)-1,'Producer completed-step count')
        if stop=='outer_step_limit':require(len(r['iterations'])==8 and len(trajectory)==9,'Eight-step finite stop')
        elif stop=='no_actual_safeguarded_improvement':require(r['iterations'][-1]['status']=='no_actual_improvement','No-improvement stop')
        elif stop=='no_eligible_fraction':require(r['iterations'][-1]['status']=='no_eligible_fraction','No eligible stop')
        elif stop=='development_target_met':require(endpoint_metrics['full_development_target'],'Declared target stop')
        else:require(stop=='hard_budget_exhausted','Declared stop only')
    cost=r['cost']
    for k in ('wall_seconds','controller_CPU_seconds','child_user_CPU_seconds','child_system_CPU_seconds'):require(math.isfinite(cost[k]) and cost[k]>=0,'Measured producer cost')
    return dict(status='completed_saved_radial14_terminal_comparison',producer_status=r['status'],producer_error=r.get('error'),producer_completed=successful,
                implementation=c['implementation'],configuration=c['campaign_configuration'],result=c['result'],declaration=c['declaration'],source_review=c['source_review'],root_release=c['root_release'],
                stopping_reason=stop,trajectory=trajectory,eligible_and_other_trial_outcomes=trials,proposal_statuses=proposals,
                all_durable_proposal_records=[{k:p.get(k) for k in ('status','error','outer_step','origin_key','solver_calls','solver_constructions','solver_status','solver_success','attempt','new_solver_cost','imported_saved_proposal')} for p in durable],endpoint_key=current,
                endpoint_is_producer_choice_only=True,candidate_selected_by_comparator=False,journal=journal,limits=cfg['limits'],
                vector_status_counts=dict(Counter(v['status'] for v in r['vectors'].values())),seed_vectors_available=len(registry),
                seed_vectors_used_for_explicit_probe_or_trial=sum(bool(r['vectors'][k]['requested_roles']) for k in registry),seed_origin_reused_separately=1,
                producer_cost=cost,parent_scheduler_accounting=scheduler,measured_process_CPU_seconds=math.fsum(cost[k] for k in ('controller_CPU_seconds','child_user_CPU_seconds','child_system_CPU_seconds')),
                scope=['Stored scalar physics/score consistency and exact acceptance decisions; no full raw scientific output or Jacobian replay.',
                       'Fresh native input/output and cell artifact hashes checked through one journal pass; inherited scientific admissions retained.',
                       'Electrical nonincrease and box faces are development safeguards, not experimental accuracy tolerances.',
                       'Consumed development only; no global convergence, heldout accuracy, periodic admission or free-energy claim.',
                       'Solver status retained literally; feasible fractions do not upgrade AlmostSolved.'],
                native_calls=0,solver_calls=0,Jacobian_replays=0,model_admission=False)


def main():
    p=argparse.ArgumentParser();p.add_argument('--config',required=True);p.add_argument('--out',required=True);a=p.parse_args()
    require(all(os.environ.get(k)=='1' for k in THREADS),'One numerical thread')
    out=Path(a.out).resolve();out.mkdir(parents=True,exist_ok=False);tick,cpu=time.monotonic(),time.process_time();c=json.loads(Path(a.config).read_text())
    try:r=run(c,out)
    except Exception as e:r=dict(status='failed_saved_radial14_terminal_comparison',error=repr(e),configuration=pin(a.config),implementation=pin(__file__),candidate_selected_by_comparator=False,native_calls=0,solver_calls=0)
    r.update(comparison_configuration=pin(a.config),checks=CHECKS,comparison_cost=dict(wall_seconds=time.monotonic()-tick,CPU_seconds=time.process_time()-cpu,threads=1))
    save(out/'RESULT.json',r)
    lines=[r['status']]
    if 'trajectory' in r:
        lines+=['Producer: '+r['producer_status']+'; stop: '+r['stopping_reason'],'step | q | P RMS kcal/mol | B RMS kcal/mol | physical25 | three electrical | signed W9']
        for x in r['trajectory']:
            m=x['point'];lines.append(f"{x['step']} | {m['q']:.10g} | {m['component_RMS_kcal_mol']['P']:.8g} | {m['component_RMS_kcal_mol']['B']:.8g} | {m['physical25_pass_count']}/25 | {sum(v['passed'] for v in m['electrical'].values())}/3 | {m['signed_W9']['passed']}")
        lines+=['Raw solver statuses: '+', '.join(x['solver_status'] for x in r['proposal_statuses']),f"New charged native attempts: {r['journal']['charged_program_attempts']}; new vectors: {r['journal']['unique_new_vectors']}; parent allocated CPU-s: {r['parent_scheduler_accounting']['parent_allocated_CPU_seconds']}",'Electrical/box safeguards are developmental. No new candidate selection, physical accuracy or free-energy admission.']
    else:lines.append(r['error'])
    (out/'REPORT.txt').write_text('\n'.join(lines)+'\n');print(json.dumps(dict(status=r['status'],result=pin(out/'RESULT.json'),error=r.get('error'))));return int(r['status'].startswith('failed'))


if __name__=='__main__':raise SystemExit(main())
