#!/usr/bin/env python3
"""Terminal-only bounded saved audit; no optimizer or native invocations."""
import argparse
import copy
from collections import Counter
import datetime
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import time
from types import SimpleNamespace
import numpy as np

CHECKS=0
VERIFIED={}
PIN_CACHE={}
JSON_CACHE={}
PROGRESS={}
THREADS=('OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','BLIS_NUM_THREADS')
def require(ok,msg):
    global CHECKS
    CHECKS+=1
    if not ok:raise ValueError(msg)
def pin(p):
    p=Path(p).resolve();st=p.stat();key=(str(p),st.st_ino,st.st_size,st.st_mtime_ns)
    if key not in PIN_CACHE:PIN_CACHE[key]={'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
    return dict(PIN_CACHE[key])
def verify(p):
    k=(p['path'],p['sha256'])
    if k not in VERIFIED:require(pin(p['path'])==p,'Changed pin '+p['path']);VERIFIED[k]=True
    return Path(p['path'])
def read(p):
    key=(p['path'],p['sha256'])
    if key not in JSON_CACHE:JSON_CACHE[key]=json.loads(verify(p).read_text())
    return JSON_CACHE[key]
def token(v):return hashlib.sha256(json.dumps(v,sort_keys=True,separators=(',',':'),allow_nan=False).encode()).hexdigest()
def module(p,name):
    s=importlib.util.spec_from_file_location(name,verify(p));m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
def canon(v):return json.loads(json.dumps(v))
def save(p,v):p.write_text(json.dumps(v,indent=2,allow_nan=False)+'\n')
def forbidden(*args,**kwargs):raise ValueError('Audit forbids optimizer/native execution')
def close(a,b,msg):
    require(math.isfinite(a) and math.isfinite(b) and abs(a-b)<=1e-10*max(1.,abs(a),abs(b)),msg+' arithmetic')
    return abs(a-b)


def independent25(v,data):
    model=v['model'];old=v['score']['old_EF'];gates={};maxdelta=0.
    for kind in ('La','Dy','paired'):
        for geo in ('water_minus','water_central','water_plus','aquo8_central','aquo9_central'):
            la,dy=data[geo+'_La'],data[geo+'_Dy'];n=len(la['geometry']['positions_A']);residual=[]
            for i in range(n):
                residual.append([(model[geo+'_Dy']['forces'][i][j]-model[geo+'_La']['forces'][i][j])-(dy['reference']['forces'][i][j]-la['reference']['forces'][i][j]) if kind=='paired' else model[geo+'_'+kind]['forces'][i][j]-data[geo+'_'+kind]['reference']['forces'][i][j] for j in range(3)])
            raw=[];xyz=la['geometry']['positions_A']
            for category in ('ion','water'):
                for body in (b for b in la['bodies'] if b['category']==category):
                    ids=body['indices'];center=[math.fsum(xyz[i][j] for i in ids)/len(ids) for j in range(3)]
                    raw.extend(math.fsum(residual[i][j] for i in ids) for j in range(3))
                    if category=='water':
                        for a,b in ((1,2),(2,0),(0,1)):raw.append(math.fsum((xyz[i][a]-center[a])*residual[i][b]-(xyz[i][b]-center[b])*residual[i][a] for i in ids))
            waters=(n-1)//3;weights=[1/6]*3+[1/(12*waters)]*(6*waters);rms=math.sqrt(math.fsum(w*x*x for w,x in zip(weights,raw)));maximum=max(abs(x) for x in raw);stored=old['per_fixture'][kind][geo]
            for a,b in zip(raw,stored['raw_body_components_kcal_mol_A']):maxdelta=max(maxdelta,close(a,b,'Independent body projection'))
            close(rms,stored['rigid_RMS_kcal_mol_A'],'Body RMS');close(maximum,stored['maximum_body_component_error_kcal_mol_A'],'Body maximum')
            require(rms<=2 and maximum<=5 and stored['passed'] is True,'Selected original body gate')
            gates['force/'+kind+'/'+geo]={'RMS':rms,'maximum':maximum,'passed':True}
    def E(x,g,m):return x[g+'_'+m]['energy']
    refs={k:r['reference'] for k,r in data.items()};errors={}
    for g in ('water_minus','water_central','water_plus','aquo9_central'):
        errors['balanced_exchange/'+g]=(E(model,g,'Dy')-E(model,g,'La')-E(model,'aquo8_central','Dy')+E(model,'aquo8_central','La'))-(E(refs,g,'Dy')-E(refs,g,'La')-E(refs,'aquo8_central','Dy')+E(refs,'aquo8_central','La'))
    for m in ('La','Dy'):
        for s in ('minus','plus'):errors['individual_radial_work/'+m+'_'+s]=(E(model,'water_'+s,m)-E(model,'water_central',m))-(E(refs,'water_'+s,m)-E(refs,'water_central',m))
    for s in ('minus','plus'):errors['paired_radial_work/'+s]=errors['individual_radial_work/Dy_'+s]-errors['individual_radial_work/La_'+s]
    for name,error in errors.items():
        group,key=name.split('/');stored=old['energy_controls'][group][key];close(error,stored['error_kcal_mol'],'Independent energy cycle');require(abs(error)<=1 and stored['passed'] is True,'Selected energy gate');gates['energy/'+name]={'error':error,'passed':True}
    require(len(gates)==25 and v['score']['original_EF_checks_passed'],'All25 physical gates')
    return {'gates':gates,'maximum_body_arithmetic_difference':maxdelta}


def stream_journal(path,journal):
    """One forward read. Preserve a partial suffix; never repair or resume it."""
    state=journal.State();previous='0'*64;sequence=0;partial=b'';intents={};returns={};cell_pins={};initial=None
    whole=hashlib.sha256()
    with Path(path).open('rb') as handle:
        for line in handle:
            whole.update(line)
            if not line.endswith(b'\n'):partial=line;break
            rec=json.loads(line)
            require(set(rec)=={'sequence','previous_sha256','event','sha256'},'Journal exact record schema')
            body={k:v for k,v in rec.items() if k!='sha256'}
            require(rec['sequence']==sequence+1 and rec['previous_sha256']==previous and journal.sha(body)==rec['sha256'],'Journal committed sequence/hash')
            event=rec['event'];kind=event['kind'];data=event['data'];state.apply(event)
            if kind=='initialize':initial=data
            elif kind=='attempt_intent':intents[data['attempt']['number']]=data
            elif kind=='attempt_return':returns[data['attempt']['number']]=data
            elif kind=='request_cell' and 'result' in data['evidence']:
                pp=data['evidence']['result'];cell_pins[(pp['path'],pp['sha256'])]=pp
            sequence+=1;previous=rec['sha256']
    require(sequence>0,'Committed journal initialization')
    summary={'journal':{'path':str(Path(path).resolve()),'sha256':whole.hexdigest()},'sequence':sequence,'last_event_sha256':previous,
             'partial_suffix_bytes':len(partial),'partial_suffix_sha256':hashlib.sha256(partial).hexdigest() if partial else None,
             'unresolved_attempt_numbers':[n for n,x in state.attempts.items() if x['status']=='invoke_boundary_entered'],
             'attempt_status_counts':{},'request_status_counts':{},'cell_status_counts':{},'automatic_retry_permitted':False}
    for x in state.ledger.get('native_attempts',[]):summary['attempt_status_counts'][x['status']]=summary['attempt_status_counts'].get(x['status'],0)+1
    for x in state.ledger['requests']:
        summary['request_status_counts'][x['status']]=summary['request_status_counts'].get(x['status'],0)+1
        for y in x['cells'].values():summary['cell_status_counts'][y['status']]=summary['cell_status_counts'].get(y['status'],0)+1
    return state.ledger,summary,initial,intents,returns,cell_pins


def journal_contract(r,cfg,journal,out,successful):
    lp=r['evaluator']['ledger'];ledger=read(lp);path=Path(lp['path']).with_name('EVENTS.jsonl')
    replay,summary,initial,intents,returns,cell_pins=stream_journal(path,journal)
    sidepin=pin(Path(lp['path']).with_name(Path(lp['path']).name+'.journal.json'));side=read(sidepin)
    require(side['policy']==cfg['persistence_policy'] and side['journal_path']==str(path),'Journal policy/path')
    require(side['ledger']==lp and side['sequence']==summary['sequence'] and side['last_event_sha256']==summary['last_event_sha256'],'Final legacy ledger committed sidecar')
    require(replay==ledger,'Single journal replay equals exact final compatibility ledger')
    require(initial['authority']=={'configuration_fingerprint':r['configuration_fingerprint'],'implementation':cfg['adapter'],'journal_implementation':cfg['journal_implementation'],'root_release':cfg['root_release']},'Journal initialization authority')
    require(initial['policy']==cfg['persistence_policy'],'Journal initialization policy')
    require(ledger['limits']==cfg['limits'] and ledger['configuration_fingerprint']==r['configuration_fingerprint'],'Journal limits/fingerprint')
    attempts=ledger.get('native_attempts',[])
    require(len(attempts)==ledger['native_program_invocations']==r['native_programs']==r['evaluator']['native_program_invocations']<=6496,'Actual invocation-boundary count/cap')
    require(len(ledger['vector_keys'])==r['new_parameter_vectors']==r['evaluator']['new_parameter_vectors']<=200,'Actual new-vector count/cap')
    require(ledger['native_returncodes_observed']==sum(bool(x.get('returncode_observed')) for x in attempts)==r['evaluator']['native_returncodes_observed'],'Observed return count')
    require(set(intents)==set(range(1,len(attempts)+1)),'Every charged attempt has exactly one durable intent')
    if successful:
        require(summary['partial_suffix_bytes']==0 and not summary['unresolved_attempt_numbers'],'Successful journal fully committed and all attempts terminal')
        require(not ledger['stopped'] and not r['evaluator']['stopped'],'Successful result has no hidden stop')
        require(set(returns)==set(intents),'Every successful intent has terminal return evidence')
        require(all(x['status']=='returned_success' and x['returncode_observed'] for x in attempts),'All actual attempts successful')
    summary.update(ledger=lp,sidecar=sidepin,new_vectors=len(ledger['vector_keys']),native_programs=len(attempts),native_returncodes_observed=ledger['native_returncodes_observed'])
    save(out/'JOURNAL_REVIEW.json',summary);PROGRESS['journal_review']=pin(out/'JOURNAL_REVIEW.json')
    return ledger,summary,intents,returns,cell_pins


def check_proposal(p,problem,fixed,core,cfg,step,expectedkeys):
    """Reconstruct secants and scalar witnesses; never instantiate a solver."""
    require(p['origin_key']==p['current_key']==core.token(problem['center'].tolist()) and p['derivative_keys']==expectedkeys,'Proposal exact origin/stencil')
    require(p['derivative_spans']==problem['spans'] and p['Jacobian_original531']==problem['Jacobian_original531'].tolist() and p['Jacobian_components']==problem['Jacobian_components'].tolist(),'Exact saved secant arrays')
    require(p['fixed_ceilings']==fixed['ceilings'] and p['fixed_component_denominators']==fixed['denominators'] and p['current_metrics']==problem['current_metrics'],'Fixed diagnostic ceilings/denominators/origin')
    require(p['old_score_ceilings_enforced'] is False and p['new_native_calls']==0,'New physical-only policy')
    model,removed,indices=core.removal.remove_joint_ceilings(core.conic.assemble(problem))
    require(p['removed_cones']==removed and p['retained_original_row_indices']==indices,'Exactly four old-score cones removed')
    settings=core.conic.frozen_settings(core.DEPENDENCIES['settings'])
    require(p['effective_settings']==core.conic.attributes(settings),'Frozen Clarabel settings')
    if step==1:
        require(p['solver_calls']==0 and p['imported_saved_proposal']==cfg['saved_initial_proposal'],'First proposal imported without repeated solve')
        require(p['imported_saved_review']==core.DEPENDENCIES['saved_initial_review'],'Imported actual arithmetic authority')
        saved=read(cfg['saved_initial_proposal']);review=read(p['imported_saved_review'])
        require(review['status']=='passed_ceiling_tradeoff_actual_arithmetic_review' and review['result']==cfg['saved_initial_proposal'] and review['certificate_checks_passed'],'Imported proposal reviewed')
        require(saved['origin_key']==p['origin_key'] and saved['derivative_keys']==expectedkeys and saved['effective_settings']==p['effective_settings'],'Imported exact origin/stencil/settings')
        original=core.removal.problem_arrays(read(saved['linearized_problem']))
        for name in ('center','lower','upper','width','z','bounds_lower','bounds_upper','spans','Jacobian_original531','Jacobian_components'):
            require(np.array_equal(problem[name],original[name]),'Imported exact problem '+name)
        require(problem['q_current']==original['q_current'],'Imported exact q origin')
        with np.load(verify(saved['conic_data']),allow_pickle=False) as z:
            for name,value in (('A_data',model['A'].data),('A_indices',model['A'].indices),('A_indptr',model['A'].indptr),('b',model['b']),('c',model['c'])):require(np.array_equal(z[name],value),'Imported exact conic '+name)
            require(tuple(z['A_shape'])==model['A'].shape,'Imported exact conic dimensions')
        for name in ('solution','info','solver_status','solver_success','independent_linear_residuals'):require(p[name]==saved[name],'Imported saved result '+name)
    else:
        require(p['solver_calls'] in (0,1) and 'imported_saved_proposal' not in p,'At most one later new solve')
        # Zero calls later can only be an explicit failure before solver.solve().
        if p['solver_calls']==0:require(p['status']=='proposal_failure_without_qualified_return' and not p['trial_slots'],'Pre-solve failure retained')
    if 'solver_status' in p:
        require(p['solver_success']==(p['solver_status']=='Solved'),'Solver success never relabeled from arithmetic feasibility')
        require(p['solution']['status']==p['solver_status'] and p['info']['status']==p['solver_status'],'Raw solution/info statuses preserved')
    if p['status']=='proposal_failure_without_qualified_return':
        require(p.get('error') and not p['trial_slots'],'Failed proposal has no eligible trial');return []
    require(p['solver_calls']==(0 if step==1 else 1),'Exactly one solve per later complete proposal')
    x=np.asarray(p['solution']['x'],dtype=float);require(x.shape==(13,) and np.all(np.isfinite(x)) and x[-1]>=0,'Finite conic return')
    require(p['returned_q']==float(x[-1]**2),'Conic sqrt(q) converted before fractional q certification')
    correspondence=core.removal.correspondence(problem,model,x,SimpleNamespace(constraints=core.constraints,near=core.near))
    stored=p['retained_correspondence']
    require(set(stored)==set(correspondence) and stored['retained_source_inequalities']==847 and stored['nonnegative_t']==correspondence['nonnegative_t'],'Retained correspondence schema')
    # Existing CORE.near is arithmetic-only; candidate inequalities stay exact.
    core.near([stored['maximum_arithmetic_difference'],stored['minimum_retained_source_margin']],[correspondence['maximum_arithmetic_difference'],correspondence['minimum_retained_source_margin']],'Retained correspondence arithmetic')
    raw=np.r_[x[:12],x[-1]**2]
    witnesses=[core.witness(problem,raw,f,p['solver_success'],SimpleNamespace(certify_fraction=core.certify_fraction),core.removal.REMOVED,problem['q_current']) for f in core.FRACTIONS]
    require(witnesses==p['witnesses'],'All exact fraction scalar witnesses, no solve')
    slots=[]
    for w in witnesses:
        cc=w['constructed_source_certificate']
        slots.append({'fraction':w['fraction'],'parameters':cc['parameters'],'key':cc['key'],'eligible':w['feasible_gain_exceeds_numerical_comparison'],
                      'certificate':dict(cc['certificate'],strict_retained_feasible=w['strict_retained_feasible'],minimum_retained_margin=w['minimum_retained_margin'],retained_scalar_margin_count=w['retained_scalar_margin_count'],actual_component_rho=w['actual_component_rho'],auxiliary_cushion=w['auxiliary_cushion'],old_score_ceilings_enforced=False)})
    require(slots==p['trial_slots'],'Exact four proposal slots')
    require(p['status']==('qualified_declared_fractions' if any(s['eligible'] for s in slots) else 'no_qualified_declared_fraction'),'Proposal qualification status')
    return slots


def actual_choice(q_current,trials):
    best=q_current;chosen=None;epsilon=1e-10*max(1.,abs(q_current))
    for i,trial in enumerate(trials):
        if trial is not None and trial['all25'] and trial['q']<best-epsilon:best=trial['q'];chosen=i
    return chosen,best


def run(c,out):
    r=read(c['result']);successful=r['status']=='bounded_physical_gate_component_fit_complete'
    require(r['status'] in ('bounded_physical_gate_component_fit_complete','failed_physical_gate_component_fit'),'Terminal campaign only')
    if successful:require(r['error'] is None and r['candidate_frozen'] is True,'Successful frozen candidate')
    cfg=read(r['configuration']);a=module(cfg['adapter'],'review_physical_adapter');base,parent,inv,seeds=a.context(cfg)
    require(cfg['implementation']==r['implementation'] and cfg['limits']==r['limits']==a.LIMITS,'Frozen implementation/limits')
    require(r['configuration_fingerprint']==a.fingerprint(cfg),'Actual configuration fingerprint')
    release=read(cfg['root_release']);source=read(cfg['source_review']);prep=read(cfg['preparation']);execution=read(cfg['execution_context'])
    require(release['status']=='released_bounded_physical_gate_component_fit' and source['status']=='passed_physical_gate_component_fit_source_review' and prep['status']=='prepared_physical_gate_component_fit_zero_native','Release/source/preparation statuses')
    for q in (release,source,prep):require(q['configuration_fingerprint']==r['configuration_fingerprint'],'Exact preparation/review/release fingerprint')
    for k in ('implementation','adapter','core','declaration','execution_context','parent_configuration','reuse_inventory','reuse_implementation','journal_implementation','persistence_policy','solver_environment','saved_initial_proposal','source_review','preparation','limits','output_directory'):
        require(release[k]==cfg[k],'Exact release '+k)
    for k in ('implementation','adapter','core'):require(source[k]==cfg[k],'Reviewed source '+k)
    require(execution['execution']['mode']=='login' and r['job_id'] is None and release['cache_directory']==cfg['output_directory']+'/cache','Truthful executor/cache')
    require(Path(c['result']['path']).parent==Path(cfg['output_directory']),'Exact declared terminal path')
    require(r['reserved_outcomes_opened'] is False and r['physical_accuracy_qualified'] is False and r['automatic_extension'] is False,'No widened scope')
    journal=module(cfg['journal_implementation'],'review_physical_journal')
    ledger,js,intents,returns,jcellpins=journal_contract(r,cfg,journal,out,successful)
    if not successful:
        require(r['error'] and r['stopping_reason']=='failure_no_retry' and r['candidate_frozen'] is False,'Failure retained, no candidate admission')
        evidence=[]
        for n,data in returns.items():
            for rr in data['receipt_delta']:
                if rr.get('output'):verify(rr['output'])
            evidence.append({'number':n,'attempt':data['attempt'],'receipts':data['receipt_delta']})
        return {'status':'reviewed_failed_physical_gate_component_fit','result':c['result'],'configuration':r['configuration'],'implementation':r['implementation'],
                'reviewer_implementation':pin(__file__),'selected':None,'selected_native_admission_passed':False,'campaign_error':r['error'],
                'journal_review':pin(out/'JOURNAL_REVIEW.json'),'terminal_attempts':evidence,'iterations':r['iterations'],'actual_cost':r['cost'],
                'scope':'Failure and invocation provenance only; no selected candidate or physical admission. Partial outputs remain immutable.',
                'new_native_calls':0,'new_optimizer_calls':0,'physical_accuracy_qualified':False}
    for p in parent['executables'].values():verify(p)
    build=read(parent['build_receipt']);qualification=read(parent['qualification_review']);verify(qualification['implementation']);verify(qualification['result'])
    require(build['status']=='complete' and build['manifest']==parent['source_manifest'],'Admitted native build')
    require(qualification['status']=='passed_local_response_native_qualification' and all(qualification[k]==parent[k] for k in ('source_manifest','build_receipt','executables')),'Same native qualification')
    require(all(qualification[k] is True for k in ('zero_baseline_parity_passed','force_checks_passed','moment_response_checks_passed','export_checks_passed')),'Original numerical qualification flags')
    h={k:module(parent[k],'review_physical_'+k) for k in base.HELPERS};h['native_helper'].invoke=forbidden
    # All file checks are memoized by exact pins; never repeat old raw admission.
    base.verify=verify;base.pin=pin
    m=read(parent['manifest']);domain=read(parent['domain']);data=read(parent['development_data']);core=module(cfg['core'],'review_physical_core')
    core.clarabel.DefaultSolver=forbidden
    core.propose=forbidden
    ctx=SimpleNamespace(h=h,m=m,config=parent,data=data,targets=read(parent['targets']),target_inventory=read(parent['target_inventory']))
    rendered={};identity_cache={}
    def render(values,cid,exports=None):
        row=m['rows'][cid];enabled=row['needs_export'] if exports is None else exports
        key=(row['source_input']['key']['sha256'],row['source_input']['parameters']['sha256'],row['metal'],tuple(base.effective(values,row['metal'],domain)),enabled)
        if key not in rendered:
            for p in row['source_input'].values():verify(p)
            rendered[key]=base.render(values,row,domain,h['resolver'],enabled)
        return rendered[key]
    ctx.render=render
    def identity(values,cid):
        row=m['rows'][cid];memo=(cid,tuple(base.effective(values,row['metal'],domain)))
        if memo not in identity_cache:
            text,meta=render(values,cid)
            q={'cell_id':cid,'effective_values':meta['physical_values_a_b_c_d_R_epsilon'],'xyz':row['source_input']['xyz'],'parameters':row['source_input']['parameters'],'key_sha256':hashlib.sha256(text.encode()).hexdigest(),'implementation':parent['implementation'],'source_manifest':parent['source_manifest'],'executables':parent['executables'],'helpers':{k:parent[k] for k in base.HELPERS}}
            identity_cache[memo]=(token(q),q,meta)
        return identity_cache[memo]
    fixed=core.make_fixed(seeds[inv['origin_key']],m,read(cfg['decomposition']),domain,read(cfg['declaration']))
    def paired(vector):
        energies=dict(fixed['static']);energies.update({k:v['energy'] for k,v in vector['cells'].items()});parents={}
        for name,channels in fixed['terms'].items():parents[name]={k:math.fsum(c*energies[cid] for cid,c in terms.items()) for k,terms in channels.items()}
        errors={};ratios={};squares={}
        for component in ('P','B'):
            errors[component]=[parents[f'aquo{n}_central_Dy'][component]-parents[f'aquo{n}_central_La'][component]-fixed['references']['W'+str(n)][component] for n in (8,9)]
            squares[component]=math.fsum(x*x for x in errors[component])/2;ratios[component]=squares[component]/fixed['denominators']['J_'+component]
        return {'parents':parents,'paired_errors_kcal_mol':errors,'J_P':squares['P'],'J_B':squares['B'],'q':max(ratios.values()),'component_RMS':{k:math.sqrt(v) for k,v in squares.items()}}
    registry=inv['registry'];loaded=dict(seeds);authority_cells={};metal_index={};loaded_authority_keys=set()
    def index_authority(key,value):
        if key in loaded_authority_keys:return
        require(token(value['parameters'])==key and value['status']=='complete' and set(value['cells'])==set(m['rows']),'Admitted historical vector structure')
        for cid,cell in value['cells'].items():authority_cells.setdefault(cell['cache_key'],{})[a.cell_digest(cell)]={'vector_key':key,'authority':registry[key]}
        loaded_authority_keys.add(key)
    for key,v in seeds.items():index_authority(key,v)
    for key,entry in sorted(registry.items()):
        for metal in ('La','Dy'):metal_index.setdefault((metal,tuple(base.effective(entry['parameters'],metal,domain))),[]).append(key)
    def historical(key):
        if key not in loaded:loaded[key]=a.R.load(registry[key]);index_authority(key,loaded[key])
        return loaded[key]
    def inherited_authority(record,values,cid):
        ck=record['cache_key'];digest=a.cell_digest(record)
        if digest in authority_cells.get(ck,{}):return authority_cells[ck][digest]
        candidates=[]
        first=record.get('first_requested_vector')
        if first is not None and token(first) in registry:candidates.append(token(first))
        candidates.extend(metal_index.get((m['rows'][cid]['metal'],tuple(base.effective(values,m['rows'][cid]['metal'],domain))),[]))
        for key in dict.fromkeys(candidates):
            historical(key)
            if digest in authority_cells.get(ck,{}):return authority_cells[ck][digest]
        raise ValueError('No exact admitted historical cell authority '+cid+' '+ck)
    vectors={};newcells={};inherited={};baseline=set()
    for key,rec in r['vectors'].items():
        if rec['status']!='complete':require(rec['status']=='budget_not_attempted','Only explicit budget-not-attempted terminal vector');continue
        v=read(rec['artifact']);require(v['status']=='complete' and token(v['parameters'])==key and v['parameters']==rec['parameters'],'Exact complete vector identity');vectors[key]=v
        require(rec['inherited']==(key in registry),'Exact inherited/new vector classification')
        if key in registry:require(v==historical(key),'Immutable full-vector import including original score')
        else:
            assembled=base.Evaluator._assemble(ctx,v['cells']);require(canon(assembled)=={k:v[k] for k in ('model','response','score','cells')},'Stored model/response/objective exact assembly')
        core.checked(v['parameters'],domain);metrics=core.assess(v,fixed);component=paired(v)
        require(component['q']==metrics['q'] and component['J_P']==metrics['J_P'] and component['J_B']==metrics['J_B'],'Independent literal paired objective')
        for cid,cell in v['cells'].items():
            ck,ident,meta=identity(v['parameters'],cid);request=cell['request_provenance']
            require(cell['status']=='complete' and cell['cache_key']==ck and cell['cache_identity']==ident and request['cache_key']==ck and request['cache_identity']==ident and request['render']==meta and request['vector']==v['parameters'],'Exact scientific/request identity')
            canonical={k:x for k,x in cell.items() if k!='request_provenance'}
            if cell.get('root_release')==cfg['root_release'] and not cell['inherited_native']:
                if ck in newcells:require(canonical=={k:x for k,x in newcells[ck].items() if k!='request_provenance'},'Same new native cell reused exactly')
                else:newcells[ck]=cell
            elif cell['inherited_native']:
                require(base.effective(v['parameters'],m['rows'][cid]['metal'],domain)==base.effective(m['initial_vector'],m['rows'][cid]['metal'],domain),'Original baseline effective parameters')
                expected=copy.deepcopy(m['baseline_cells'][cid]);expected.update(cache_identity=ident,request_render=meta,cache_key=ck)
                require(a.cell_digest(cell)==a.cell_digest(expected),'Original baseline exact cell authority');baseline.add(ck)
            else:
                variant=(ck,a.cell_digest(cell))
                if variant not in inherited:inherited[variant]=inherited_authority(cell,v['parameters'],cid)
    require(set(seeds)<=set(vectors) and len(seeds)==27 and r['origin_key']==inv['origin_key'],'Complete27 imports and exact start')
    require(r['import_counts']['origin']==1 and r['import_counts']['initial_derivatives']==24 and r['import_counts']['additional_pilot']==2 and r['import_counts']['available_saved_vectors']==len(registry),'Separate import counts')
    current=vectors[inv['origin_key']];require(r['initial_metrics']==core.assess(current,fixed),'Exact initial metrics')
    iterations=[];real_calls=imported_calls=0;requested=[];trial_gate_reports=[]
    for step,it in enumerate(r['iterations'],1):
        require(step<=8 and it['outer_step']==step and it['origin_key']==token(current['parameters']),'Contiguous incumbent origins')
        points=core.probe_points(current['parameters'],domain);expectedkeys=[token(x) for x in points]
        if step==1:require(expectedkeys==inv['derivative_keys'],'Exact24 initial imported stencil')
        require(it['derivative_keys']==expectedkeys[:len(it['derivative_keys'])],'Ordered exact derivative prefix')
        requested.extend((step,'probe',k) for k in it['derivative_keys'])
        if 'proposal' not in it:
            require(it['status']=='hard_budget_exhausted' and step==len(r['iterations']),'Only terminal hard-budget prefix may omit proposal');break
        require(it['derivative_keys']==expectedkeys,'All24 derivatives before proposal')
        p=read(it['proposal']);problem=core.linearize(current,[vectors[k] for k in expectedkeys],fixed)
        slots=copy.deepcopy(check_proposal(p,problem,fixed,core,cfg,step,expectedkeys));real_calls+=p['solver_calls'];imported_calls+=int('imported_saved_proposal' in p)
        require(len(it['trial_slots'])==len(slots),'Runtime complete fraction denominator')
        for stored,slot in zip(it['trial_slots'],slots):require({k:stored[k] for k in slot if k!='certificate'}=={k:v for k,v in slot.items() if k!='certificate'},'Runtime exact fraction/key/eligibility')
        if it['status']=='no_eligible_fraction':
            require(not any(s['eligible'] for s in slots) and 'selection' not in it and step==len(r['iterations']),'No eligible fraction stop');iterations.append({'step':step,'status':it['status'],'proposal_status':p['status'],'solver_status':p.get('solver_status'),'solver_success':p.get('solver_success')});break
        require(len(slots)==4 and [s['fraction'] for s in slots]==[1.,.5,.25,.125],'Four fixed fractions')
        for i,slot in enumerate(slots):
            runtime=it['trial_slots'][i]
            if not slot['eligible']:require('artifact' not in runtime and 'execution_status' not in runtime,'No ineligible native trial');continue
            if runtime.get('execution_status')=='complete':
                require(runtime['artifact']==r['vectors'][slot['key']]['artifact'] and runtime['native_admitted'] is True,'Explicit eligible completion')
                require(runtime['exact_saved_vector']==(slot['key'] in registry),'Exact full-vector reuse declaration')
                if runtime['exact_saved_vector']:require(runtime['required_new_programs']==0,'Imported full vector never rerun')
                else:require(runtime['required_new_programs'] in (0,28,56),'Finite one-/two-metal native count declaration')
                requested.append((step,'trial',slot['key']));slot.update(execution_status='complete',native_admitted=True,vector=vectors[slot['key']])
        if it['status']=='hard_budget_exhausted':require('selection' not in it and step==len(r['iterations']),'No partial budget-batch selection');break
        selection=core.select(current,slots,fixed);require(selection==it['selection'],'Complete-batch exact trial/current selection')
        q0=paired(current)['q'];choices=[None]*4
        for i,s in enumerate(slots):
            if not s['eligible']:continue
            metrics=core.assess(s['vector'],fixed);q=paired(s['vector'])['q']
            require(metrics['physical_constraints_pass']==metrics['all25_original_physical_gates'] and metrics['old_score_ceilings_enforced'] is False,'Selection only original25 physical gates')
            trial_gate_reports.append({'step':step,'fraction':s['fraction'],'key':s['key'],'q':q,'all25':metrics['all25_original_physical_gates'],'old_ceiling_flags_descriptive':metrics['fixed_ceiling_pass'],'electrical_RMS':metrics['electrical_residual_RMS']})
            choices[i]={'all25':metrics['all25_original_physical_gates'],'q':q}
        chosen,best=actual_choice(q0,choices)
        require(chosen==selection['selected_trial_index'],'Independent actual q/current-epsilon tie selection')
        if chosen is None:require(it['status']=='no_actual_improvement' and step==len(r['iterations']),'No improvement stops campaign');iterations.append({'step':step,'status':it['status'],'solver_status':p.get('solver_status'),'solver_success':p.get('solver_success')});break
        require(it['status']=='complete' and it['selected_key']==slots[chosen]['key'],'Selected eligible trial identity')
        current=slots[chosen]['vector'];iterations.append({'step':step,'origin':it['origin_key'],'selected':it['selected_key'],'q_before':q0,'q_after':best,'solver_status':p.get('solver_status'),'solver_success':p.get('solver_success'),'solver_calls':p['solver_calls']})
    require(real_calls==r['real_proposal_calls']<=7 and imported_calls==r['imported_proposal_calls']<=1 and real_calls+imported_calls==r['logical_proposal_calls']<=8,'Imported versus new solver accounting')
    selectedkey=token(current['parameters']);require(selectedkey==r['selected_key']==r['current_key'] and r['selected']==r['vectors'][selectedkey]['artifact'],'Selected lineage')
    require(read(r['selected'])==current and core.assess(current,fixed)==r['selected_metrics'] and r['selected_metrics']['physical_constraints_pass'],'Final exact selected metrics and25')
    require(r['completed_outer_steps']==sum(it['status']=='complete' for it in r['iterations']),'Completed-step count')
    stop=r['stopping_reason']
    if stop=='outer_step_limit':require(len(r['iterations'])==8 and all(it['status']=='complete' for it in r['iterations']),'Exact outer cap')
    elif stop=='development_target_met':require(r['selected_metrics']['full_development_target'],'Declared development stop')
    else:require(stop in ('no_eligible_fraction','no_actual_admitted_improvement','hard_budget_exhausted'),'Declared terminal stop')
    for step,role,key in requested:require({'outer_step':step,'role':role} in r['vectors'][key]['requested_roles'],'Explicit probe/trial role provenance')
    expected_roles={}
    for step,role,key in requested:expected_roles.setdefault(key,[]).append((step,role))
    for key in vectors:require(Counter((x['outer_step'],x['role']) for x in r['vectors'][key]['requested_roles'])==Counter(expected_roles.get(key,[])),'No omitted or extraneous completed-vector roles')
    PROGRESS['selection_phases_complete']=True
    require(set(ledger['vector_keys'])==set(vectors)-set(registry),'New vectors versus exact admitted full-vector imports')
    require(set(ledger['vector_keys'])<={key for step,role,key in requested},'Every new vector has declared probe/trial role')
    requests=ledger['requests'];require(len(requests)==len(ledger['vector_keys']) and [q['vector_key'] for q in requests]==ledger['vector_keys'],'Unique ordered new requests')
    for q in requests:
        v=vectors[q['vector_key']];require(q['status']=='complete' and q['values']==v['parameters'] and set(q['cells'])==set(m['rows']),'Complete46 new request')
        role=q['role'];require(role['parameters']==v['parameters'] and (role['outer_step'],role['role'],q['vector_key']) in requested,'Request role identity')
        require(q['objective']==v['score']['objective'],'Request saved objective')
        require(all(x['status']=='complete' and x['cache_key']==v['cells'][cid]['cache_key'] for cid,x in q['cells'].items()),'Complete request/cell joins')
    cache=Path(release['cache_directory']);files={p.parent.name:p for p in (cache/'cells').glob('*/RESULT.json')}
    require(set(files)==set(newcells),'Complete unique new-cache denominator')
    receipts=[];events=[];native_rows={};attempts=ledger.get('native_attempts',[]);attempt_index={(x['cache_key'],x['program'],x['folder']):x for x in attempts}
    require(len(attempt_index)==len(attempts),'Unique native intent identity')
    request_attempt_counts=[0]*len(requests)
    for ck,record in newcells.items():
        resultpin=pin(files[ck]);disk=read(resultpin);require({k:v for k,v in record.items() if k!='request_provenance'}==disk,'Exact new-cache record')
        require((resultpin['path'],resultpin['sha256']) in jcellpins,'New cell committed evidence pin')
        cid=record['cell_id'];row=m['rows'][cid]
        require(record['implementation']==parent['implementation'] and record['root_release']==cfg['root_release'] and record['job_id'] is None,'New producer/root/null job identity')
        require(record['executor']['implementation']==cfg['adapter'] and record['executor']['driver']==cfg['implementation'] and record['executor']['mode']=='login' and record['executor']['job_id'] is None,'Truthful native login identity')
        base.Evaluator._saved_inputs(ctx,record,row,record['first_requested_vector'],files[ck].parent)
        require(len(record['receipts'])==1+int(row['needs_gradient']),'Unique cell program denominator')
        for program,receipt in zip(('analyze','testgrad'),record['receipts']):
            attempt=attempt_index[(ck,program,receipt['cwd'])];number=attempt['number'];intent=intents[number];returned=returns[number]
            require(attempt['status']=='returned_success' and attempt['returncode_observed'],'Successful unique native attempt')
            require(intent['program_inputs']==returned['program_inputs']==record['program_inputs'][program],'Durable prelaunch input pin identity')
            require(returned['receipt_delta']==[receipt] and returned['attempt']==attempt,'Actual terminal receipt equals journal return')
            ri=intent['request_index'];require(requests[ri]['values']==record['first_requested_vector'] and requests[ri]['cells'][cid]['cache_key']==ck,'Native invocation request lineage');request_attempt_counts[ri]+=1
            require(receipt['executor']==record['executor'] and receipt['timeout_seconds']==300 and receipt['returncode']==0 and receipt['elapsed_seconds']>=0,'Actual unchanged invocation policy')
            start=datetime.datetime.fromisoformat(receipt['started_utc']);end=datetime.datetime.fromisoformat(receipt['ended_utc']);require(end>=start,'Receipt chronological order');events.extend([(start,1),(end,-1)]);receipts.append(receipt)
        native_rows[ck]={'cell_id':cid,'result':resultpin,'programs':len(record['receipts'])}
    require(len(receipts)==len(attempts)==ledger['native_program_invocations']==ledger['native_returncodes_observed'],'Receipts/attempts/observed counts agree')
    for q,n in zip(requests,request_attempt_counts):require(q['required_new_programs']==n,'Whole-vector predicted versus actual new-program count')
    overlap=maximum=0
    for _,delta in sorted(events,key=lambda x:(x[0],x[1])):overlap+=delta;maximum=max(maximum,overlap)
    require(overlap==0 and maximum<=12,'Actual maximum native concurrency')
    PROGRESS['native_cache_receipts_complete']=True
    replay=module(c['derived_replay_helper'],'review_derived_arithmetic');rawrows={}
    for cid,record in current['cells'].items():
        row=m['rows'][cid];ck=record['cache_key']
        if ck not in newcells:
            authority={'original_baseline_review':parent['baseline_review']} if record['inherited_native'] else inherited_authority(record,current['parameters'],cid)
            rawrows[cid]={'status':'inherited_exact_admitted_native_authority','cache_key':ck,'authority':authority};continue
        folder=Path(record['program_inputs']['analyze']['xyz']['path']).parent.parent
        # Inputs/receipt hashes already checked once per new cache cell above.
        parsed=canon(base.Evaluator._parse(ctx,record,row));derived={'export_qualification','chi_a0_cubed','chi_solve_residual'}
        require({k:record[k] for k in parsed if k not in derived}=={k:v for k,v in parsed.items() if k not in derived},'Exact selected raw energy/moment/u-v/analytic-force replay '+cid)
        q=None;chi=[]
        if 'export_qualification' in parsed:
            q=replay.compare_qualification(record['export_qualification'],parsed['export_qualification'],row,h['scout_helper'])
            for report in (record['export_qualification'],parsed['export_qualification']):
                for xp in report['exports']:
                    lam=xp['spectrum']['A'];scaled=xp['spectrum']['normalized_A'];require(0<lam[0] and scaled[0]>1e-6*scaled[-1]>0,'Original stored/replayed SPD gates')
            dim=3*row['atoms'];matrix=np.asarray(parsed['export_qualification']['exports'][1]['A']);z=np.tile([0.,0.,1.],row['atoms']);x=np.linalg.solve(matrix,z)
            bound=h['scout_helper'].arithmetic(dim,matrix,x)*max(1.,float(np.max(abs(matrix))),float(np.max(abs(x))))
            require(0<=record['chi_solve_residual']<=bound and 0<=parsed['chi_solve_residual']<=bound,'Both conditioned chi solve gates')
            for k in ('chi_a0_cubed','chi_solve_residual'):replay.derived(record[k],parsed[k],k,dim,h['scout_helper'].arithmetic,chi)
        rawrows[cid]={'status':'passed_selected_raw_native_replay','cache_key':ck,'input':record['program_inputs'],'outputs':[v['output'] for v in record['receipts']],'derived_export_replay':q,'derived_chi_replay':chi}
    require(len(rawrows)==46 and sum(row['needs_gradient'] for row in m['rows'].values())==10,'All selected46/ten analytic-force rows covered')
    gates=independent25(current,data['rows']);component=paired(current)
    electrical={k:math.sqrt(current['score']['block_objectives'][v]) for k,v in (('pair_Debye','pair'),('C_Debye','C'),('chi_a0_cubed','chi'))}
    require(electrical==r['selected_metrics']['electrical_residual_RMS'],'Absolute-unit electrical RMS')
    oldblocks={k:current['score']['block_objectives'][v] for k,v in core.SCORE_KEYS.items()}
    descriptive={k:{'selected':v,'origin24':fixed['ceilings'][k],'difference':v-fixed['ceilings'][k],'nonincrease':v<=fixed['ceilings'][k]} for k,v in oldblocks.items()}
    require({k:v['nonincrease'] for k,v in descriptive.items()}==r['selected_metrics']['fixed_ceiling_pass'] and r['selected_metrics']['electrical_physical_acceptance_established'] is False,'Old ceilings descriptive, no electrical acceptance')
    save(out/'SELECTED_RAW_REPLAY.json',{'status':'passed_selected_native_admission','selected':r['selected'],'rows':rawrows,'original25':gates,'literal_components':component,'old_four_blocks_descriptive':descriptive,'electrical_RMS':electrical,'new_native_calls':0})
    require(r['cost']['scheduler_allocated_CPU_seconds']==0 and r['cost']['execution_mode']=='login' and r['cost']['scheduler_allocation'] is None,'Actual login cost, no fabricated Slurm allocation')
    require(r['cost']['numerical_threads']==1 and r['cost']['native_workers']==12,'Declared resource envelope')
    for k in ('wall_seconds','controller_CPU_seconds','child_user_CPU_seconds','child_system_CPU_seconds'):require(math.isfinite(r['cost'][k]) and r['cost'][k]>=0,'Finite measured cost')
    return {'status':'passed_physical_gate_component_fit_actual_review','result':c['result'],'configuration':r['configuration'],'implementation':r['implementation'],
            'driver_implementation':r['implementation'],'reviewer_implementation':pin(__file__),'selected':r['selected'],'selected_native_admission_passed':True,'selected_key':selectedkey,
            'selected_metrics':r['selected_metrics'],'literal_components':component,'old_four_blocks_descriptive':descriptive,'electrical_RMS':electrical,
            'root_release':cfg['root_release'],'source_review':cfg['source_review'],'immutable_reuse_inventory':cfg['reuse_inventory'],'journal_review':pin(out/'JOURNAL_REVIEW.json'),
            'counts':{'all_complete_vectors':len(vectors),'new_vectors':r['new_parameter_vectors'],'initial_imported_vectors':27,'used_complete_imported_vectors':len(set(vectors)&set(registry)),
                      'historical_vector_bodies_read_once':len(loaded),'inherited_cell_variants':len(inherited),'original_baseline_cells':len(baseline),'new_native_cells':len(newcells),
                      'native_programs':len(receipts),'maximum_native_overlap':maximum,'outer_steps':len(r['iterations']),'completed_outer_steps':r['completed_outer_steps'],
                      'logical_proposals':r['logical_proposal_calls'],'new_solver_calls':real_calls,'imported_proposals':imported_calls,'selected_cells':46,'selected_EF_cells':10},
            'selected_raw_review':pin(out/'SELECTED_RAW_REPLAY.json'),'native_cache':native_rows,'trajectory':iterations,'executed_trial_gates':trial_gate_reports,
            'stopping_reason':stop,'actual_cost':r['cost'],'new_native_calls':0,'new_optimizer_calls':0,'new_reference_calls':0,'physical_accuracy_qualified':False,'reserved_outcomes_opened':False,
            'limits':['Consumed development only; no global capacity, electrical accuracy or physical model qualification.',
                      'Imported historic authorities are reused; their old vectors, ledger and raw native outputs are not re-admitted. Each needed historical body is loaded at most once.',
                      'Journal and new-cache receipts are audited once; only final-selected new cells receive raw scientific parser replay.',
                      'This reviewer authored the driver/adapter and evaluator; numerical CORE has separate source review. Explicit reconstructed arithmetic is not an independent QM oracle.',
                      'Raw inputs, native energy/moment/force arrays and stored scores remain exact. Only derived LAPACK export/chi replay uses the previously qualified arithmetic bounds.',
                      'All solver statuses remain literal. A strict fraction witness does not relabel AlmostSolved or establish nonlinear/global convergence.',
                      'No consumed-trimer result or fresh-validation output is opened. Follow-up requires its separately bound release.']}


def main():
    p=argparse.ArgumentParser();p.add_argument('--result',type=Path,required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--derived-replay-helper',type=Path,required=True);args=p.parse_args()
    require(all(os.environ.get(k)=='1' for k in THREADS),'Exactly one numerical thread')
    if not args.result.exists():print(json.dumps({'status':'pending_terminal_result'}));return 0
    terminal=json.loads(args.result.read_text())
    if terminal.get('status')=='started':print(json.dumps({'status':'pending_terminal_result'}));return 0
    args.out.mkdir(parents=True,exist_ok=False);tick,cpu=time.monotonic(),time.process_time()
    c={'reviewer_implementation':pin(__file__),'result':pin(args.result),'derived_replay_helper':pin(args.derived_replay_helper)};save(args.out/'CONFIG_RUN.json',c)
    try:r=run(c,args.out)
    except Exception as exc:
        r={'status':'unavailable_physical_gate_component_fit_actual_review','error':repr(exc),'result':c['result'],'configuration':terminal.get('configuration'),
           'implementation':terminal.get('implementation'),'reviewer_implementation':pin(__file__),'selected':terminal.get('selected'),'selected_native_admission_passed':False,
           'completed_audit_phases':PROGRESS,'new_native_calls':0,'new_optimizer_calls':0}
    r.update(review_configuration=pin(args.out/'CONFIG_RUN.json'),checks=CHECKS,review_cost={'wall_seconds':time.monotonic()-tick,'CPU_seconds':time.process_time()-cpu,'threads':1})
    save(args.out/'REVIEW.json',r)
    explanation=r.get('error') or (r.get('scope') if r['status']=='reviewed_failed_physical_gate_component_fit' else 'Terminal journal, explicit fraction selection, unique new receipts and final-selected native/25-gate/component arithmetic checked within stated scope. Old score ceilings descriptive; no physical accuracy claim.')
    (args.out/'REPORT.txt').write_text(r['status']+'\n'+explanation+'\n')
    print(json.dumps({'review':pin(args.out/'REVIEW.json'),'status':r['status'],'error':r.get('error')}))
    return int(r['status']!='passed_physical_gate_component_fit_actual_review')

if __name__=='__main__':raise SystemExit(main())
