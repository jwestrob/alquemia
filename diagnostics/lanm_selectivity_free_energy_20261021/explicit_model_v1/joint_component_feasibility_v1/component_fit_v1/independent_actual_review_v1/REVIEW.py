#!/usr/bin/env python3
"""Terminal-only bounded saved audit; no optimizer or native invocations."""
import argparse
import copy
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
THREADS=('OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS')
def require(ok,msg):
    global CHECKS
    CHECKS+=1
    if not ok:raise ValueError(msg)
def pin(p):
    p=Path(p).resolve();return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def verify(p):
    k=(p['path'],p['sha256'])
    if k not in VERIFIED:require(pin(p['path'])==p,'Changed pin '+p['path']);VERIFIED[k]=True
    return Path(p['path'])
def read(p):return json.loads(verify(p).read_text())
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


def run(c,out):
    r=read(c['result']);require(r['status']=='bounded_component_fit_complete' and r['error'] is None and r['candidate_frozen'],'Successful terminal frozen candidate only')
    cfg=read(r['configuration']);a=module(cfg['adapter'],'review_component_adapter');base,parent,inv,seeds=a.context(cfg)
    require(cfg['implementation']==r['implementation'] and cfg['limits']==r['limits']==a.LIMITS,'Frozen implementation/limits')
    require(r['configuration_fingerprint']==a.fingerprint(cfg),'Actual scientific fingerprint')
    release=read(cfg['root_release']);source=read(cfg['source_review']);prep=read(cfg['preparation']);execution=read(cfg['execution_context'])
    require(release['status']=='released_bounded_component_fit' and source['status']=='passed_component_fit_source_review' and prep['status']=='prepared_component_fit_zero_native','Actual release/source/preparation statuses')
    for q in (release,source,prep):require(q['configuration_fingerprint']==r['configuration_fingerprint'],'Exact source/preparation/release fingerprint')
    for k in ('implementation','adapter','core','declaration','execution_context','parent_configuration','reuse_inventory','source_review','preparation','limits','output_directory'):require(release[k]==cfg[k],'Exact release '+k)
    for k in ('implementation','adapter','core'):require(source[k]==cfg[k],'Independent reviewed source '+k)
    require(execution['execution']['mode']=='login' and r['job_id'] is None and release['cache_directory']==cfg['output_directory']+'/cache','Truthful execution/cache')
    require(Path(c['result']['path']).parent==Path(cfg['output_directory']),'Exact declared terminal output')
    for p in parent['executables'].values():verify(p)
    build=read(parent['build_receipt']);qualification=read(parent['qualification_review']);verify(qualification['implementation']);verify(qualification['result'])
    require(build['status']=='complete' and build['manifest']==parent['source_manifest'],'Admitted native build')
    require(qualification['status']=='passed_local_response_native_qualification' and all(qualification[k]==parent[k] for k in ('source_manifest','build_receipt','executables')),'Unchanged native qualification')
    require(all(qualification[k] is True for k in ('zero_baseline_parity_passed','force_checks_passed','moment_response_checks_passed','export_checks_passed')),'Original numerical admission flags')
    h={k:module(parent[k],'review_component_'+k) for k in base.HELPERS};h['native_helper'].invoke=forbidden
    m=read(parent['manifest']);domain=read(parent['domain']);data=read(parent['development_data']);core=module(cfg['core'],'review_component_core');core.minimize=forbidden
    ctx=SimpleNamespace(h=h,m=m,config=parent,data=data,targets=read(parent['targets']),target_inventory=read(parent['target_inventory']))
    rendered={}
    def render(values,cid,exports=None):
        row=m['rows'][cid];enabled=row['needs_export'] if exports is None else exports
        k=(row['source_input']['key']['sha256'],row['source_input']['parameters']['sha256'],row['metal'],tuple(base.effective(values,row['metal'],domain)),enabled)
        if k not in rendered:
            for p in row['source_input'].values():verify(p)
            rendered[k]=base.render(values,row,domain,h['resolver'],enabled)
        return rendered[k]
    ctx.render=render
    def identity(values,cid):
        row=m['rows'][cid];text,meta=render(values,cid)
        q={'cell_id':cid,'effective_values':meta['physical_values_a_b_c_d_R_epsilon'],'xyz':row['source_input']['xyz'],'parameters':row['source_input']['parameters'],'key_sha256':hashlib.sha256(text.encode()).hexdigest(),'implementation':parent['implementation'],'source_manifest':parent['source_manifest'],'executables':parent['executables'],'helpers':{k:parent[k] for k in base.HELPERS}}
        return token(q),q,meta
    vectors={};newcells={};inherited=set();decomposition=read(cfg['decomposition']);fixed=core.make_fixed(seeds[inv['origin_key']],m,decomposition,domain,read(cfg['declaration']))
    def independent_q(vector):
        energies={k:v['energy_kcal_mol'] for k,v in m['static'].items()};energies.update({k:v['energy'] for k,v in vector['cells'].items()})
        errors={'P':[],'B':[]}
        for n in (8,9):
            for component in ('P','B'):
                totals={metal:math.fsum(coefficient*energies[cid] for cid,coefficient in decomposition['parents'][f'aquo{n}_central_{metal}']['model_terms'][component].items()) for metal in ('La','Dy')}
                errors[component].append(totals['Dy']-totals['La']-decomposition['paired_Dy_minus_La']['W'+str(n)]['reference_total_kcal_mol'][component])
        return max(math.fsum(x*x for x in errors[k])/2/fixed['denominators']['J_'+k] for k in ('P','B'))
    for key,rec in r['vectors'].items():
        if rec['status']!='complete':require(rec['status']=='budget_not_attempted','Only an explicit unattempted budget slot may be incomplete');continue
        v=read(rec['artifact']);require(v['status']=='complete' and token(v['parameters'])==key and v['parameters']==rec['parameters'],'Exact complete vector identity');vectors[key]=v
        if key in seeds:require(v==seeds[key] and rec['inherited'] is True,'Exact immutable inherited vector')
        else:
            require(rec['inherited'] is False,'New vector classification')
            assembled=base.Evaluator._assemble(ctx,v['cells']);require(canon(assembled)=={k:v[k] for k in ('model','response','score','cells')},'Full stored objective/response assembly')
        for cid,cell in v['cells'].items():
            ck,ident,meta=identity(v['parameters'],cid);request=cell['request_provenance']
            require(cell['cache_key']==ck and cell['cache_identity']==ident and request['cache_key']==ck and request['cache_identity']==ident and request['render']==meta and request['vector']==v['parameters'],'Exact vector/cell scientific identity')
            if ck in inv['cells']:
                require(a.cell_digest(cell)==inv['cells'][ck]['record_sha256'],'Exact inherited cell bytes');inherited.add(ck)
            elif cell['inherited_native']:
                require(base.effective(v['parameters'],m['rows'][cid]['metal'],domain)==base.effective(m['initial_vector'],m['rows'][cid]['metal'],domain),'Exact original baseline reuse')
                expected=copy.deepcopy(m['baseline_cells'][cid]);expected.update(cache_identity=ident,request_render=meta,cache_key=ck)
                require(a.cell_digest(cell)==a.cell_digest(expected),'Exact baseline authority')
            elif ck in newcells:require(a.cell_digest(cell)==a.cell_digest(newcells[ck]),'Same new-cache record reused exactly')
            else:newcells[ck]=cell
    require(set(seeds)<=set(vectors) and r['origin_key']==inv['origin_key'],'Complete26 import denominator')
    current=vectors[inv['origin_key']];iterations=[];proposal_calls=0;requested=[]
    for step,it in enumerate(r['iterations'],1):
        require(step<=8 and it['outer_step']==step and it['origin_key']==token(current['parameters']),'Contiguous actual incumbent origins')
        points=core.probe_points(current['parameters'],domain);expectedkeys=[token(x) for x in points]
        require(it['derivative_keys']==expectedkeys[:len(it['derivative_keys'])],'Exact ordered derivative prefix')
        requested.extend((step,'probe',k) for k in it['derivative_keys'])
        if 'proposal' not in it:
            require(it['status']=='hard_budget_exhausted' and step==len(r['iterations']),'Only final hard budget may omit proposal');break
        require(it['derivative_keys']==expectedkeys,'Complete24 stencil before proposal')
        probes=[vectors[k] for k in expectedkeys];p=read(it['proposal']);problem=core.linearize(current,probes,fixed);proposal_calls+=1
        require(p['solver_calls']==1 and p['current_key']==it['origin_key'] and p['derivative_keys']==expectedkeys and p['derivative_spans']==problem['spans'],'One declared proposal with exact stencil')
        require(p['fixed_ceilings']==fixed['ceilings'] and p['fixed_component_denominators']==fixed['denominators'] and p['current_metrics']==core.assess(current,fixed),'Fixed target normalization and actual origin metrics')
        require(p['Jacobian_original531']==problem['Jacobian_original531'].tolist() and p['Jacobian_components']==problem['Jacobian_components'].tolist(),'Exact saved residual/component secants')
        slots=copy.deepcopy(p['trial_slots'])
        require(len(it['trial_slots'])==len(slots),'Runtime/proposal fraction denominator')
        for stored,slot in zip(it['trial_slots'],slots):require({k:stored[k] for k in slot if k!='certificate'}=={k:v for k,v in slot.items() if k!='certificate'},'Runtime exact fraction/parameter/key/eligibility')
        if slots:
            require(len(slots)==4 and [s['fraction'] for s in slots]==[1.,.5,.25,.125],'Exact four fractions')
            solver=p['solver'];require(solver['method']=='SLSQP' and solver['options']==core.OPTIONS,'Original solver/options')
            certificates=[core.certify_fraction(problem,np.asarray(solver['returned']),f,solver['success']) for f in core.FRACTIONS]
            require(slots==certificates,'Exact scalar-certified fraction identities; no optimizer rerun')
        if it['status']=='no_eligible_fraction':
            require(not any(s['eligible'] for s in slots) and 'selection' not in it and step==len(r['iterations']),'Stop on no eligible fraction');iterations.append({'step':step,'status':it['status']});break
        for i,slot in enumerate(slots):
            if not slot['eligible']:continue
            runtime=it['trial_slots'][i]
            if runtime.get('execution_status')=='complete':
                require(runtime['artifact']==r['vectors'][slot['key']]['artifact'] and runtime['native_admitted'] is True,'Explicit eligible trial completion')
                requested.append((step,'trial',slot['key']));slot.update(execution_status='complete',native_admitted=True,vector=vectors[slot['key']])
        if it['status']=='hard_budget_exhausted':require('selection' not in it and step==len(r['iterations']),'No selection from budget-incomplete batch');break
        selection=core.select(current,slots,fixed);require(selection==it['selection'],'Exact complete-batch trial/current selection')
        # Independently enforce the declared current-based tie and actual gates.
        q0=independent_q(current);require(q0==core.assess(current,fixed)['q'],'Independent paired minimax origin');best=q0;chosen=None;epsilon=1e-10*max(1.,abs(q0))
        for i,s in enumerate(slots):
            if not s['eligible']:continue
            measured=core.assess(s['vector'],fixed);require(independent_q(s['vector'])==measured['q'],'Independent paired minimax trial')
            if measured['physical_constraints_pass'] and measured['q']<best-epsilon:best=measured['q'];chosen=i
        require(chosen==selection['selected_trial_index'],'Independent explicit tie selection')
        if chosen is None:require(it['status']=='no_actual_improvement' and step==len(r['iterations']),'No-improvement stop');break
        require(it['status']=='complete' and it['selected_key']==slots[chosen]['key'],'Actual selected trial identity');current=slots[chosen]['vector']
        iterations.append({'step':step,'origin':it['origin_key'],'selected':it['selected_key'],'q_before':q0,'q_after':best,'solver_success':p['solver']['success']})
    require(proposal_calls==r['real_proposal_calls']<=8,'Actual proposal count')
    selectedkey=token(current['parameters']);require(selectedkey==r['selected_key']==r['current_key'] and r['selected']==r['vectors'][selectedkey]['artifact'],'Final selected lineage')
    require(read(r['selected'])==current and core.assess(current,fixed)==r['selected_metrics'] and r['selected_metrics']['physical_constraints_pass'],'Selected exact data/gates/ceilings')
    require(r['completed_outer_steps']==sum(it['status']=='complete' for it in r['iterations']),'Completed-step accounting')
    stop=r['stopping_reason']
    if stop=='outer_step_limit':require(len(r['iterations'])==8 and all(it['status']=='complete' for it in r['iterations']),'Exact outer-step cap')
    elif stop=='development_target_met':require(r['selected_metrics']['full_development_target'],'Actual declared development stop')
    else:require(stop in ('no_eligible_fraction','no_actual_admitted_improvement','hard_budget_exhausted'),'Declared stopping reason')
    for step,role,key in requested:require({'outer_step':step,'role':role} in r['vectors'][key]['requested_roles'],'Explicit derivative/trial role lineage')
    ledger=read(r['evaluator']['ledger']);require(not ledger['stopped'] and not r['evaluator']['stopped'],'No hidden native failure')
    require(len(ledger['vector_keys'])==r['new_parameter_vectors']==r['evaluator']['new_parameter_vectors']<=200 and len(set(ledger['vector_keys']))==len(ledger['vector_keys']),'Exact distinct new-vector count')
    require(set(ledger['vector_keys'])==set(vectors)-set(seeds),'New vectors versus immutable imports')
    require(set(ledger['vector_keys'])<={key for step,role,key in requested},'Every new vector has declared probe/trial provenance')
    requests=ledger['requests'];require(len(requests)==len(ledger['vector_keys']) and {q['vector_key'] for q in requests}==set(ledger['vector_keys']),'One new request per vector')
    for q in requests:
        v=vectors[q['vector_key']];require(q['status']=='complete' and q['values']==v['parameters'] and set(q['cells'])==set(m['rows']),'Complete new request46')
        role=q['role'];require(role['parameters']==v['parameters'] and (role['outer_step'],role['role'],q['vector_key']) in requested,'Native request declared role')
        require(all(x['status']=='complete' and x['cache_key']==v['cells'][cid]['cache_key'] for cid,x in q['cells'].items()),'Request/cell joins')
    cache=Path(release['cache_directory']);files={p.parent.name:p for p in (cache/'cells').glob('*/RESULT.json')};require(set(files)==set(newcells),'Exact new native-cache denominator')
    receipts=[];events=[];native_rows={};attempt_index={(x['cache_key'],x['program'],x['folder']):x for x in ledger['native_attempts']}
    for ck,record in newcells.items():
        disk=read(pin(files[ck]));require(a.cell_digest(record)==token(disk),'Exact new vector/cache bytes');cid=record['cell_id'];row=m['rows'][cid]
        require(record['implementation']==parent['implementation'] and record['root_release']==cfg['root_release'] and record['job_id'] is None,'Actual new producer/root/null job')
        require(record['executor']['implementation']==cfg['adapter'] and record['executor']['driver']==cfg['implementation'] and record['executor']['mode']=='login' and record['executor']['job_id'] is None,'Actual login execution identity')
        base.Evaluator._saved_inputs(ctx,record,row,record['first_requested_vector'],files[ck].parent)
        require(len(record['receipts'])==1+int(row['needs_gradient']),'Native program denominator per unique cell')
        for program,receipt in zip(('analyze','testgrad'),record['receipts']):
            attempt=attempt_index[(ck,program,receipt['cwd'])];require(attempt['status']=='returned_success' and attempt['returncode_observed'],'Actual unique native attempt')
            require(receipt['executor']==record['executor'] and receipt['timeout_seconds']==300 and receipt['returncode']==0 and receipt['elapsed_seconds']>=0,'Unchanged native execution receipt')
            start=datetime.datetime.fromisoformat(receipt['started_utc']);end=datetime.datetime.fromisoformat(receipt['ended_utc']);require(end>=start,'Native receipt time order');events.extend([(start,1),(end,-1)]);receipts.append(receipt)
        native_rows[ck]={'cell_id':cid,'result':pin(files[ck]),'programs':len(record['receipts'])}
    attempts=ledger['native_attempts'];require(len(attempt_index)==len(attempts)==len(receipts)==ledger['native_program_invocations']==ledger['native_returncodes_observed']==r['native_programs']<=6496,'Exact native accounting')
    require([x['number'] for x in attempts]==list(range(1,len(attempts)+1)),'Contiguous invocation boundary counts')
    overlap=maximum=0
    for _,delta in sorted(events,key=lambda x:(x[0],x[1])):overlap+=delta;maximum=max(maximum,overlap)
    require(overlap==0 and maximum<=12,'Actual native concurrency bound')
    # Reuse the already reviewed strict raw/derived arithmetic distinction.
    replay=module(c['derived_replay_helper'],'review_derived_arithmetic');rawrows={}
    for cid,record in current['cells'].items():
        row=m['rows'][cid]
        if record['inherited_native']:
            rawrows[cid]={'status':'inherited_admitted_original_baseline','authority':parent['baseline_review']};continue
        folder=Path(record['program_inputs']['analyze']['xyz']['path']).parent.parent
        base.Evaluator._saved_inputs(ctx,record,row,current['parameters'],folder)
        parsed=canon(base.Evaluator._parse(ctx,record,row));derived={'export_qualification','chi_a0_cubed','chi_solve_residual'}
        require({k:record[k] for k in parsed if k not in derived}=={k:v for k,v in parsed.items() if k not in derived},'Exact selected raw energy/moment/u-v/force replay '+cid)
        q=None;chi=[]
        if 'export_qualification' in parsed:
            q=replay.compare_qualification(record['export_qualification'],parsed['export_qualification'],row,h['scout_helper'])
            for report in (record['export_qualification'],parsed['export_qualification']):
                for xp in report['exports']:
                    lam=xp['spectrum']['A'];scaled=xp['spectrum']['normalized_A'];require(0<lam[0] and scaled[0]>1e-6*scaled[-1]>0,'Original stored/replayed SPD criteria')
            dim=3*row['atoms'];matrix=np.asarray(parsed['export_qualification']['exports'][1]['A']);z=np.tile([0.,0.,1.],row['atoms']);x=np.linalg.solve(matrix,z)
            bound=h['scout_helper'].arithmetic(dim,matrix,x)*max(1.,float(np.max(abs(matrix))),float(np.max(abs(x))))
            require(0<=record['chi_solve_residual']<=bound and 0<=parsed['chi_solve_residual']<=bound,'Both original chi solve gates')
            for k in ('chi_a0_cubed','chi_solve_residual'):replay.derived(record[k],parsed[k],k,dim,h['scout_helper'].arithmetic,chi)
        rawrows[cid]={'status':'passed_selected_raw_native_replay','cache_key':record['cache_key'],'input':record['program_inputs'],'outputs':[v['output'] for v in record['receipts']],'derived_export_replay':q,'derived_chi_replay':chi}
    require(len(rawrows)==46,'All selected46 covered');gates=independent25(current,data['rows'])
    save(out/'SELECTED_RAW_REPLAY.json',{'status':'passed_selected_raw_native_replay','selected':r['selected'],'rows':rawrows,'original25':gates,'new_native_calls':0})
    require(r['cost']['scheduler_allocated_CPU_seconds']==0 and r['cost']['execution_mode']=='login' and r['cost']['scheduler_allocation'] is None,'Truthful final cost envelope')
    for k in ('wall_seconds','controller_CPU_seconds','child_user_CPU_seconds','child_system_CPU_seconds'):require(math.isfinite(r['cost'][k]) and r['cost'][k]>=0,'Finite actual cost')
    return {'status':'passed_component_fit_actual_review','result':c['result'],'configuration':r['configuration'],'implementation':r['implementation'],'driver_implementation':r['implementation'],'reviewer_implementation':pin(__file__),'selected':r['selected'],
            'selected_native_admission_passed':True,'selected_key':selectedkey,'selected_metrics':r['selected_metrics'],'root_release':cfg['root_release'],'source_review':cfg['source_review'],'immutable_reuse_inventory':cfg['reuse_inventory'],
            'counts':{'all_vectors':len(vectors),'new_vectors':r['new_parameter_vectors'],'imported_vectors':26,'unique_inherited_cells_used':len(inherited),'new_native_cells':len(newcells),'native_programs':len(receipts),'maximum_native_overlap':maximum,'outer_steps':len(r['iterations']),'completed_outer_steps':r['completed_outer_steps'],'proposal_calls':proposal_calls,'selected_cells':46,'selected_EF_cells':10},
            'selected_raw_review':pin(out/'SELECTED_RAW_REPLAY.json'),'ledger':r['evaluator']['ledger'],'native_cache':native_rows,'trajectory':iterations,'stopping_reason':stop,'actual_cost':r['cost'],
            'new_native_calls':0,'new_optimizer_calls':0,'new_reference_calls':0,'physical_accuracy_qualified':False,'reserved_outcomes_opened':False,
            'limits':['One bounded consumed-development campaign; no global model-capacity or physical-accuracy claim.','Original failed solver statuses are retained; certified fractions are finite feasibility witnesses.','This reviewer authored the driver/adapter and original scientific evaluator/objective; their frozen helpers are reused and authorship is not presented as independent oracle evidence. CORE arithmetic/selection and final25 body/energy gates are explicitly reconstructed; the separate source review remains authoritative.','Raw native inputs/tokens/energies/moments/forces stay exact. Only derived export/chi replay uses the previously qualified arithmetic allowance; authoritative stored scores remain unchanged.','No postselection regression or fresh validation outcome is read.']}


def main():
    p=argparse.ArgumentParser();p.add_argument('--result',required=True,type=Path);p.add_argument('--out',required=True,type=Path);p.add_argument('--derived-replay-helper',required=True,type=Path);a=p.parse_args()
    require(all(os.environ.get(k)=='1' for k in THREADS),'One numerical thread')
    if not a.result.exists():print(json.dumps({'status':'pending_terminal_result'}));return 0
    terminal=json.loads(a.result.read_text())
    if terminal.get('status')=='started':print(json.dumps({'status':'pending_terminal_result'}));return 0
    a.out.mkdir(parents=True,exist_ok=False);tick,cpu=time.monotonic(),time.process_time()
    c={'reviewer_implementation':pin(__file__),'result':pin(a.result),'derived_replay_helper':pin(a.derived_replay_helper)};save(a.out/'CONFIG_RUN.json',c)
    try:r=run(c,a.out)
    except Exception as exc:r={'status':'unavailable_component_fit_actual_review','error':repr(exc),'result':c['result'],'configuration':terminal.get('configuration'),'implementation':terminal.get('implementation'),'reviewer_implementation':pin(__file__),'selected':terminal.get('selected'),'selected_native_admission_passed':False,'new_native_calls':0,'new_optimizer_calls':0}
    r.update(review_configuration=pin(a.out/'CONFIG_RUN.json'),checks=CHECKS,review_cost={'wall_seconds':time.monotonic()-tick,'CPU_seconds':time.process_time()-cpu,'threads':1})
    save(a.out/'REVIEW.json',r);(a.out/'REPORT.txt').write_text(r['status']+'\n'+(r.get('error') or 'Terminal proposal/selection/cache receipts and selected46 raw admission checked. No new molecular/optimizer calls; physical qualification remains false.')+'\n')
    print(json.dumps({'review':pin(a.out/'REVIEW.json'),'status':r['status'],'error':r.get('error')}));return int(r['status']!='passed_component_fit_actual_review')
if __name__=='__main__':raise SystemExit(main())
