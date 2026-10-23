#!/usr/bin/env python3
"""Terminal-only saved audit. Never instantiate an evaluator or run a solver/engine."""
import argparse
import copy
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import time
from types import SimpleNamespace
import numpy as np

THREADS = ('OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS')
LIMITS = {'outer_steps':24,'new_vectors':672,'native_programs':24000,'workers':12}
TERMINAL = ('bounded_joint_response_fit_complete','failed_or_unavailable')
PARENTS = ('aquo8_central_La','aquo8_central_Dy','aquo9_central_La','aquo9_central_Dy')
STATES = ('La_dimer','La_trimer','Dy_dimer','Dy_trimer')
CHECKS = 0
VERIFIED = {}


def require(ok, message):
    global CHECKS
    CHECKS += 1
    if not ok: raise ValueError(message)


def token(value):
    return hashlib.sha256(json.dumps(value,sort_keys=True,separators=(',',':'),allow_nan=False).encode()).hexdigest()


def pin(path):
    p=Path(path).resolve()
    return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}


def verify(q):
    require(isinstance(q,dict) and set(q)=={'path','sha256'},'Exact file pin required')
    key=(q['path'],q['sha256'])
    if key not in VERIFIED:
        require(pin(q['path'])==q,'Changed file: '+q['path']);VERIFIED[key]=True
    return Path(q['path'])


def read(q): return json.loads(verify(q).read_text())
def fingerprint(c): return token({k:v for k,v in c.items() if k not in ('root_release','preparation_result')})


def module(q, name):
    spec=importlib.util.spec_from_file_location(name,verify(q));m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m


def close(a,b,message):
    a,b=np.asarray(a,dtype=float),np.asarray(b,dtype=float)
    require(a.shape==b.shape and np.isfinite(a).all() and np.isfinite(b).all(),message+' shape/finite')
    # Saved arithmetic roundoff only; no scientific acceptance tolerance changes.
    scale=max(1.,float(np.max(abs(a),initial=0)),float(np.max(abs(b),initial=0)))
    require(float(np.max(abs(a-b),initial=0))<=256*np.finfo(float).eps*max(1,a.size)*scale,message)


def terminal_guard(r):
    require(r.get('stage')=='run' and r.get('status') in TERMINAL and isinstance(r.get('cost'),dict),
            'Refuse nonterminal result before reading any vectors or ledger')


def choose(origin, keys, scores):
    selected=origin
    for k in keys:
        if scores[k]['objective']<scores[selected]['objective']: selected=k
    return selected


def probes(center,lower,upper):
    width=upper-lower;z=(center-lower)/width;out=[]
    for j in range(12):
        for sign in (-1,1):
            t=np.clip(z[j]+sign*.005,0,1);v=center.copy()
            v[j]=np.clip(center[j]+(t-z[j])*width[j],lower[j],upper[j]);out.append(v.tolist())
    return out


def score_audit(v,targets,data,h):
    s=v['score'];r=v['response'];old=s['old_EF']
    require(len(old['weighted_residuals'])==413 and old['available_energy_controls']==10 and old['available_force_controls']==15,'Original 413/25 denominator')
    require(old['block_weights']=={'energy':.25,'La':.25,'Dy':.25,'paired':.25},'Original EF weights')
    require(s['original_EF_checks_passed']==old['development_checks_passed'],'Original EF pass flag')
    residual=[x/2 for x in old['weighted_residuals']];parts={'pair':[],'C':[],'chi':[]}
    require(set(r['pair_dipoles_Debye'])=={x['id'] for x in targets['pair_rows']},'All 34 pair predictions')
    require(set(r['C_dipoles_Debye'])==set(PARENTS) and set(r['chi_states_a0_cubed'])==set(STATES),'All four C and chi predictions')
    for parent in PARENTS:
        rows=[x for x in targets['pair_rows'] if x['parent']==parent]
        require(len(rows)==(8 if parent.startswith('aquo8') else 9),'Literal pair denominator')
        for row in rows:
            parts['pair'].extend((x-y)/math.sqrt(12*len(rows)) for x,y in zip(r['pair_dipoles_Debye'][row['id']],row['reference_Debye']))
        parts['C'].extend((x-y)/math.sqrt(12) for x,y in zip(r['C_dipoles_Debye'][parent],targets['C_Debye'][parent]))
    parts['chi']=[(r['chi_states_a0_cubed'][k]-targets['chi_a0_cubed'][k])/2 for k in STATES]
    blocks={'old_EF':math.fsum(x*x for x in old['weighted_residuals'])}
    for k,values in parts.items(): blocks[k]=math.fsum(x*x for x in values);residual.extend(x/2 for x in values)
    require(len(residual)==531,'Joint residual denominator')
    close(residual,s['weighted_residuals'],'Independent 531 weights/order')
    for k,x in blocks.items():close(x,s['block_objectives'][k],'Block '+k)
    close(sum(blocks.values())/4,s['objective'],'Four equal objective blocks')
    close(math.fsum(x*x for x in residual),s['objective'],'Residual norm objective')
    require(s['response_physical_acceptance'] is None and s['joint_physical_acceptance'] is None and not s['production_released'],'No invented physical response gate')
    # This helper was authored by this reviewer: exact replay, not an independent oracle.
    replay=h['objective'].evaluate(data,v['model'],r,targets,h['old_objective'],h['body_helper'])
    require(replay==s,'Exact frozen objective/EF replay')


def audit(cfg):
    require(cfg['implementation']==pin(__file__),'Reviewer source differs')
    require(isinstance(cfg.get('fit_result'),dict),'Root must bind an exact terminal result pin')
    require(cfg['fit_result']['path']==cfg['expected_result_path'],'Different fit endpoint')
    r=read(cfg['fit_result']);terminal_guard(r)
    # No vector/ledger read occurs before terminal_guard above.
    dc=read(cfg['driver_config']);ec=read(dc['evaluator_config']);release=read(dc['root_release'])
    require(r['configuration']==cfg['driver_config'] and r['implementation']==dc['implementation'],'Exact driver provenance')
    require(r['scientific_fingerprint']==fingerprint(dc)==release['scientific_fingerprint'],'Driver scientific fingerprint')
    require(r['limits']==dc['limits']==release['limits']==LIMITS,'Original 24/672/24000/12 caps')
    require(release['status']=='released_bounded_joint_response_fit' and release['evaluator_config_sha256']==fingerprint(ec),'Evaluator release fingerprint')
    require(not r['validation_outcomes_opened'] and not r['automatic_extension'] and not r['physical_accuracy_qualified'],'Scope boundaries')
    for k in ('implementation','protocol','root_declaration','parameter_domain','directions','evaluator','objective_ready'):verify(dc[k])
    for q in cfg['source_reviews'].values():verify(q)
    verify(cfg['data_contract'])
    dr=read(cfg['source_reviews']['driver']);er=read(cfg['source_reviews']['evaluator'])
    require(dr['status']=='passed_bounded_joint_response_driver_review' and dr['implementation']==dc['implementation'] and dr['directions']==dc['directions'] and dr['protocol']==dc['protocol'],'Independent driver source review joins')
    require(er['status']=='passed_joint_evaluator_source_review' and er['implementation']==ec['implementation'],'Independent evaluator source review join')
    for k in ('manifest','domain','fit_protocol','source_manifest','build_receipt','executables','qualification_review'):
        require(release[k]==ec[k],'Root evaluator binding '+k)
    require(release['reserved_outcomes_authorized'] is False,'Reserved outcomes not released')
    for q in ec['executables'].values():verify(q)
    build=read(ec['build_receipt']);qual=read(ec['qualification_review'])
    require(build['status']=='complete' and build['manifest']==ec['source_manifest'] and all(build['executables'][k]==v for k,v in ec['executables'].items()),'Built executable identity')
    require(qual['status']=='passed_local_response_native_qualification','Actual numerical qualification')
    for k in ('source_manifest','build_receipt','executables'):require(qual[k]==ec[k],'Qualification '+k)
    require(all(qual[k] is True for k in ('zero_baseline_parity_passed','force_checks_passed','moment_response_checks_passed','export_checks_passed')),'Four qualification gates')
    verify(qual['implementation']);verify(qual['result']);verify(ec['source_manifest'])
    e=module(dc['evaluator'],'fit_review_evaluator');h={k:module(ec[k],'fit_review_'+k) for k in e.HELPERS}
    m=read(ec['manifest']);domain=read(ec['domain']);targets=read(ec['targets']);data=read(ec['development_data']);inventory=read(ec['target_inventory'])
    require(m['sources']=={k:ec[k] for k in e.SOURCES} and m['implementation']==ec['implementation']==dc['evaluator'],'Prepared source joins')
    for q in m['sources'].values():verify(q)
    require(h['objective'].targets_from_inventory(inventory)==targets,'Frozen reference estimator/targets')
    parameters=domain['parameters'];initial=[p['initial'] for p in parameters]
    lower=np.array([p['lower'] for p in parameters]);upper=np.array([p['upper'] for p in parameters]);width=upper-lower
    require(len(parameters)==12 and np.all(width>0) and m['initial_vector']==initial,'Exact initial/domain')
    require(len(m['rows'])==46 and sum(x['needs_gradient'] for x in m['rows'].values())==10 and sum(x['needs_export'] for x in m['rows'].values())==4 and len(m['static'])==48,'46/10/4 and 46 metal-free plus two ions')
    prep=read(dc['preparation_result']);require(prep['status']=='prepared_joint_response_fit' and prep['scientific_fingerprint']==fingerprint(dc) and prep['evaluator']['native_program_invocations']==0,'Zero-call driver preparation')
    baseline=read(prep['vectors'][prep['initial_key']]['artifact'])
    records=r['vectors'];require(len(records)<=673 and r['new_parameter_vectors']==max(0,len(records)-1),'All vectors count including failures')
    require(r.get('initial_key')==token(initial) and records[r['initial_key']]['role']=='initial','Exact single initial vector')
    vectors={};scores={};unique={};rendered={};failures=[];cell_links={}
    def render(values,cid,exports=None):
        row=m['rows'][cid];enabled=row['needs_export'] if exports is None else exports
        key=(row['source_input']['key']['sha256'],row['source_input']['parameters']['sha256'],row['metal'],tuple(e.effective(values,row['metal'],domain)),enabled)
        if key not in rendered:rendered[key]=e.render(values,row,domain,h['resolver'],enabled)
        return rendered[key]
    def identity(values,cid):
        row=m['rows'][cid];key,meta=render(values,cid)
        q={'cell_id':cid,'effective_values':meta['physical_values_a_b_c_d_R_epsilon'],'xyz':row['source_input']['xyz'],'parameters':row['source_input']['parameters'],'key_sha256':hashlib.sha256(key.encode()).hexdigest(),'implementation':ec['implementation'],'source_manifest':ec['source_manifest'],'executables':ec['executables'],'helpers':{k:ec[k] for k in e.HELPERS}}
        return token(q),q,meta
    selected_raw=None
    for key,rec in records.items():
        values=rec['parameters'];require(token(values)==key and e.checked(values,domain)==values,'Vector identity/bounds')
        require('artifact' in rec,'Terminal vector artifact missing; keep review unavailable')
        v=read(rec['artifact']);require(v['parameters']==values and v['status']==rec['status'],'Vector metadata exact')
        vectors[key]={'parameters':values,'status':v['status']}
        if v['status']!='complete':failures.append({'key':key,**rec});continue
        score_audit(v,targets,data,h);scores[key]=v['score']
        for name in ('objective','block_objectives','original_EF_checks_passed'):require(rec[name]==v['score'][name],'Vector index score '+name)
        require(set(v['cells'])==set(m['rows']) and set(v['model'])==set(data['rows']),'Full vector fixture denominator')
        values_mu={cid:c['mu_Debye'] for cid,c in m['static'].items()};links={}
        for cid,cell in v['cells'].items():
            ck,ident,meta=identity(values,cid);row=m['rows'][cid]
            require(cell['status']=='complete' and cell['cache_key']==ck and cell['cache_identity']==ident,'Complete cell identity')
            require(cell['request_provenance']=={'vector':values,'cache_key':ck,'cache_identity':ident,'render':meta},'Cell requesting vector provenance')
            clean={k:x for k,x in cell.items() if k!='request_provenance'};signature=token(clean)
            if ck in unique:require(unique[ck]['signature']==signature,'Reused cell changed')
            else:unique[ck]={'signature':signature,'cell_id':cid,'values':values,'inherited':cell['inherited_native']}
            if cell['inherited_native']:
                require(e.effective(values,row['metal'],domain)==e.effective(initial,row['metal'],domain),'Baseline reuse only unchanged metal')
                require({k:x for k,x in clean.items() if k not in ('cache_identity','cache_key','request_render')}==m['baseline_cells'][cid],'Exact inherited baseline record')
                require(cell['request_render']==meta,'Inherited requested render')
            values_mu[cid]=cell['mu_Debye'];links[cid]={'cache_key':ck,'inherited_native':cell['inherited_native']}
            if cid in data['rows']:require(v['model'][cid]==cell['model'],'EF source cell')
        pairs={x['id']:[math.fsum(c*values_mu[cid][j] for cid,c in x['model_terms'].items()) for j in range(3)] for x in inventory['groups']['pair_and_collective_dipoles']['pair_cycles']}
        C={p:[values_mu[p][j]-values_mu['ion_'+p.rsplit('_',1)[1]][j]-values_mu[p+'__shell'][j] for j in range(3)] for p in PARENTS}
        chi={state:v['cells'][alias['fixture']]['chi_a0_cubed'] for state,alias in m['field_aliases'].items()}
        require(v['response']=={'pair_dipoles_Debye':pairs,'C_dipoles_Debye':C,'chi_states_a0_cubed':chi},'Exact cycle/alias assembly')
        cell_links[key]=links
        if key==r['initial_key']:
            require(all(v[k]==baseline[k] for k in ('model','response','score')),'Exact saved initial arithmetic replay')
            require(all(c['inherited_native'] for c in v['cells'].values()),'Initial native calls zero; runtime cache identities may bind newly qualified executable')
        if key==r.get('selected_key'):selected_raw=v
    expected={r['initial_key']:('initial',None)};current=r['initial_key'];trajectory=[]
    def register(points,role,step):
        keys=[token(p) for p in points]
        for key,p in zip(keys,points):
            if key not in expected:expected[key]=(role,step)
            require(key in records and records[key]['parameters']==p,'Exact scheduled '+role+' vector')
        return keys
    require(len(r['iterations'])<=24,'Outer step cap')
    for i,it in enumerate(r['iterations'],1):
        require(it['outer_step']==i and it['origin_key']==current,'Sequential accepted origins')
        center=np.asarray(vectors[current]['parameters']);z=(center-lower)/width
        keys=register(probes(center,lower,upper),'derivative',i)
        if it['status']=='derivatives':require(r['status']=='failed_or_unavailable' and i==len(r['iterations']),'Only final failed derivative stage');break
        require(it['derivative_keys']==keys and all(k in scores for k in keys),'Full ordered derivative batch')
        spans=[(vectors[keys[2*j+1]]['parameters'][j]-vectors[keys[2*j]]['parameters'][j])/width[j] for j in range(12)]
        require(all(x>0 for x in spans),'Positive clipped derivative spans')
        J=np.column_stack([(np.asarray(scores[keys[2*j+1]]['weighted_residuals'])-scores[keys[2*j]]['weighted_residuals'])/spans[j] for j in range(12)])
        close(spans,it['derivative_spans'],'Derivative spans');close(J,it['Jacobian_normalized'],'Saved Jacobian arithmetic')
        rr=np.asarray(scores[current]['weighted_residuals']);diag=it['linear_directions'];lo=np.maximum(-.025,-z);hi=np.minimum(.025,1-z)
        close(lo,diag['trust_lower'],'Trust lower');close(hi,diag['trust_upper'],'Trust upper')
        full=np.asarray(diag['unrestricted_GN_not_evaluated']);factor=1.
        # SVD reconstructs the saved linear equation only; no new residual,
        # candidate point, nonlinear optimizer or BVLS call is made.
        U,singular,Vt=np.linalg.svd(J,full_matrices=False);active=singular>singular[0]*1e-6
        reconstructed=Vt[active].T@((U[:,active].T@(-rr))/singular[active])
        close(singular,diag['singular_values'],'Recorded GN spectrum')
        require(int(sum(active))==diag['rank'],'Relative-rcond GN rank')
        # Different LAPACK factorizations can differ by condition-amplified roundoff.
        condition=float(singular[0]/singular[active][-1]) if active.any() else 1.
        bound=512*np.finfo(float).eps*max(J.shape)*condition*max(1.,float(np.linalg.norm(full)),float(np.linalg.norm(reconstructed)))
        require(np.isfinite(full).all() and float(np.linalg.norm(full-reconstructed))<=bound,'Saved truncated GN direction')
        for x,l,u in zip(full,lo,hi):
            if x>0:factor=min(factor,u/x)
            elif x<0:factor=min(factor,l/x)
        close(full*factor,diag['scaled_GN'],'Scaled recorded GN direction')
        bounded=np.asarray(diag['bounded_LS']);require(np.isfinite(bounded).all() and np.all(bounded>=lo-1e-12) and np.all(bounded<=hi+1e-12),'Bounded LS feasibility')
        require(diag['bounded_solver']=='scipy.optimize.lsq_linear/bvls' and diag['status']>0 and (diag['iterations'] is None or diag['iterations']<=1000),'Declared successful BVLS')
        gradient=J.T@(J@bounded+rr);kkt=gradient.copy();kkt[bounded<=lo+1e-12]=np.minimum(gradient[bounded<=lo+1e-12],0);kkt[bounded>=hi-1e-12]=np.maximum(gradient[bounded>=hi-1e-12],0)
        close(gradient,diag['raw_gradient'],'BVLS gradient');close(max(abs(kkt)),diag['KKT_residual_max'],'Recorded KKT residual')
        points=[]
        for name in ('scaled_GN','bounded_LS'):
            for fraction in (.5,1.):
                p=np.clip(center+fraction*np.asarray(diag[name])*width,lower,upper);delta=(p-center)/width
                require(max(abs(delta))<=.025+1e-12,'Declared trust region')
                points.append({'method':name,'fraction':fraction,'parameters':p.tolist(),'predicted_objective':float(np.linalg.norm(rr+J@delta)**2)})
        require(len(it['trial_points'])==4,'Four fixed trials')
        for a,b in zip(points,it['trial_points']):
            require(all(a[k]==b[k] for k in ('method','fraction','parameters')),'Trial generation/order');close(a['predicted_objective'],b['predicted_objective'],'Predicted score')
        trial_keys=register([p['parameters'] for p in points],'trial',i)
        if it['status']=='trials':require(r['status']=='failed_or_unavailable' and i==len(r['iterations']),'Only final failed trial stage');break
        require(it['status']=='complete' and it['trial_keys']==trial_keys and all(k in scores for k in trial_keys),'Full completed four-trial batch')
        selected=choose(current,trial_keys,scores)
        require(it['selected_key']==selected and it['actual_improvement']==(selected!=current),'Strict current-wins-ties selection; derivative-only points excluded')
        trajectory.append({'outer_step':i,'origin_key':current,'selected_key':selected,'objective_before':scores[current]['objective'],'objective_after':scores[selected]['objective'],'actual_improvement':selected!=current,'KKT_residual_max':diag['KKT_residual_max']})
        require(selected!=current or i==len(r['iterations']),'Stop after no improvement');current=selected
    require(set(expected)==set(records),'No undeclared vector or omitted scheduled vector')
    for k,(role,step) in expected.items():require(records[k]['role']==role and records[k].get('outer_step')==step,'First-request role/step retained')
    require(r['current_key']==current,'Final accepted current key')
    if r['status']=='bounded_joint_response_fit_complete':
        require(not failures and r['error'] is None and r['candidate_frozen'] is True,'Completed candidate with no hidden failure')
        require(r['selected_key']==current and r['selected']==records[current]['artifact'] and r['completed_outer_steps']==len(trajectory)>0,'Frozen selected authority')
        require(r['last_Jacobian_origin_key']==r['iterations'][-1]['origin_key'],'Last Jacobian origin is distinct from selected point')
        reason='no_actual_improvement_at_complete_trial_batch' if not trajectory[-1]['actual_improvement'] else 'outer_step_limit'
        require(r['stopping_reason']==reason and (reason!='outer_step_limit' or len(trajectory)==24),'Original finite stop')
        require(r['selected_objective']==scores[current]['objective'] and r['selected_blocks']==scores[current]['block_objectives'] and r['original_EF_checks_passed']==scores[current]['original_EF_checks_passed'],'Selected score flags')
        contacts=[{'name':p['name'],'lower':x==p['lower'],'upper':x==p['upper']} for p,x in zip(parameters,vectors[current]['parameters'])]
        require(r['bound_contacts']==contacts,'Exact bound contacts')
    else:require(r['candidate_frozen'] is False and bool(r['error']),'Failed terminal outcome preserved')
    snapshot=r['evaluator'];ledger=read(snapshot['ledger'])
    require(snapshot['config_sha256']==ledger['config_sha256']==fingerprint(ec) and snapshot['implementation']==ec['implementation'],'Terminal ledger identity')
    require(snapshot['requests']==ledger['requests'] and snapshot['native_program_invocations']==ledger['native_program_invocations'],'Final ledger/snapshot exact')
    require(snapshot['distinct_parameter_vectors']==len(ledger['vector_keys'])<=673 and snapshot['native_program_invocations']<=24000 and snapshot['workers']==12,'Evaluator caps')
    require(len(ledger['vector_keys'])==len(set(ledger['vector_keys'])) and set(ledger['vector_keys'])<=set(records),'No undeclared evaluator vector')
    require(snapshot['reserved_outcomes_opened'] is False,'No reserved outcomes')
    requests={}
    for request in ledger['requests']:
        key=request['vector_key'];require(token(request['values'])==key and key not in requests and key in records,'Exactly one evaluator request per attempted vector');requests[key]=request
        require(set(request['cells'])==set(m['rows']),'46 ledger slots preserved')
        if request['status']=='complete':
            require(key in scores and request['objective']==scores[key]['objective'],'Completed evaluator score')
            require(request['cells']=={cid:{'status':'complete',**v} for cid,v in cell_links[key].items()},'All cell ledger joins')
    require(set(requests)==set(ledger['vector_keys']) and set(scores)<=set(requests),'Attempted versus scheduled requests')
    cache=Path(release['cache_directory']);require(snapshot['ledger']['path']==str(cache/'LEDGER.json'),'Exact cache directory')
    receipts=[];native_cells=[];raw_selected=0;ctx=SimpleNamespace(h=h,m=m,config=ec)
    for path in sorted((cache/'cells').glob('*/RESULT.json')):
        cp=pin(path);cell=read(cp);ck=path.parent.name;cid=cell['cell_id'];require(cid in m['rows'],'Declared cache cell')
        row=m['rows'][cid];values=cell['first_requested_vector'];expected_key,ident,meta=identity(values,cid)
        require(ck==expected_key and cell['cache_key']==ck and cell['cache_identity']==ident and cell['render']==meta,'Cache artifact identity')
        require(token(values) in requests and cell['implementation']==ec['implementation'] and cell['root_release']==dc['root_release'],'Actual producer vector/release')
        if ck in unique:require(token(cell)==unique[ck]['signature'],'Vector contains exact saved native cache result')
        programs=['analyze']+(['testgrad'] if row['needs_gradient'] else [])
        require(len(cell['receipts'])<=len(programs) and set(cell['program_inputs'])<=set(programs),'Finite program slots')
        for program,receipt in zip(programs,cell['receipts']):
            folder=path.parent/program;inp=cell['program_inputs'][program]
            for name,file in (('xyz','system.xyz'),('parameters','amoeba09.prm')):
                require(verify(inp[name])==folder/file and inp[name]['sha256']==row['source_input'][name]['sha256'],'Literal native geometry/parameters')
            key_text,_=render(values,cid,row['needs_export'] if program=='analyze' else False)
            require(verify(inp['key'])==folder/'system.key' and verify(inp['key']).read_text()==key_text,'Exact permitted native key')
            argv=[ec['executables'][program]['path'],'system.xyz']+(['EM'] if program=='analyze' else ['Y','N'])
            require(receipt['argv']==argv and receipt['cwd']==str(folder),'Native argv/cwd')
            require(set(receipt['threads'])==set(THREADS) and all(receipt['threads'][k]=='1' for k in THREADS),'Single-thread native receipts')
            require(verify(receipt['output'])==folder/(program+'.out'),'Native output pin')
            receipts.append({'cache_key':ck,'cell_id':cid,'program':program,**receipt})
        if cell['status']=='complete':require(len(cell['receipts'])==len(programs) and all(q.get('returncode')==0 for q in cell['receipts']),'Completed native denominator')
        native_cells.append({'cache_key':ck,'cell_id':cid,'result':cp,'status':cell['status'],'error':cell.get('error'),'programs':len(cell['receipts'])})
        if selected_raw is not None and selected_raw['cells'][cid]['cache_key']==ck:
            parsed=e.Evaluator._parse(ctx,cell,row)  # Saved parser only; no Evaluator construction.
            require(all(cell[k]==v for k,v in parsed.items()),'Selected native raw replay');raw_selected+=1
    require({x['cache_key'] for x in native_cells}>={k for k,v in unique.items() if not v['inherited']},'All noninherited cells have native authority')
    require(len(receipts)==snapshot['native_program_invocations'],'Every attempted native program counted once, including failures')
    if r['status']=='bounded_joint_response_fit_complete':require(all(x['status']=='complete' for x in native_cells) and not snapshot['stopped'] and not ledger['stopped'],'No concealed native failure')
    for x in r['cost'].values():require(isinstance(x,(float,int)) and math.isfinite(x) and x>=0,'Finite actual driver cost')
    summary={'status':'passed_terminal_joint_fit_review','fit_status':r['status'],'fit_result':cfg['fit_result'],'implementation':pin(__file__),'driver_configuration':cfg['driver_config'],
        'selected':r.get('selected'),'selected_parameters':vectors[current]['parameters'] if r.get('candidate_frozen') else None,
        'original_EF_checks_passed':r.get('original_EF_checks_passed'),'initial_objective':scores[r['initial_key']]['objective'],'initial_blocks':scores[r['initial_key']]['block_objectives'],
        'selected_objective':r.get('selected_objective'),'selected_blocks':r.get('selected_blocks'),'stopping_reason':r.get('stopping_reason',r.get('error')),
        'bound_contacts':r.get('bound_contacts'),'trajectory':trajectory,'all_vectors':len(records),'new_vectors':len(records)-1,'complete_vectors':len(scores),'failed_vectors':failures,
        'native_programs':len(receipts),'native_cells':native_cells,'selected_new_raw_cells_replayed':raw_selected,'selected_inherited_cells':46-raw_selected if selected_raw else None,
        'cost':{'driver':r['cost'],'sum_native_elapsed_seconds':math.fsum(x['elapsed_seconds'] for x in receipts),'native_receipt_child_CPU_not_summed':'Process-global child CPU deltas overlap across concurrent workers; use driver child CPU.',
                'scheduler_accounting':cfg.get('scheduler_accounting'),'scheduler_allocation_cost_available':cfg.get('scheduler_accounting') is not None},
        'new_native_calls':0,'new_parameter_evaluations':0,'new_optimization_calls':0,'physical_accuracy_qualified':False,'reserved_outcomes_opened':False,
        'independence':'Driver selection/provenance and response normalization independently reconstructed. Reviewer authored evaluator/objective; their frozen pure parsers/scorer are reused and not labeled independent oracles.',
        'limits':['Consumed development only; original EF gates remain separate.','Response scales are engineering weights, not physical error limits.','Recorded directions, trust bounds, scaling and KKT arithmetic checked; no BVLS or new optimization rerun.','No trimer regression or fresh validation outcomes opened.','Selected raw cells replayed; other native outputs are pin/receipt checked with exact stored scores.']}
    if cfg.get('scheduler_accounting'):summary['cost']['scheduler_accounting_record']=read(cfg['scheduler_accounting'])
    return summary


def main():
    p=argparse.ArgumentParser();p.add_argument('--config',required=True,type=Path);p.add_argument('--out',required=True,type=Path);args=p.parse_args()
    require(all(os.environ.get(k)=='1' for k in THREADS),'One local numerical thread')
    cfg=json.loads(args.config.read_text());require(cfg.get('fit_result') is not None,'Unbound source-only configuration; root must bind terminal result')
    tick,cpu=time.monotonic(),time.process_time();args.out.mkdir(parents=True,exist_ok=False)
    try:result=audit(cfg)
    except Exception as exc:result={'status':'unavailable_terminal_joint_fit_review','error':repr(exc),'implementation':pin(__file__),'fit_result':cfg.get('fit_result'),'new_native_calls':0,'new_parameter_evaluations':0}
    result.update(configuration=pin(args.config),checks=CHECKS,review_cost={'wall_seconds':time.monotonic()-tick,'CPU_seconds':time.process_time()-cpu,'threads':1})
    (args.out/'REVIEW.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    text='# Terminal joint-fit review\n\n'+result['status']+'\n\n'
    if 'error' in result:text+=result['error']+'\n'
    else:text+=f"Fit status: {result['fit_status']}. Original 25 EF gates pass: {result['original_EF_checks_passed']}.\n\nVectors {result['all_vectors']}; native programs {result['native_programs']}; stop {result['stopping_reason']}.\n\n"+result['independence']+'\n\n'+'\n'.join('- '+x for x in result['limits'])+'\n'
    (args.out/'REPORT.md').write_text(text)
    print(json.dumps({'status':result['status'],'review':pin(args.out/'REVIEW.json')}))
    return 0 if result['status']=='passed_terminal_joint_fit_review' else 1


if __name__=='__main__':raise SystemExit(main())
