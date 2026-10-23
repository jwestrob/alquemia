#!/usr/bin/env python3
"""Terminal saved-only audit. No molecular parser rerun, solve, or matrix factor."""
import argparse
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import time
from types import SimpleNamespace
import numpy as np
from scipy import sparse

HERE=Path(__file__).resolve().parent
CHECKS=0
CACHE={}
PINS={}

def require(ok,msg):
    global CHECKS
    CHECKS+=1
    if not ok:raise ValueError(msg)
def pin(path):
    p=Path(path).resolve();st=p.stat();k=(str(p),st.st_size,st.st_mtime_ns)
    if k not in PINS:PINS[k]={'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
    return PINS[k]
def verify(p):require(pin(p['path'])==p,'Immutable pin '+p['path']);return Path(p['path'])
def read(p):
    k=(p['path'],p['sha256'])
    if k not in CACHE:CACHE[k]=json.loads(verify(p).read_text())
    return CACHE[k]
def token(v):return hashlib.sha256(json.dumps(v,sort_keys=True,separators=(',',':'),allow_nan=False).encode()).hexdigest()
def module(p,name):
    s=importlib.util.spec_from_file_location(name,verify(p));m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
def forbidden(*a,**k):raise ValueError('Saved review forbids native/solver calls')
def near(a,b,label):
    a,b=np.asarray(a,dtype=float),np.asarray(b,dtype=float)
    require(a.shape==b.shape and np.isfinite(a).all() and np.isfinite(b).all(),label+' finite shape')
    require(np.all(abs(a-b)<=1e-10*np.maximum(1.,np.maximum(abs(a),abs(b)))),label+' inherited arithmetic agreement')
def compare(a,b,label):
    if isinstance(a,dict):
        require(isinstance(b,dict) and set(a)==set(b),label+' keys')
        for k in a:compare(a[k],b[k],label+'/'+k)
    elif isinstance(a,list):
        require(isinstance(b,list) and len(a)==len(b),label+' length')
        for i,(x,y) in enumerate(zip(a,b)):compare(x,y,label+'/'+str(i))
    elif isinstance(a,bool) or a is None or isinstance(a,str):require(a==b,label+' exact flag/text')
    elif isinstance(a,(float,int)):near(a,b,label)
    else:raise TypeError(label)

def arrays(p):
    for k in ('center','lower','upper','width','z','bounds_lower','bounds_upper','Jacobian_original531','Jacobian_components'):p[k]=np.asarray(p[k])
    for b in p['blocks'].values():b['r']=np.asarray(b['r']);b['J']=np.asarray(b['J'])
    p['energy']['r']=np.asarray(p['energy']['r']);p['energy']['J']=np.asarray(p['energy']['J'])
    for b in p['forces'].values():
        for k in ('r','J','weights'):b[k]=np.asarray(b[k])
    for b in p['electrical_unweighted'].values():b['r']=np.asarray(b['r']);b['J']=np.asarray(b['J'])
    return p

def verify_secants(core,p,center,vectors,cfg,sources):
    """Saved physical quantities/cycles only; no evaluator or new linearization."""
    parent=read(cfg['parent_configuration']);manifest=read(parent['manifest']);domain=read(parent['domain'])
    decomposition=read(sources['decomposition']);declaration=read(cfg['campaign_declaration'])
    require(len(manifest['rows'])==46 and len(manifest['static'])==48,'Original46/48 denominator')
    require(set(center['cells'])==set(manifest['rows']),'Exact original model cell IDs')
    lower=np.array([x['lower'] for x in domain['parameters']]);upper=np.array([x['upper'] for x in domain['parameters']])
    near(p['lower'],lower,'Physical lower box');near(p['upper'],upper,'Physical upper box')
    near(p['width'],upper-lower,'Physical widths');near(p['z'],(p['center']-lower)/(upper-lower),'Normalized center')
    require([v['parameters'] for v in vectors]==core.probe_points(center['parameters'],domain),'Exact clipped24 stencil')
    spans=[(vectors[2*j+1]['parameters'][j]-vectors[2*j]['parameters'][j])/p['width'][j] for j in range(12)]
    require(spans==p['spans'],'Actual normalized spans including asymmetric Dy_v')
    def secant(values):return np.column_stack([(values[2*j+1]-values[2*j])/spans[j] for j in range(12)])
    cp=core.physical(center['score']);vp=[core.physical(v['score']) for v in vectors]
    require(cp['all25'] and cp['energy_names']==p['energy']['names'] and all(x['energy_names']==cp['energy_names'] for x in vp),'Original25 physical naming and center pass')
    near(p['energy']['r'],cp['energy'],'Ten central physical energy errors');near(p['energy']['J'],secant([x['energy'] for x in vp]),'Ten physical energy secants')
    require(set(p['forces'])==set(cp['force']) and len(p['forces'])==15,'All15 physical force groups')
    for name,b in p['forces'].items():
        near(b['r'],cp['force'][name]['raw'],'Physical force center '+name);near(b['weights'],cp['force'][name]['weights'],'Physical body weights '+name)
        require(all(np.array_equal(x['force'][name]['weights'],b['weights']) for x in vp),'Unchanged body weights')
        near(b['J'],secant([x['force'][name]['raw'] for x in vp]),'Physical force secants '+name)
    static={k:v['energy_kcal_mol'] for k,v in manifest['static'].items()}
    def component(v):
        energies=dict(static,**{k:x['energy'] for k,x in v['cells'].items()});parents={}
        for name,row in decomposition['parents'].items():
            terms=row['model_terms'];n=int(name[4]);m=name.rsplit('_',1)[1]
            require(terms['C']=={name:1,'ion_'+m:-1,name+'__shell':-1} and terms['P']['ion_'+m]==-n,'Literal full-shell/pair cycle')
            b=dict(terms['C'])
            for k,a in terms['P'].items():b[k]=b.get(k,0)-a
            require(b==terms['B'],'Literal collective C-P cycle')
            parents[name]={k:math.fsum(float(a)*energies[cid] for cid,a in cycle.items()) for k,cycle in terms.items()}
        return np.asarray([parents[f'aquo{n}_central_Dy'][k]-parents[f'aquo{n}_central_La'][k]-decomposition['paired_Dy_minus_La']['W'+str(n)]['reference_total_kcal_mol'][k] for n in (8,9) for k in ('P','B')])
    ce=component(center);ve=[component(v) for v in vectors];jc=secant(ve)
    near(p['Jacobian_components'],jc,'Literal paired P/B secants')
    for name,ix,den in (('P_energy',[0,2],'J_P'),('B_energy',[1,3],'J_B')):
        b=p['blocks'][name];near(b['r'],ce[ix]/math.sqrt(2),'Component center '+name);near(b['J'],jc[ix]/math.sqrt(2),'Component normalized secants '+name)
        require(b['denominator']==declaration['objective']['fixed_denominators'][den],'Original origin24 component scale')
    for name,(a,b) in core.SLICES.items():
        block=p['blocks'][name];near(block['r'],center['score']['weighted_residuals'][a:b],'Original slice '+name);near(block['J'],p['Jacobian_original531'][a:b],'Original slice Jacobian '+name)
        require(block['denominator']==math.fsum(float(x*x) for x in block['r']),'Current center norm denominator')
    return {'derivative_physical_pass_count':sum(x['all25'] for x in vp),'derivative_physical_pass_required':False,'paired_component_order':['P8','B8','P9','B9'],'center_paired_component_errors_kcal_mol':ce.tolist(),'asymmetric_Dy_v_span':spans[7]}

def verify_matrix(p,case,model,layout):
    """Literal coefficient/sign correspondence, no call to conic assembly."""
    seen=[];offset=0;eye=np.eye(13)
    for row in layout:
        name=row['name'];require(row['start']==offset,'Contiguous cone layout');offset=row['stop'];seen.append(name)
        if row['kind']=='SOC':
            if name.startswith('block/'):
                key=name.split('/',1)[1];block=p['blocks'][key];scale=1/math.sqrt(block['denominator']);component=key in ('P_energy','B_energy')
                require(component or (case=='B' and key in ('pair_mu','C_mu','chi')),'Retained block family')
                top=-eye[12] if component else np.zeros(13);rhs=[0. if component else 1.]
            else:
                require(case!='D' and name.startswith('force/RMS/'),'Retained force cone')
                block=p['forces'][name[len('force/RMS/'):]];scale=np.sqrt(block['weights'])/2;top=np.zeros(13);rhs=[1.]
            J=block['J']*(np.asarray(scale)[...,None] if np.ndim(scale) else scale)
            expected=np.vstack([top,np.column_stack([-J,np.zeros(len(J))])]);rhs.extend(block['r']*scale)
        else:
            require(row['stop']==row['start']+1,'Scalar side')
            parts=name.split('/');kind,side=parts[:2];sign=1 if side=='upper' else -1
            if kind=='energy':
                require(case!='D','Physical energy removed only in D');i=p['energy']['names'].index('/'.join(parts[2:]));expected=np.r_[sign*p['energy']['J'][i],0.][None,:];rhs=[1-sign*p['energy']['r'][i]]
            elif kind=='force':
                require(case!='D','Physical components removed only in D');i=int(parts[-1]);b=p['forces']['/'.join(parts[2:-1])];expected=np.r_[sign*b['J'][i]/5,0.][None,:];rhs=[1-sign*b['r'][i]/5]
            elif kind in ('trust','box'):
                i=int(parts[-1]);expected=(sign*eye[i])[None,:]
                if kind=='trust':rhs=[.025]
                else:require(case!='C','Box removed only in C');rhs=[1-p['z'][i] if side=='upper' else p['z'][i]]
            else:
                require(kind=='epigraph','Known scalar family');expected=(sign*eye[12])[None,:];rhs=[math.sqrt(p['q_current']) if side=='upper' else 0.]
        near(model['A'][row['start']:row['stop']].toarray(),expected,case+'/'+name+' matrix')
        near(model['b'][row['start']:row['stop']],rhs,case+'/'+name+' rhs')
    require(offset==model['A'].shape[0] and len(seen)==len(set(seen)),'Complete unique row coverage')
    require(sum(r['kind']=='SOC' for r in layout)=={'A':17,'B':20,'C':17,'D':2}[case],'Exact case SOC count')
    require(sum(n.startswith('trust/') for n in seen)==24 and sum(n.startswith('box/') for n in seen)==(0 if case=='C' else 24),'Trust/domain faces')
    require(sum(n.startswith('energy/') for n in seen)==(0 if case=='D' else 20),'All energy inequalities')
    require(sum(n.startswith('force/') for n in seen)==(0 if case=='D' else 825),'All force RMS/component inequalities')

def run(c,out):
    tick,cpu=time.monotonic(),time.process_time();out=Path(out);out.mkdir(parents=True,exist_ok=False)
    report={'status':'started','implementation':pin(__file__),'configuration':c,'new_native_calls':0,'new_solver_calls':0,'error':None,
        'scope':'Saved journal/immutable cell and vector provenance, receipt hashes, actual secants/cone coefficients, retained-case certificates. Native scientific parsers and full objective/native evaluation are not rerun; their admitted producer is reused.'}
    try:
        require(c['implementation']==pin(__file__) and str(out.resolve())==c['output_directory'],'Exact reviewer/output binding')
        require(all(os.environ.get(k)=='1' for k in ('OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS')),'One review thread')
        r=read(c['result']);cfg=read(r['configuration']);ready=read(c['runner_ready'])
        require(cfg['implementation']==ready['implementation'] and cfg['solver']==ready['solver'],'Frozen runner/solver')
        require(r['configuration_fingerprint']==ready['configuration_fingerprint']==token({k:v for k,v in cfg.items() if k not in ('root_release','source_review','preparation')}),'Exact runtime fingerprint')
        runner=module(cfg['implementation'],'review_final_runner');runner.check_release(cfg)
        old=module(c['saved_reviewer'],'review_inherited_journal');journal=module(cfg['journal_implementation'],'review_final_journal')
        inv=read(cfg['inventory']);inp=read(inv['inputs']);sources=inp['sources'];center=read(sources['selected'])
        report.update(result=c['result'],producer=r['configuration'],root_release=cfg['root_release'],source_review=cfg['source_review'],center=sources['selected'],inventory=cfg['inventory'])
        successful=r['status']=='completed_final_local_attribution_panel' and r['error'] is None
        if 'evaluator' not in r:
            require(not successful,'Successful result needs evaluator');report.update(status='reviewed_failed_panel_before_evaluator',producer_status=r['status'],producer_error=r['error']);return report
        ledger=read(r['evaluator']['ledger']);lp=Path(r['evaluator']['ledger']['path']);jp=lp.with_name('EVENTS.jsonl')
        replay,js,initial,intents,returns,cellpins=old.stream_journal(jp,journal)
        require(replay==ledger,'Inherited journal replay equals complete ledger')
        side=read(pin(lp.with_name(lp.name+'.journal.json')))
        require(side['ledger']==r['evaluator']['ledger'] and side['sequence']==js['sequence'] and side['last_event_sha256']==js['last_event_sha256'],'Final journal sidecar')
        require(initial['authority']=={'configuration_fingerprint':r['configuration_fingerprint'],'implementation':cfg['implementation'],'root_release':cfg['root_release']},'Exact journal initialization authority')
        require(ledger['limits']==cfg['limits'] and ledger['native_program_invocations']==r['native_programs']<=672,'Actual charged cap')
        require(len(ledger['requests'])==len(ledger['vector_keys'])<=24 and set(intents)==set(range(1,r['native_programs']+1)),'Exact request/intent accounting')
        report['journal']=js
        if not successful:
            report.update(status='reviewed_failed_or_incomplete_final_local_panel',producer_status=r['status'],producer_error=r['error'],actual_native_programs=r['native_programs'],solver_calls=r['solver_calls'],cost=r['cost']);return report
        require(not ledger['stopped'] and not js['partial_suffix_bytes'] and not js['unresolved_attempt_numbers'],'No hidden failure/partial journal')
        require(r['native_programs']==672 and len(r['vectors'])==24 and len(ledger['requests'])==24 and len(returns)==672,'Exact24/672 completion')
        require(all(a['status']=='returned_success' and a['returncode_observed'] for a in ledger['native_attempts']),'All actual attempts returned successfully')
        expected=[x['key'] for x in inv['probe_rows']];require(ledger['vector_keys']==expected and list(r['vectors'])==expected,'Exact frozen probe ordering')
        require(all(q['status']=='complete' and len(q['cells'])==46 and all(v['status']=='complete' for v in q['cells'].values()) for q in ledger['requests']),'All46 cells per vector')
        vectors=[];new_cells={};imported=0
        for row in inv['probe_rows']:
            meta=r['vectors'][row['key']];require(meta['status']=='complete' and meta['parameters']==row['parameters'],'Exact complete vector metadata')
            v=read(meta['artifact']);require(v['status']=='complete' and v['parameters']==row['parameters'] and token(v['parameters'])==row['key'],'Exact saved vector identity');require(set(v['cells'])==set(center['cells']),'Full46 denominator');vectors.append(v)
            for cid,rec in v['cells'].items():
                require(rec['status']=='complete' and rec['request_provenance']['vector']==row['parameters'],'Complete cell/request identity')
                # The inventory's literal cell lists are authoritative, not names.
                metal=next(m for m,x in row['metal_cache_matches'].items() if cid in x['cells'])
                require(rec['cache_key']==rec['request_provenance']['cache_key'] and rec['cache_identity']==rec['request_provenance']['cache_identity'],'Scientific cache joins')
                if row['metal_cache_matches'][metal]['matching_saved_vector'] is not None:
                    require({k:x for k,x in rec.items() if k!='request_provenance'}=={k:x for k,x in center['cells'][cid].items() if k!='request_provenance'},'Exact unchanged-metal record, no LAPACK replay');imported+=1
                else:
                    key=rec['cache_key'];require(key not in new_cells,'Unique new cell executed once');new_cells[key]=rec
                    body=read(pin(lp.parent/'cells'/key/'RESULT.json'))
                    require(body=={k:x for k,x in rec.items() if k!='request_provenance'},'Immutable final new-cell body')
                    require(rec['root_release']==cfg['root_release'] and rec['job_id'] is None,'Truthful new-cell release/login')
                    for inputs in rec['program_inputs'].values():
                        for p in inputs.values():verify(p)
        require(len(new_cells)==552 and imported==552,'Exactly552 new and552 unchanged-metal cell uses')
        observed=[]
        for num,returned in returns.items():
            attempt=returned['attempt'];rec=new_cells[attempt['cache_key']];receipts=returned['receipt_delta'];require(len(receipts)==1 and receipts[0]['returncode']==0,'One successful invocation receipt')
            receipt=receipts[0];require(receipt in rec['receipts'] and returned['program_inputs']==rec['program_inputs'][attempt['program']],'Receipt/cell/program/input join')
            verify(receipt['output']);require(receipt['executor']['job_id'] is None and receipt['executor']['mode']=='login','Truthful native receipt');observed.append((attempt['cache_key'],attempt['program']))
        require(len(set(observed))==672 and sum(len(x['receipts']) for x in new_cells.values())==672,'672 unique programs, no step doublecount')
        report['native_accounting']={'vectors':24,'native_programs':672,'new_cells':552,'unchanged_metal_cell_uses':552,'unique_selected_cells_reused':46,'raw_native_parser_repeated':False}
        core=module(cfg['core'],'review_final_core');solver=module(cfg['solver'],'review_final_solver');solver.clarabel.DefaultSolver=forbidden
        affine=read(r['affine_result']);require(affine['status']=='completed_four_fixed_affine_diagnostics' and affine['error'] is None and affine['solver_calls']==affine['solver_constructions']==r['solver_calls']==4,'Exactly four completed solve attempts')
        require(affine['center_key']==inv['center_key'] and affine['derivative_keys']==expected and affine['settings']==cfg['settings'],'One exact origin/stencil/settings for all cases')
        require(affine['effective_settings']==read(cfg['settings']),'All stored effective settings frozen')
        p=arrays(read(affine['linearization']));near(p['center'],center['parameters'],'Final linearization center');require(p['spans']==[x['actual_normalized_span'] for x in inv['spans']],'Actual asymmetric secant spans')
        near(p['Jacobian_original531'],np.column_stack([(np.asarray(vectors[2*j+1]['score']['weighted_residuals'])-np.asarray(vectors[2*j]['score']['weighted_residuals']))/p['spans'][j] for j in range(12)]),'Actual531 saved secants')
        report['linearization_review']=verify_secants(core,p,center,vectors,cfg,sources)
        require(p['raw_center_q']==r['center_metrics']['q'] and p['q_current']==p['raw_center_q']+1e-8,'Frozen auxiliary center cushion')
        for channel,b in p['electrical_unweighted'].items():
            near(b['r'],[center['response'][channel][k] for k in b['names']],channel+' center')
            J=np.stack([(np.asarray([vectors[2*j+1]['response'][channel][k] for k in b['names']])-np.asarray([vectors[2*j]['response'][channel][k] for k in b['names']]))/p['spans'][j] for j in range(12)],axis=-1)
            near(b['J'],J,channel+' signed secants')
        result_cases={}
        for case in ('A','B','C','D'):
            row=affine['cases'][case];require(row['status']=='complete','Complete case '+case);attempt=read(row['attempt'])
            require(attempt['case']==case and attempt['maximum_attempts']==1 and attempt['configuration_fingerprint']==r['configuration_fingerprint'],'Exact one-attempt case receipt')
            require(row['solver_success']==(row['solver_status']=='Solved') and row['solution']['status']==row['info']['status']==row['solver_status'],'Raw solver status preserved '+case)
            with np.load(verify(row['conic_data']),allow_pickle=False) as z:
                model={'A':sparse.csc_matrix((z['A_data'],z['A_indices'],z['A_indptr']),shape=tuple(z['A_shape'])),'b':z['b'].copy(),'c':z['c'].copy(),'layout':read(row['row_layout'])}
            verify_matrix(p,case,model,model['layout']);near(model['c'],np.r_[np.zeros(12),1.],'Minimize t objective')
            residual=core.conic.diagnostics(model,SimpleNamespace(**row['solution']));compare(row['independent_residuals'],core.conic.plain(residual),'Residuals/'+case)
            z=np.asarray(row['solution']['z'])
            for item in model['layout']:
                a,b=item['start'],item['stop'];dual=row['dual_by_named_constraint'][item['name']]
                near(dual['vector'],z[a:b],'Exact named dual vector');near(dual['minus_b_dot_z'],-model['b'][a:b]@z[a:b],'Named dual contribution')
            x=np.asarray(row['solution']['x']);require(x.shape==(13,) and np.isfinite(x).all() and x[-1]>=0,'Finite nonnegative-t return '+case)
            require(len(row['predictions'])==2,'Full and half predictions only')
            for frac,saved in zip((1.,.5),row['predictions']):
                calculated=solver.certificate(core,p,np.r_[x[:12],x[-1]**2],frac,case,row['solver_success'])
                compare(saved,calculated,'Certificate/'+case+'/'+str(frac))
            full=row['predictions'][0]
            result_cases[case]={'solver_status':row['solver_status'],'solver_success':row['solver_success'],'actual_component_q':full['actual_component_q'],
                'q_gain_from_exact_final':full['q_gain_from_exact_final'],'P_B_RMS_kcal_mol':full['P_B_RMS_kcal_mol'],
                'all25_predicted':full['all25_predicted'],'box_pass':full['literal_original_box_pass'],
                'retained_exact_sign_feasible':full['retained_exact_sign_feasible'],'retained_numerically_feasible_1e_8':full['retained_numerically_feasible_1e_8'],
                'minimum_retained_margin':full['minimum_retained_margin'],'raw_epigraph_margin':full['raw_epigraph_margin'],
                'electrical_RMS':full['electrical_RMS'],'signed_electrical_errors':full['signed_electrical_errors'],
                'chi_conditional_increments_a0_cubed':full['chi_conditional_increments_a0_cubed'],'W9_C_Dy_minus_La_Debye':full['W9_C_Dy_minus_La_Debye'],
                'physical_gate_failures':{k:v for k,v in full['original25'].items() if not v['passed']},
                'negative_retained_margins':{k:v for k,v in full['all_retained_scalar_margins'].items() if v<0},
                'half_step':{k:row['predictions'][1][k] for k in ('actual_component_q','q_gain_from_exact_final','retained_exact_sign_feasible','minimum_retained_margin','all25_predicted','literal_original_box_pass')},
                'primal_residual_inf':residual['primal_residual_inf'],'dual_residual_inf':residual['dual_residual_inf'],'signed_gap_t':residual['signed_gap_t'],
                'cost':row['cost'],'no_candidate_admission':True}
        report.update(status='passed_final_local_attribution_saved_review',affine_result=r['affine_result'],linearization=affine['linearization'],cases=result_cases,
            actual_cost=r['cost'],measured_total_CPU_seconds=math.fsum(r['cost'][k] for k in ('controller_CPU_seconds','child_user_CPU_seconds','child_system_CPU_seconds')),
            interpretation='Verified saved finite local predictions and all actual derivative receipts. Solver status/removed constraints remain explicit; no nonlinear candidate or physical accuracy admission.')
        return report
    except Exception as exc:report.update(status='failed_final_local_saved_review',error=repr(exc));return report
    finally:
        report['checks']=CHECKS;report['inherited_journal_checks']=getattr(locals().get('old'), 'CHECKS',0)
        report['review_cost']={'wall_seconds':time.monotonic()-tick,'CPU_seconds':time.process_time()-cpu,'threads':1}
        (out/'REVIEW.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')

if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('--config',required=True);a.add_argument('--out',required=True);args=a.parse_args()
    result=run(json.loads(Path(args.config).read_text()),args.out)
    print(json.dumps({'status':result['status'],'review':pin(Path(args.out)/'REVIEW.json'),'error':result['error']}))
    raise SystemExit(bool(result['error']))
