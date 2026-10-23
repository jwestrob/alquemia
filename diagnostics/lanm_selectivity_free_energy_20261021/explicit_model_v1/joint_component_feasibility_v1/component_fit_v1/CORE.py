"""Pure bounded campaign mathematics. Only propose() calls one numerical solver.

No file I/O, molecular evaluator, scheduler, native subprocess or reference API.
The driver owns immutable input/native admission and the campaign's finite budget.
"""
import hashlib
import json
import math
import numpy as np
from scipy.optimize import minimize

FRACTIONS=(1.,.5,.25,.125)
GEOMETRIES=('water_minus','water_central','water_plus','aquo8_central','aquo9_central')
SLICES={'old_EF':(0,413),'pair_mu':(413,515),'C_mu':(515,527),'chi':(527,531)}
SCORE_KEYS={'old_EF':'old_EF','pair_mu':'pair','C_mu':'C','chi':'chi'}
COMPONENT_ORDER=('P8','B8','P9','B9')
OPTIONS={'maxiter':200,'ftol':1e-12,'disp':False}


def require(ok,msg):
    if not ok:raise ValueError(msg)


def token(v):return hashlib.sha256(json.dumps(v,sort_keys=True,separators=(',',':'),allow_nan=False).encode()).hexdigest()


def near(a,b,msg):
    a,b=np.asarray(a,dtype=float),np.asarray(b,dtype=float)
    require(a.shape==b.shape and np.all(np.isfinite(a)) and np.all(np.isfinite(b)),msg+' shape/finite')
    require(np.all(abs(a-b)<=1e-10*np.maximum(1.,np.maximum(abs(a),abs(b)))),msg+' arithmetic')


def checked(values,domain):
    require(len(values)==12 and len(domain['parameters'])==12,'Twelve box parameters')
    require(all(isinstance(v,(int,float)) and not isinstance(v,bool) and math.isfinite(v) and p['lower']<=v<=p['upper'] for v,p in zip(values,domain['parameters'])),'Finite in-box vector')
    return list(values)


def probe_points(values,domain):
    center=np.asarray(checked(values,domain));lower=np.array([p['lower'] for p in domain['parameters']]);upper=np.array([p['upper'] for p in domain['parameters']]);width=upper-lower;z=(center-lower)/width
    result=[]
    for j in range(12):
        for sign in (-1,1):
            target=z.copy();target[j]=np.clip(target[j]+sign*.005,0,1)
            value=center.copy();value[j]=np.clip(center[j]+(target[j]-z[j])*width[j],lower[j],upper[j]);result.append(value.tolist())
    return result


def physical(score):
    old=score['old_EF'];energy=[];energy_names=[];force={};reconstruct=[];gate_flags=[]
    for group,rows in old['energy_controls'].items():
        for name,row in rows.items():
            error=row['error_kcal_mol'];require(math.isfinite(error),'Finite physical energy error')
            energy.append(error);energy_names.append(group+'/'+name);gate_flags.append(abs(error)<=1)
            require(row['passed']==gate_flags[-1],'Unchanged energy gate')
    require(len(energy)==10,'Ten energy gates')
    errors=old['energy_errors_kcal_mol']
    for group in ('balanced_exchange','individual_radial_work'):
        reconstruct.extend(v*math.sqrt(.5/4)*.5 for v in errors[group].values())
    for kind in ('La','Dy','paired'):
        for geo in GEOMETRIES:
            row=old['per_fixture'][kind][geo];raw=np.asarray(row['raw_body_components_kcal_mol_A'],dtype=float)
            n=8 if geo=='aquo8_central' else 9 if geo=='aquo9_central' else 1
            require(raw.shape==(3+6*n,) and np.all(np.isfinite(raw)),'Complete finite body components')
            weights=np.r_[np.full(3,1/6),np.full(6*n,1/(12*n))]
            rms=math.sqrt(math.fsum(float(w*x*x) for w,x in zip(weights,raw)));maximum=max(abs(float(x)) for x in raw)
            near([rms,maximum],[row['rigid_RMS_kcal_mol_A'],row['maximum_body_component_error_kcal_mol_A']],'Original physical force replay')
            # Preserve the actual native scorer's exact-threshold gate byte; the
            # independent reconstruction is an arithmetic consistency check.
            expected=row['rigid_RMS_kcal_mol_A']<=2 and row['maximum_body_component_error_kcal_mol_A']<=5
            require(row['passed']==expected,'Unchanged force gate');gate_flags.append(expected)
            force[kind+'/'+geo]={'raw':raw,'weights':weights,'RMS':row['rigid_RMS_kcal_mol_A'],'maximum':row['maximum_body_component_error_kcal_mol_A'],'passed':expected}
            reconstruct.extend((raw*np.sqrt(weights)/(4*math.sqrt(5))).tolist())
    near(reconstruct,old['weighted_residuals'],'Complete original413 reconstruction')
    near(np.asarray(reconstruct)/2,score['weighted_residuals'][:413],'Original531 EF slice')
    for sign in ('minus','plus'):near(errors['paired_radial_work'][sign],errors['individual_radial_work']['Dy_'+sign]-errors['individual_radial_work']['La_'+sign],'Paired radial closure')
    require(len(gate_flags)==25 and all(gate_flags)==score['original_EF_checks_passed']==old['development_checks_passed'],'Original25 summary')
    return {'energy':np.asarray(energy),'energy_names':energy_names,'force':force,'all25':all(gate_flags)}


def components(vector,fixed):
    energies=dict(fixed['static'],**{k:v['energy'] for k,v in vector['cells'].items()})
    require(all(math.isfinite(e) for e in energies.values()),'Finite saved cell energies')
    parents={}
    for name,channels in fixed['terms'].items():
        parents[name]={k:math.fsum(c*energies[cid] for cid,c in terms.items()) for k,terms in channels.items()}
        near(parents[name]['C'],parents[name]['P']+parents[name]['B'],'C=P+B')
    error=[]
    for n in (8,9):
        for k in ('P','B'):
            paired=parents[f'aquo{n}_central_Dy'][k]-parents[f'aquo{n}_central_La'][k]
            error.append(paired-fixed['references']['W'+str(n)][k])
    return np.asarray(error),parents


def measure(vector,fixed):
    require(vector['status']=='complete' and set(vector['cells'])==set(fixed['row_ids']) and all(v['status']=='complete' for v in vector['cells'].values()),'Complete native-admitted46-vector contract')
    checked(vector['parameters'],fixed['domain']);score=vector['score'];r=np.asarray(score['weighted_residuals'],dtype=float)
    require(r.shape==(531,) and np.all(np.isfinite(r)),'Finite531 residuals')
    require({n:score['residual_layout'][SCORE_KEYS[n]] for n in SLICES}=={n:list(ix) for n,ix in SLICES.items()},'Exact531 layout')
    for name,(a,b) in SLICES.items():near(4*float(r[a:b]@r[a:b]),score['block_objectives'][SCORE_KEYS[name]],'Authoritative original block score')
    p=physical(score);e,parents=components(vector,fixed)
    return {'r':r,'physical':p,'components':e,'parents':parents}


def make_fixed(origin,manifest,decomposition,domain,declaration):
    require(declaration['limits']=={'outer_steps':8,'new_parameter_vectors':200,'new_native_programs':6496,'reused_origin_vectors':1,'reused_initial_derivative_vectors':24,'workers':12,'numerical_threads_per_process':1,'aggregate_login_threads':16,'native_timeout_seconds':300},'Declared one finite campaign')
    require(declaration['proposal']['fractions']==list(FRACTIONS) and declaration['proposal']['maxiter']==200 and declaration['proposal']['ftol']==1e-12 and declaration['proposal']['normalized_trust_radius']==.025,'Frozen proposal settings')
    require(token(origin['parameters'])==declaration['start']['key'],'Fixed origin24 identity')
    require(len(manifest['rows'])==46 and len(manifest['static'])==48,'Original46+48 denominator')
    terms={k:v['model_terms'] for k,v in decomposition['parents'].items()}
    require(set(terms)=={f'aquo{n}_central_{m}' for n in (8,9) for m in ('La','Dy')},'All four literal parents')
    for name,channels in terms.items():
        n=int(name[4]);metal=name.rsplit('_',1)[1]
        require(channels['C']=={name:1,'ion_'+metal:-1,name+'__shell':-1} and channels['P']['ion_'+metal]==-n,'Literal full-shell/N-fold ion cycles')
        b=dict(channels['C'])
        for k,c in channels['P'].items():b[k]=b.get(k,0)-c
        require(b==channels['B'],'Expanded B=C-P cycle')
        require(decomposition['parents'][name]['reference_electronic_plus_D4_kcal_mol'] is None,'No invented absolute-parent total target')
    refs={k:v['reference_total_kcal_mol'] for k,v in decomposition['paired_Dy_minus_La'].items()}
    for n,row in decomposition['paired_Dy_minus_La'].items():
        for k in ('C','P','B'):near(refs[n][k],row['reference_electronic_kcal_mol'][k]+row['reference_D4_kcal_mol'][k],'Separate electronic plus D4 authority')
    fixed={'domain':domain,'row_ids':list(manifest['rows']),'static':{k:v['energy_kcal_mol'] for k,v in manifest['static'].items()},'terms':terms,'references':refs,
           'ceilings':dict(declaration['fixed_old_block_ceilings']),'denominators':dict(declaration['objective']['fixed_denominators']),
           'origin_key':declaration['start']['key'],'origin_parameters':list(origin['parameters']),'limits':dict(declaration['limits'])}
    require(manifest['initial_vector']==[p['initial'] for p in domain['parameters']],'Original initial-model identity')
    baseline=[decomposition['paired_Dy_minus_La']['W'+str(n)]['errors']['baseline']['model_minus_reference_kcal_mol'][k] for n in (8,9) for k in ('P','B')]
    fixed['original_initial_component_RMS']={k:math.sqrt(math.fsum(baseline[i]**2 for i in ix)/2) for k,ix in (('P',(0,2)),('B',(1,3)))}
    m=measure(origin,fixed);e=m['components']
    near([math.fsum(float(x*x) for x in e[[0,2]])/2,math.fsum(float(x*x) for x in e[[1,3]])/2],[fixed['denominators']['J_P'],fixed['denominators']['J_B']],'Frozen origin component denominators')
    require(all(v>1e-24 for v in (*fixed['denominators'].values(),*fixed['ceilings'].values())),'Nondegenerate fixed normalization')
    require(m['physical']['all25'],'Qualified fixed origin')
    for k,v in fixed['ceilings'].items():require(origin['score']['block_objectives'][SCORE_KEYS[k]]==v,'Exact fixed original score ceiling')
    return fixed


def assess(vector,fixed):
    m=measure(vector,fixed);e=m['components'];score=vector['score']
    JP=math.fsum(float(x*x) for x in e[[0,2]])/2;JB=math.fsum(float(x*x) for x in e[[1,3]])/2
    ratio={'P':JP/fixed['denominators']['J_P'],'B':JB/fixed['denominators']['J_B']}
    old={k:score['block_objectives'][SCORE_KEYS[k]] for k in SLICES};ceilings={k:old[k]<=fixed['ceilings'][k] for k in SLICES}
    rms={'P':math.sqrt(JP),'B':math.sqrt(JB)}
    baseline_improvement={k:rms[k]<v for k,v in fixed['original_initial_component_RMS'].items()}
    gates={'energy/'+name:{'error_kcal_mol':float(v),'passed':bool(abs(v)<=1)} for name,v in zip(m['physical']['energy_names'],m['physical']['energy'])}
    for name,row in m['physical']['force'].items():
        gates['force/'+name]={'rigid_RMS_kcal_mol_A':row['RMS'],'maximum_body_component_error_kcal_mol_A':row['maximum'],'passed':row['passed']}
    return {'key':token(vector['parameters']),'q':max(ratio.values()),'J_P':JP,'J_B':JB,'component_RMS':rms,
            'component_errors_kcal_mol':dict(zip(COMPONENT_ORDER,e.tolist())),'component_squared_error_ratios':ratio,
            'original_block_scores':old,'fixed_ceiling_pass':ceilings,'all25_original_physical_gates':m['physical']['all25'],
            'original25':gates,
            'physical_constraints_pass':m['physical']['all25'] and all(ceilings.values()),
            'full_development_target':math.sqrt(JP)<=1 and math.sqrt(JB)<=1 and m['physical']['all25'] and all(ceilings.values()),
            'original_initial_component_RMS':fixed['original_initial_component_RMS'],
            'improved_over_original_initial_component_RMS':baseline_improvement,
            'substantive_recovery_milestone':max(ratio.values())<=.25 and all(baseline_improvement.values()) and m['physical']['all25'] and all(ceilings.values()),
            'q_at_or_below_quarter':max(ratio.values())<=.25,'native_admission_is_driver_responsibility':True}


def linearize(current,probes,fixed):
    require(len(probes)==24,'Exactly24 derivative probes')
    expected=probe_points(current['parameters'],fixed['domain'])
    center=np.asarray(current['parameters']);lower=np.asarray([p['lower'] for p in fixed['domain']['parameters']]);upper=np.asarray([p['upper'] for p in fixed['domain']['parameters']]);width=upper-lower;z=(center-lower)/width
    c=measure(current,fixed);metrics=assess(current,fixed);require(metrics['physical_constraints_pass'],'Current accepted incumbent')
    measured=[]
    for vector,point in zip(probes,expected):
        require(vector['parameters']==point,'Exact original normalized probe order/coordinates');measured.append(measure(vector,fixed))
    cols=[];extra=[];energy=[];forces={k:[] for k in c['physical']['force']};spans=[]
    for j in range(12):
        a,b=measured[2*j:2*j+2];span=(probes[2*j+1]['parameters'][j]-probes[2*j]['parameters'][j])/width[j]
        require(span>0 and math.isfinite(span),'Actual nonzero normalized span');spans.append(float(span))
        cols.append((b['r']-a['r'])/span);extra.append((b['components']-a['components'])/span)
        require(a['physical']['energy_names']==b['physical']['energy_names']==c['physical']['energy_names'],'Physical energy gate order')
        energy.append((b['physical']['energy']-a['physical']['energy'])/span)
        for name in forces:
            for item in (a,b):require(np.array_equal(item['physical']['force'][name]['weights'],c['physical']['force'][name]['weights']),'Fixed geometry/body weights')
            forces[name].append((b['physical']['force'][name]['raw']-a['physical']['force'][name]['raw'])/span)
    J=np.column_stack(cols);jc=np.column_stack(extra)
    blocks={k:{'r':c['r'][a:b],'J':J[a:b],'denominator':fixed['ceilings'][k]/4} for k,(a,b) in SLICES.items()}
    for name,ix,den in (('P_energy',[0,2],'J_P'),('B_energy',[1,3],'J_B')):blocks[name]={'r':c['components'][ix]/math.sqrt(2),'J':jc[ix]/math.sqrt(2),'denominator':fixed['denominators'][den]}
    return {'center':center,'lower':lower,'upper':upper,'width':width,'z':z,'q_current':metrics['q'],'blocks':blocks,
            'energy':{'r':c['physical']['energy'],'J':np.column_stack(energy),'names':c['physical']['energy_names']},
            'forces':{k:{'r':v['raw'],'J':np.column_stack(forces[k]),'weights':v['weights']} for k,v in c['physical']['force'].items()},
            'bounds_lower':np.r_[np.maximum(-.025,-z),0.],'bounds_upper':np.r_[np.minimum(.025,1-z),metrics['q']],
            'spans':spans,'Jacobian_original531':J,'Jacobian_components':jc,'current_metrics':metrics}


def constraints(problem,x):
    d,q=x[:-1],x[-1];values=[];rows=[];names=[]
    def add(name,value,gradient):names.append(name);values.append(float(value));rows.append(gradient)
    for name,b in problem['blocks'].items():
        v=b['r']+b['J']@d;component=name in ('P_energy','B_energy');den=b['denominator']
        add('block/'+name,(q if component else 1.)-float(v@v)/den,np.r_[-2*v@b['J']/den,1. if component else 0.])
    e=problem['energy'];v=e['r']+e['J']@d
    for i,name in enumerate(e['names']):
        add('energy/upper/'+name,1-v[i],np.r_[-e['J'][i],0.]);add('energy/lower/'+name,1+v[i],np.r_[e['J'][i],0.])
    for name,b in problem['forces'].items():
        v=b['r']+b['J']@d;w=b['weights']
        add('force/RMS/'+name,1-float(np.dot(w,v*v))/4,np.r_[-(w*v)@b['J']/2,0.])
        for i in range(len(v)):
            add('force/upper/'+name+'/'+str(i),1-v[i]/5,np.r_[-b['J'][i]/5,0.]);add('force/lower/'+name+'/'+str(i),1+v[i]/5,np.r_[b['J'][i]/5,0.])
    return np.asarray(values),np.asarray(rows),names


def certify_fraction(problem,returned,fraction,solver_success):
    """Independent scalar fsum arithmetic; never calls constraints() or optimizer."""
    require(fraction in FRACTIONS,'Only four declared fractions')
    raw=[float(x) for x in returned];require(len(raw)==13 and all(math.isfinite(x) for x in raw),'Finite returned solver point')
    d=[fraction*x for x in raw[:12]];q=(1-fraction)*problem['q_current']+fraction*raw[-1]
    params=[float(x+w*y) for x,w,y in zip(problem['center'],problem['width'],d)]
    margins={};blocks={};gates={}
    def affine(r,j):return [float(v)+math.fsum(float(a)*b for a,b in zip(row,d)) for v,row in zip(r,j)]
    for name,b in problem['blocks'].items():
        v=affine(b['r'],b['J']);sq=math.fsum(x*x for x in v);ratio=sq/b['denominator'];component=name in ('P_energy','B_energy')
        margins['block/'+name]=(q if component else 1.)-ratio
        blocks[name]={'predicted_squared_norm':sq,'fixed_denominator':b['denominator'],'ratio':ratio}
    e=problem['energy'];values=affine(e['r'],e['J'])
    for name,v in zip(e['names'],values):
        margins['energy/upper/'+name]=1-v;margins['energy/lower/'+name]=1+v
        gates['energy/'+name]={'error_kcal_mol':v,'passed':abs(v)<=1}
    for name,b in problem['forces'].items():
        values=affine(b['r'],b['J']);rms2=math.fsum(float(w)*v*v for w,v in zip(b['weights'],values));maximum=max(abs(v) for v in values)
        margins['force/RMS/'+name]=1-rms2/4
        for i,v in enumerate(values):margins['force/upper/'+name+'/'+str(i)]=1-v/5;margins['force/lower/'+name+'/'+str(i)]=1+v/5
        gates['force/'+name]={'rigid_RMS_kcal_mol_A':math.sqrt(rms2),'maximum_body_component_error_kcal_mol_A':maximum,'passed':rms2<=4 and maximum<=5}
    require(len(margins)==851 and len(gates)==25,'Exact independent851/25 certificate')
    box=[y-float(lo) for y,lo in zip(d,problem['bounds_lower'][:12])]+[float(hi)-y for y,hi in zip(d,problem['bounds_upper'][:12])]
    physical_box=all(float(lo)<=p<=float(hi) for p,lo,hi in zip(params,problem['lower'],problem['upper']))
    valid_box=physical_box and min(box)>=0 and 0<=q<=problem['q_current']
    all_gates=all(v['passed'] for v in gates.values());old_ceiling=all(blocks[k]['ratio']<=1 for k in SLICES)
    residual_min=min(margins.values());epigraph=min(margins['block/'+k] for k in ('P_energy','B_energy'))
    improved=q<problem['q_current']
    # Successful status also needs an exact nonrelaxed physical/ceiling check.
    # The failed-status branch admits only strict positive residual margins.
    witness=(residual_min>=0 if solver_success else residual_min>0)
    eligible=valid_box and all_gates and old_ceiling and improved and witness
    return {'fraction':fraction,'parameters':params,'key':token(params),'eligible':eligible,
            'certificate':{'method':'independent_scalar_fsum','solver_success':bool(solver_success),'failed_solver_status_remains_failed':not solver_success,
                           'q_fraction':q,'q_current':problem['q_current'],'strict_q_improvement':improved,'in_box_and_trust':valid_box,
                           'minimum_box_margin':min(box),'minimum_residual_constraint_margin':residual_min,'minimum_epigraph_margin':epigraph,
                           'all_residual_margins_strictly_positive':residual_min>0,'original25':gates,'all25_predicted':all_gates,
                           'fixed_old_ceiling_pass':old_ceiling,'blocks':blocks,'residual_constraints':margins,
                           'interpretation':'Unevaluated secant feasibility only; no optimizer convergence or nonlinear/native acceptance inferred.'}}


def propose(current,probes,fixed):
    p=linearize(current,probes,fixed);x0=np.r_[np.zeros(12),p['q_current']]
    zero,_,_=constraints(p,x0);require(np.min(zero)>=-1e-10,'Accepted current must be numerically feasible at zero step')
    result={'current_key':token(current['parameters']),'current_metrics':p['current_metrics'],'derivative_spans':p['spans'],
            'derivative_keys':[token(v['parameters']) for v in probes],'Jacobian_original531':p['Jacobian_original531'].tolist(),'Jacobian_components':p['Jacobian_components'].tolist(),
            'fixed_ceilings':fixed['ceilings'],'fixed_component_denominators':fixed['denominators'],'solver_calls':1,'new_native_calls':0}
    try:
        opt=minimize(lambda x:float(x[-1]),x0,jac=lambda x:np.r_[np.zeros(12),1.],method='SLSQP',bounds=list(zip(p['bounds_lower'],p['bounds_upper'])),
                     constraints={'type':'ineq','fun':lambda x:constraints(p,x)[0],'jac':lambda x:constraints(p,x)[1]},options=OPTIONS)
        require(opt.x.shape==(13,) and np.all(np.isfinite(opt.x)),'Finite solver return')
        result['solver']={'method':'SLSQP','options':OPTIONS,'success':bool(opt.success),'status':int(opt.status),'message':str(opt.message),'iterations':int(opt.nit),'function_evaluations':int(opt.nfev),'Jacobian_evaluations':int(opt.njev),'returned':opt.x.tolist(),'multipliers':np.asarray(opt.get('multipliers',[])).tolist()}
        result['trial_slots']=[certify_fraction(p,opt.x,f,bool(opt.success)) for f in FRACTIONS]
        result['status']='qualified_declared_fractions' if any(s['eligible'] for s in result['trial_slots']) else 'no_qualified_declared_fraction'
    except Exception as exc:
        result.update(status='proposal_solver_failure_without_qualified_return',error=repr(exc),trial_slots=[])
    return result


def select(current,trial_slots,fixed):
    require(len(trial_slots)==4 and [r['fraction'] for r in trial_slots]==list(FRACTIONS),'Exact declared trial ordering')
    initial=assess(current,fixed);require(initial['physical_constraints_pass'],'Current accepted constraints')
    best=initial;selected=None;tol=1e-10*max(1.,abs(initial['q']));reports=[]
    for i,slot in enumerate(trial_slots):
        if not slot['eligible']:
            require('vector' not in slot,'Ineligible fractions cannot be evaluated candidates');reports.append({'fraction':slot['fraction'],'status':'not_qualified'});continue
        require(slot.get('execution_status')=='complete' and slot.get('native_admitted') is True,'Never select from incomplete/failed native batch')
        vector=slot['vector'];require(vector['parameters']==slot['parameters'] and token(vector['parameters'])==slot['key'],'Exact certified candidate vector')
        metrics=assess(vector,fixed);reports.append({'fraction':slot['fraction'],'status':'admitted','metrics':metrics})
        if metrics['physical_constraints_pass'] and metrics['q']<best['q']-tol:best=metrics;selected=i
    return {'status':'improved_actual_candidate' if selected is not None else 'no_actual_admitted_improvement',
            'selected_trial_index':selected,'selected_key':best['key'],'selected_metrics':best,'current_metrics':initial,
            'strict_q_improvement_tolerance':tol,'trials':reports,'current_wins_ties':True,'earlier_fraction_wins_ties':True,
            'derivative_probes_are_candidates':False,'global_convergence_claim':False}
