"""Four same-linearization descriptive SOCPs. No candidate evaluation API."""
import copy
import json
import math
import os
from pathlib import Path
import time
import numpy as np
from scipy import sparse
import clarabel

CASES=('A','B','C','D')
CUSHION=1e-8

def require(ok,msg):
    if not ok:raise ValueError(msg)

def normalized_problem(core,center,probes,fixed):
    p=core.linearize(center,probes,fixed)
    p['raw_center_q']=p['q_current'];p['q_current']=p['raw_center_q']+CUSHION
    p['bounds_upper'][-1]=p['q_current']
    # All descriptive norm ratios now refer to the final center; only B retains
    # pair/C/chi as constraints, while old_EF remains descriptive in every case.
    for name in core.SLICES:
        b=p['blocks'][name];b['denominator']=math.fsum(float(v*v) for v in b['r'])
        require(b['denominator']>1e-24,'Nondegenerate current block norm')
    extra={}
    for channel in ('chi_states_a0_cubed','C_dipoles_Debye','pair_dipoles_Debye'):
        names=list(center['response'][channel]);arr=np.asarray([center['response'][channel][k] for k in names],dtype=float)
        columns=[]
        for j,span in enumerate(p['spans']):
            a,b=probes[2*j:2*j+2]
            av=np.asarray([a['response'][channel][k] for k in names]);bv=np.asarray([b['response'][channel][k] for k in names])
            columns.append((bv-av)/span)
        extra[channel]={'names':names,'r':arr,'J':np.stack(columns,axis=-1)}
    p['electrical_unweighted']=extra
    p['mapping_S']=fixed['domain']['S']
    return p

def retain(name,case):
    if name=='block/old_EF':return False
    if name in ('block/pair_mu','block/C_mu','block/chi'):return case=='B'
    if name.startswith(('force/','energy/')):return case!='D'
    return True

def assemble(core,p,case):
    require(case in CASES,'Frozen four cases')
    full=core.conic.assemble(p);matrices=[];rhs=[];layout=[];cones=[];offset=0
    linear=[]
    for item in full['layout']:
        name=item['name']
        if name.startswith('bound/') or not retain(name,case):continue
        ids=slice(item['start'],item['stop']);a=full['A'][ids];b=full['b'][ids]
        if item['kind']=='SOC':
            matrices.append(a);rhs.extend(b);size=len(b);cones.append(clarabel.SecondOrderConeT(size))
            layout.append(dict(item,start=offset,stop=offset+size));offset+=size
        else:linear.append((name,a,b))
    def side(name,constant,gradient):
        linear.append((name,sparse.csc_matrix(-np.asarray(gradient,dtype=float)[None,:]),np.array([constant],dtype=float)))
    eye=np.eye(13)
    for j in range(12):
        side('trust/lower/'+str(j),.025,eye[j]);side('trust/upper/'+str(j),.025,-eye[j])
        if case!='C':
            side('box/lower/'+str(j),float(p['z'][j]),eye[j]);side('box/upper/'+str(j),float(1-p['z'][j]),-eye[j])
    side('epigraph/lower',0.,eye[12]);side('epigraph/upper',math.sqrt(p['q_current']),-eye[12])
    linear_count=0
    for name,a,b in linear:
        matrices.append(a);rhs.extend(b);layout.append({'name':name,'kind':'nonnegative','start':offset,'stop':offset+1});offset+=1;linear_count+=1
    cones.append(clarabel.NonnegativeConeT(linear_count))
    require(len(cones)-1=={'A':17,'B':20,'C':17,'D':2}[case],'Exact retained SOC families')
    model={'P':full['P'],'c':full['c'],'A':sparse.vstack(matrices,format='csc'),'b':np.asarray(rhs),'cones':cones,'layout':layout}
    require(model['A'].shape==(offset,13) and np.isfinite(model['A'].data).all() and np.isfinite(model['b']).all(),'Finite13-column problem')
    return model

def affine(r,J,d):
    r=np.asarray(r);J=np.asarray(J)
    flat=[float(v)+math.fsum(float(a)*float(b) for a,b in zip(row,d)) for v,row in zip(r.reshape(-1),J.reshape(-1,12))]
    return np.asarray(flat).reshape(r.shape)

def certificate(core,p,returned,fraction,case,solver_success):
    raw=core.certify_fraction(p,returned,fraction,solver_success)
    cert=raw['certificate'];rho=max(cert['blocks'][k]['ratio'] for k in ('P_energy','B_energy'))
    d=np.asarray(returned[:12],dtype=float)*fraction;qaux=rho+CUSHION
    constructed=core.certify_fraction(p,np.r_[d,qaux],1.,solver_success)
    require(constructed['parameters']==raw['parameters'],'Auxiliary construction leaves physical fraction exact')
    c=constructed['certificate'];retained={k:v for k,v in c['residual_constraints'].items() if retain(k,case)}
    bounds={}
    for j,x in enumerate(d):
        bounds['trust/lower/'+str(j)]=float(x+.025);bounds['trust/upper/'+str(j)]=float(.025-x)
        if case!='C':
            # Literal physical-parameter checks accompany normalized inequalities.
            bounds['box/lower/'+str(j)]=float(x+p['z'][j]);bounds['box/upper/'+str(j)]=float(1-p['z'][j]-x)
    bounds.update({'epigraph/lower':qaux,'epigraph/upper':p['q_current']-qaux})
    all_margins={**retained,**bounds};minimum=min(all_margins.values())
    in_physical_box=all(lo<=x<=hi for x,lo,hi in zip(raw['parameters'],p['lower'],p['upper']))
    current_scores={name:4*b['predicted_squared_norm'] for name,b in c['blocks'].items() if name in core.SLICES}
    signed={}
    for channel,b in p['electrical_unweighted'].items():
        vals=affine(b['r'],b['J'],d)
        signed[channel]={k:float(v) if np.ndim(v)==0 else v.tolist() for k,v in zip(b['names'],vals)}
    chi=signed['chi_states_a0_cubed'];increments={m:chi[m+'_trimer']-chi[m+'_dimer'] for m in ('La','Dy')}
    increments['Dy_minus_La']=increments['Dy']-increments['La']
    C=signed['C_dipoles_Debye'];contrast=[a-b for a,b in zip(C['aquo9_central_Dy'],C['aquo9_central_La'])]
    targets=p['electrical_targets']
    errors={'chi_a0_cubed':{k:v-targets['chi_a0_cubed'][k] for k,v in chi.items()},
        'C_Debye':{k:[a-b for a,b in zip(v,targets['C_Debye'][k])] for k,v in C.items()}}
    reference_increments={m:targets['chi_a0_cubed'][m+'_trimer']-targets['chi_a0_cubed'][m+'_dimer'] for m in ('La','Dy')}
    reference_increments['Dy_minus_La']=reference_increments['Dy']-reference_increments['La']
    mapped={}
    for metal,i in [('La',0),('Dy',4)]:
        a,b,u,v=raw['parameters'][i:i+4]
        mapped[metal]={'a':a,'b':b,'c':p['mapping_S']*(u+v)/4,'d':p['mapping_S']*(u-v)/4}
    components={}
    for name,label in [('P_energy','P'),('B_energy','B')]:
        vals=math.sqrt(2)*affine(p['blocks'][name]['r'],p['blocks'][name]['J'],d)
        for n,v in zip((8,9),vals):components[label+str(n)]=float(v)
    physical_ok=bool(c['all25_predicted']);retained_exact=minimum>=0 and (case=='C' or in_physical_box)
    return {'fraction':fraction,'parameters':raw['parameters'],'key':raw['key'],'normalized_direction':d.tolist(),
        'raw_certificate':raw,'constructed_auxiliary_certificate':constructed,'raw_q':cert['q_fraction'],
        'raw_epigraph_margin':cert['minimum_epigraph_margin'],'actual_component_q':rho,'auxiliary_q':qaux,'cushion':CUSHION,
        'physical_fraction_changed':False,'all_retained_scalar_margins':all_margins,'minimum_retained_margin':minimum,
        'retained_exact_sign_feasible':retained_exact,'retained_strict_positive':minimum>0 and (case=='C' or in_physical_box),
        'retained_numerically_feasible_1e_8':minimum>=-1e-8,'literal_original_box_pass':bool(in_physical_box),
        'all25_predicted':physical_ok,'original25':c['original25'],'component_errors_kcal_mol':components,
        'P_B_RMS_kcal_mol':{k:math.sqrt(math.fsum(components[k+str(n)]**2 for n in (8,9))/2) for k in ('P','B')},
        'original_block_scores':current_scores,'electrical_RMS':{'pair_Debye':math.sqrt(current_scores['pair_mu']),
            'C_Debye':math.sqrt(current_scores['C_mu']),'chi_a0_cubed':math.sqrt(current_scores['chi'])},
        'signed_electrical_models':signed,'signed_electrical_errors':errors,
        'chi_conditional_increments_a0_cubed':{k:{'model':v,'reference':reference_increments[k],'error':v-reference_increments[k]} for k,v in increments.items()},
        'W9_C_Dy_minus_La_Debye':{'model':contrast,'reference':[a-b for a,b in zip(targets['C_Debye']['aquo9_central_Dy'],targets['C_Debye']['aquo9_central_La'])]},
        'mapped_native_coefficients':mapped,'box_distances':{'lower':(np.asarray(raw['parameters'])-p['lower']).tolist(),'upper':(p['upper']-np.asarray(raw['parameters'])).tolist()},
        'q_gain_from_exact_final':p['raw_center_q']-rho,'model_admitted':False,'native_trial_authorized':False,
        'physical_gate_relaxation_case':case=='D','parameter_box_relaxation_case':case=='C'}

def run_cases(core,center,probes,fixed,cfg,out):
    out=Path(out);out.mkdir(parents=True,exist_ok=False);base=core.conic
    tick,cpu=time.monotonic(),time.process_time()
    result={'status':'started','solver_calls':0,'solver_constructions':0,'cases':{k:{'status':'not_attempted'} for k in CASES},
        'error':None,'native_calls':0,'candidate_evaluations':0,'selected_candidate':None,'configuration_fingerprint':cfg['configuration_fingerprint_for_solver']}
    save=lambda:base.save(out/'RESULT.json',result)
    try:
        p=normalized_problem(core,center,probes,fixed)
        p['electrical_targets']=base.read(cfg['targets'])
        base.save(out/'LINEARIZATION.json',base.plain(p));result['linearization']=base.pin(out/'LINEARIZATION.json')
        result.update(center_key=core.token(center['parameters']),derivative_keys=[core.token(v['parameters']) for v in probes],
            spans=p['spans'],raw_center_q=p['raw_center_q'],auxiliary_center_q=p['q_current'],settings=cfg['settings'],targets=cfg['targets'])
        settings=base.frozen_settings(cfg['settings']);result['effective_settings']=base.attributes(settings);save()
        for case in CASES:
            folder=out/case;folder.mkdir();row=result['cases'][case];row.update(status='preparing');save()
            model=assemble(core,p,case)
            np.savez(folder/'CONIC_DATA.npz',A_data=model['A'].data,A_indices=model['A'].indices,A_indptr=model['A'].indptr,A_shape=model['A'].shape,b=model['b'],c=model['c'])
            base.save(folder/'ROW_LAYOUT.json',model['layout']);row.update(conic_data=base.pin(folder/'CONIC_DATA.npz'),row_layout=base.pin(folder/'ROW_LAYOUT.json'))
            row['baseline_certificate']=certificate(core,p,np.r_[np.zeros(12),p['q_current']],1.,case,False)
            require(row['baseline_certificate']['retained_numerically_feasible_1e_8'],'Center numerical baseline feasible')
            with (folder/'ATTEMPT.json').open('x') as h:
                json.dump({'configuration_fingerprint':cfg['configuration_fingerprint_for_solver'],'case':case,'pid':os.getpid(),'maximum_attempts':1},h);h.flush();os.fsync(h.fileno())
            row['status']='attempt_reserved';row['attempt']=base.pin(folder/'ATTEMPT.json');result['solver_constructions']+=1;save()
            start,cp=time.monotonic(),time.process_time()
            solver=clarabel.DefaultSolver(model['P'],model['c'],model['A'],model['b'],model['cones'],settings)
            require(base.attributes(solver.get_settings())==base.attributes(settings),'Every effective setting unchanged')
            solver.print_to_buffer();result['solver_calls']+=1;row['status']='solve_entered';save()
            solution=solver.solve();info=solver.get_info();(folder/'SOLVER.log').write_text(solver.get_print_buffer())
            details=base.attributes(info);details['linsolver']=base.attributes(info.linsolver)
            row.update(status='returned',solution=base.plain(base.attributes(solution)),info=base.plain(details),
                solver_status=str(solution.status),solver_success=str(solution.status)=='Solved',
                cost={'wall_seconds':time.monotonic()-start,'CPU_seconds':time.process_time()-cp,'threads':1});save()
            row['independent_residuals']=base.plain(base.diagnostics(model,solution))
            x=np.asarray(solution.x,dtype=float)
            if x.shape==(13,) and np.isfinite(x).all() and x[-1]>=0:
                row['predictions']=[certificate(core,p,np.r_[x[:12],x[-1]**2],f,case,row['solver_success']) for f in (1.,.5)]
                z=np.asarray(solution.z)
                row['dual_by_named_constraint']={a['name']:{'vector':z[a['start']:a['stop']].tolist(),
                     'minus_b_dot_z':float(-model['b'][a['start']:a['stop']]@z[a['start']:a['stop']])} for a in model['layout']}
            else:row['predictions']=[];row['prediction_rejection']='Nonfinite/malformed/negative-t solution; no repair'
            row['status']='complete';save()
        result['status']='completed_four_fixed_affine_diagnostics'
    except Exception as exc:
        result.update(status='failed_four_case_diagnostic_no_retry',error=repr(exc));raise
    finally:
        result['cost']={'wall_seconds':time.monotonic()-tick,'CPU_seconds':time.process_time()-cpu,'threads':1};save()
    return {'result':base.pin(out/'RESULT.json'),'status':result['status'],'solver_calls':result['solver_calls']}
