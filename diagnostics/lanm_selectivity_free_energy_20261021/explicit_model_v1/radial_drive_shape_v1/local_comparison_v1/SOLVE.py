"""Exactly two matched fixed-center comparisons; no candidate evaluation API."""
import json,math,os,time
from pathlib import Path
import numpy as np
from scipy import sparse
import clarabel
import MATH as M
CASES=('original12','radial14')

def assemble(core,attr,p,case):
    M.require(case in CASES,'Two declared comparisons')
    model=attr.assemble(core,p,'B');s=p['signed_W9_z'];offset=len(model['b']);J=np.r_[s['J_Debye'],0.]
    extra=sparse.csc_matrix(np.stack([J,-J]));model['A']=sparse.vstack([model['A'],extra],format='csc')
    model['b']=np.r_[model['b'],s['absolute_ceiling_Debye']-s['r_Debye'],s['absolute_ceiling_Debye']+s['r_Debye']]
    for i,name in enumerate(('upper','lower')):model['layout'].append({'name':'signed_W9_z/'+name,'kind':'nonnegative','start':offset+i,'stop':offset+i+1})
    model['cones'][-1]=clarabel.NonnegativeConeT(sum(v['kind']=='nonnegative' for v in model['layout']))
    columns=list(range(12))+[14] if case=='original12' else list(range(15))
    model['A']=model['A'][:,columns];model['P']=model['P'][columns,:][:,columns];model['c']=model['c'][columns];model['columns']=columns
    M.require(len(model['cones'])==21 and model['A'].shape[1]==len(columns),'Same20 SOC families; only radial columns differ')
    return model

def expand(x,case):
    x=np.asarray(x,dtype=float);M.require(x.shape==((13,) if case=='original12' else (15,)),'Case dimension')
    return np.r_[x[:12],0.,0.,x[-1]] if case=='original12' else x

def certificate(core,attr,p,point,fraction,solver_success):
    report=attr.certificate(core,p,point,fraction,'B',solver_success);d=report['normalized_direction'];s=p['signed_W9_z']
    residual=s['r_Debye']+math.fsum(float(a)*float(b) for a,b in zip(s['J_Debye'],d))
    extra={'signed_W9_z/upper':s['absolute_ceiling_Debye']-residual,'signed_W9_z/lower':s['absolute_ceiling_Debye']+residual}
    report['all_retained_scalar_margins'].update(extra);minimum=min(report['all_retained_scalar_margins'].values())
    report.update(minimum_retained_margin=minimum,retained_exact_sign_feasible=minimum>=0 and report['literal_original_box_pass'],retained_strict_positive=minimum>0 and report['literal_original_box_pass'],retained_numerically_feasible_1e_8=minimum>=-1e-8,
                  signed_W9_z_residual_Debye=residual,signed_W9_z_absolute_ceiling_Debye=s['absolute_ceiling_Debye'],signed_W9_safeguard_pass=abs(residual)<=s['absolute_ceiling_Debye'],radial_coefficients={'La':report['parameters'][12],'Dy':report['parameters'][13]},model_admitted=False,native_trial_authorized=False)
    return report

def run(core,attr,center,probes,fixed,cfg,out):
    out=Path(out);out.mkdir(parents=True,exist_ok=False);base=core.conic;tick,cpu=time.monotonic(),time.process_time()
    result={'status':'started','solver_calls':0,'solver_constructions':0,'cases':{k:{'status':'not_attempted'} for k in CASES},'error':None,'native_calls':0,'candidate_evaluations':0,'selected_candidate':None,'configuration_fingerprint':cfg['fingerprint']}
    save=lambda:base.save(out/'RESULT.json',result)
    try:
        p=M.problem(core,attr,center,probes,fixed,base.read(cfg['targets']));M.require(p['spans']==cfg['inherited_actual_spans']+[.01,.01],'Exact inherited and new actual normalized spans');base.save(out/'LINEARIZATION.json',base.plain(p));result['linearization']=base.pin(out/'LINEARIZATION.json')
        result.update(center_key=core.token(center['parameters']),derivative_keys=[core.token(v['parameters']) for v in probes],spans=p['spans'],raw_center_q=p['raw_center_q'],auxiliary_center_q=p['q_current'],settings=cfg['settings'],targets=cfg['targets'],dimension_only_edits={'core':core.dimension_edits,'attribution':attr.dimension_edits})
        result['column_description']={'original531_L2':np.linalg.norm(p['Jacobian_original531'],axis=0).tolist(),'component_L2':np.linalg.norm(p['Jacobian_components'],axis=0).tolist(),'full531_singular_values':np.linalg.svd(p['Jacobian_original531'],compute_uv=False).tolist(),'rank_cutoff':None,'global_identifiability_claim':False}
        settings=base.frozen_settings(cfg['settings']);result['effective_settings']=base.attributes(settings);save()
        for case in CASES:
            folder=out/case;folder.mkdir();row=result['cases'][case];row['status']='preparing';save();model=assemble(core,attr,p,case)
            np.savez(folder/'CONIC_DATA.npz',A_data=model['A'].data,A_indices=model['A'].indices,A_indptr=model['A'].indptr,A_shape=model['A'].shape,b=model['b'],c=model['c'],columns=model['columns']);base.save(folder/'ROW_LAYOUT.json',model['layout'])
            row.update(conic_data=base.pin(folder/'CONIC_DATA.npz'),row_layout=base.pin(folder/'ROW_LAYOUT.json'),baseline_certificate=certificate(core,attr,p,np.r_[np.zeros(14),p['q_current']],1.,False))
            M.require(row['baseline_certificate']['retained_numerically_feasible_1e_8'],'Center feasible under unchanged physical/development safeguards')
            with (folder/'ATTEMPT.json').open('x') as h:json.dump({'case':case,'configuration_fingerprint':cfg['fingerprint'],'pid':os.getpid(),'maximum_attempts':1},h);h.flush();os.fsync(h.fileno())
            row['status']='attempt_reserved';row['attempt']=base.pin(folder/'ATTEMPT.json');result['solver_constructions']+=1;save();t,c=time.monotonic(),time.process_time()
            solver=clarabel.DefaultSolver(model['P'],model['c'],model['A'],model['b'],model['cones'],settings)
            M.require(base.attributes(solver.get_settings())==base.attributes(settings),'All effective settings unchanged')
            solver.print_to_buffer();result['solver_calls']+=1;row['status']='solve_entered';save();solution=solver.solve();info=solver.get_info();(folder/'SOLVER.log').write_text(solver.get_print_buffer())
            details=base.attributes(info);details['linsolver']=base.attributes(info.linsolver)
            row.update(status='returned',solution=base.plain(base.attributes(solution)),info=base.plain(details),solver_status=str(solution.status),solver_success=str(solution.status)=='Solved',cost={'wall_seconds':time.monotonic()-t,'CPU_seconds':time.process_time()-c,'threads':1});save()
            row['independent_residuals']=base.plain(base.diagnostics(model,solution));x=np.asarray(solution.x,dtype=float)
            if x.shape==(len(model['columns']),) and np.isfinite(x).all() and x[-1]>=0:
                full=expand(x,case);row['expanded_normalized_solution_t']=full.tolist();point=np.r_[full[:14],full[-1]**2]
                row['predictions']=[certificate(core,attr,p,point,f,row['solver_success']) for f in (1.,.5)]
                z=np.asarray(solution.z);row['dual_by_named_constraint']={a['name']:{'vector':z[a['start']:a['stop']].tolist(),'minus_b_dot_z':float(-model['b'][a['start']:a['stop']]@z[a['start']:a['stop']])} for a in model['layout']}
            else:row['predictions']=[];row['prediction_rejection']='Nonfinite/malformed/negative t; no fallback'
            row['status']='complete';save()
        result['status']='completed_two_matched_radial_affine_comparisons'
    except Exception as e:result.update(status='failed_radial_affine_comparison_no_retry',error=repr(e));raise
    finally:result['cost']={'wall_seconds':time.monotonic()-tick,'CPU_seconds':time.process_time()-cpu,'threads':1};save()
    return {'result':base.pin(out/'RESULT.json'),'status':result['status'],'solver_calls':result['solver_calls']}
