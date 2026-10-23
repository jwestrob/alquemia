#!/usr/bin/env python3
"""Independent saved joins/arithmetic only; no native/evaluator/optimizer imports."""
import hashlib, json, math, os, time
from pathlib import Path
import numpy as np

BASE=Path(__file__).resolve().parents[5]
D=BASE/'diagnostics/lanm_selectivity_free_energy_20261021/explicit_model_v1/joint_component_feasibility_v1'
W=BASE/'workspaces/lanm_selectivity_free_energy_20261021/explicit_model_v1/joint_component_feasibility_v1'
OUT=Path(__file__).resolve().parent
CHECKS=[]; PINS={}; MAXDIFF={}
SLICES={'old_EF':(0,413),'pair_mu':(413,515),'C_mu':(515,527),'chi':(527,531)}
SCOREMAP={'old_EF':'old_EF','pair_mu':'pair','C_mu':'C','chi':'chi'}
CHANNELS=('C','P','B'); ORDER=('P8','B8','P9','B9')
MODES={'components_only':[],'protect_EF':['old_EF'],'protect_all_original':list(SLICES)}

def check(ok,name):
    if not ok: raise ValueError(name)
    CHECKS.append(name)
def pin(p):
    p=Path(p).resolve(); return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def read(p,label):
    if isinstance(p,dict):
        check(pin(p['path'])==p,'pin '+label); PINS[label]=p; p=p['path']
    else: PINS[label]=pin(p)
    return json.loads(Path(p).read_text())
def close(a,b,name):
    a,b=np.asarray(a,float),np.asarray(b,float)
    check(a.shape==b.shape and np.isfinite(a).all() and np.isfinite(b).all(),name+' finite shape')
    err=float(np.max(np.abs(a-b),initial=0));MAXDIFF[name]=err
    check(bool(np.all(np.abs(a-b)<=1e-10*np.maximum(1,np.maximum(abs(a),abs(b))))),name+' declared arithmetic')
def sub(a,b):return {k:a[k]-b[k] for k in CHANNELS}
def scores(r,e):
    s={k:4*math.fsum(float(x*x) for x in r[a:b]) for k,(a,b) in SLICES.items()}
    s.update(P_energy=(e[0]**2+e[2]**2)/2,B_energy=(e[1]**2+e[3]**2)/2);return s
def flags(a,b):
    ni={k:a[k]<=b[k]+1e-10*max(1,abs(a[k]),abs(b[k])) for k in a}
    st={k:a[k]<b[k]-1e-10*max(1,abs(a[k]),abs(b[k])) for k in a}
    return {'strict_descent_by_block':st,'nonincrease_by_block':ni,'both_component_scores_decrease':st['P_energy'] and st['B_energy'],'old_EF_nonincrease':ni['old_EF'],'all_original_four_nonincrease':all(ni[k] for k in SLICES)}

def review():
    t,c=time.monotonic(),time.process_time()
    check(all(os.environ.get(k)=='1' for k in ('OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS')),'one numerical thread')
    result=read(W/'actual_v1/RESULT.json','diagnostic');src=read(result['inputs'],'inputs')
    check(result['status']=='complete_saved_joint_component_feasibility' and result['error'] is None,'terminal diagnostic')
    check(src==result['sources'],'exact source inventory')
    for name,p in [('implementation',result['implementation']),('protocol',src['protocol']),('LP',src['lp'])]:
        check(pin(p['path'])==p,'pin '+name);PINS[name]=p
    release=read(result['root_release'],'root release')
    check(release['implementation']==result['implementation'] and release['inputs']==result['inputs'] and release['protocol']==src['protocol'] and release['lp']==src['lp'],'executed release joins')
    data={}
    for k,p in src['artifacts'].items():
        if p['path'].endswith('.json'): data[k]=read(p,k)
        else: check(pin(p['path'])==p,'pin '+k);PINS[k]=p
    fit,dec,em,common,pair,targets=(data[k] for k in ('fit_result','decomposition_result','evaluator_manifest','common_basis','matched_pair','targets'))
    panel=read(common['panel'],'common basis panel')
    check(dec['sources']['selected']==fit['selected']==src['artifacts']['selected'],'selected fit/decomposition identity')
    for k in ('common_basis','matched_pair','targets','evaluator_manifest'):check(dec['sources'][k]==src['artifacts'][k],'decomposition source '+k)
    check(common['pair_collection']==src['artifacts']['matched_pair'],'common basis D4 authority')
    check(em['sources']['target_inventory']==src['artifacts']['targets'],'moment/energy literal target identity')
    check(data['evaluator_configuration']['manifest']==src['artifacts']['evaluator_manifest'],'executed static authority')
    check(data['root_fit_review']['independent_review']==src['artifacts']['full_fit_review'],'prior raw fit audit retained')
    check(common['new_complete']==68 and common['reused_complete']==12 and pair['reference_complete']==36,'reference 80 and original D4 denominator')
    parents={};d4={};HA=627.509474
    for parent,cycles in result['literal_cycles'].items():
        n=int(parent[4]);metal=parent.rsplit('_',1)[1]
        pairs=[x for x in targets['groups']['pair_and_collective_dipoles']['pair_cycles'] if x['parent']==parent]
        check(len(pairs)==n and {p['water_index'] for p in pairs}==set(range(1,n+1)),'literal waters '+parent)
        C={parent:1,'ion_'+metal:-1,parent+'__shell':-1};P={}
        for row in pairs:
            for key,coef in row['model_terms'].items():P[key]=P.get(key,0)+coef
        B=dict(C)
        for key,coef in P.items():B[key]=B.get(key,0)-coef
        check(cycles=={'C':C,'P':P,'B':B}==dec['parents'][parent]['model_terms'] and P['ion_'+metal]==-n,'N-fold/expanded cycles '+parent)
        check(dec['parents'][parent]['reference_electronic_plus_D4_kcal_mol'] is None,'no absolute parent target '+parent)
        check(em['rows'][parent]['positions_A']==panel['geometries'][parent]['positions_A'],'parent physical frame '+parent)
        for row in pairs:
            i=row['water_index'];rr=next(x for x in panel['cells'] if x['parent_geometry_id']==parent and x['role']=='metal_water' and x['water_index']==i)
            check(em['rows'][parent+f'__dimer_{i:02d}']['positions_A']==[em['rows'][parent]['positions_A'][j] for j in rr['real_parent_atom_indices']],'literal dimer frame '+parent+str(i))
        e=common['decomposition']['per_parent_electronic'][parent]
        parents[parent]={k:e[v]*HA for k,v in zip(CHANNELS,('cluster_interaction_hartree','pair_sum_hartree','cooperative_remainder_hartree'))}
        d4[parent]={'complex':common['reused_terms'][parent]['D4_energy_hartree'],'ion':pair['reference_rows']['pbe0__'+metal+'__bare_ion']['D4_energy_hartree'],'pairs':[pair['reference_rows'][f'pbe0__aquo{n}__{metal}__water_{i:02d}']['D4_energy_hartree'] for i in range(1,n+1)]}
    target_report={}
    for n in (8,9):
        la,dy=[f'aquo{n}_central_{m}' for m in ('La','Dy')]
        check(panel['geometries'][la]['positions_A']==panel['geometries'][dy]['positions_A'],'matched dispersion geometry '+str(n))
        el=sub(parents[dy],parents[la]);a,b=d4[dy],d4[la]
        dc=math.fsum([a['complex'],-b['complex'],-a['ion'],b['ion']])*HA
        dp=math.fsum(a['pairs']+[-x for x in b['pairs']]+[-n*a['ion'],n*b['ion']])*HA
        disp={'C':dc,'P':dp,'B':dc-dp};total={k:el[k]+disp[k] for k in CHANNELS}
        for k in CHANNELS:
            close(total[k],result['reference_targets'][str(n)][k],'independent paired target '+str(n)+k)
            close(el[k],dec['paired_Dy_minus_La']['W'+str(n)]['reference_electronic_kcal_mol'][k],'electronic join '+str(n)+k)
            close(disp[k],dec['paired_Dy_minus_La']['W'+str(n)]['reference_D4_kcal_mol'][k],'D4 once '+str(n)+k)
        target_report[str(n)]={'electronic':el,'D4':disp,'total':total}
    selected=data['selected'];energies={k:x['energy_kcal_mol'] for k,x in em['static'].items()};energies.update({k:x['energy'] for k,x in selected['cells'].items()})
    sv=result['vectors'][result['selected_key']];selected_replay={}
    for parent,cycles in result['literal_cycles'].items():
        v={k:math.fsum(coef*energies[key] for key,coef in terms.items()) for k,terms in cycles.items()};selected_replay[parent]=v
        for k in CHANNELS:
            close(v[k],sv['parents_kcal_mol'][parent][k],'selected literal '+parent+k)
            close(v[k],dec['parents'][parent]['selected_kcal_mol'][k],'admitted decomposition '+parent+k)
    selected_errors={}
    for n in (8,9):
        v=sub(selected_replay[f'aquo{n}_central_Dy'],selected_replay[f'aquo{n}_central_La']);err=sub(v,target_report[str(n)]['total']);selected_errors[str(n)]=err
        for k in CHANNELS:close(err[k],sv['paired']['W'+str(n)]['error_kcal_mol'][k],'selected paired '+str(n)+k)
    check(set(result['vectors'])==set(fit['vectors']) and len(result['vectors'])==673,'all 673 stored vector identities')
    for key,v in result['vectors'].items():
        old=fit['vectors'][key]
        check(v['artifact']==old['artifact'] and v['parameters']==old['parameters'] and v['role']==old['role'],'stored vector join '+key)
        check(len(v['original25'])==25 and sum(k.startswith('energy/') for k in v['original25'])==10,'25 denominator '+key)
        check(v['original25_pass']==old['original_EF_checks_passed']==all(g['passed'] for g in v['original25'].values()),'retained gates '+key)
        e=[v['component_errors_kcal_mol'][k] for k in ORDER]
        close([v['scores']['P_energy'],v['scores']['B_energy']],[(e[0]**2+e[2]**2)/2,(e[1]**2+e[3]**2)/2],'component scores '+key)
        for k,orig in SCOREMAP.items():close(v['scores'][k],old['block_objectives'][orig],'four old scores '+key+k)
    initial=read(fit['vectors'][fit['initial_key']]['artifact'],'initial artifact')
    residual={fit['initial_key']:np.array(initial['score']['weighted_residuals'])}
    for tr in result['trials']:residual[tr['trial_key']]=np.array(tr['actual_original_residuals531'])
    check(len(result['origins'])==24 and len(result['trials'])==96,'24 origins 96 trials')
    width=np.array([p['upper']-p['lower'] for p in data['domain']['parameters']]);low=np.array([p['lower'] for p in data['domain']['parameters']])
    origin_report=[];lp_stats={k:[] for k in MODES};trial_summary={'both_components':0,'also_four_blocks':0,'also_25_gates':0,'also_origin_25_gates':0};trial_rows=[]
    for o,it in zip(result['origins'],fit['iterations']):
        key=o['key'];v=result['vectors'][key];r=residual[key];J=np.array(o['Jacobian_original531']);JC=np.array(o['Jacobian_components']);e=np.array([v['component_errors_kcal_mol'][k] for k in ORDER]);z=(np.array(v['parameters'])-low)/width
        check(key==it['origin_key'] and o['outer_step']==it['outer_step'] and o['derivative_keys']==it['derivative_keys'] and o['derivative_spans']==it['derivative_spans'],'exact origin and probes '+str(o['outer_step']))
        check(o['Jacobian_original531']==it['Jacobian_normalized'],'existing Jacobian exact '+str(o['outer_step']))
        jc=np.column_stack([(np.array([result['vectors'][b]['component_errors_kcal_mol'][k] for k in ORDER])-np.array([result['vectors'][a]['component_errors_kcal_mol'][k] for k in ORDER]))/span for a,b,span in zip(it['derivative_keys'][::2],it['derivative_keys'][1::2],it['derivative_spans'])])
        close(jc,JC,'component derivative '+str(o['outer_step']))
        blocks={k:(r[a:b],J[a:b]) for k,(a,b) in SLICES.items()};blocks.update(P_energy=(e[[0,2]]/math.sqrt(2),JC[[0,2]]/math.sqrt(2)),B_energy=(e[[1,3]]/math.sqrt(2),JC[[1,3]]/math.sqrt(2)))
        g={k:rv@jj/(rv@rv) for k,(rv,jj) in blocks.items()}
        check(o['original25_pass']==v['original25_pass'],'origin failed gates retained '+str(o['outer_step']))
        check(o['LP']['status']=='complete_tangent_screen' and set(o['LP']['modes'])==set(MODES),'three LP modes '+str(o['outer_step']))
        for mode,m in o['LP']['modes'].items():
            dx=np.array(m['normalized_direction']);lo=np.maximum(-.025,-z);hi=np.minimum(.025,1-z);x=np.r_[dx,m['t']]
            names=['P_energy','B_energy']+MODES[mode];A=np.array([np.r_[g[k],float(k in ('P_energy','B_energy'))] for k in names]);y=np.array(m['inequality_marginals']);lm=np.array(m['lower_marginals']);um=np.array(m['upper_marginals']);obj=np.r_[np.zeros(12),-1.]
            check(m['protected_blocks']==MODES[mode] and m['constraint_names']==names,'LP protection '+str(o['outer_step'])+mode)
            cert={'primal_violation':max(0,float(np.max(A@x)),float(np.max(lo-dx)),float(np.max(dx-hi)),-m['t']),'primal_minus_dual_gap':-m['t']-float(lo@lm[:-1]+hi@um[:-1]),'stationarity_max':float(np.max(abs(obj-A.T@y-lm-um))),'dual_sign_violation':max(0,float(np.max(y)),float(np.max(-lm)),float(np.max(um)),abs(float(um[-1]))),'complementarity_max':max(float(np.max(abs(y*(-A@x)))),float(np.max(abs(lm[:-1]*(dx-lo)))),float(np.max(abs(um[:-1]*(hi-dx)))),abs(float(lm[-1]*m['t'])))}
            check(all(abs(value)<=1e-7 for value in cert.values()),'independent saved LP certificate '+str(o['outer_step'])+mode)
            ratios={}
            for k,(rv,jj) in blocks.items():
                pred=m['block_predictions'][k];vals=[np.linalg.norm(rv),g[k]@dx,np.linalg.norm(rv+.5*(jj@dx)),np.linalg.norm(rv+jj@dx)]
                close(vals,[pred[q] for q in ('origin_norm','fractional_norm_directional_derivative','half_step_linearized_norm','full_step_linearized_norm')],'saved LP prediction '+str(o['outer_step'])+mode+k)
                ratios[k]=float(vals[-1]/vals[0])
            lp_stats[mode].append({'step':o['outer_step'],'t':m['t'],'full_norm_ratios':ratios,'certificate_max':max(abs(x) for x in cert.values())})
        origin_report.append({'step':o['outer_step'],'key':key,'all25':o['original25_pass'],'failed_gates':[k for k,g0 in v['original25'].items() if not g0['passed']]})
    origins={o['outer_step']:o for o in result['origins']}
    for tr in result['trials']:
        o=origins[tr['outer_step']];it=fit['iterations'][tr['outer_step']-1];v=result['vectors'][tr['trial_key']];orig=result['vectors'][o['key']]
        check(tr['trial_key']==it['trial_keys'][tr['trial_index']] and tr['origin_key']==o['key'],'trial origin identity '+str(len(trial_rows)))
        computed=flags(v['scores'],orig['scores']);check(computed==tr['flags'],'trial classification '+str(len(trial_rows)))
        dz=(np.array(v['parameters'])-np.array(orig['parameters']))/width;pr=residual[o['key']]+np.array(o['Jacobian_original531'])@dz;pe=np.array([orig['component_errors_kcal_mol'][k] for k in ORDER])+np.array(o['Jacobian_components'])@dz
        close(pr,tr['predicted_original_residuals531'],'trial predicted old '+str(len(trial_rows)));close(pe,tr['predicted_component_errors_kcal_mol'],'trial predicted components '+str(len(trial_rows)))
        check(tr['trial_original25_pass']==v['original25_pass'] and tr['origin_original25_pass']==orig['original25_pass'],'trial gate retention '+str(len(trial_rows)))
        both=computed['both_component_scores_decrease'];four=both and computed['all_original_four_nonincrease'];passes=four and v['original25_pass'];oldpasses=passes and orig['original25_pass']
        for k,b in zip(trial_summary,(both,four,passes,oldpasses)):trial_summary[k]+=int(b)
        trial_rows.append({'step':tr['outer_step'],'index':tr['trial_index'],'both':both,'four':four,'trial25':v['original25_pass'],'origin25':orig['original25_pass']})
    check(result['counts']=={'vectors':673,'origins':24,'saved_trials':96,'expected_nondegenerate_LPs':72,'LPs':72},'reported denominator')
    check(result['last_Jacobian_origin']==fit['last_Jacobian_origin_key'] and result['selected_key']==fit['selected_key'] and result['last_origin_is_selected'] is False,'no selected point Jacobian substitution')
    check(result['cost']['native_programs']==result['new_native_calls']==result['new_parameter_evaluations']==0,'zero new native/evaluation count')
    out={'status':'passed_saved_joint_component_feasibility_actual_review','result':PINS['diagnostic'],'implementation':pin(__file__),'sources':PINS,'checks_passed':len(CHECKS),'counts':result['counts'],'reference_targets_kcal_mol':target_report,'selected_errors_kcal_mol':selected_errors,'selected_coordination_error_kcal_mol':sub(selected_errors['9'],selected_errors['8']),'origins':origin_report,'trial_counts':trial_summary,'trial_classification':trial_rows,'LP_saved_certificate_and_prediction_review':lp_stats,'original_diagnostic_cost':result['cost'],'scope':{'direct_authority_hashes_checked':True,'all673_fit_index_and_score_joins_checked':True,'selected_literal_saved_energy_replay':True,'all24_component_derivatives_and_saved_LP_certificates_checked':True,'all96_trial_classifications_checked':True,'native_readmission':False,'LPs_resolved':False,'all673_original_vector_bodies_reloaded':False,'no_candidate_selected':True},'limits':['Consumed development results only; failed origins retained.','Paired common-basis electronic plus original D4 targets only; no absolute-parent total target or raw-cluster CC substitution.','Twenty-five E/F gates do not bound pair/collective errors. Four original blocks remain unmodified.','LPs describe tangent constraints; finite full-step residuals can grow. No LP direction was evaluated.','No selected-point Jacobian, global family-capacity conclusion, physical accuracy, regression recovery or fresh-validation release.'],'cost':{'wall_seconds':time.monotonic()-t,'CPU_seconds':time.process_time()-c,'numerical_threads':1,'native_calls':0,'LP_solves':0,'new_parameter_evaluations':0}}
    with (OUT/'REVIEW.json').open('x') as f:json.dump(out,f,indent=2,allow_nan=False);f.write('\n')
    print(json.dumps({k:out[k] for k in ('status','checks_passed','counts','trial_counts','cost')},indent=2))

if __name__=='__main__':review()
