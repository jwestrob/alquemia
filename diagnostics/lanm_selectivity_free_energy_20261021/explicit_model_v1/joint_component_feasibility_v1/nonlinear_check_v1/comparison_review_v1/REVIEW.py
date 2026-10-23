#!/usr/bin/env python3
"""Saved component/prediction/scorer audit; no native parser or engine calls."""
import hashlib, importlib.util, json, math, os, time
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[6]
D=ROOT/'diagnostics/lanm_selectivity_free_energy_20261021/explicit_model_v1/joint_component_feasibility_v1/nonlinear_check_v1'
W=ROOT/'workspaces/lanm_selectivity_free_energy_20261021/explicit_model_v1/joint_component_feasibility_v1/nonlinear_check_v1'
PINS={};CHECKS=[];MAX_ERROR=0.
def check(x,s):
    if not x:raise ValueError(s)
    CHECKS.append(s)
def pin(p):
    p=Path(p).resolve();return dict(path=str(p),sha256=hashlib.sha256(p.read_bytes()).hexdigest())
def read(p,name):
    if isinstance(p,dict):check(pin(p['path'])==p,'pin '+name);PINS[name]=p;p=p['path']
    else:PINS[name]=pin(p)
    return json.loads(Path(p).read_text())
def module(p,name):
    check(pin(p['path'])==p,'source '+name);PINS[name]=p
    s=importlib.util.spec_from_file_location(name,p['path']);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
def close(a,b,s):
    global MAX_ERROR
    if isinstance(a,dict):
        check(set(a)==set(b),s+' keys')
        for k in a:close(a[k],b[k],s+'/'+k)
    elif isinstance(a,(tuple,list)):
        check(len(a)==len(b),s+' length')
        for i,(x,y) in enumerate(zip(a,b)):close(x,y,s+'/'+str(i))
    elif isinstance(a,(int,float)) and not isinstance(a,bool):
        check(isinstance(b,(int,float)) and not isinstance(b,bool) and math.isfinite(a) and math.isfinite(b),s+' finite')
        error=abs(a-b);MAX_ERROR=max(MAX_ERROR,error)
        check(error<=1e-10*max(1,abs(a),abs(b)),s+' existing arithmetic')
    else:check(a==b,s+' exact')

def main():
    t,c=time.monotonic(),time.process_time()
    check(all(os.environ.get(x)=='1' for x in ('OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS')),'one thread')
    n=read(W/'run_v1/execution/RESULT.json','native result');s=read(W/'comparison_v1/RESULT.json','comparison')
    cfg=read(n['configuration'],'run config');decl=read(s['declaration'],'declaration');cert=read(s['half_review'],'half certificate');failed=read(s['failed_optimizer'],'failed proposal');saved=read(s['saved_diagnostic'],'saved diagnostic')
    parent=read(cfg['parent_configuration'],'parent configuration');manifest=read(parent['manifest'],'manifest');target=read(parent['target_inventory'],'target inventory');dec=read(cfg['decomposition'],'decomposition');origin=read(cfg['origin'],'origin')
    check(n['status']=='complete_one_fixed_half_native_pilot' and n['error'] is None,'terminal native result claim')
    check(s['status']=='complete_consumed_native_prediction_comparison' and s['native_result']==PINS['native result'],'summary exact native join')
    for name,p in [('summarizer',s['implementation']),('runner',n['implementation'])]:check(pin(p['path'])==p,'source pin '+name);PINS[name]=p
    check(n['implementation']==cfg['implementation'] and cfg['declaration']==PINS['declaration'] and cfg['component_diagnostic']==PINS['saved diagnostic'],'implementation/declaration/diagnostic binding')
    check(n['candidate']==cfg['candidate'] and cfg['candidate']['parameters']==cert['physical_half_parameters']==decl['physical_parameters'],'sole physical half vector')
    check(cfg['candidate']['key']==cert['half_parameter_token']==decl['parameter_token'] and cfg['candidate']['fraction']==.5,'exact half token')
    check(cfg['origin']==decl['origin_artifact']==cert['origin_artifact'] and decl['normalized_direction']==cert['normalized_half_direction'],'half origin/direction')
    check(failed['status']=='failed_numerical_local_proposal' and failed['solver']['success'] is False and s['original_optimizer_passed'] is False,'failed optimizer preserved')
    check(dec['sources']['evaluator_manifest']==parent['manifest'] and dec['sources']['targets']==parent['target_inventory'],'fixed component/static authority')
    check(saved['component_order']==['P8','B8','P9','B9'],'fixed component ordering')
    group=target['groups']['pair_and_collective_dipoles'];check(len(group['pair_cycles'])==34 and len(group['parents'])==4,'literal34/fourparent denominator')
    terms={}
    for p in group['parents']:
        N=int(p[4]);metal=p.rsplit('_',1)[1];C={p:1,'ion_'+metal:-1,p+'__shell':-1};P={}
        pairs=[r for r in group['pair_cycles'] if r['parent']==p]
        check(len(pairs)==N and {r['water_index'] for r in pairs}==set(range(1,N+1)),'water inventory '+p)
        for row in pairs:
            for k,a in row['model_terms'].items():P[k]=P.get(k,0)+a
        B=dict(C)
        for k,a in P.items():B[k]=B.get(k,0)-a
        terms[p]={'C':C,'P':P,'B':B}
        check(P['ion_'+metal]==-N and terms[p]==dec['parents'][p]['model_terms']==saved['literal_cycles'][p],'N-fold literal coefficients '+p)
    comp={}
    for name,answer in [('origin',origin),('candidate',n['answer'])]:
        check(set(answer['cells'])==set(manifest['rows']) and len(answer['cells'])==46,'complete saved46 '+name)
        energies={k:v['energy_kcal_mol'] for k,v in manifest['static'].items()};energies.update({k:v['energy'] for k,v in answer['cells'].items()})
        pvals={p:{k:math.fsum(a*energies[x] for x,a in tt.items()) for k,tt in channels.items()} for p,channels in terms.items()}
        paired={}
        for size in (8,9):
            ref=dec['paired_Dy_minus_La']['W'+str(size)]['reference_total_kcal_mol'];check(ref==saved['reference_targets'][str(size)],'unchanged paired target '+str(size)+name)
            model={k:pvals[f'aquo{size}_central_Dy'][k]-pvals[f'aquo{size}_central_La'][k] for k in ('C','P','B')}
            paired['W'+str(size)]={'model':model,'reference':ref,'error':{k:model[k]-ref[k] for k in model}}
        comp[name]={'parents':pvals,'paired':paired,'J_P':math.fsum(x['error']['P']**2 for x in paired.values())/2,'J_B':math.fsum(x['error']['B']**2 for x in paired.values())/2}
        close(comp[name],n['components'][name],'literal component replay '+name)
    ans=n['answer'];mu={k:v['mu_Debye'] for k,v in manifest['static'].items()};mu.update({k:v['mu_Debye'] for k,v in ans['cells'].items()})
    response={'pair_dipoles_Debye':{r['id']:[math.fsum(a*mu[k][j] for k,a in r['model_terms'].items()) for j in range(3)] for r in group['pair_cycles']},'C_dipoles_Debye':{p:[mu[p][j]-mu['ion_'+p.rsplit('_',1)[1]][j]-mu[p+'__shell'][j] for j in range(3)] for p in group['parents']},'chi_states_a0_cubed':{k:ans['cells'][v['fixture']]['chi_a0_cubed'] for k,v in manifest['field_aliases'].items()}}
    close(response,ans['response'],'response cycle assembly')
    check(ans['model']=={k:ans['cells'][k]['model'] for k in ans['model']} and len(ans['model'])==10,'ten E/F models exact cell join')
    objective=module(parent['objective'],'frozen_objective');oldobj=module(parent['old_objective'],'frozen_old_objective');body=module(parent['body_helper'],'frozen_body_helper')
    data=read(parent['development_data'],'development data');targets=read(parent['targets'],'objective targets')
    score=objective.evaluate(data,ans['model'],response,targets,oldobj,body);close(score,ans['score'],'frozen objective arithmetic')
    o=next(x for x in saved['origins'] if x['key']==decl['origin']);dx=np.array(cert['normalized_half_direction']);J=np.array(o['Jacobian_original531']);JC=np.array(o['Jacobian_components'])
    check(J.shape==(531,12) and JC.shape==(4,12),'saved derivative dimensions')
    pr=np.array(origin['score']['weighted_residuals'])+J@dx;errors=[comp['origin']['paired']['W'+k[1]]['error'][k[0]] for k in saved['component_order']];pe=np.array(errors)+JC@dx
    channels={}
    for i,k in enumerate(saved['component_order']):
        actual=comp['candidate']['paired']['W'+k[1]]['error'][k[0]]
        channels[k]={'origin_error_kcal_mol':errors[i],'predicted_error_kcal_mol':float(pe[i]),'actual_error_kcal_mol':actual,'actual_minus_predicted_kcal_mol':actual-float(pe[i])}
    close(channels,s['signed_component_errors'],'summary signed errors/prediction')
    blocks={}
    for k,sk,a,b in [('old_EF','old_EF',0,413),('pair_mu','pair',413,515),('C_mu','C',515,527),('chi','chi',527,531)]:
        old=origin['score']['block_objectives'][sk];new=score['block_objectives'][sk];pred=4*math.fsum(float(x*x) for x in pr[a:b])
        close(pred,failed['half']['blocks'][k]['predicted_score'],'independent old-block prediction '+k)
        blocks[k]={'origin_squared_norm_unweighted':old,'predicted_squared_norm_unweighted':pred,'actual_squared_norm_unweighted':new,'actual_nonincrease_exact':new<=old,'change':new-old}
        # Flags use exact stored values; replay arithmetic cannot change classification.
        check((ans['score']['block_objectives'][sk]<=old)==s['blocks'][k]['actual_nonincrease_exact'],'exact unrelaxed block flag '+k)
    for k,term,indices in [('P_energy','P',[0,2]),('B_energy','B',[1,3])]:
        old,new=comp['origin']['J_'+term],comp['candidate']['J_'+term];pred=math.sqrt(math.fsum(float(pe[j]**2) for j in indices)/2)
        close(pred,failed['half']['blocks'][k]['predicted_norm'],'independent component RMS prediction '+k)
        blocks[k]={'origin_RMS_kcal_mol':math.sqrt(old),'predicted_RMS_kcal_mol':pred,'actual_RMS_kcal_mol':math.sqrt(new),'actual_squared_norm':new,'actual_strict_improvement':new<old,'fractional_RMS_improvement':1-math.sqrt(new/old)}
    close(blocks,s['blocks'],'six summary blocks')
    gates={}
    for g,rows in score['old_EF']['energy_controls'].items():
        for k,row in rows.items():check(row['passed']==(abs(row['error_kcal_mol'])<=1),'unchanged energy gate '+g+k);gates['energy/'+g+'/'+k]=row
    for g,rows in score['old_EF']['per_fixture'].items():
        for k,row in rows.items():check(row['passed']==(row['rigid_RMS_kcal_mol_A']<=2 and row['maximum_body_component_error_kcal_mol_A']<=5),'unchanged force gate '+g+k);gates['force/'+g+'/'+k]={x:row[x] for x in ('rigid_RMS_kcal_mol_A','maximum_body_component_error_kcal_mol_A','passed')}
    check(len(gates)==25 and sum(k.startswith('energy/') for k in gates)==10,'10 energy15 force gates');close(gates,s['all25'],'summary25 gates')
    flags={'all25_pass':all(v['passed'] for v in gates.values()),'both_components_improve':all(blocks[k]['actual_strict_improvement'] for k in ('P_energy','B_energy')),'all_old_blocks_nonincrease_exact':all(s['blocks'][k]['actual_nonincrease_exact'] for k in ('old_EF','pair_mu','C_mu','chi'))}
    flags['declared_joint_native_outcome_passes']=all(flags.values())
    check(all(s[k]==v for k,v in flags.items()),'independent joint flags')
    check(s['native_cost']==n['cost'] and n['native_programs']==56,'reported cost/56 count join only')
    check(s['physical_model_qualified'] is False and s['fresh_validation'] is False,'physical limitations retained')
    out={'status':'passed_fixed_half_component_comparison_review','implementation':pin(__file__),'native_result':PINS['native result'],'comparison_result':PINS['comparison'],'sources':PINS,'checks_passed':len(CHECKS),'maximum_saved_arithmetic_difference':MAX_ERROR,'counts':{'pairs':34,'parents':4,'variable_cells':46,'static_terms':len(manifest['static']),'objective_residuals':len(score['weighted_residuals']),'old_energy_gates':10,'old_force_gates':15},'flags':flags,'blocks':blocks,'signed_component_errors':channels,'maximum_absolute_component_prediction_error_kcal_mol':max(abs(v['actual_minus_predicted_kcal_mol']) for v in channels.values()),'native_reported_cost':n['cost'],'scope':{'saved_literal_energy_and_response_cycles':True,'frozen_original_scorer_arithmetic':True,'native_raw_inputs_receipts_binaries_admission_reviewed_here':False,'native_admission_review_owner':'Planck','native_calls':0,'new_parameter_evaluations':0},'limitations':['Only this fixed consumed half point; failed SLSQP remains failed.','Pair/collective RMS errors remain25.59/21.39kcal/mol; no interaction-model qualification.','No new physical error tolerance, reference method qualification, trimer regression recovery or fresh validation.'],'cost':{'wall_seconds':time.monotonic()-t,'CPU_seconds':time.process_time()-c,'threads':1,'native_calls':0}}
    dest=W/'comparison_review_v1';dest.mkdir(exist_ok=False)
    with (dest/'REVIEW.json').open('x') as f:json.dump(out,f,indent=2,allow_nan=False);f.write('\n')
    print(json.dumps({k:out[k] for k in ('status','checks_passed','maximum_saved_arithmetic_difference','flags','maximum_absolute_component_prediction_error_kcal_mol','cost')},indent=2));print(json.dumps(pin(dest/'REVIEW.json')))

if __name__=='__main__':main()
