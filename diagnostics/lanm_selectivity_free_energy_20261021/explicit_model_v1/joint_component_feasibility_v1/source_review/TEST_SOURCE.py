"""Narrow synthetic arithmetic/absence controls; no saved-vector diagnostic or LP."""
import copy
import importlib.util
import json
from pathlib import Path
import tempfile
import time

HERE=Path(__file__).resolve().parent
s=importlib.util.spec_from_file_location('saved_components',HERE/'CALCULATE.py')
m=importlib.util.module_from_spec(s);s.loader.exec_module(m)

def controls():
    tick,cpu=time.monotonic(),time.process_time(); checks={}
    group={'parents':dict.fromkeys(m.PARENTS),'pair_cycles':[]}; dec={'parents':{},'paired_Dy_minus_La':{}}
    energy={}; expected={}
    for parent in m.PARENTS:
        n=int(parent[4]);metal=parent.rsplit('_',1)[1];ion='ion_'+metal
        energy[ion]=2 if metal=='La' else 3;energy[parent+'__shell']=7.;energy[parent]=100+n+(metal=='Dy')
        C={parent:1,ion:-1,parent+'__shell':-1};P={}
        for i in range(1,n+1):
            d,w=parent+f'__dimer_{i:02d}',parent+f'__water_{i:02d}'
            coeff={d:1,ion:-1,w:-1};energy[d]=10+i;energy[w]=1.
            group['pair_cycles'].append({'parent':parent,'water_index':i,'model_terms':coeff})
            for k,v in coeff.items():P[k]=P.get(k,0)+v
        B=dict(C)
        for k,v in P.items():B[k]=B.get(k,0)-v
        dec['parents'][parent]={'model_terms':{'C':C,'P':P,'B':B},'reference_electronic_plus_D4_kcal_mol':None}
    for n in (8,9):dec['paired_Dy_minus_La']['W'+str(n)]={'reference_total_kcal_mol':dict(C=0,P=0,B=0),'reference_electronic_kcal_mol':dict(C=0,P=0,B=0),'reference_D4_kcal_mol':dict(C=0,P=0,B=0)}
    terms,refs=m.cycles(dec,{'groups':{'pair_and_collective_dipoles':group}})
    parents,paired,e,closure=m.component_values(energy,terms,refs)
    checks['N_fold_ions'] = all(terms[p]['P']['ion_'+p.rsplit('_',1)[1]]==-int(p[4]) for p in m.PARENTS)
    checks['independent_signed_components'] = e.tolist()==[-8.,8.,-9.,9.]
    checks['paired_C_cancellation'] = all(x['model_kcal_mol']['C']==0 for x in paired.values())
    checks['cycle_closure']=closure==0
    r=m.np.arange(531,dtype=float)/100
    scores=m.block_scores(r,e)
    checks['component_mean_squares']=scores['P_energy']==72.5 and scores['B_energy']==72.5
    jac=m.np.tile(m.np.arange(12,dtype=float),(531,1));jc=m.np.tile(m.np.arange(12,dtype=float),(4,1))
    b=m.blocks(r,e,jac,jc)
    checks['six_block_dimensions']=[len(x['residual']) for x in b.values()]==[413,102,12,4,2,2]
    checks['shared_extra_scaling']=m.np.array_equal(b['P_energy']['jacobian'],jc[[0,2]]/m.math.sqrt(2.))
    checks['original_four_factor']=scores['old_EF']==4*float(r[:413]@r[:413])
    def rejected(fn):
        try:fn()
        except (ValueError,KeyError):return True
        return False
    missing=dict(energy);missing.pop(next(k for k in missing if '__water_' in k))
    checks['missing_literal_term_rejected']=rejected(lambda:m.component_values(missing,terms,refs))
    bad=copy.deepcopy(dec);bad['parents'][m.PARENTS[0]]['model_terms']['P']['ion_La']=-1
    checks['changed_ion_coefficient_rejected']=rejected(lambda:m.cycles(bad,{'groups':{'pair_and_collective_dipoles':group}}))
    with tempfile.TemporaryDirectory() as d:
        p=Path(d)/'a.json';p.write_text('{"x":1}');record=m.pin(p);p.write_text('{"x":2}')
        checks['changed_hash_rejected']=rejected(lambda:m.read(record))
    flags=m.comparisons(dict(old_EF=1,pair_mu=1,C_mu=1,chi=1,P_energy=.9,B_energy=.8),dict(old_EF=1,pair_mu=1,C_mu=1,chi=1,P_energy=1,B_energy=1))
    checks['joint_descent_and_nonincrease_flags']=flags['both_component_scores_decrease'] and flags['all_original_four_nonincrease']
    m.require(all(checks.values()),str(checks))
    return {'status':'passed_narrow_saved_component_controls','implementation':m.pin(HERE/'CALCULATE.py'),'tests':m.pin(__file__),'checks':checks,'count':len(checks),'actual_diagnostic_executed':False,'native_calls':0,'cost':{'wall_seconds':time.monotonic()-tick,'CPU_seconds':time.process_time()-cpu,'threads':1}}

if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser();p.add_argument('--out',required=True);a=p.parse_args();m.save(a.out,controls());print(a.out)
