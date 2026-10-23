#!/usr/bin/env python3
"""Focused saved fixtures, rendering and interface controls; no native execution."""
import argparse
import copy
from decimal import Decimal, localcontext
import importlib.util
import json
import math
import os
from pathlib import Path
import re
import time
import numpy as np

HERE=Path(__file__).resolve().parent


def main():
    p=argparse.ArgumentParser();p.add_argument('--preparation',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    tick,cpu=time.monotonic(),time.process_time();a.out.mkdir(parents=True,exist_ok=False)
    spec=importlib.util.spec_from_file_location('joint_test_evaluator',HERE/'EVALUATE.py');e=importlib.util.module_from_spec(spec);spec.loader.exec_module(e)
    cfg=json.loads((a.preparation/'CONFIG_EVALUATOR_TEMPLATE.json').read_text());prep=json.loads((a.preparation/'RESULT.json').read_text())
    e.require(prep['status']=='prepared_joint_evaluator_zero_calls' and prep['new_native_calls']==0,'Actual saved preparation')
    initial=e.read(prep['initial']);ev=e.Evaluator(cfg,a.out/'test_cache');m=ev.m;checks={}
    def check(name,ok):e.require(bool(ok),'Control failed: '+name);checks[name]=True
    def reject(name,fn):
        try:fn()
        except (ValueError,AssertionError,KeyError):checks[name]=True;return
        raise AssertionError('Invalid case accepted: '+name)
    def forbidden(*args,**kwargs):raise AssertionError('Native invocation forbidden in saved controls')
    ev.h['native_helper'].invoke=forbidden
    check('exact46_10_4_48_denominator',len(m['rows'])==46 and len(m['static'])==48 and sum(x['needs_gradient'] for x in m['rows'].values())==10 and sum(x['needs_export'] for x in m['rows'].values())==4)
    check('23_rows_per_metal',all(sum(r['metal']==metal for r in m['rows'].values())==23 for metal in ('La','Dy')))
    check('literal_selected10_model_and_producer',initial['model']==ev.data['baseline_model'])
    check('unchanged_old_EF_score',initial['score']['old_EF']==ev.data['baseline_score'])
    check('531_residuals_four_blocks',len(initial['score']['weighted_residuals'])==531 and set(initial['score']['block_objectives'])=={'old_EF','pair','C','chi'})
    check('exact_frozen_reference_targets',ev.targets==ev.h['objective'].targets_from_inventory(ev.target_inventory))
    check('exact_initial_domain',ev.initial_vector==[x['initial'] for x in ev.domain['parameters']])
    check('baseline46_all_saved',set(initial['cells'])==set(m['rows']) and all(x['inherited_native'] for x in initial['cells'].values()) and prep['snapshot']['native_program_invocations']==0)
    # Inspect every prepared key: only selected VDWPR and declared reporting/zero
    # response lines differ; every other physical baseline line remains literal.
    for cid,row in m['rows'].items():
        old=e.verify(row['source_input']['key']).read_text().splitlines();new=e.verify(row['prepared_initial_input']['key']).read_text().splitlines()
        def keep(lines):
            kept=[]
            for s in lines:
                t=s.split('#',1)[0].split()
                if not t:continue
                if t[0].lower() in ('polar-print','save-uinduce','response-matrix-export','polar-local-response'):continue
                if t[0].lower()=='vdwpr' and sorted(map(int,t[1:3]))==[34,500 if row['metal']=='La' else 501]:continue
                kept.append(s)
            return kept
        check('unchanged_physical_key_'+cid,keep(old)==keep(new))
    v=ev.initial_vector[:];v[0]=.5;v[2]=1.;v[3]=-1.;v[8]=ev.domain['parameters'][8]['upper']
    la='aquo9_central_La__dimer_01';dy='aquo9_central_Dy__dimer_01'
    key,meta=ev.render(v,la);physical=e.effective(v,'La',ev.domain)
    check('physical_c_d_rotation',physical[:4]==[.5,0.,0.,ev.domain['S']/2])
    check('changed_La_pair_fixed_POLPAIR',meta['resolved_pair']['thole']==ev.domain['fixed_native_POLPAIR']['La'] and meta['resolved_pair']['Rmin_pair_A']==v[8])
    check('unchanged_Dy_identity',ev.cache_identity(v,dy)[:2]==ev.cache_identity(ev.initial_vector,dy)[:2])
    check('changed_La_identity',ev.cache_identity(v,la)[0]!=ev.cache_identity(ev.initial_vector,la)[0])
    cached=len(ev.rendered);same=ev.render(v,'aquo9_central_La__dimer_02',exports=True)
    check('literal_key_memoization_preserves_bytes',same==(key,meta) and len(ev.rendered)==cached)
    same[1]['POLPAIR_unchanged']=False
    check('memoized_render_cannot_be_mutated_by_caller',ev.render(v,la)[1]['POLPAIR_unchanged'] is True)
    reused=ev._cell(v,dy);check('unchanged_metal_true_producer',reused['inherited_native'] and reused['producer_result']==cfg['baseline_result'] and reused['producer_cell_id']==dy)
    check('testgrad_disables_export','RESPONSE-MATRIX-EXPORT' not in ev.render(v,la,exports=False)[0])
    for name,bad in [('short',v[:-1]),('NaN',[float('nan')]+v[1:]),('boolean',[True]+v[1:]),('out_of_box',[1.0001]+v[1:])]:reject('reject_'+name,lambda bad=bad:e.checked(bad,ev.domain))
    reject('unbound_native_gate',ev._gate)
    saved=ev.ledger['native_program_invocations'];ev.ledger['native_program_invocations']=e.MAX_PROGRAMS;ev._gate=lambda:None
    reject('program_cap_before_invoke',lambda:ev._invoke('analyze',a.out,{'receipts':[]}));ev.ledger['native_program_invocations']=saved
    old=ev.ledger['vector_keys'];ev.ledger['vector_keys']=['reserved'+str(i) for i in range(e.MAX_VECTORS)]
    reject('vector_cap_before_cells',lambda:ev.evaluate(v));ev.ledger['vector_keys']=old
    # Synthetic interface trace exercises the full-precision conjugate mean,
    # missing/duplicate/nonfinite slots, and native spectral rejection.
    row=m['rows'][la];n=row['atoms'];tokens=[];u=np.full((n,3),.01);w=np.full((n,3),.02)
    for solve in range(2):
        tokens.append('LOCAL_RESPONSE_PARAMETERS 500 36 '+' '.join(map(str,physical[:4])))
        for oi,t in enumerate(row['types'],1):
            if t==36:
                r=np.linalg.norm(np.array(row['positions_A'][0])-row['positions_A'][oi-1]);kappa=row['fixed_POLPAIR']/math.sqrt(row['alpha_A3'][0]*row['alpha_A3'][oi-1]);tokens.append(f'LOCAL_RESPONSE_PAIR 1 {oi} {row["fixed_POLPAIR"]} {kappa} {r}')
        tokens.extend(['LOCAL_RESPONSE_ENVIRONMENT 1 .5 0 0 0','LOCAL_RESPONSE_SPECTRUM .5 2 0 0 0'])
        tokens += ['LOCAL_RESPONSE_INDUCED '+str(i+1)+' '+' '.join(map(str,list(u[i])+list(w[i]))) for i in range(n)]
        tokens.append('LOCAL_RESPONSE_DUAL 0 .01')
    text='\n'.join(tokens);q=e.local_readback(text,row,physical,m['units'])
    check('conjugate_mean_uses_both_channels',np.allclose(q['conjugate_induced_Debye'],.015*m['units']['debye_per_eA'],rtol=0,atol=1e-15))
    check('zero_native_bypass',e.local_readback('',row,[0.]*4,m['units']) is None)
    reject('missing_native_slot',lambda:e.local_readback('\n'.join(tokens[:-1]),row,physical,m['units']))
    reject('duplicate_native_slot',lambda:e.local_readback(text+'\n'+tokens[-1],row,physical,m['units']))
    reject('nonfinite_native_readback',lambda:e.local_readback(text.replace('SPECTRUM .5','SPECTRUM NaN'),row,physical,m['units']))
    reject('unstable_native_spectrum',lambda:e.local_readback(text.replace('SPECTRUM .5','SPECTRUM -1'),row,physical,m['units']))
    reject('wrong_native_physical_coefficient',lambda:e.local_readback(text.replace('PARAMETERS 500 36 0.5','PARAMETERS 500 36 0.6'),row,physical,m['units']))
    reject('changed_native_damping_length',lambda:e.local_readback(text.replace(str(kappa),str(kappa*1.01)),row,physical,m['units']))
    # The reference estimator is unchanged. Independently reconstruct old model
    # centered energy curvature from pinned native tokens at h=.001.
    inv=e.read(cfg['fixture_inventory']);fm=e.read(inv['artifacts']['field_manifest']);fr=e.read(inv['artifacts']['field_model'])
    comparison=e.read(ev.target_inventory['artifacts']['field_comparison']);proof={};H=m['units']['hartree_kcal_mol'];h=.001
    for state in m['field_aliases']:
        moments=[fm['rows'][state]['baseline'],fr['rows'][state+'__plus_h']['moments'],fr['rows'][state+'__minus_h']['moments']]
        energies=[];bounds=[];pins=[]
        for moment in moments:
            path=e.verify(moment['analyze_output']);parsed=ev.h['moment_adapter'].parse_energy(path)
            check('field_raw_energy_'+state+'_'+str(len(energies)),parsed['energy_kcal_mol']==moment['energy_kcal_mol'])
            token=re.findall(r'Total Potential Energy\s*:\s*(\S+)',path.read_text());e.require(len(token)==1,'One energy token')
            energies.append(Decimal(token[0]));bounds.append(Decimal('.5')*Decimal(10)**energies[-1].as_tuple().exponent);pins.append(moment['analyze_output'])
        with localcontext() as ctx:
            ctx.prec=40;denom=Decimal(str(H))*Decimal(str(h))**2
            curved=float(-(energies[1]+energies[2]-2*energies[0])/denom);printing=float((bounds[1]+bounds[2]+2*bounds[0])/denom)
        saved=comparison['model_response']['states'][state]['steps']['h'];analytic=initial['response']['chi_states_a0_cubed'][state]
        ar=16*np.finfo(float).eps*max(abs(float(x)) for x in energies)/(H*h*h)
        check('saved_field_curvature_reconstruction_'+state,abs(curved-saved['chi_E_a0_cubed'])<=ar)
        check('analytic_within_original_energy_print_bound_'+state,abs(analytic-curved)<=printing)
        check('original_print_bound_retained_'+state,printing<=saved['energy_curvature_print_bound_a0_cubed'])
        proof[state]={'raw_energy_sources_zero_plus_minus':pins,'raw_energy_tokens_kcal_mol':[str(x) for x in energies],
                      'field_au':h,'reconstructed_energy_curvature_a0_cubed':curved,'original_reported_energy_curvature_a0_cubed':saved['chi_E_a0_cubed'],
                      'analytic_matrix_curvature_a0_cubed':analytic,'difference_a0_cubed':analytic-curved,'raw_energy_print_bound_a0_cubed':printing,
                      'original_retained_print_bound_a0_cubed':saved['energy_curvature_print_bound_a0_cubed'],'passed':True}
    old_response=copy.deepcopy(initial['response'])
    for s in old_response['chi_states_a0_cubed']:old_response['chi_states_a0_cubed'][s]=comparison['model_response']['states'][s]['steps']['h']['chi_E_a0_cubed']
    old_score=ev.h['objective'].evaluate(ev.data,initial['model'],old_response,ev.targets,ev.h['old_objective'],ev.h['body_helper'])
    check('old_EF_and_response_targets_unchanged',old_score['old_EF']==initial['score']['old_EF'] and ev.targets==ev.h['objective'].targets_from_inventory(ev.target_inventory))
    e.save(a.out/'SUSCEPTIBILITY_CONVENTION.json',{'status':'passed_saved_model_curvature_convention','states':proof,'reference_estimator_unchanged':'h=.001 energy curvature',
           'model_evaluation':'Exact inverse-operator linear susceptibility; h=.001 centered energy curvature apart from solver/print error',
           'original_printed_field_score':old_score['objective'],'analytic_model_score':initial['score']['objective'],'difference':initial['score']['objective']-old_score['objective'],
           'source_comparison':ev.target_inventory['artifacts']['field_comparison'],'independent_matrix_review':cfg['baseline_review'],'new_native_calls':0})
    e.save(a.out/'RESULT.json',{'status':'passed_saved_evaluator_controls','implementation':cfg['implementation'],'preparation':e.pin(a.preparation/'RESULT.json'),'test_implementation':e.pin(__file__),
           'checks':checks,'count':len(checks),'curvature_convention':e.pin(a.out/'SUSCEPTIBILITY_CONVENTION.json'),'new_native_calls':0,'new_parameter_molecular_points':0,
           'cost':{'wall_seconds':time.monotonic()-tick,'CPU_seconds':time.process_time()-cpu,'threads':1}})
    print(json.dumps({'result':e.pin(a.out/'RESULT.json'),'checks':len(checks),'new_native_calls':0}))


if __name__=='__main__':main()
