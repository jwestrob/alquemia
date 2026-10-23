#!/usr/bin/env python3
"""Pinned admitted-value extraction only; no evaluator or numerical solver import."""
import hashlib
import json
import math
from pathlib import Path
import time

HERE = Path(__file__).resolve().parent

def pin(p):
    return {'path': str(p.resolve()), 'sha256': hashlib.sha256(p.read_bytes()).hexdigest()}

def main():
    tick, cpu = time.monotonic(), time.process_time()
    inp = json.loads((HERE / 'INPUTS.json').read_text())
    src = inp['sources']
    for p in src.values():
        assert pin(Path(p['path'])) == p
    def read(name):
        return json.loads(Path(src[name]['path']).read_text())
    fit, audit, vec = read('fit_result'), read('fit_review'), read('selected')
    tri, tria = read('trimer_result'), read('trimer_review')
    context, prior, domain = read('electrical_context'), read('prior_values'), read('domain')
    target, prop = read('targets'), read('last_proposal')
    assert fit['status'] == 'bounded_physical_gate_component_fit_complete'
    assert audit['status'] == 'passed_physical_gate_component_fit_actual_review'
    assert audit['result'] == src['fit_result'] and audit['selected'] == src['selected'] == fit['selected']
    assert tria['status'] == 'passed_consumed_trimer_arithmetic_review'
    assert tria['result'] == src['trimer_result'] and tria['fit_result'] == src['fit_result']
    assert tria['fit_review'] == src['fit_review'] and tri['status'] == 'consumed_trimer_regression_complete'
    m, origin = fit['selected_metrics'], fit['initial_metrics']
    assert m == audit['selected_metrics'] and len(vec['parameters']) == 12
    assert m['all25_original_physical_gates'] and audit['selected_native_admission_passed']
    assert target['chi_a0_cubed'] == context['reference_chi_states_a0_cubed']
    chi = {k: {'model': v, 'reference': target['chi_a0_cubed'][k], 'error': v-target['chi_a0_cubed'][k]}
           for k,v in vec['response']['chi_states_a0_cubed'].items()}
    increments = {}
    for metal in ('La', 'Dy'):
        increments[metal] = {k: chi[metal+'_trimer'][k]-chi[metal+'_dimer'][k] for k in ('model','reference','error')}
    increments['Dy_minus_La'] = {k: increments['Dy'][k]-increments['La'][k] for k in ('model','reference','error')}
    cd = {k: {'model': v, 'reference': target['C_Debye'][k], 'error': [a-b for a,b in zip(v,target['C_Debye'][k])]}
          for k,v in vec['response']['C_dipoles_Debye'].items()}
    cd_pair = {str(n): {k: [a-b for a,b in zip(cd[f'aquo{n}_central_Dy'][k],cd[f'aquo{n}_central_La'][k])]
                       for k in ('model','reference','error')} for n in (8,9)}
    parameters = []
    for p,v in zip(domain['parameters'],vec['parameters']):
        width = p['upper']-p['lower']
        parameters.append(dict(p, value=v, lower_distance=v-p['lower'], upper_distance=p['upper']-v,
                               normalized_nearest_edge_distance=min(v-p['lower'],p['upper']-v)/width))
    mapped = {}
    for metal,i in [('La',0),('Dy',4)]:
        a,b,u,v=vec['parameters'][i:i+4]
        mapped[metal]={'a':a,'b':b,'u_box_coordinate':u,'v_box_coordinate':v,
                       'c_native':domain['S']*(u+v)/4,'d_native':domain['S']*(u-v)/4}
    margins = {}
    for k,g in m['original25'].items():
        if k.startswith('energy/'):
            margins[k]={'error_kcal_mol':g['error_kcal_mol'],'absolute_error_margin_to_1':1-abs(g['error_kcal_mol'])}
        else:
            margins[k]={**g,'RMS_margin_to_2':2-g['rigid_RMS_kcal_mol_A'],
                        'component_margin_to_5':5-g['maximum_body_component_error_kcal_mol_A']}
    cost = fit['cost']
    total = math.fsum(cost[k] for k in ('controller_CPU_seconds','child_user_CPU_seconds','child_system_CPU_seconds'))
    values = {'status':'admitted_terminal_saved_synthesis_not_physical_accuracy_admission',
        'implementation':pin(Path(__file__)), 'inputs':pin(HERE/'INPUTS.json'), 'sources':src,
        'selected_key':fit['selected_key'], 'selected_metrics':m, 'origin24_metrics':origin,
        'previous_component_terminal':prior['prior_terminal'],
        'original_prefit_component_RMS':m['original_initial_component_RMS'],
        'RMS_fractional_reductions':{base:{k:1-m['component_RMS'][k]/ref[k] for k in ('P','B')}
             for base,ref in [('origin24',origin['component_RMS']),('original_prefit',m['original_initial_component_RMS'])]},
        'derived_C_paired_errors_kcal_mol':{str(n):m['component_errors_kcal_mol']['P'+str(n)]+m['component_errors_kcal_mol']['B'+str(n)] for n in (8,9)},
        'chi_states_a0_cubed':chi, 'conditional_chi_increments_a0_cubed':increments,
        'C_dipoles_Debye':cd, 'C_dipoles_Dy_minus_La_Debye':cd_pair,
        'electrical_RMS_reference_scale':context['reference_RMS_with_original_weights'],
        'parameters_and_box_distances':parameters, 'mapped_native_response_coefficients':mapped,
        'physical_gate_margins':margins, 'trimer_cycles':tria['cycles'],
        'trimer_paired_Dy_minus_La':tria['paired_Dy_minus_La'], 'trimer_comparison_counts':tria['comparison_counts'],
        'fit_counts':audit['counts'], 'fit_stopping_reason':fit['stopping_reason'], 'fit_trajectory':audit['trajectory'],
        'last_Jacobian_origin':prop['origin_key'], 'last_Jacobian_at_selected':prop['origin_key']==fit['selected_key'],
        'fit_cost_literal_field':cost, 'fit_total_measured_CPU_seconds':total,
        'new_native_calls':0,'new_solver_calls':0,'new_reference_calls':0,'heldout_reads':0,
        'scope':'Simple arithmetic on already admitted artifacts; no raw native readmission or new proposal.'}
    assert not values['last_Jacobian_at_selected']
    assert len(chi)==4 and len(cd)==4 and len(tria['cycles'])==8 and len(tria['paired_Dy_minus_La'])==4
    values['extraction_cost']={'wall_seconds':time.monotonic()-tick,'CPU_seconds':time.process_time()-cpu,'threads':1}
    with (HERE/'VALUES.json').open('x') as h:
        json.dump(values,h,indent=2,allow_nan=False);h.write('\n')
    print(json.dumps({'values':pin(HERE/'VALUES.json'),'cost':values['extraction_cost'],
        'chi_increments':increments,'native_coefficients':mapped,'RMS_reductions':values['RMS_fractional_reductions'],
        'derived_C_errors':values['derived_C_paired_errors_kcal_mol'],'measured_fit_CPU_seconds':total},indent=2))

if __name__ == '__main__':
    main()
