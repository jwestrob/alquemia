"""Report all five declared multisite descriptors without borrowing PQQ bands."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import shutil
from affordable_common import InvalidArtifact, read_json, record, verify, write_new
from mace_file_checks import cached_file_checks
from mace_hybrid import EV_TO_KCAL
from mace_omol_locality import endpoint
from mace_omol_multisite import audit
from mace_omol_ablation_run import PROTOCOL

ORDER = ('PARV_4CPV_CD','PARV_4CPV_EF','AEQ_1SL8_EF1','AEQ_1SL8_EF3','AEQ_1SL8_EF4')


@cached_file_checks
def report(reports, ggr_report, agreement, output):
    if len(reports) != len(ORDER): raise InvalidArtifact('all five ordered site reports required')
    scores = {}; checks = []; reference = None; model_rows = []
    for name, path in zip(ORDER, reports):
        r = read_json(path); p = read_json(verify(r['preparation']))
        if p['case_id'] != name or r['protocol_id'] != PROTOCOL:
            raise InvalidArtifact('wrong site order or protocol')
        audit(verify(r['preparation']))
        if r['calibrated_class'] is not None: raise InvalidArtifact('no multisite biological band supplied')
        if reference is None: reference = r['factorization_reference']
        if r['factorization_reference'] != reference: raise InvalidArtifact('incompatible descriptor reference')
        ref = read_json(verify(reference)); energies = {}
        for metal in ('La','Ca'):
            row = r['rows'][metal+'_bound_primary']; t, arrays = endpoint(row)
            if t['metal'] != metal or t['case_id'] != name or t['preparation'] != r['preparation']:
                raise InvalidArtifact('receipt identity differs')
            energies[metal] = row['energy_eV']; model_rows.append(row)
        value = (energies['Ca']-energies['La']-ref['Ca_minus_La_disconnected_atom_model_eV'])*EV_TO_KCAL
        if value != r['R_mask_model_kcal']: raise InvalidArtifact('score algebra differs from actual endpoints')
        checks.append({'name': name+'_endpoint_accounting', 'pass': r['numerical_gate_pass']})
        scores[name] = {'R_mask_model_kcal': value, 'bound_model_eV': energies,
                        'report': record(path), 'preparation': r['preparation'],
                        'evidence': p['evidence'], 'physical_atom_count': len(p['physical_atoms']),
                        'background_metal_count': len(p['background_metals']),
                        'explicit_water_count': len(p['explicit_waters']),
                        'charge_e': {m:p['endpoints'][m]['charge'] for m in ('La','Ca')},
                        'calibrated_class': None}
    for group in (ORDER[:2], ORDER[2:]):
        base = scores[group[0]]['bound_model_eV']['Ca']
        for name in group[1:]:
            error = (scores[name]['bound_model_eV']['Ca']-base)*EV_TO_KCAL
            checks.append({'name': name+'_all_Ca_permutation', 'error_model_kcal': error, 'pass': abs(error)<=.01})
    ggr = read_json(ggr_report)
    if ggr['protocol_id'] != PROTOCOL or not ggr['numerical_gate_pass']:
        raise InvalidArtifact('actual compatible GGR reference required')
    margins = []
    for name in ORDER[:2]:
        for negative in ggr['GGR_order']:
            value = scores[name]['R_mask_model_kcal']-ggr['scores'][negative]['reported_R_mask_model_kcal']
            margins.append({'positive_case': name, 'negative_case': negative,
                            'margin_model_kcal': value, 'required_min_model_kcal': .02, 'pass': value>.02})
    result = {'status': 'complete', 'protocol_id': PROTOCOL, 'agreement': record(agreement),
              'ggr_report': record(ggr_report), 'site_order': ORDER, 'scores': scores,
              'checks': checks, 'numerical_gate_pass': all(c['pass'] for c in checks),
              'parvalbumin_supporting_contrasts': margins,
              'parvalbumin_supporting_all_case_gate_pass': all(c['pass'] for c in margins),
              'aequorin_ordered_vector_model_kcal': [scores[n]['R_mask_model_kcal'] for n in ORDER[2:]],
              'aequorin_site_resolved_direction': None, 'aequorin_assay_reproduction_claimed': False,
              'parvalbumin_sub_kcal_site_ordering_target': False, 'biological_groups_added': 2,
              'gold_same_assay_site_resolved_groups_added': 0, 'broad_affinity_validated': False,
              'baseline_changed': False, 'production_promotion': False, 'new_report_model_calls': 0,
              'model_cost': {'new_forwards': len(model_rows),
                             'summed_evaluation_seconds': sum(r['evaluation_seconds'] for r in model_rows),
                             'summed_worker_wall_seconds': sum(r['wall_seconds'] for r in model_rows),
                             'peak_cuda_allocated_bytes': max(r['peak_cuda_allocated_bytes'] for r in model_rows),
                             'peak_worker_host_RSS_KiB': max(r['peak_host_RSS_KiB'] for r in model_rows)}}
    out = Path(output).resolve(); out.mkdir(parents=True,exist_ok=False)
    shutil.copyfile(__file__,out/Path(__file__).name)
    result['report_implementation'] = record(out/Path(__file__).name)
    write_new(out/'result.json',result)
    return result


if __name__=='__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--reports',required=True,nargs=5)
    for key in ('ggr-report','agreement','output'): p.add_argument('--'+key,required=True)
    r=report(**vars(p.parse_args()))
    print(json.dumps({k:r[k] for k in ('numerical_gate_pass','parvalbumin_supporting_all_case_gate_pass',
                                      'aequorin_ordered_vector_model_kcal','model_cost')},indent=2))
