"""Summarize actual joint candidates and matched common-pool scores; no calculations."""
import argparse
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'scripts'))
from affordable_common import read_json, record, verify, write_new
from nikasha_pool import JOINT_PROTOCOL, choose_rows, score
from nikasha_pool_compare import old_bands


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for key in ('collection', 'proposals', 'accounting', 'output'):
        p.add_argument('--' + key, required=True)
    a = p.parse_args()
    data = read_json(a.collection); proposals = read_json(a.proposals)
    m = read_json(verify(data['manifest'])); pm = read_json(verify(proposals['manifest']))
    assert data['protocol_id'] == JOINT_PROTOCOL
    assert m['adaptive_proposals'] == record(a.proposals)
    parent = read_json(verify(m['source_manifest']))
    bands = old_bands(read_json(verify(parent['frozen_comparison'])))
    def decision(r):
        if r is None: return 'unavailable'
        return 'Ca_supported' if r <= bands['Ca_max'] else 'La_supported' if r >= bands['La_min'] else 'inconclusive'
    cases = []
    for c in data['cases']:
        if c['status'] == 'prepared':
            assert c['pool'] == choose_rows(c['matrix'], [q['id'] for q in c['candidates']])
        old = c['prior_pool']['operational']['composite_R_model_kcal_mol']
        new = c['pool']['operational']
        new = new['composite_R_model_kcal_mol'] if new else None
        origin = score(c['matrix']['Ca']['origin']['components'], c['matrix']['La']['origin']['components'])
        cases.append({'case_id': c['case_id'], 'expected_class': c['old_result']['expected_class'],
                      'label_scope': c['old_result']['label_scope'], 'status': c['pool']['status'],
                      'reason': c['reason'], 'origin_R': origin['composite_R_model_kcal_mol'],
                      'angular_R': old, 'joint_R': new, 'delta_R_from_angular': new-old if new is not None else None,
                      'old_band_angular_decision': decision(old), 'old_band_joint_decision': decision(new),
                      'joint_rows': c['pool']['rows']})
    endpoints = []; fresh_requests = cache_reuses = native_calls = 0
    for e in proposals['endpoints']:
        r = read_json(verify(e['proposal_receipt']))
        for request in r['requests']:
            cache_reuses += request['reused']
            if not request['reused']:
                fresh_requests += 1
                point = read_json(verify(request['result']))
                if point.get('MACE'):
                    native = read_json(verify(point['MACE']))
                    native_calls += bool(native.get('model_call_started'))
        endpoints.append({'task_id': e['task_id'], 'status': e['status'], 'reason': e['reason'],
                          'receipt': e['proposal_receipt'], 'physical_geometry': r.get('final_geometry'),
                          'optimizer': r['optimizer'], 'residual': r.get('residual'),
                          'MACE_work_kcal_mol': r.get('MACE_proposal_work_kcal_mol'),
                          'infeasible_completed_requests': r['infeasible_completed_MACE_requests'],
                          'maximum_trial_heavy_extent_A': r['maximum_trial_heavy_extent_A']})
    lines = Path(a.accounting).read_text().splitlines()
    header = lines[0].split('|'); jobs = []
    for line in lines[1:]:
        row = dict(zip(header, line.split('|')))
        if row['JobID'].isdigit(): jobs.append(row)
    assert {r['JobID'] for r in jobs} == {'1209901', '1209902', '1209903'}
    assert all(r['State'] == 'COMPLETED' for r in jobs)
    result = {'protocol_id': JOINT_PROTOCOL, 'collection': record(a.collection),
              'proposals': record(a.proposals), 'accounting': record(a.accounting),
              'implementation': record(__file__), 'old_reference': parent['frozen_comparison'],
              'old_bands': bands, 'new_calibration': None, 'cases': cases, 'endpoints': endpoints,
              'denominator': 4, 'available': data['available'],
              'fresh_proposal_requests': fresh_requests, 'proposal_cache_reuses': cache_reuses,
              'proposal_MACE_calls': native_calls, 'cross_MACE_calls': data['MACE_calls_started'],
              'GFN2_calls': data['required_new_GFN2_calls'],
              'allocated_core_seconds': sum(int(r['ElapsedRaw'])*int(r['AllocCPUS']) for r in jobs),
              'allocated_GPU_seconds': sum(int(r['ElapsedRaw']) for r in jobs if r['JobID'] in ('1209901', '1209902')),
              'production_changed': False, 'new_DFT_calls': 0,
              'interpretation': 'proposal-Hamiltonian stationary candidates in a bounded subspace; no composite minimum, affinity or biological PLM label'}
    write_new(a.output, result)
    print({k: result[k] for k in ('available', 'proposal_MACE_calls', 'cross_MACE_calls', 'GFN2_calls', 'allocated_core_seconds', 'allocated_GPU_seconds')})


if __name__ == '__main__': main()
