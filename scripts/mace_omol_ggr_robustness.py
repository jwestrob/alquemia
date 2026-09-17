"""Compare the frozen three-structure GGR test and explain actual saved readouts."""
from __future__ import annotations
import argparse
import csv
import json
import math
from pathlib import Path
import shutil
import numpy as np
from affordable_common import InvalidArtifact, read_json, record, verify, write_new
from mace_file_checks import cached_file_checks
from mace_hybrid import EV_TO_KCAL
from mace_omol_ablation_run import PROTOCOL
from mace_omol_locality import endpoint, group_rows

GGR = ('GGR_1GLG', 'GGR_2FW0', 'GGR_2FVY')
ALPHA = ('ALPHA_1F6S', 'ALPHA_6IP9')


def seconds(value):
    parts = list(map(float, value.split(':')))
    return sum(x * 60**i for i, x in enumerate(reversed(parts)))


@cached_file_checks
def report(development_report, ggr_2fw0_report, ggr_2fvy_report,
           agreement, diagnosis_plan, sacct, output):
    old = read_json(development_report)
    if old['protocol_id'] != PROTOCOL or not old['numerical_gate_pass']:
        raise InvalidArtifact('qualified actual masked development report required')
    sources = {name: old['scores'][name] for name in (*GGR[:1], *ALPHA)}
    new = [read_json(p) for p in (ggr_2fw0_report, ggr_2fvy_report)]
    for name, r in zip(GGR[1:], new):
        if r['protocol_id'] != PROTOCOL or r['status'] != 'complete' or not r['numerical_gate_pass']:
            raise InvalidArtifact('complete numerically valid GGR result required')
        if read_json(verify(r['preparation']))['case_id'] != name:
            raise InvalidArtifact('GGR source/report identity mismatch')
        sources[name] = {'R_mask_model_kcal': r['R_mask_model_kcal'],
                         'preparation': r['preparation'],
                         'endpoints': {m: {'bound': r['rows'][m+'_bound_primary']} for m in ('La', 'Ca')}}
    factorization = read_json(verify(new[0]['factorization_reference']))
    if new[0]['factorization_reference'] != new[1]['factorization_reference']:
        raise InvalidArtifact('incompatible disconnected references')
    reference = factorization['Ca_minus_La_disconnected_atom_model_eV']
    scores = {}; checks = []; model_rows = []
    for name, s in sources.items():
        p = read_json(verify(s['preparation'])); arrays = {}; energies = {}
        for metal in ('La', 'Ca'):
            r = s['endpoints'][metal]['bound']; t, a = endpoint(r)
            if t['case_id'] != name or t['metal'] != metal or t['position'] != 'bound':
                raise InvalidArtifact('bound endpoint identity mismatch')
            if t['preparation'] != s['preparation']:
                raise InvalidArtifact('physical mapping does not belong to endpoint')
            arrays[metal] = a['node_energy_eV'] + a['embedding_energy_eV']
            energies[metal] = r['energy_eV']
            if name in GGR[1:]: model_rows.append(r)
        values = arrays['Ca'] - arrays['La']; values[-1] -= reference
        values *= EV_TO_KCAL
        value = (energies['Ca']-energies['La']-reference)*EV_TO_KCAL
        for label, error in [('readout_closure', math.fsum(map(float, values))-value),
                             ('published_score', value-s['R_mask_model_kcal'])]:
            checks.append({'name': name+'_'+label, 'error_model_kcal': error, 'pass': abs(error) <= .01})
        atoms = p['physical_atoms']; pos = np.array([a['xyz_A'] for a in atoms])
        distances = np.linalg.norm(pos-pos[-1], axis=1)
        groups, detail = group_rows(atoms, values, distances, len(atoms)-1)
        # Names in the shared bookkeeping helper are historical. Values here are
        # explicitly descriptor units, not quantum energies or atomic observables.
        scores[name] = {'reported_R_mask_model_kcal': s['R_mask_model_kcal'],
                        'two_call_replay_model_kcal': value, 'preparation': s['preparation'],
                        'physical_atom_count': len(atoms), 'protein_charge_e': p['protein_charge_e'],
                        'explicit_water_count': len(p['explicit_waters']), 'groups': groups,
                        'nearby_heteroatoms_3p5_A': [
                            {'id': a['id'], 'resname': a.get('resname'), 'element': a['element'],
                             'distance_A': float(d)} for a, d in zip(atoms, distances)
                            if a['element'] in ('O', 'N', 'S') and d <= 3.5],
                        'node_contributions': detail}
    margins = [{'positive_case': a, 'negative_case': g,
                'margin_model_kcal': sources[a]['R_mask_model_kcal']-sources[g]['R_mask_model_kcal'],
                'required_min_model_kcal': .02} for a in ALPHA for g in GGR]
    for row in margins: row['pass'] = row['margin_model_kcal'] > row['required_min_model_kcal']
    raw = list(csv.DictReader(Path(sacct).read_text().splitlines(), delimiter='|'))
    jobs = [r for r in raw if r['JobID'] in ('1200845', '1200846')]
    if len(jobs) != 2 or any(r['State'] != 'COMPLETED' for r in jobs):
        raise InvalidArtifact('complete actual job accounting required')
    cost = {'new_successful_forwards': len(model_rows), 'failed_forward_attempts': 0,
            'GPU_allocation_seconds': sum(int(r['ElapsedRaw']) for r in jobs),
            'allocated_core_seconds': sum(int(r['CPUTimeRAW']) for r in jobs),
            'reported_actual_CPU_seconds': sum(seconds(r['TotalCPU']) for r in jobs),
            'summed_model_evaluation_seconds': sum(r['evaluation_seconds'] for r in model_rows),
            'summed_worker_wall_seconds': sum(r['wall_seconds'] for r in model_rows),
            'peak_cuda_allocated_bytes': max(r['peak_cuda_allocated_bytes'] for r in model_rows),
            'peak_worker_host_RSS_KiB': max(r['peak_host_RSS_KiB'] for r in model_rows),
            'scope': 'two inference allocations only; local preparation, tests and reports separately receipted'}
    if any(len(r['attempts']) != 2 or any(not a['accepted'] for a in r['attempts']) for r in new):
        raise InvalidArtifact('unexpected extra/failed attempts require explicit accounting')
    result = {'status': 'complete', 'protocol_id': PROTOCOL,
              'source_reports': [record(p) for p in (development_report, ggr_2fw0_report, ggr_2fvy_report)],
              'agreement': record(agreement), 'diagnosis_plan': record(diagnosis_plan), 'sacct': record(sacct),
              'GGR_order': GGR, 'alpha_order': ALPHA, 'scores': scores, 'margins': margins,
              'robustness_gate_pass': all(r['pass'] for r in margins),
              'passed_margins': sum(r['pass'] for r in margins), 'total_margins': len(margins),
              'GGR_range_model_kcal': max(sources[g]['R_mask_model_kcal'] for g in GGR)-min(sources[g]['R_mask_model_kcal'] for g in GGR),
              'checks': checks, 'numerical_gate_pass': all(c['pass'] for c in checks), 'cost': cost,
              'biological_groups': 2, 'independent_assay_conditions': False,
              'broad_affinity_validated': False, 'unique_physical_cause_established': False,
              'baseline_changed': False, 'production_promotion': False, 'new_report_model_calls': 0,
              'readout_units': 'kcal_equivalent_model_units; inherited helper kcal_mol keys are model bookkeeping'}
    out = Path(output).resolve(); out.mkdir(parents=True, exist_ok=False)
    shutil.copyfile(__file__, out/Path(__file__).name)
    result['report_implementation'] = record(out/Path(__file__).name)
    write_new(out/'result.json', result)
    return result


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    for key in ('development-report', 'ggr-2fw0-report', 'ggr-2fvy-report', 'agreement', 'diagnosis-plan', 'sacct', 'output'):
        p.add_argument('--'+key, required=True)
    r = report(**vars(p.parse_args()))
    print(json.dumps({k: r[k] for k in ('robustness_gate_pass', 'margins', 'GGR_range_model_kcal', 'numerical_gate_pass', 'cost')}, indent=2))
