#!/usr/bin/env python3
"""Finite saved-only fourteen-state termination-sensitivity comparison."""
import argparse
from decimal import Decimal
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import re
import time

HERE = Path(__file__).resolve().parent
LABELS = ('original_pre_response_baseline', 'feca12_physical_gate_terminal',
          'frozen_radial14_terminal5d992')
CYCLES = ('aquo9_central_La__far', 'aquo9_central_Dy__far')
THREADS = ('OMP_NUM_THREADS', 'OMP_THREAD_LIMIT', 'OPENBLAS_NUM_THREADS',
           'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS', 'BLIS_NUM_THREADS')


def require(ok, message):
    if not ok:
        raise ValueError(message)


def pin(path):
    p = Path(path).resolve()
    return {'path': str(p), 'sha256': hashlib.sha256(p.read_bytes()).hexdigest()}


def verify(p):
    require(isinstance(p, dict) and pin(p['path']) == p, 'Missing/changed pin: '+str(p))
    return Path(p['path'])


def read(p):
    return json.loads(verify(p).read_text())


def load(p):
    spec = importlib.util.spec_from_file_location('precision_cycle_helpers', verify(p))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def finite(x):
    require(isinstance(x, (int, float)) and not isinstance(x, bool) and math.isfinite(x), 'Nonfinite scalar')
    return x


def vector(x):
    require(isinstance(x, list) and len(x) == 3, 'Three signed components required')
    return [finite(v) for v in x]


def sub(a, b):
    return [x-y for x, y in zip(vector(a), vector(b))]


def norm(x):
    return math.sqrt(math.fsum(v*v for v in vector(x)))


def denominator(expected, records):
    missing = [cid for cid in expected if cid not in records or records[cid].get('status') != 'complete']
    require(len(expected) == len(set(expected)) == 14 and not missing and set(records) == set(expected),
            'Incomplete fourteen-state denominator: '+repr(missing))


def phase(record, ids, manifest, implementation, phase_name, count):
    require(record['status'] == ('reference_precision_scout_complete' if phase_name == 'scout'
                                else 'reference_precision_all14_complete') and record['error'] is None,
            'Reference phase is not successfully terminal')
    require(record['manifest'] == manifest and record['implementation'] == implementation and
            record['phase'] == phase_name and record['native_attempts_this_phase'] == count and
            record['native_attempts_total'] == (2 if phase_name == 'scout' else 14) and
            record['scheduler_job_id'] is None, 'Reference phase identity/count differs')
    require(set(record['cells']) == set(ids), 'Reference phase inventory differs')


def current_gate(c):
    limits = {'Energy change': Decimal('1e-11'), 'MAX-Density change': Decimal('1e-8'),
              'RMS-Density change': Decimal('1e-9'), 'Orbital Gradient': Decimal('1e-8'),
              'Orbital Rotation': Decimal('1e-8'), 'DIIS Error': Decimal('1e-8')}
    required = {'Energy change', 'MAX-Density change', 'RMS-Density change'}
    require(c['algorithm'] in ('SOSCF', 'DIIS'), 'Unknown terminal solver mode')
    required |= {'Orbital Gradient', 'Orbital Rotation'} if c['algorithm'] == 'SOSCF' else {'DIIS Error'}
    require(c['status'] == 'passed_current_scf_residuals' and c['all_required_passed'] is True and
            c['mode'] == 'All-Criteria' and set(c['required_checks']) == required and
            all(v is True for v in c['required_checks'].values()), 'Current residual admission differs')
    for name in required:
        f = c['summary'][name]
        value, tolerance = (Decimal(f[k].replace('D', 'E').replace('d', 'e'))
                            for k in ('value_token', 'tolerance_token'))
        require(value.is_finite() and tolerance == limits[name] and abs(value) <= tolerance and
                (name == 'Energy change' or value >= 0), 'Failed current residual: '+name)


def energy_printing(output, electronic):
    # Printing only: do not rerun the native/property admission or SCF parser.
    text = verify(output).read_text()
    hits = re.findall(r'^\s*FINAL SINGLE POINT ENERGY\s+([-+0-9.eEdD]+)\s*$', text, re.M)
    require(len(hits) == 1, 'Missing/duplicate final-energy printing token')
    value = Decimal(hits[0].replace('D', 'E').replace('d', 'e'))
    require(value.is_finite(), 'Nonfinite printed energy')
    half = Decimal(10)**value.as_tuple().exponent/2
    require(abs(float(value)-finite(electronic)) <= float(half)+2*math.ulp(electronic),
            'Electronic energy is not the same no-D4 final printing stage')
    return {'token': hits[0], 'half_ulp_hartree': float(half),
            'scope': 'Literal printed rounding allowance only; excludes SCF/grid/method error.'}


def moment_gate(m, cell, output, property_file):
    require(m['output'] == output and m['property_file'] == property_file and
            m['charge'] == cell['charge'] and m['labels'] == cell['symbols'] and
            m['ghost_mask'] == cell['ghost_mask'] and m['native_electrons'] == cell['expected_electrons'],
            'Admitted common-origin moment identity differs')
    for name in ('total_dipole_au', 'electronic_dipole_au', 'nuclear_dipole_au',
                 'precise_origin_Bohr', 'printed_origin_Bohr'):
        vector(m[name])
    require(finite(m['print_bounds_au']['total_dipole_au']) >= 0, 'Negative print allowance')


def endpoints(cfg, sources, m, inv, d, progress):
    ids = [r['cell_id'] for r in m['rows']]
    expected = [r['cell_id'] for r in inv['rows']]
    require(ids == expected and len(ids) == 14 and m['inventory'] == sources['inventory'] and
            m['declaration'] == sources['declaration'] and m['implementation'] == sources['runtime'] and
            m['residual_parser'] == sources['residual_parser'], 'Frozen input/source identity differs')
    scout, full = read(cfg['scout_result']), read(cfg['all14_result'])
    progress['missing_states'] = [cid for cid in ids if full.get('cells', {}).get(cid, {}).get('status') != 'complete']
    phase(scout, ids, sources['runtime_manifest'], sources['runtime'], 'scout', 2)
    phase(full, ids, sources['runtime_manifest'], sources['runtime'], 'remaining', 12)
    denominator(ids, full['cells'])
    for p, result_pin, status in ((cfg['scout_root_review'], cfg['scout_result'], 'admitted_reference_precision_scout'),
                                 (cfg['all14_root_review'], cfg['all14_result'], 'admitted_reference_precision_all14')):
        r = read(p)
        require(r['status'] == status and r['result'] == result_pin and
                r['manifest'] == sources['runtime_manifest'], 'Actual reference admission differs')
    phase_cfgs, phase_releases = {}, {}
    for name, record in (('scout', scout), ('remaining', full)):
        pc, rel = read(record['configuration']), read(record['release'])
        require(pc['status'] == 'bound_reference_precision_login_configuration' and pc['phase'] == name and
                pc['implementation'] == sources['runtime'] and pc['manifest'] == sources['runtime_manifest'] and
                pc['campaign_root'] == m['campaign_root'] and
                rel['status'] == 'approved_reference_precision_login_phase' and
                rel['configuration'] == record['configuration'] and rel['implementation'] == sources['runtime'] and
                rel['manifest'] == sources['runtime_manifest'] and rel['phase'] == name and
                rel['maximum_total_attempts'] == 14, 'Phase configuration/release differs')
        require(rel['cell_ids'] == d['staging']['scout_order' if name == 'scout' else 'continuation_order'] and
                rel['maximum_phase_attempts'] == (2 if name == 'scout' else 12), 'Phase release denominator differs')
        phase_cfgs[name], phase_releases[name] = pc, rel
    require(phase_cfgs['remaining']['scout_result'] == cfg['scout_result'] and
            phase_cfgs['remaining']['scout_root_review'] == cfg['scout_root_review'] and
            phase_releases['remaining']['scout_result'] == cfg['scout_result'] and
            phase_releases['remaining']['scout_root_review'] == cfg['scout_root_review'], 'Scout reuse chain differs')
    terms, ledger = {}, []
    for row in m['rows']:
        cid, old = row['cell_id'], row['old']
        progress['current_cell'] = cid
        require(old == next(x for x in inv['rows'] if x['cell_id'] == cid), 'Old reference row differs')
        which = 'scout' if cid in d['staging']['scout_order'] else 'remaining'
        producer = scout if which == 'scout' else full
        if which == 'scout':
            require(full['cells'][cid] == scout['cells'][cid], 'Scout result was replaced')
        rp = full['cells'][cid]['result']; r = read(rp)
        require(r['status'] == 'complete' and r['error'] is None and r['cell_id'] == cid and
                r['cell'] == old['cell'] and r['old_reference'] == old['result'] and
                r['manifest'] == sources['runtime_manifest'] and r['implementation'] == sources['runtime'] and
                r['configuration'] == producer['configuration'] and r['release'] == producer['release'] and
                r['native_invocation_attempts'] == 1, 'Endpoint identity/status differs: '+cid)
        folder = Path(m['campaign_root'])/'cells'/cid
        require(rp == pin(folder/'RESULT.json') and r['input'] == pin(folder/'endpoint.inp') and
                r['input']['sha256'] == row['input']['sha256'], 'Endpoint deck/path differs')
        start = read(pin(folder/'START.json'))
        require(start['status'] == 'started' and all(start[k] == r[k] for k in
                ('cell_id', 'implementation', 'manifest', 'configuration', 'release', 'cell',
                 'old_reference', 'input', 'native_invocation_attempts')), 'Durable attempt intent differs')
        receipt = read(r['native_receipt'])
        require(receipt['returncode'] == 0 and receipt['output'] == r['output'] and
                receipt['stderr'] == r['stderr'] and receipt['scheduler_job_id'] is None and
                receipt['argv'] == [m['orca']['path'], 'endpoint.inp', '--bind-to none'] and
                receipt['cwd'] == str(folder) and receipt['host'].split('.')[0] == 'biotite' and
                all(receipt['threads'].get(k) == '1' for k in THREADS), 'Native receipt differs')
        for p in (r['output'], r['stderr'], r['property_file'], receipt['process_samples'], old['result'], old['property_file']):
            verify(p)
        parsed = read(r['native_readback']); current = read(r['current_residuals_file'])
        require(parsed == r['parsed'] and current == r['current_residuals'] and parsed['status'] == 'complete' and
                parsed['checks'] and all(v is True for v in parsed['checks'].values()) and
                parsed['D4_energy_hartree'] == 0, 'Saved native readback differs')
        current_gate(current)
        moment_gate(r['moment'], old['cell'], r['output'], r['property_file'])
        moment_gate(old['old_moment'], old['cell'], old['output'], old['property_file'])
        old_e = finite(old['old_native_metadata']['electronic_energy_hartree'])
        new_e = finite(parsed['electronic_energy_hartree'])
        old_grid, new_grid = old['old_native_metadata']['final_COSX_correction_hartree'], parsed['final_COSX_correction_hartree']
        require(len(old_grid) == len(new_grid) == 1, 'Final grid correction count differs')
        finite(old_grid[0]); finite(new_grid[0])
        terms[cid] = {'result': rp, 'coefficient': old['coefficient'], 'cycle_id': old['cycle_id'],
            'old': {'electronic_energy_hartree': old_e, 'moment': old['old_moment'],
                    'native_metadata': old['old_native_metadata'], 'terminal_trace': old['SCF_trace'],
                    'energy_printing': energy_printing(old['output'], old_e)},
            'new': {'electronic_energy_hartree': new_e, 'moment': r['moment'], 'current_residuals': current,
                    'native_metadata': parsed, 'energy_printing': energy_printing(r['output'], new_e)},
            'electronic_energy_shift_hartree': new_e-old_e,
            'common_origin_dipole_shift_au': sub(r['moment']['total_dipole_au'], old['old_moment']['total_dipole_au']),
            'gridx3_correction_shift_hartree': new_grid[0]-old_grid[0]}
        ledger.append({'cell_id': cid, 'phase': which, 'result': rp, 'start': pin(folder/'START.json'),
                       'native_receipt': r['native_receipt'], 'cost_and_process_observations': receipt})
        progress['checked_states'].append(cid)
    require(len({x['native_receipt']['path'] for x in ledger}) == 14, 'Duplicate attempt receipt')
    return terms, ledger, {'scout': scout['cost'], 'remaining': full['cost']}


def cycle_values(terms, coefficients, side, helper):
    require(len(coefficients) == 7 and sorted(coefficients.values()) == [-1, -1, -1, 1, 1, 1, 1],
            'Signed seven-term formula differs')
    rows = {k: terms[k][side]['moment'] for k in coefficients}
    require(math.fsum(c*rows[k]['charge'] for k, c in coefficients.items()) == 0, 'Nonneutral cycle')
    return {'E3_kcal_mol': helper.HA*math.fsum(c*terms[k][side]['electronic_energy_hartree'] for k,c in coefficients.items()),
            'mu3_Debye': [helper.CONVERSION*v for v in helper.combine(rows, coefficients, 'total_dipole_au')],
            'dipole_print_bound_Debye_per_component': helper.CONVERSION*math.fsum(
                abs(c)*rows[k]['print_bounds_au']['total_dipole_au'] for k,c in coefficients.items()),
            'energy_print_bound_kcal_mol': helper.HA*math.fsum(abs(c)*terms[k][side]['energy_printing']['half_ulp_hartree']
                                                             for k,c in coefficients.items())}


def shift(old, new):
    delta = sub(new['mu3_Debye'], old['mu3_Debye'])
    return {'E3_kcal_mol': new['E3_kcal_mol']-old['E3_kcal_mol'], 'mu3_Debye': delta,
            'mu3_vector_shift_norm_Debye': norm(delta)}


def model_error(model, reference):
    delta = sub(model['mu3_Debye'], reference['mu3_Debye'])
    return {'E3_kcal_mol': finite(model['E3_kcal_mol'])-reference['E3_kcal_mol'],
            'mu3_Debye': delta, 'mu3_error_norm_Debye': norm(delta)}


def arithmetic(terms, inv, model_result, helper):
    denominator([r['cell_id'] for r in inv['rows']], {k: {'status': 'complete'} for k in terms})
    cycles = {}
    for cid in CYCLES:
        authority = inv['cycles'][cid]; coefficients = authority['reference_terms']
        old, new = (cycle_values(terms, coefficients, side, helper) for side in ('old', 'new'))
        require(old['E3_kcal_mol'] == authority['reference_electronic_E3_kcal_mol'] and
                old['mu3_Debye'] == authority['reference_mu3_Debye'] and
                old['dipole_print_bound_Debye_per_component'] == authority['reference_print_bound_Debye_per_component'],
                'Original cycle authority no longer replays exactly')
        stored = model_result['comparison']['cycles'][cid]
        require(stored['status'] == 'complete' and stored['reference_electronic_E3_kcal_mol'] == old['E3_kcal_mol'] and
                stored['reference_mu3_Debye'] == old['mu3_Debye'] and set(stored['models']) == set(LABELS),
                'Fixed model/reference join differs')
        models = {label: {'unchanged_model': stored['models'][label],
                         'error_against_old': model_error(stored['models'][label], old),
                         'error_against_tightened': model_error(stored['models'][label], new)} for label in LABELS}
        cycles[cid] = {'coefficients': coefficients, 'old': old, 'tightened': new,
                       'measured_protocol_shift': shift(old, new), 'fixed_models': models}
    paired = {}
    la, dy = (cycles[k] for k in CYCLES)
    for side in ('old', 'tightened'):
        paired[side] = {'E3_kcal_mol': dy[side]['E3_kcal_mol']-la[side]['E3_kcal_mol'],
                        'mu3_Debye': sub(dy[side]['mu3_Debye'], la[side]['mu3_Debye']),
                        'dipole_print_bound_Debye_per_component': sum(c[side]['dipole_print_bound_Debye_per_component'] for c in (la, dy)),
                        'energy_print_bound_kcal_mol': sum(c[side]['energy_print_bound_kcal_mol'] for c in (la, dy))}
    paired['measured_protocol_shift'] = shift(paired['old'], paired['tightened'])
    paired['fixed_models'] = {}
    for label in LABELS:
        a, b = (c['fixed_models'][label]['unchanged_model'] for c in (la, dy))
        model = {'E3_kcal_mol': b['E3_kcal_mol']-a['E3_kcal_mol'], 'mu3_Debye': sub(b['mu3_Debye'], a['mu3_Debye'])}
        paired['fixed_models'][label] = {'unchanged_model_Dy_minus_La': model,
            'error_against_old': model_error(model, paired['old']),
            'error_against_tightened': model_error(model, paired['tightened'])}
    return cycles, paired


def compare(cfg, progress):
    require(cfg['implementation'] == pin(__file__), 'Comparator source differs')
    sources = read(cfg['inputs']); d = read(sources['declaration']); inv = read(sources['inventory'])
    require(d['inventory'] == sources['inventory'] and d['comparison']['term_count'] == 14 and
            d['comparison']['cycle_count'] == 2 and set(inv['cycles']) == set(CYCLES), 'Comparison declaration differs')
    m = read(sources['runtime_manifest'])
    sr = read(sources['runtime_source_review'])
    require(sr['status'] == 'admitted_reference_precision_runtime_source' and sr['implementation'] == sources['runtime'] and
            sr['manifest'] == sources['runtime_manifest'], 'Runtime source admission differs')
    ir = read(sr['independent_review'])
    require(ir['status'] == 'passed_reference_precision_runtime_source_review' and
            ir['implementation'] == sources['runtime'] and ir['manifest'] == sources['runtime_manifest'] and
            ir['residual_parser'] == sources['residual_parser'], 'Independent runtime source join differs')
    helper = load(sources['cycle_helper'])
    require(helper.HA == 627.509474 and helper.CONVERSION == 0.529177210544*4.803206802, 'Inherited units differ')
    models, model_review = read(sources['model_result']), read(sources['model_root_review'])
    require(model_review['result'] == sources['model_result'] and model_review['status'] == sources['model_root_status'] and
            models['error'] is None, 'Fixed model actual admission differs')
    terms, ledger, phase_costs = endpoints(cfg, sources, m, inv, d, progress)
    for term in terms.values():
        term['electronic_energy_shift_kcal_mol'] = helper.HA*term['electronic_energy_shift_hartree']
        term['common_origin_dipole_shift_Debye'] = [helper.CONVERSION*v for v in term['common_origin_dipole_shift_au']]
        term['common_origin_dipole_shift_norm_Debye'] = norm(term['common_origin_dipole_shift_Debye'])
    cycles, paired = arithmetic(terms, inv, models, helper)
    return {'status': 'complete_saved_reference_precision_comparison', 'terms': terms, 'cycles': cycles,
            'paired_Dy_minus_La': paired, 'attempt_ledger': ledger, 'missing_states': [],
            'reference_attempts': 14, 'reference_phase_costs': phase_costs,
            'reference_cost_summary': {
                'native_wall_seconds_sum': math.fsum(x['cost_and_process_observations']['wall_seconds'] for x in ledger),
                'native_child_CPU_seconds_sum': math.fsum(x['cost_and_process_observations']['child_user_CPU_seconds']+
                    x['cost_and_process_observations']['child_system_CPU_seconds'] for x in ledger),
                'max_observed_RSS_pages_sum': max(x['cost_and_process_observations']['process_observations']['peak_observed_RSS_pages_sum'] for x in ledger),
                'note': 'Login observations; no allocated CPU-seconds. Phase wall contains native wall and is not added to it. RSS is sampled or controller highwater, not additive.'},
            'printing_scope': 'Literal rounding allowances and inherited dipole print bounds are separate from measured protocol shifts; neither is an SCF error bar.',
            'limitations': [
                'One stricter termination protocol measures sensitivity, not a rigorous error bound or method truth.',
                'PBE0 electronic energies, full-parent ghost basis and no D4 remain unchanged; no grid/basis/correlation sensitivity is tested.',
                'Current residual admission does not prove final-grid self-consistency or electronic-root continuity.',
                'Restricted effective multiplicity1 remains distinct from physical La1/Dy6 parent states.',
                'Two consumed W9far cycles do not establish full-shell B, forces, general model adequacy or global family infeasibility.',
                'All model values remain fixed; no fitting, candidate selection, target replacement or new physical threshold.']}


def run(config_path, release_path, output):
    out = Path(output).resolve(); out.mkdir(parents=True, exist_ok=False)
    tick, cpu = time.monotonic(), time.process_time()
    progress = {'missing_states': [], 'checked_states': []}
    result = {'status': 'failed_saved_reference_precision_comparison', 'error': None}
    try:
        require(all(os.environ.get(k) == '1' for k in THREADS), 'Single-thread environment required')
        cfg_pin, release_pin = pin(config_path), pin(release_path)
        cfg, rel = read(cfg_pin), read(release_pin)
        sources = read(cfg['inputs']); inv = read(sources['inventory'])
        progress['missing_states'] = [r['cell_id'] for r in inv['rows']]
        require(rel['status'] == 'approved_saved_reference_precision_comparison' and rel['configuration'] == cfg_pin and
                rel['implementation'] == pin(__file__) and rel['inputs'] == cfg['inputs'], 'Saved comparison release differs')
        result = compare(cfg, progress)
        result.update(configuration=cfg_pin, release=release_pin, inputs=cfg['inputs'])
    except Exception as exc:
        result.update(status='failed_saved_reference_precision_comparison', error=repr(exc), admission_progress=progress)
    result.update(implementation=pin(__file__), new_native_reference_solver_calls=0, raw_native_admission_replays=0,
                  physical_acceptance_threshold=None, reference_replacement=False,
                  cost={'wall_seconds': time.monotonic()-tick, 'CPU_seconds': time.process_time()-cpu, 'threads': 1})
    (out/'RESULT.json').write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    return result


if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('--config', required=True); ap.add_argument('--release', required=True); ap.add_argument('--out', required=True)
    args = ap.parse_args(); r = run(args.config, args.release, args.out)
    print(json.dumps({'status': r['status'], 'error': r.get('error')}))
    raise SystemExit(0 if r['status'] == 'complete_saved_reference_precision_comparison' else 1)
