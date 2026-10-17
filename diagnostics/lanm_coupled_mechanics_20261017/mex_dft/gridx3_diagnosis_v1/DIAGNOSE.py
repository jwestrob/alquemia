#!/usr/bin/env python3
"""Saved-only diagnosis of the failed Mex La O97 derivative; Slurm only."""
import argparse
from decimal import Decimal, localcontext
import hashlib
import json
import math
import os
from pathlib import Path
import re
import resource
import sys
import time
import traceback

PROTOCOL = 'mex_gridx3_saved_diagnosis_v1'
COLLECTOR_SHA = '334bb6b4374538a5f293f1e6175340309566de92e42af0588485bca5dfe20953'
QUALIFY_SHA = '74b56bf12de0543251c59dcf7a6bb1bb3b84afdedc078ac1ff3045471dec7f48'
REVIEW_SHA = '9b1166b9c906b2b5a4c776279b88d9c57b8be135b1c7b4b23f98c9c54a7f4af9'
HA, BOHR = 627.509474, .529177210903
DEC_HA = Decimal('627.509474')
H = Decimal('0.001')
SID, PI, FI = 'A/97//O', 1019, 832
PAIR = ('MM_O97_La_plus', 'MM_O97_La_minus')
NUMBER = r'([-+]?\d+(?:\.\d*)?(?:[Ee][-+]?\d+)?)'


def require(ok, message):
    if not ok:
        raise ValueError(message)


def pin(path):
    p = Path(path).resolve()
    with p.open('rb') as stream:
        sha = hashlib.file_digest(stream, 'sha256').hexdigest()
    return {'path': str(p), 'sha256': sha}


def verify(record):
    require(isinstance(record, dict) and set(record) == {'path', 'sha256'}, 'Invalid pin')
    require(pin(record['path'])['sha256'] == record['sha256'], 'Changed pin: ' + record['path'])
    return Path(record['path'])


def read(record):
    return json.loads(verify(record).read_text())


def write(path, value):
    with Path(path).open('x') as stream:
        json.dump(value, stream, indent=2, allow_nan=False)
        stream.write('\n')
    return pin(path)


def half_quantum(value):
    return Decimal(1).scaleb(Decimal(value).as_tuple().exponent) / 2


def extract(output):
    text = verify(output).read_text()
    require('ORCA TERMINATED NORMALLY' in text and 'DFTD4 V3.4.0' in text
            and 'Recomputing exchange energy using gridx3' in text, 'Expected native ending/method absent')
    lines = text.splitlines()
    converged = [(i, re.search(r'SCF CONVERGED AFTER\s+(\d+) CYCLES', line))
                 for i, line in enumerate(lines) if 'SCF CONVERGED AFTER' in line]
    require(len(converged) == 1 and converged[0][1], 'Single native SCF ending required')
    ci, cm = converged[0]
    cycles = [(i, re.match(r'^\s*(\d+)\s+' + NUMBER + r'\s+[-+0-9.eE]+\s+', line))
              for i, line in enumerate(lines[:ci])]
    cycles = [(i, match) for i, match in cycles if match]
    require(cycles and cycles[-1][1].group(1) == cm.group(1), 'Last SCF cycle identity differs')
    i, match = cycles[-1]
    fields = {'last_cycle': dict(line=i+1, text=lines[i], value=match.group(2))}
    patterns = {
        'cosx_printed': r'^\s*Exchange energy change after final integration\s*:\s*' + NUMBER + r'\s+Eh\s*$',
        'post_grid': r'^\s*Total Energy\s*:\s*' + NUMBER + r'\s+Eh\s+',
        'post_grid_short': r'^\s*Total energy after final integration\s*:\s*' + NUMBER + r'\s+Eh\s*$',
        'exchange_old': r'^\s*Old exchange energy\s*:\s*' + NUMBER + r'\s+Eh\s*$',
        'exchange_new': r'^\s*New exchange energy\s*:\s*' + NUMBER + r'\s+Eh\s*$',
        'D4_printed': r'^\s*Dispersion correction\s+' + NUMBER + r'\s*$',
        'final': r'^\s*FINAL SINGLE POINT ENERGY\s+' + NUMBER + r'\s*$'}
    for name, pattern in patterns.items():
        found = [(i, re.search(pattern, line)) for i, line in enumerate(lines)]
        found = [(i, match) for i, match in found if match]
        require(len(found) == 1, 'Expected one printed field: ' + name)
        i, match = found[0]
        require(i > ci, 'Energy label precedes SCF convergence: ' + name)
        fields[name] = dict(line=i+1, text=lines[i], value=match.group(1))
    require(fields['cosx_printed']['line'] < fields['post_grid']['line']
            < fields['D4_printed']['line'] < fields['final']['line'], 'Native energy ordering differs')
    for item in fields.values():
        require(Decimal(item['value']).is_finite(), 'Nonfinite printed energy')
        item['half_quantum_hartree'] = str(half_quantum(item['value']))
    return dict(output=output, converged_cycle=int(cm.group(1)), fields=fields)


def scalar_closure(fields, coefficients):
    value = sum((Decimal(fields[k]['value']) * n for k, n in coefficients.items()), Decimal(0))
    bound = sum((half_quantum(fields[k]['value']) * abs(n) for k, n in coefficients.items()), Decimal(0))
    return dict(residual_hartree=str(value), printed_rounding_bound_hartree=str(bound),
                within_printed_rounding=abs(value) <= bound)


def fd(fields, coefficients):
    values, bounds = {}, {}
    for side in ('plus', 'minus'):
        f = fields[side]['fields']
        values[side] = sum((Decimal(f[k]['value']) * n for k, n in coefficients.items()), Decimal(0))
        bounds[side] = sum((half_quantum(f[k]['value']) * abs(n) for k, n in coefficients.items()), Decimal(0))
    native = (values['plus'] - values['minus']) / (2*H)
    return dict(coefficients=coefficients, native_hartree_per_A=str(native),
                kcal_mol_A=float(native*DEC_HA),
                printed_half_quantum_bound_kcal_mol_A=float(sum(bounds.values())*DEC_HA/(2*H)))


def receipt_evidence(row, task, cfg):
    receipt = read(row['receipt'])
    require(receipt['returncode'] == 0 and receipt['normal_termination'] and receipt['scf_converged']
            and receipt['orca_executable'] == cfg['orca'] and receipt['task_id'] == task['task_id'],
            'Native receipt identity/admission differs')
    require(receipt['artifacts']['output'] == row['output']
            and receipt['artifacts']['template_input'] == task['input']
            and receipt['artifacts']['xyz'] == task['xyz'], 'Native input/output receipt differs')
    for record in receipt['artifacts'].values():
        verify(record)
    sidecar = read(receipt['artifacts']['runtime_sidecar'])
    return dict(receipt=row['receipt'], memory_policy=sidecar['memory_policy'],
                started_at_utc=receipt['started_at_utc'], finished_at_utc=receipt['finished_at_utc'])


def diagnose(cfg, result):
    require(cfg['protocol'] == PROTOCOL and cfg['implementation']['sha256'] == pin(__file__)['sha256'],
            'Diagnostic protocol/source differs')
    verify(cfg['implementation']); verify(cfg['specification']); verify(cfg['python'])
    require(pin(sys.executable)['sha256'] == cfg['python']['sha256'], 'Interpreter differs')
    require(cfg['collector_result']['sha256'] == COLLECTOR_SHA and cfg['source_review']['sha256'] == REVIEW_SHA,
            'Frozen failed collector/source review differs')
    col = read(cfg['collector_result']); ccfg = read(col['config']); review = read(cfg['source_review'])
    require(review['authoritative_failure'] == cfg['collector_result'], 'Source review parent differs')
    for record in review['reviewed_source_pins']:
        verify(record)
    verify(col['implementation'])
    require(col['implementation']['sha256'] == QUALIFY_SHA and col['slurm_job_id'] == '13271'
            and col['stage'] == 'collect' and col['scope'] == 'scout_pair'
            and col['status'] == 'failed_derivative_gate' and col['scout_pair_passed'] is False
            and col['selected_tasks'] == list(PAIR) and len(col['checks']) == 1,
            'Original failed scout-pair identity differs')
    require(ccfg['settings'] == col['settings'] and col['settings']['h_A'] == float(H)
            and col['settings']['absolute_kcal_mol_A'] == .005 and col['settings']['relative'] == .0001,
            'Original numerical settings differ')
    prep = read(col['preparation_result']); origins = read(prep['origins']); origin = origins['La']
    require(col['preparation_result'] == ccfg['preparation_result'] and prep['status'] == 'prepared'
            and prep['slurm_job_id'] == '12976' and prep['new_molecular_calls'] == 0,
            'Original preparation identity differs')
    require(origin['status'] == 'complete' and read(origin['result']) == {k:v for k,v in origin.items() if k != 'result'},
            'Archived assembled origin differs')
    origin_parent = read(ccfg['origin_result'])
    require(origin_parent['status'] == 'complete' and origin_parent['slurm_job_id'] == '12849'
            and origin in origin_parent['rows'], 'Actual central-origin parent differs')
    fixture = read(prep['directions'])['MM_O97']; mapping = read(origin['physical_mapping'])
    require(fixture['id'] == SID and fixture['index'] == PI and len(fixture['direction']) == 3
            and origin['source_ids'][PI] == SID and len(origin['source_ids']) == 2049,
            'Frozen direction/physical identity differs')
    central = origin['reference_evidence']['reference_row']; original_task = origin['reference_evidence']['reference_task']
    for key in ('xyz', 'pointcharges', 'input', 'core_mapping', 'boundary_mapping'):
        verify(original_task[key])
    for key in ('output', 'receipt', 'engrad', 'pointcharge_gradient'):
        verify(central[key])
    require(original_task['electronic_state'] == central['electronic_state'] == origin['electronic_state'],
            'Central state differs')
    parsed = {'central': extract(central['output'])}
    resources = {'central': receipt_evidence(central, original_task, ccfg)}
    rows = {}
    original_pc = verify(original_task['pointcharges']).read_text().splitlines()
    require(len(original_pc) == 1841 and int(original_pc[0]) == 1840, 'Original field inventory differs')
    for side, tid in zip(('plus', 'minus'), PAIR):
        row = col['rows'][tid]; parent = read(row['result']); task = row['DFT']['task']; endpoint_cfg = read(parent['config'])
        require(row['status'] == parent['status'] == 'complete' and parent['DFT'] == row['DFT']
                and parent['native_energies'] == row['native_energies'] and parent['task_id'] == tid
                and row['result'] == ccfg['endpoint_results'][tid] and read(parent['task']) == task,
                'Endpoint/collector exact replay differs')
        for key in ('settings', 'origin_result', 'orca', 'basis_assets', 'implementation', 'preparation_result'):
            require(endpoint_cfg[key] == ccfg[key], 'Endpoint scientific binding differs: ' + key)
        require(task['origin_task'] == original_task and task['electronic_state'] == origin['electronic_state']
                and task['signed_h_A'] == (float(H) if side == 'plus' else -float(H))
                and task['fixture'] == 'MM_O97' and task['metal'] == 'La', 'Fixed endpoint/state differs')
        require(verify(task['xyz']).read_bytes() == verify(original_task['xyz']).read_bytes()
                and verify(task['core_mapping']).read_bytes() == verify(original_task['core_mapping']).read_bytes(),
                'Core/caps differ under MM-only displacement')
        require(verify(task['input']).read_text() == verify(original_task['input']).read_text().replace(
            ' VeryTightSCF EnGrad\n', ' VeryTightSCF\n', 1), 'Only declared EnGrad removal permitted')
        field = verify(task['pointcharges']).read_text().splitlines()
        require(len(field) == len(original_pc) and field[0] == original_pc[0]
                and [i for i,(a,b) in enumerate(zip(field[1:],original_pc[1:])) if a != b] == [FI]
                and all(a.split()[0] == b.split()[0] for a,b in zip(field[1:],original_pc[1:])),
                'External field change differs from sole displaced row832')
        parsed[side] = extract(row['DFT']['output'])
        require(parsed[side]['fields']['final']['value'] == row['DFT']['printed_energy'], 'Raw final energy string differs')
        resources[side] = receipt_evidence(row['DFT'], task, ccfg)
        rows[side] = row
    require(float(parsed['central']['fields']['final']['value']) == origin['native_energies']['PBE0_embedded_hartree'],
            'Central raw final energy differs')
    require(len({p['fields']['D4_printed']['value'] for p in parsed.values()}) == 1,
            'Printed fixed-core D4 values differ')
    for item in parsed.values():
        f = item['fields']
        item['printed_closures'] = {
            'post_minus_last_cycle_minus_printed_COSX': scalar_closure(f, {'post_grid':1,'last_cycle':-1,'cosx_printed':-1}),
            'final_minus_post_minus_printed_D4': scalar_closure(f, {'final':1,'post_grid':-1,'D4_printed':-1}),
            'new_minus_old_exchange_minus_printed_COSX': scalar_closure(f, {'exchange_new':1,'exchange_old':-1,'cosx_printed':-1})}
    result.update(raw_energy_fields=parsed, native_runtime_evidence=resources,
                  input_identity=dict(core_bytes_identical=True, core_mapping_bytes_identical=True,
                                      only_changed_PC_row=FI, charges_unchanged=True, printed_D4_identical=True),
                  preparation_result=col['preparation_result'], origins=prep['origins'], directions=prep['directions'],
                  original_collector_implementation=col['implementation'], original_collector_config=col['config'],
                  source_review=cfg['source_review'], original_failure=cfg['collector_result'])

    decomposition = {name:fd(parsed, terms) for name, terms in {
        'original_final':{'final':1}, 'last_cycle':{'last_cycle':1},
        'printed_COSX':{'cosx_printed':1}, 'net_post_SCF_shift':{'post_grid':1,'last_cycle':-1},
        'printed_D4':{'D4_printed':1}, 'final_minus_post':{'final':1,'post_grid':-1},
        'post_SCF_printed_COSX_remainder':{'post_grid':1,'last_cycle':-1,'cosx_printed':-1},
        'final_printed_D4_remainder':{'final':1,'post_grid':-1,'D4_printed':-1}}.items()}
    old = col['checks'][0]
    require(old['fixture'] == 'MM_O97' and old['metal'] == 'La' and old['passed'] is False, 'Failed comparison differs')
    reported = old['checks']['PBE0_embedded']; final_fd = decomposition['original_final']['kcal_mol_A']
    require(final_fd == reported['FD_kcal_mol_A']
            and final_fd-reported['analytic_kcal_mol_A'] == reported['residual_kcal_mol_A'],
            'Original final FD/residual failed exact replay')
    classical_fd = {
        'MM_ledger': (rows['plus']['native_energies']['MM_kcal_mol']-rows['minus']['native_energies']['MM_kcal_mol'])/.002,
        'fixed_FF_charge_ALPB': (rows['plus']['native_energies']['ALPB_hartree']-rows['minus']['native_energies']['ALPB_hartree'])*HA/.002}
    require(all(v == old['checks'][k]['FD_kcal_mol_A'] for k,v in classical_fd.items())
            and math.fsum([final_fd,*classical_fd.values()]) == old['checks']['total']['FD_kcal_mol_A'],
            'Original MM/ALPB/total FD differs')
    result.update(original_checks=old, original_FD_replay_exact=True, derivative_decomposition=decomposition,
                  diagnostic_residuals_kcal_mol_A={
                      'final_minus_printed_COSX_minus_analytic': final_fd-decomposition['printed_COSX']['kcal_mol_A']-reported['analytic_kcal_mol_A'],
                      'final_minus_net_post_SCF_shift_minus_analytic': final_fd-decomposition['net_post_SCF_shift']['kcal_mol_A']-reported['analytic_kcal_mol_A']},
                  decomposition_caveat='Net post-SCF shift is post-grid Total Energy minus last SCF cycle, not assumed to be pure COSX. Printed closure remainders remain explicit. These are diagnostic energy definitions only.')

    # This physical atom is reached once by a PC row and has no QM/cap role.
    require([a['id'] for a in mapping['physical_atoms']] == origin['source_ids']
            and len(mapping['field']) == 1840 and len(mapping['core']) == 209
            and [i for i,a in enumerate(mapping['field']) if a['id'] == SID] == [FI]
            and mapping['field'][FI]['source_index'] == PI and mapping['field'][FI]['charge_e'] == -.5679
            and not any(a.get('source_id') == SID for a in mapping['core'])
            and not any(PI in (a['retained_index'],a['omitted_index']) for a in mapping['caps']),
            'Unique point-charge-only physical mapping differs')
    require(mapping['xyz']['sha256'] == original_task['xyz']['sha256']
            and mapping['pointcharges']['sha256'] == original_task['pointcharges']['sha256'], 'Central core/field map differs')
    lines = verify(central['pointcharge_gradient']).read_text().splitlines()
    require(len(lines) == 1841 and int(lines[0]) == 1840, 'Native PC-gradient inventory differs')
    for line in lines[1:]:
        require(len(line.split()) == 3 and all(math.isfinite(float(v)) for v in line.split()), 'Invalid native PC gradient')
    tokens = lines[FI+1].split(); vector = [float(v)*(HA/BOHR) for v in tokens]
    stored = origin['component_gradients_kcal_mol_A']['PBE0_embedded'][PI]
    direction = fixture['direction']; projected = math.fsum(v*d for v,d in zip(vector,direction))
    roundoff = 64*sys.float_info.epsilon*max(1., math.fsum(abs(v*d) for v,d in zip(vector,direction)))
    vector_residual = [a-b for a,b in zip(vector,stored)]
    projection_residual = projected-reported['analytic_kcal_mol_A']
    matches = all(abs(x) <= roundoff for x in vector_residual) and abs(projection_residual) <= roundoff
    pc_quantum = sum((half_quantum(v)*abs(Decimal(str(d))) for v,d in zip(tokens,direction)), Decimal(0))*DEC_HA/Decimal(str(BOHR))
    result['direct_native_PC_replay'] = dict(
        mapping=origin['physical_mapping'], native_pcgrad=central['pointcharge_gradient'],
        source_id=SID, physical_index=PI, field_row=FI, file_line=FI+2, raw_line=lines[FI+1],
        native_gradient_Hartree_bohr=tokens, conversion='one multiplication by 627.509474/0.529177210903; no charge factor or sign reversal',
        physical_gradient_kcal_mol_A=vector, assembled_PBE0_gradient_kcal_mol_A=stored,
        component_residual_kcal_mol_A=vector_residual, direction=direction,
        projected_gradient_kcal_mol_A=projected, reported_analytic_kcal_mol_A=reported['analytic_kcal_mol_A'],
        projection_residual_kcal_mol_A=projection_residual, roundoff_only_replay_bound=roundoff,
        replay_matches=matches, projected_force_kcal_mol_A=-projected,
        printed_pcgrad_half_quantum_bound_kcal_mol_A=float(pc_quantum))
    result.update(status='complete' if matches else 'replay_mismatch',
                  qualification_status='failed_derivative_gate', scout_pair_passed=False,
                  unrun_endpoint_slots={k:v['status'] for k,v in col['rows'].items() if k not in PAIR},
                  remaining_six_singlepoints_authorized=False, response_references_authorized=False,
                  physical_state_qualified=False, global_force_qualified=False,
                  printed_rounding_is_SCF_error_bound=False, new_energy_or_gradient_definition_admitted=False,
                  interpretation='Saved decomposition only. Original final-energy PBE0 and total derivative gates remain failed; no promotion or numerical recovery follows.')


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--config', required=True, type=Path)
    ap.add_argument('--output', required=True, type=Path)
    args = ap.parse_args()
    require(os.environ.get('SLURM_JOB_ID','').isdigit(), 'Saved diagnosis requires Slurm')
    out = args.output.resolve(); out.mkdir(parents=True, exist_ok=False); tick = time.monotonic()
    result = dict(protocol=PROTOCOL, status='execution_failure', implementation=pin(__file__),
                  config=pin(args.config), slurm_job_id=os.environ['SLURM_JOB_ID'], command=sys.argv,
                  new_molecular_calls=0, new_DFT_calls=0, new_MACE_calls=0, new_GFN_calls=0,
                  new_MM_calls=0, new_ALPB_calls=0, optimization_steps=0)
    try:
        with localcontext() as context:
            context.prec = 50
            diagnose(read(result['config']), result)
    except Exception as exc:
        result.update(status='execution_failure', error=repr(exc), traceback=traceback.format_exc())
    result['application_wall_seconds'] = time.monotonic()-tick
    usage = resource.getrusage(resource.RUSAGE_SELF)
    result['measured_parent_CPU_seconds'] = usage.ru_utime+usage.ru_stime
    result['parent_peak_RSS_KiB'] = usage.ru_maxrss
    result['scheduler_cost_scope'] = 'Actual allocation receipt supplied by root; application time is not scheduler elapsed.'
    write(out/'RESULT.json', result)
    return 0 if result['status'] == 'complete' else 1


if __name__ == '__main__':
    raise SystemExit(main())
