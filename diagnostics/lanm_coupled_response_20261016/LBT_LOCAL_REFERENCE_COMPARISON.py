#!/usr/bin/env python3
"""Saved-only four-state MACE/GFN2/PBE0 source-response comparison."""
import argparse
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import sys
import time
import traceback

IDS = ('7CCO_La', '7CCN_La', '7CCO_Dy', '7CCN_Dy')
HA = 627.509474
EV = 23.06054783061903


def require(ok, message):
    if not ok:
        raise ValueError(message)


def pin(path):
    p = Path(path).resolve()
    return {'path': str(p), 'sha256': hashlib.sha256(p.read_bytes()).hexdigest()}


def verify(record):
    p = Path(record['path'])
    require(pin(p)['sha256'] == record['sha256'], 'Changed pin: ' + str(p))
    return p


def read(record):
    return json.loads(verify(record).read_text())


def module(name, record):
    spec = importlib.util.spec_from_file_location(name, verify(record))
    obj = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(obj)
    return obj


def write(path, value):
    with Path(path).open('x') as handle:
        json.dump(value, handle, indent=2, allow_nan=False)
        handle.write('\n')


def reference_cell(tid, result_pin, cfg, refcfg, reference, model, collector, residual):
    row = {'task_id': tid, 'status': 'unavailable', 'energy_hartree': None, 'result': result_pin}
    if result_pin is None:
        return dict(row, reason='No completed reference result pinned for this declared state')
    try:
        result = read(result_pin); actual_cfg = read(result['config'])
        require(result['status'] == 'complete' and result['preflight'] == 'passed'
                and result['stage'] in ('scout', 'remaining') and tid in result['selected_tasks'],
                'Run not complete or cell not an owned selected task')
        require((result['stage'] == 'scout') == (tid == '7CCO_La'), 'Scout/remaining ownership differs')
        for key in ('plan', 'compact_manifest', 'preparation_result', 'MACE_result', 'tasks', 'python', 'health_monitor'):
            require(actual_cfg[key] == refcfg[key], 'Reference scientific inputs differ: ' + key)
        verify(result['implementation']); verify(actual_cfg['implementation'])
        require(result['implementation']['sha256'] == actual_cfg['implementation']['sha256']
                == cfg['reference_adapter']['sha256'], 'Reference implementation differs')
        manifest = read(result['manifest'])
        require({t['task_id'] for t in manifest['tasks']} == set(result['selected_tasks'])
                and manifest['inputs_config'] == result['config'], 'Executed task/config set differs')
        if result['stage'] == 'remaining':
            require(actual_cfg['scout_result'] == result['scout_result'] == cfg['reference_results']['7CCO_La'],
                    'Remaining scout lineage differs')
            reference.scout_gate(actual_cfg, model, collector, residual)
        native = reference.admit(verify(result['manifest']), model, collector, residual)[tid]
        require(native == result['rows'][tid] and native['status'] == 'available', 'Raw native replay differs')
        require(native['task']['source_task'] == refcfg['tasks'][tid]
                and native['task']['xyz']['sha256'] == refcfg['tasks'][tid]['xyz']['sha256'], 'Reference source/state differs')
        row.update(status='available', energy_hartree=native['energy_hartree'],
                   electronic_state=native['electronic_state'], SCF_evidence=native['SCF_evidence'],
                   receipt=native['receipt'], output=native['output'], engrad=native['engrad'],
                   molecular_wall_seconds=native['molecular_wall_seconds'],
                   runtime_memory_policy=native['runtime_memory_policy'])
    except Exception as exc:
        row.update(reason=repr(exc))
    return row


def strains(energies, factor):
    native = {}
    for metal in ('La', 'Dy'):
        values = [energies['7CCN_' + metal], energies['7CCO_' + metal]]
        native[metal] = None if None in values else math.fsum([values[0], -values[1]])
    differential = None if None in energies.values() else math.fsum([
        energies['7CCN_Dy'], -energies['7CCO_Dy'], -energies['7CCN_La'], energies['7CCO_La']])
    return {'definition': 'S_M=E_M(7CCN)-E_M(7CCO)',
            'La_kcal_mol': None if native['La'] is None else native['La'] * factor,
            'Dy_kcal_mol': None if native['Dy'] is None else native['Dy'] * factor,
            'Dy_minus_La_kcal_mol': None if differential is None else differential * factor,
            'same_metal_native_differences': native, 'balanced_native_difference': differential}


def compare(cfg, result):
    verify(cfg['plan']); require(pin(__file__)['sha256'] == cfg['implementation']['sha256'], 'Comparator differs')
    require(set(cfg['reference_results']) == set(IDS), 'Four reference slots required')
    refcfg = read(cfg['reference_config'])
    require(refcfg['plan'] == cfg['plan'] and refcfg['MACE_result'] == cfg['MACE_result']
            and refcfg['implementation']['sha256'] == cfg['reference_adapter']['sha256'], 'Reference parent differs')
    reference = module('local_reference_saved_parent', cfg['reference_adapter'])
    compact, model, runner, collector, residual, runtime, states = reference.inputs(refcfg)
    gfn_result = read(cfg['GFN_result']); gcfg = read(gfn_result['config'])
    verify(gfn_result['implementation']); verify(gcfg['implementation'])
    require(gfn_result['status'] == 'complete' and gfn_result['slurm_job_id'] == '12825'
            and gfn_result['preflight'] == 'passed' and gfn_result['plan'] == cfg['plan']
            and gcfg['MACE_result'] == gfn_result['MACE_result'] == cfg['MACE_result']
            and gfn_result['implementation']['sha256'] == gcfg['implementation']['sha256'], 'Actual12825 parent differs')
    gfn = module('local_reference_saved_gfn', gfn_result['implementation'])
    world = gfn.inputs(gcfg)
    require(world['preparation'] == refcfg['preparation_result'], 'GFN/MACE/reference source differs')
    require(len(gfn_result['matrix']) == 4 and
            {r['source_id'] + '_' + r['metal'] for r in gfn_result['matrix']} == set(IDS), 'GFN matrix denominator differs')
    matrix = []
    for cell in gfn_result['matrix']:
        source, metal = cell['source_id'], cell['metal']; tid = source + '_' + metal
        endpoint = world['endpoints'][(source, metal)]
        require(cell['endpoint'] == endpoint and cell['GFN_status'] == 'accepted'
                and cell['GFN_reused'] == (source == '7CCO'), 'GFN matrix source/status differs')
        require({k: v for k, v in endpoint.items() if k != 'cache_key'} == refcfg['tasks'][tid], 'Reference/GFN input differs')
        saved = read(cell['GFN_result']); task = saved['task']
        for artifact in saved['artifacts'].values():
            verify(artifact)
        require(read(task['task_record']) == {k: v for k, v in task.items() if k != 'task_record'},
                'GFN task record differs')
        gfn.task_check(task, endpoint, saved['expected'])
        parsed = gfn.admit(Path(task['directory']), task, read(saved['execution']), saved['expected'],
                           world['backend'], world['tolerances'])
        require(parsed['energy_hartree'] == saved['energy_hartree'] == cell['GFN_energy_hartree']
                and parsed['charges_e'] == saved['charges_e']
                and parsed['standalone_observation'] == saved['standalone_observation'], 'GFN raw replay differs')
        mace = world['MACE'][(source, metal)]
        require(mace['energy_eV'] == cell['MACE_energy_eV']
                and mace['input_state_check'] == cell['MACE_native_state'], 'MACE matrix differs')
        dft = reference_cell(tid, cfg['reference_results'][tid], cfg, refcfg, reference, model, collector, residual)
        matrix.append({'task_id': tid, 'source_id': source, 'metal': metal, 'xyz': endpoint['xyz'],
                       'charge': -1, 'physical_multiplicity': endpoint['multiplicity'],
                       'MACE_energy_eV': mace['energy_eV'], 'GFN2_energy_hartree': saved['energy_hartree'],
                       'GFN2_result': cell['GFN_result'], 'GFN2_effective_multiplicity': 1, 'PBE0': dft})
    summaries = {}
    for method, factor in (('MACE', EV), ('GFN2', HA), ('PBE0', HA)):
        energies = {r['task_id']: (r['PBE0']['energy_hartree'] if method == 'PBE0'
                    else r[method + ('_energy_eV' if method == 'MACE' else '_energy_hartree')]) for r in matrix}
        summaries[method] = strains(energies, factor)
    for method, oldkey in (('MACE', 'MACE'), ('GFN2', 'GFN2')):
        require(all(summaries[method][key] == gfn_result['source_strains'][oldkey][key]
                    for key in ('La_kcal_mol', 'Dy_kcal_mol', 'Dy_minus_La_kcal_mol')), '12825 source strains do not replay')
    errors = {}
    for method in ('MACE', 'GFN2'):
        errors[method] = {}
        for key in ('La_kcal_mol', 'Dy_kcal_mol', 'Dy_minus_La_kcal_mol'):
            value, ref = summaries[method][key], summaries['PBE0'][key]
            delta = None if value is None or ref is None else math.fsum([value, -ref])
            errors[method][key] = {'signed_method_minus_PBE0': delta,
                'absolute_difference': None if delta is None else abs(delta),
                'within_1kcal_interpretive_resolution': None if delta is None else abs(delta) <= 1.0}
    result.update(matrix=matrix, source_strains=summaries, comparison_to_PBE0=errors,
        status='complete' if all(r['PBE0']['status'] == 'available' for r in matrix) else 'unavailable',
        available_PBE0_cells=sum(r['PBE0']['status'] == 'available' for r in matrix),
        complete_MACE_GFN_replay=True, interpretation_resolution_kcal_mol=1.0,
        interpretation='Conditional48atom vacuum source response; reference is approximate. The1kcal scale is not a biological threshold/error bar. Native differences precede unit conversion; no raw unlike-metal affinity.',
        production_default_changed=False, PQQ_changed=False, physical_Dy_state_qualified=False,
        force_consistency_qualified=False, automatic_search_extension=False)


def report(result):
    lines = ['# Local LBT48 electronic reference comparison', '', 'Status: ' + result['status'] + '.', '',
             'Same-metal source strain S = E(7CCN) − E(7CCO), kcal/mol. No new molecular calls.', '',
             '| Method | La S | Dy S | S(Dy) − S(La) |', '|---|---:|---:|---:|']
    fmt = lambda value: 'unavailable' if value is None else f'{value:.6f}'
    for method, row in result.get('source_strains', {}).items():
        lines.append('| ' + method + ' | ' + ' | '.join(fmt(row[k]) for k in ('La_kcal_mol','Dy_kcal_mol','Dy_minus_La_kcal_mol')) + ' |')
    lines += ['', '| Method minus PBE0 | La strain error | Dy strain error | Differential error |',
              '|---|---:|---:|---:|']
    for method, row in result.get('comparison_to_PBE0', {}).items():
        lines.append('| ' + method + ' | ' + ' | '.join(fmt(row[k]['signed_method_minus_PBE0'])
                     for k in ('La_kcal_mol','Dy_kcal_mol','Dy_minus_La_kcal_mol')) + ' |')
    lines += ['', '| Reference state | Availability | All printed residuals pass |', '|---|---|---|']
    for row in result.get('matrix', []):
        dft = row['PBE0']
        lines.append(f"| {row['task_id']} | {dft['status']} | {dft.get('SCF_evidence', {}).get('all_printed_residuals_pass', 'unavailable')} |")
    lines += ['', 'PBE0 admits native convergence plus the existing state/gradient gates. Printed residual failures remain visible in RESULT.json; producing a gradient does not qualify its derivatives or electronic stability.', '',
              'MACE uses physical La1/Dy6 metadata. GFN2 and the frozen-f PBE0 Dy model use effective closed valence states; spin–orbit physics is unqualified. This finite comparison does not change a default, establish affinity, or select a new search.', '',
              'Exact parents, all four availability records, raw receipts and continuous method-minus-PBE0 differences are retained in RESULT.json.']
    return '\n'.join(lines) + '\n'


def run(args):
    require(os.environ.get('SLURM_JOB_ID', '').isdigit(), 'Slurm required for saved scientific comparison')
    require(int(os.environ.get('SLURM_CPUS_PER_TASK', '0')) == 1, 'One shared CPU required')
    output = Path(args.output).resolve(); output.mkdir(parents=True, exist_ok=False)
    cfg = json.loads(Path(args.config).read_text()); start = time.monotonic(); cpu = time.process_time()
    result = {'status': 'execution_failure', 'protocol_id': 'lbt3_local_reference_saved_comparison_v1',
              'config': pin(args.config), 'implementation': pin(__file__), 'plan': cfg['plan'],
              'slurm_job_id': os.environ['SLURM_JOB_ID'], 'declared_states': list(IDS),
              'MACE_result': cfg['MACE_result'], 'GFN_result': cfg['GFN_result'],
              'reference_results': cfg['reference_results'], 'new_molecular_calls': 0, 'new_geometries': 0}
    try:
        compare(cfg, result)
    except Exception as exc:
        result.update(status='execution_failure', failure=repr(exc), traceback=traceback.format_exc())
    result.update(elapsed_seconds=time.monotonic()-start, process_CPU_seconds=time.process_time()-cpu)
    write(output / 'RESULT.json', result)
    (output / 'REPORT.md').write_text(report(result))
    return 0 if result['status'] == 'complete' else 1


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', required=True); parser.add_argument('--output', required=True)
    raise SystemExit(run(parser.parse_args()))
