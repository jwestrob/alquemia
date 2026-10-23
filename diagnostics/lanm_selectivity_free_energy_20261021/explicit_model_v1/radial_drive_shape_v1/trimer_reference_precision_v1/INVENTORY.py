#!/usr/bin/env python3
"""Fourteen exact saved W9far reference states; metadata only, no native imports."""
import hashlib
import json
import os
from pathlib import Path
import re
import socket
import time

HERE = Path(__file__).resolve().parent
ROOT = next(p for p in HERE.parents if (p / 'diagnostics').is_dir())
OUT = ROOT / 'workspaces/lanm_selectivity_free_energy_20261021/explicit_model_v1/radial_drive_shape_v1/trimer_reference_precision_v1/inventory_v1'
CHECKS = []
PINS = {}


def require(test, label):
    if not test:
        raise ValueError(label)
    CHECKS.append(label)


def pin(path):
    p = Path(path).resolve()
    return {'path': str(p), 'sha256': hashlib.sha256(p.read_bytes()).hexdigest()}


def verify(p):
    key = (p['path'], p['sha256'])
    if key not in PINS:
        require(pin(p['path']) == p, 'pin:' + p['path'])
        PINS[key] = True
    return Path(p['path'])


def read(p):
    return json.loads(verify(p).read_text())


def save(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')


def trace(text):
    lines = text.splitlines()
    evidence = []
    iterations = []
    in_soscf = False
    for n, line in enumerate(lines, 1):
        if 'Initializing SOSCF' in line:
            in_soscf = True
        tokens = line.split()
        if in_soscf and len(tokens) == 7 and tokens[0].isdigit() and 'e' in tokens[2].lower():
            iterations.append({'line_number': n, 'text': line, 'iteration': tokens[0],
                               'delta_E': tokens[2], 'RMSDP': tokens[3],
                               'MaxDP': tokens[4], 'MaxGrad': tokens[5]})
        if any(k in line for k in ('Convergence check mode', 'signals convergence',
                                    'SCF CONVERGED AFTER', 'Initializing SOSCF',
                                    'Energy change', 'MAX-Density change', 'RMS-Density change',
                                    'DIIS Error', 'Orbital Gradient', 'Orbital Rotation',
                                    'Recomputing exchange energy using gridx3',
                                    '# of basis functions in Aux-J')):
            evidence.append({'line_number': n, 'text': line})
        if 'SCF CONVERGED AFTER' in line:
            in_soscf = False
    return {'lines': evidence, 'first_SOSCF_iteration': iterations[0] if iterations else None,
            'last_SOSCF_iteration': iterations[-1] if iterations else None,
            'scope': 'Literal trace extraction only; no energy/moment or SCF replay.'}


def main():
    start, cpu = time.monotonic(), time.process_time()
    OUT.mkdir(parents=True, exist_ok=False)
    result = {'status': 'started', 'implementation': pin(__file__), 'scope': pin(HERE / 'SCOPE.json'),
              'new_reference_calls': 0, 'new_model_calls': 0, 'new_solver_calls': 0,
              'scientific_cycle_recomputations': 0, 'error': None}
    try:
        sources = json.loads((HERE / 'SOURCE_INPUTS.json').read_text())
        for p in sources.values():
            verify(p)
        authority = read(sources['root_decision'])
        require(authority['source_only'] is True and authority['new_numerical_execution_release'] is None,
                'Source-only root authority')
        original = read(sources['original_manifest'])
        comparison = read(sources['original_comparison'])
        consumed = read(sources['consumed_inventory'])
        panel = read(original['sources']['common_panel'])
        for p in (original['orca'], original['sources']['common_adapter']):
            verify(p)
        for assets in original['reference_assets'].values():
            for p in assets.values():
                verify(p)
        cycles, rows = {}, []
        for metal in ('La', 'Dy'):
            cycle_id = f'aquo9_central_{metal}__far'
            cycle = consumed['cycles'][cycle_id]
            require(cycle == comparison['cycles'][cycle_id], cycle_id + ': unchanged literal cycle')
            terms = cycle['reference_terms']
            expected = {
                f'pbe0__aquo9__{metal}__far__metal_two_water': 1,
                f'pbe0__aquo9__{metal}__far__two_water': -1,
                f'pbe0__aquo9__{metal}__water_01__metal_water': -1,
                f'pbe0__aquo9__{metal}__water_03__metal_water': -1,
                f'pbe0__aquo9__{metal}__ion_ghost_shell': 1,
                f'pbe0__aquo9__{metal}__water_01__water_only': 1,
                f'pbe0__aquo9__{metal}__water_03__water_only': 1}
            require(terms == expected and cycle['selection']['water_indices'] == [1, 3],
                    cycle_id + ': exact seven terms')
            cycles[cycle_id] = {k: v for k, v in cycle.items() if k.startswith('reference_') or k == 'selection'}
            for cid, coefficient in terms.items():
                fresh_historical = cid in comparison['inputs']['reference_results']
                record = None if fresh_historical else original['reused_reference'][cid]['record']
                rp = comparison['inputs']['reference_results'][cid] if fresh_historical else record['result']
                r = read(rp)
                parsed = r['parsed'] if fresh_historical else r
                cell = r['cell']
                moment = comparison['quantum_moment_rows'][cid]
                require(r['status'] == 'complete' and all(parsed['checks'].values()), cid + ': inherited admitted status')
                require(r['cell_id'] == cell['cell_id'] == cid and cell['metal'] == metal and
                        cell['method'] == 'pbe0' and cell['d4'] is False and
                        len(cell['symbols']) == len(cell['positions_A']) == len(cell['ghost_mask']) == 28 and
                        cell['expected_basis_functions'] == 434 and cell['multiplicity'] == 1,
                        cid + ': literal parent/state/basis identity')
                require(cell['real_parent_atom_indices'] == [i for i, g in enumerate(cell['ghost_mask']) if not g],
                        cid + ': ghost mask identity')
                require(moment['output'] == r['output'] and moment['ghost_mask'] == cell['ghost_mask'] and
                        moment['labels'] == cell['symbols'] and moment['charge'] == cell['charge'],
                        cid + ': inherited common-origin moment joins')
                for p in (r['input'], r['output'], r['stderr'], cell['source_geometry'],
                          cell['template_input'], moment['property_file'], r['adapter'], r['release']):
                    verify(p)
                deck = verify(r['input']).read_text()
                scf = re.findall(r'(?im)^%scf\s*\n(.*?)^end\s*$', deck, re.S)
                require(scf == [' Guess PModel\n MaxIter 300\n'], cid + ': unchanged original SCF block')
                require(deck.startswith('! PBE0 def2-TZVP def2/J RIJCOSX DefGrid3 VeryTightSCF NoAutostart\n') and
                        '%maxcore 3000\n%pal nprocs 4 end\n' in deck,
                        cid + ': actual method/resource header')
                require(parsed['native_scf_mode'] == ['Total+1el-Energy'] and
                        parsed['D4_energy_hartree'] == 0.0, cid + ': original stopping mode and no D4')
                body = re.sub(r'(?im)^%scf\s*\n.*?^end\s*$', '%scf\n<SCF_BLOCK>\nend', deck, count=1, flags=re.S)
                raw = verify(r['output']).read_text()
                aux = re.findall(r'# of basis functions in Aux-J\s*\.\.\.\s*(\d+)', raw)
                require(aux == [str(990 if metal == 'La' else 1008)], cid + ': literal Aux-J count')
                evidence = trace(raw)
                require(evidence['first_SOSCF_iteration'] is not None and evidence['last_SOSCF_iteration'] is not None,
                        cid + ': saved SOSCF trace present')
                sp = r.get('start') or (record or {}).get('start')
                start_record = read(sp) if sp else None
                metadata_keys = ('checks', 'energy_hartree', 'electronic_energy_hartree', 'D4_energy_hartree',
                                 'printed_scf_criteria', 'native_scf_mode', 'native_scf_stop_signals',
                                 'final_COSX_correction_hartree', 'total_energy_after_final_integration_hartree',
                                 'energy_stage_metadata', 'physical_state_qualified', 'method_accuracy_qualified')
                rows.append({'cell_id': cid, 'cycle_id': cycle_id, 'coefficient': coefficient,
                             'historical_producer': 'trimer16' if fresh_historical else 'previous_fragment',
                             'result': rp, 'input': r['input'], 'output': r['output'], 'stderr': r['stderr'],
                             'property_file': moment['property_file'], 'adapter': r['adapter'], 'release': r['release'],
                             'start': sp, 'cell': cell, 'old_native_metadata': {k: parsed[k] for k in metadata_keys},
                             'old_moment': moment, 'input_header_and_scf': deck[:deck.index('%basis')],
                             'input_with_SCF_block_masked_sha256': hashlib.sha256(body.encode()).hexdigest(),
                             'raw_Aux_J_count': aux, 'SCF_trace': evidence,
                             'old_cost': {k: r[k] for k in ('elapsed_seconds', 'allocated_cpu_seconds',
                                          'allocated_elapsed_CPU_seconds_proxy', 'job_id', 'ranks') if k in r},
                             'old_runtime': None if start_record is None else {k: start_record[k] for k in
                                          ('command', 'ranks', 'memory_MiB', 'maxcore_MiB', 'job_id', 'host', 'threads') if k in start_record},
                             'stricter_protocol_reuse_eligible': False,
                             'reuse_reason': 'Actual old protocol is mode2/VeryTightSCF; no admitted identical stricter-protocol result is bound.'})
        require(len(rows) == len({r['cell_id'] for r in rows}) == 14, 'Exactly fourteen distinct states')
        require(len({r['input']['sha256'] for r in rows}) == 14, 'Fourteen distinct literal input decks')
        require(sum(r['historical_producer'] == 'trimer16' for r in rows) == 4, 'Four trimer/pair plus ten earlier fragment results')
        inv = {'status': 'inventoried_exact_W9far_reference_precision_source_only',
               'implementation': result['implementation'], 'sources': sources, 'cycles': cycles, 'rows': rows,
               'orca': original['orca'], 'reference_assets': original['reference_assets'],
               'parent_panel': original['sources']['common_panel'],
               'parent_runtime_and_parser_authorities': {k: v for k, v in panel.items() if k in
                   ('adapter', 'configuration', 'declaration', 'cp_panel', 'cp_final_admission', 'pair_final_admission')},
               'counts': {'cycles': 2, 'distinct_reference_states': 14, 'exact_stricter_protocol_reuse_bound': 0,
                          'proposed_maximum_new_states': 14, 'new_reference_model_solver_calls': 0},
               'reuse_scope': 'Exact fourteen-state lineage and its admitted inputs examined; no global claim that a tighter calculation cannot exist elsewhere. No result meets a yet-unfrozen new protocol by inference.',
               'old_reference_preserved': True, 'actual_stricter_protocol': None, 'execution_release': None}
        save(OUT / 'INVENTORY.json', inv)
        result.update(status='passed_reference_precision_metadata_inventory', inventory=pin(OUT / 'INVENTORY.json'))
    except Exception as exc:
        result.update(status='failed_reference_precision_metadata_inventory', error=repr(exc))
    result.update(checks=CHECKS, check_count=len(CHECKS), verified_pin_count=len(PINS),
                  cost={'wall_seconds': time.monotonic()-start, 'CPU_seconds': time.process_time()-cpu,
                        'host': socket.gethostname(), 'PID': os.getpid(), 'numerical_threads': 1,
                        'scope': 'stdlib metadata/hash/text inspection only'})
    save(OUT / 'RESULT.json', result)
    print(json.dumps({k: result[k] for k in ('status', 'error', 'check_count', 'verified_pin_count', 'cost')}))
    return 0 if result['error'] is None else 1


if __name__ == '__main__':
    raise SystemExit(main())
