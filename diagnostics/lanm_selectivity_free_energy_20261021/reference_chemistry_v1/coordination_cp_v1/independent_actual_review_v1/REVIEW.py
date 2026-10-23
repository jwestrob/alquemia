#!/usr/bin/env python3
"""Read-only native replay and independent raw-token CP cycle arithmetic."""
import datetime
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
ROOT = HERE.parents[4]
WS = ROOT/'workspaces/lanm_selectivity_free_energy_20261021/reference_chemistry_v1/coordination_cp_v1'
OUT = WS/'independent_actual_review_v1'
HA = Decimal('627.509474')
GIDS = [f'aquo{n}_central_{m}' for n in (8, 9) for m in ('La', 'Dy')]
ROLES = ('ion_ghost_shell', 'shell_ghost_ion')
IDS = [f'pbe0__aquo{n}__{m}__{role}' for n in (8, 9) for m in ('La', 'Dy') for role in ROLES]
SIGNS = dict(zip(GIDS, (1, -1, -1, 1)))
F = r'[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][-+]?\d+)?'
checks = []


def check(ok, label):
    if not ok:
        raise AssertionError(label)
    checks.append(label)


def pin(path):
    p = Path(path).resolve()
    return {'path': str(p), 'sha256': hashlib.sha256(p.read_bytes()).hexdigest()}


def verify(rec):
    check(pin(rec['path']) == rec, 'pin:' + rec['path'])
    return Path(rec['path'])


def read(rec):
    return json.loads(verify(rec).read_text())


def load_module(rec):
    spec = importlib.util.spec_from_file_location('independent_actual_cp_native', verify(rec))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def terms(raw, d4_expected):
    text = raw.split('****END OF INPUT****')[-1]
    patterns = {'electronic': r'^Total Energy\s*:\s*(' + F + r')\s+Eh',
                'D4': r'^Dispersion correction\s+(' + F + r')',
                'final': r'FINAL SINGLE POINT ENERGY\s+(' + F + r')'}
    strings = {k: re.findall(v, text, re.M) for k, v in patterns.items()}
    check(len(strings['electronic']) == len(strings['final']) == 1, 'unique raw electronic/final terms')
    check(len(strings['D4']) == (1 if d4_expected else 0), 'raw D4 scope')
    result = {k: Decimal(v[0].replace('D', 'E').replace('d', 'e')) if v else Decimal(0) for k, v in strings.items()}
    check(all(v.is_finite() for v in result.values()), 'raw finite terms')
    check(abs(result['electronic'] + result['D4'] - result['final']) < Decimal('1e-8'), 'raw electronic/D4 closure')
    return result, strings


def delta(values):
    return sum((SIGNS[g] * values[g] for g in GIDS), Decimal(0))


def save(path, value):
    with Path(path).open('x') as h:
        json.dump(value, h, indent=2, allow_nan=False)
        h.write('\n')


def main():
    tick, cpu = time.monotonic(), time.process_time()
    check(all(os.environ.get(k) == '1' for k in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS')), 'one-thread environment')
    OUT.mkdir(exist_ok=False)
    result_path = WS/'collection/run_v1/result/RESULT.json'
    result = json.loads(result_path.read_text())
    for k in ('implementation', 'config', 'panel', 'root_inputs', 'consumed_input_data'):
        verify(result[k])
    check(result['status'] == 'coordination_cp_collection_complete' and result['fragment_complete'] == 8 and result['original_complete'] == 4, 'complete exact 8+4 result')
    check(set(result['fragments']) == set(IDS) and set(result['original_complexes']) == set(GIDS), 'exact slot inventories')
    check(result['new_native_calls'] == 0 and not result['fit_targets_changed'] and not result['model_selected'], 'collector zero-call no target or model change')
    panel = read(result['panel'])
    adapter = load_module(panel['adapter'])
    loaded = adapter.load(result['panel']['path'])
    check(loaded == panel, 'frozen native load and all helper/asset checks')
    data = read(result['consumed_input_data'])
    original_panel = adapter.b.load_panel(verify(panel['parent_panel']))
    scheduler = {}
    for kind, path in [('native', WS/'run_v1/completion_event.json'), ('collector', WS/'collection/run_v1/completion_event.json')]:
        raw = path.read_bytes()
        snapshot = OUT/(kind + '_completion_snapshot.json')
        snapshot.write_bytes(raw)
        scheduler[kind] = {'source_path': str(path.resolve()), 'snapshot': pin(snapshot), 'data': json.loads(raw)}
    fragments, complexes, raw_f, raw_c = {}, {}, {}, {}
    native_allocated = 0.0
    for c in panel['cells']:
        cid = c['cell_id']; published = result['fragments'][cid]
        native = read(published['result']); start = read(published['start'])
        for k in ('input', 'output', 'stderr', 'release'):
            verify(native[k])
        release = read(native['release'])
        check(native['status'] == 'complete' and native['returncode'] == 0 and all(native['checks'].values()), cid + ':native admission')
        check(native['cell'] == start['cell'] == c and native['panel'] == start['panel'] == result['panel'] and native['adapter'] == panel['adapter'], cid + ':physical source join')
        check(start['input'] == native['input'] == published['input'] and start['release'] == native['release'] == published['release'], cid + ':input/release join')
        check(start['job_id'] == native['job_id'] == published['job_id'], cid + ':job join')
        check(release['status'] == 'approved_coordination_cp_eight_cells' and release['implementation'] == panel['adapter'] and release['panel'] == result['panel'], cid + ':exact root release')
        check(release['cell_ids'] == IDS and release['maximum_native_calls'] == 8 and release['ranks'] == 4, cid + ':released count/resources')
        for k in ('preparation', 'independent_review'):
            verify(release[k])
        check(Path(published['result']['path']).resolve() == Path(release['execution_root']).resolve()/'tasks'/cid/'execution/RESULT.json', cid + ':one-attempt result path')
        check(start['ranks'] == 4 and start['maxcore_MiB'] == math.floor(.75 * start['memory_MiB']/4), cid + ':actual memory/ranks')
        check(start['maximum_native_calls'] == native['native_program_invocations_attempted'] == 1, cid + ':one program')
        check(Path(native['input']['path']).read_text() == adapter.deck(panel, c, 4, start['maxcore_MiB']), cid + ':exact actual ghost deck')
        text = Path(native['output']['path']).read_text()
        parsed = adapter.parse(text, c, native['returncode'])
        check(all(native[k] == v for k, v in parsed.items()) and parsed['status'] == 'complete', cid + ':saved native parser replay')
        for k in ('checks', 'electronic_energy_hartree', 'D4_energy_hartree', 'printed_scf_criteria', 'native_scf_mode', 'native_scf_stop_signals', 'energy_stage_metadata', 'final_COSX_correction_hartree', 'total_energy_after_final_integration_hartree'):
            check(published[k] == parsed[k], cid + ':collector ' + k)
        dec, tokens = terms(text, False)
        check(abs(float(dec['electronic']) - published['electronic_energy_hartree']) < 1e-12, cid + ':independent electronic token')
        raw_f[cid] = dec
        record = scheduler['native']['data']['records'][native['job_id']]
        check(record['state'] == 'COMPLETED' and record['exit_code'] == '0:0' and int(record['allocated_cpus']) == 4, cid + ':scheduler terminal')
        check(abs(native['allocated_cpu_seconds'] - 4 * native['elapsed_seconds']) < 1e-8, cid + ':native allocation proxy')
        native_allocated += native['allocated_cpu_seconds']
        fragments[cid] = {'result': published['result'], 'raw_tokens': tokens, 'job_id': native['job_id'], 'native_allocated_CPU_seconds': native['allocated_cpu_seconds'],
                          'scheduler': record, 'native_stopping': {k: parsed[k] for k in ('printed_scf_criteria', 'native_scf_mode', 'native_scf_stop_signals', 'final_COSX_correction_hartree', 'total_energy_after_final_integration_hartree')},
                          'false_printed_scf_criteria': [v['name'] for v in parsed['printed_scf_criteria'] if not v['passes']]}
    for gid in GIDS:
        published = result['original_complexes'][gid]; native = read(published['result'])
        cell = next(c for c in original_panel['cells'] if c['cell_id'] == 'ecp_pbe0__' + gid)
        check(native['panel'] == panel['parent_panel'] and native['adapter'] == panel['base_helper'], gid + ':original producer')
        for k in ('output', 'input', 'engrad', 'geometry'):
            verify(native[k])
        parsed = adapter.b.collect_one(original_panel, cell, Path(published['result']['path']).parent, native['returncode'])
        check(parsed['status'] == 'complete' and all(native[k] == v for k, v in parsed.items()), gid + ':original raw replay')
        check(native['energy_hartree'] * float(HA) == data['rows'][gid]['reference']['energy'] and native['forces_kcal_mol_A'] == data['rows'][gid]['reference']['forces'], gid + ':consumed source')
        text = Path(native['output']['path']).read_text(); dec, tokens = terms(text, True); raw_c[gid] = dec
        for k, rawkey in [('electronic_energy_hartree', 'electronic'), ('D4_energy_hartree', 'D4'), ('final_energy_hartree', 'final')]:
            check(float(dec[rawkey]) == published[k] == panel['saved_complexes'][gid][k], gid + ':original ' + k)
        check(published['source_geometry']['sha256'] == published['consumed_geometry']['sha256'], gid + ':same geometry bytes')
        verify(published['source_geometry']); verify(published['consumed_geometry'])
        complexes[gid] = {'result': published['result'], 'raw_tokens': tokens, 'energy_stages': published['energy_stages'], 'printed_scf_criteria': native['printed_scf_criteria'], 'reused_allocated_CPU_seconds': native['allocated_cpu_seconds']}
    for n in (8, 9):
        la, dy = (panel['geometries'][f'aquo{n}_central_{m}'] for m in ('La', 'Dy'))
        check(la['positions_A'] == dy['positions_A'] and la['symbols'][1:] == dy['symbols'][1:] == ['O', 'H', 'H'] * n, f'W{n}:own-shell exact cancellation')
    for m in ('La', 'Dy'):
        gs = [panel['geometries'][f'aquo{n}_central_{m}'] for n in (8, 9)]
        inp = [panel['saved_complexes'][f'aquo{n}_central_{m}']['input'] for n in (8, 9)]
        check(inp[0]['sha256'] == inp[1]['sha256'] and gs[0]['positions_A'][0] == gs[1]['positions_A'][0] == [0., 0., 0.], m + ':own-ion same recipe and origin')
    check(all(result['cancellation_checks'].values()), 'all published cancellation checks')
    G = {g: sum((raw_f[f'pbe0__aquo{g[4]}__{g.rsplit("_", 1)[1]}__{role}']['electronic'] for role in ROLES), Decimal(0)) for g in GIDS}
    de = delta({g: raw_c[g]['electronic'] for g in GIDS}); dd = delta({g: raw_c[g]['D4'] for g in GIDS}); dt = delta({g: raw_c[g]['final'] for g in GIDS}); dg = delta(G)
    h = dict(original_electronic=de, original_D4_unchanged=dd, original_total=dt, delta_G_electronic=dg, CP_shift=-dg, CP_electronic=de-dg, CP_total=dt-dg)
    model = delta({g: Decimal(str(data['baseline_model'][g]['energy'])) for g in GIDS})
    kcal = {k: v*HA for k, v in h.items()}; kcal.update(unchanged_model=model, model_minus_original=model-dt*HA, model_minus_CP=model-(dt-dg)*HA)
    for k, v in h.items():
        check(abs(float(v) - result['cycle']['hartree'][k]) < 1e-12, 'independent Decimal Eh cycle:' + k)
    for k, v in kcal.items():
        check(abs(float(v) - result['cycle']['kcal_mol'][k]) < 2e-9, 'independent Decimal kcal cycle:' + k)
    check(kcal['CP_shift'] > 0 and abs(kcal['model_minus_CP']) > abs(kcal['model_minus_original']), 'CP increases absolute unchanged-model error')
    sched_native = sum(int(v['elapsed_seconds'])*int(v['allocated_cpus']) for v in scheduler['native']['data']['records'].values())
    sched_collect = sum(int(v['elapsed_seconds'])*int(v['allocated_cpus']) for v in scheduler['collector']['data']['records'].values())
    check(set(scheduler['native']['data']['records']) == {v['job_id'] for v in fragments.values()}, 'scheduler exact native eight jobs')
    check(abs(native_allocated - result['cost']['admitted_new_fragment_allocated_CPU_seconds']) < 1e-8, 'native allocation sum')
    review = {'status': 'passed_independent_actual_replay_and_cycle', 'created_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
              'declaration': pin(HERE/'DECLARATION.json'), 'review_source': pin(__file__), 'collection': pin(result_path), 'checks_passed': len(checks), 'checks': checks,
              'fragments': fragments, 'original_complexes': complexes, 'independent_Decimal_hartree': {k: str(v) for k, v in h.items()},
              'independent_Decimal_kcal_mol': {k: str(v) for k, v in kcal.items()}, 'printed_D4_rounding_closure_kcal_mol': str((dt-de-dd)*HA),
              'cost': {'new_native_programs': 8, 'native_elapsed_times_four_ranks_CPU_seconds': native_allocated, 'scheduler_native_allocated_CPU_seconds': sched_native,
                       'scheduler_collector_allocated_CPU_seconds': sched_collect, 'scheduler_total_new_allocated_CPU_seconds': sched_native+sched_collect,
                       'reused_original_CPU_seconds_not_new': result['cost']['reused_original_allocated_CPU_seconds'],
                       'collector_actual_parent_CPU_seconds': result['cost']['collector_CPU_seconds'], 'review_wall_seconds': time.monotonic()-tick, 'review_CPU_seconds': time.process_time()-cpu,
                       'review_native_calls': 0, 'local_threads': 1, 'actual_native_CPU_utilization_measured': False},
              'scheduler_evidence': {k: {'source_path': v['source_path'], 'snapshot': v['snapshot']} for k, v in scheduler.items()},
              'conclusion': 'The positive metal-whole-shell CP shift makes the already negative model-minus-reference error more negative and larger in magnitude. This particular BSSE correction does not explain the failed model contrast.',
              'limitations': ['Effective frozen-f singlets and provisional PBE0-D4 reference remain unchanged.', 'Two-fragment CP does not remove intra-shell water-water BSSE or establish the basis limit, functional accuracy, physical electronic state, force consistency, equilibrium hydration or exchange free energy.', 'Some printed SCF criteria fail while the native gradient criterion signals convergence; retain all raw stopping and final-integration metadata under the original admission protocol.', 'W8 to W9 rearranges the original eight waters as well as adding one; this is not pure insertion.', 'No target update, refit or new scientific release follows from this review.']}
    save(OUT/'REVIEW.json', review)
    print(json.dumps({'review': pin(OUT/'REVIEW.json'), 'checks_passed': len(checks), 'kcal_mol': {k: float(v) for k,v in kcal.items()}, 'cost': review['cost']}))


if __name__ == '__main__':
    main()
