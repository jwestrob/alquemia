#!/usr/bin/env python3
"""Saved-only replay of job17670; independently reconstruct the four checks."""
import datetime as dt
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[4]
WORK = ROOT/'workspaces/lanm_selectivity_free_energy_20261021/aqueous_qualification_design_v1/coordination_screen_v1'
CELLS = ('aquo9_central_La', 'aquo9_central_Dy')


def pin(path):
    p = Path(path).resolve()
    return {'path': str(p), 'sha256': hashlib.sha256(p.read_bytes()).hexdigest()}


def independent_body_error(xyz, errors):
    # One ion and nine rigid waters. Ion and water categories carry equal
    # weight; each water contributes three net-force and three torque entries.
    bodies = []
    for ids in [[0]]+[[i, i+1, i+2] for i in range(1, 28, 3)]:
        center = [sum(xyz[i][c] for i in ids)/len(ids) for c in range(3)]
        net = [sum(errors[i][c] for i in ids) for c in range(3)]
        torque = [0., 0., 0.]
        for i in ids:
            dx, dy, dz = [xyz[i][c]-center[c] for c in range(3)]
            fx, fy, fz = errors[i]
            torque[0] += dy*fz-dz*fy
            torque[1] += dz*fx-dx*fz
            torque[2] += dx*fy-dy*fx
        bodies.append({'indices': ids, 'net_force_error': net, 'torque_error_over_1A': torque})
    ion = bodies[0]['net_force_error']
    water = [v for b in bodies[1:] for k in ('net_force_error', 'torque_error_over_1A') for v in b[k]]
    rms = math.sqrt(sum(v*v for v in ion)/6+sum(v*v for v in water)/108)
    maximum = max(abs(v) for v in ion+water)
    return {'RMS_kcal_mol_A': rms, 'maximum_component_kcal_mol_A': maximum,
            'passed': rms <= 2 and maximum <= 5, 'body_errors': bodies,
            'component_count': len(ion+water)}


def main():
    start, cpu = time.monotonic(), time.process_time()
    out = WORK/'actual_review_v1'; out.mkdir(parents=True, exist_ok=False)
    checks = []; pins = {}
    def check(ok, name):
        if not ok:
            raise ValueError(name)
        checks.append(name)
    def verify(p):
        check(pin(p['path']) == p, 'Pin '+p['path']); pins[p['path']] = p
        return Path(p['path'])
    def read(p):
        return json.loads(verify(p).read_text())
    def module(p, name):
        spec = importlib.util.spec_from_file_location(name, verify(p))
        m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
        return m
    rp = WORK/'run_v1/execution/RESULT.json'; result_pin = pin(rp); result = read(result_pin)
    bound = read(result['bound']); cfg = read(bound['config']); release = read(result['release'])
    verify(release['root_review']); verify(cfg['specification']); verify(bound['candidate_result'])
    check(result['implementation'] == bound['implementation'] == cfg['implementation'] == release['implementation'], 'Exact released implementation')
    verify(result['implementation'])
    check(release['bound'] == result['bound'] and release['specification'] == cfg['specification'], 'Exact release binding')
    check(release['status'] == 'approved_frozen_W9_development_screen' and release['maximum_native_programs'] == 4 and release['new_FD_energies'] == 0, 'Finite approved scope')
    check(result['status'] == 'complete' and result['job_id'] == '17670' and result['error'] is None, 'Actual terminal native completion')
    check(result['native_program_invocations'] == 4 and result['new_FD_energies'] == 0 and set(result['rows']) == set(CELLS), 'Complete finite inventory')
    check(not result['physical_accuracy_qualified'], 'Physical qualification remains false')
    native = module(bound['native_helper'], 'review_native')
    for p in bound['executables'].values():
        verify(p)
    verify(bound['numerical_helper']); verify(bound['body_helper'])
    ccfg = read(bound['candidate_config'])
    for name in ('water_result', 'transfer_result'):
        verify(ccfg[name])
    models = {}; anchors = {}; new_evidence = {}; native_wall = 0.; receipts_count = 0
    for gid in CELLS:
        row, actual = bound['rows'][gid], result['rows'][gid]; folder = rp.parent/gid
        check(actual['status'] == 'complete' and actual['input'] == row['files'], 'Complete actual row '+gid)
        for key, filename in (('xyz', 'system.xyz'), ('key', 'system.key'), ('parameters', 'amoeba09.prm')):
            verify(row['files'][key]); verify(row['source_files'][key])
            check(pin(folder/filename)['sha256'] == row['files'][key]['sha256'] == row['source_files'][key]['sha256'], 'Actual input bytes '+gid+' '+key)
        parsed = native.parse(folder/'analyze.out', folder/'testgrad.out', row['input'])
        check(parsed == actual['model']['native'], 'Exact raw new native replay '+gid)
        check(parsed['total_energy_kcal_mol'] == actual['model']['energy'] and parsed['force_kcal_mol_A'] == actual['model']['forces'], 'Full28 saved E/F replay '+gid)
        check(len(actual['receipts']) == 2, 'Two native receipts '+gid)
        for rec, program, args in zip(actual['receipts'], ('analyze', 'testgrad'), (['system.xyz', 'E'], ['system.xyz', 'Y', 'N'])):
            verify(rec['output']); receipts_count += 1; native_wall += rec['elapsed_seconds']
            check(rec['argv'] == [bound['executables'][program]['path']]+args and rec['cwd'] == str(folder), 'Exact native command '+gid+' '+program)
            check(rec['returncode'] == 0 and set(rec['threads'].values()) == {'1'}, 'Native success/one thread '+gid+' '+program)
        models[gid] = actual['model']
        new_evidence[gid] = {'energy_kcal_mol': parsed['total_energy_kcal_mol'], 'outputs': {k: parsed[k] for k in ('analyze_output', 'gradient_output')}, 'input': row['files']}
    check(receipts_count == 4, 'All four new programs replayed')
    for metal in ('La', 'Dy'):
        a = bound['anchor_model'][metal]; pv = a['provenance']; authority = read(pv['native_authority'])
        for p in pv['effective_input'].values():
            verify(p)
        for k in ('analyze_output', 'gradient_output'):
            verify(a['native'][k])
        parsed = native.parse(a['native']['analyze_output']['path'], a['native']['gradient_output']['path'], {'atoms': 25, 'derivative_test': False})
        check(parsed == a['native'] and parsed['total_energy_kcal_mol'] == a['energy'] and parsed['force_kcal_mol_A'] == a['forces'], 'Exact saved W8 raw E/F replay '+metal)
        check(authority['status'] == 'complete' and all(authority[k] == a[k] for k in ('energy', 'forces', 'native')), 'Actual W8 producer authority '+metal)
        check(pv['requested_vector'] == bound['candidate_parameters'], 'Frozen candidate W8 lineage '+metal)
        for k in ('key', 'parameters'):
            check(pv['effective_input'][k]['sha256'] == bound['rows']['aquo9_central_'+metal]['files'][k]['sha256'], 'Identical W8/W9 frozen Hamiltonian '+metal+' '+k)
        anchors[metal] = {'energy_kcal_mol': a['energy'], 'native_authority': pv['native_authority'], 'raw_outputs': {k: parsed[k] for k in ('analyze_output', 'gradient_output')}}
    refs = {}; ref_evidence = {}; panel = refmod = None
    for gid, p in cfg['references'].items():
        raw = read(p)
        for k in ('adapter', 'panel', 'input', 'geometry', 'output', 'engrad'):
            verify(raw[k])
        if refmod is None:
            refmod = module(raw['adapter'], 'review_reference'); panel = refmod.load_panel(verify(raw['panel']))
        cell = next(c for c in panel['cells'] if c['cell_id'] == 'ecp_pbe0__'+gid)
        replay = refmod.collect_one(panel, cell, Path(p['path']).parent, 0)
        check(all(raw[k] == v for k, v in replay.items()) and raw['status'] == 'complete' and all(raw['checks'].values()), 'Exact full reference raw replay '+gid)
        n = 9 if gid in CELLS else 8; metal = cell['metal']
        check(cell['recipe'] == 'ecp_pbe0' and cell['charge'] == 3 and cell['multiplicity'] == 1 and cell['physical_multiplicity'] == (1 if metal == 'La' else 6) and cell['ecp_core'] == (46 if metal == 'La' else 55) and cell['explicit_electrons'] == 8+10*n, 'Preserved effective/physical states '+gid)
        check(len(raw['forces_kcal_mol_A']) == 1+3*n and all(math.isfinite(v) for f in raw['forces_kcal_mol_A'] for v in f), 'Full finite reference force inventory '+gid)
        target = bound['rows'][gid]['reference'] if n == 9 else bound['anchor_reference'][metal]
        check(raw['energy_hartree']*refmod.HA == target['energy'] and raw['forces_kcal_mol_A'] == target['forces'] and p == target['raw_result'], 'Bound reference E/F lineage '+gid)
        refs[gid] = target
        ref_evidence[gid] = {k: raw[k] for k in ('energy_hartree', 'job_id', 'checks', 'input', 'geometry', 'output', 'engrad', 'panel', 'adapter', 'physical_state_qualified', 'method_accuracy_qualified')}
    check(len(refs) == 4, 'All four saved ECP references replayed')
    for metal in ('La', 'Dy'):
        a, b = [ref_evidence['aquo'+str(n)+'_central_'+metal] for n in (8, 9)]
        check(a['input']['sha256'] == b['input']['sha256'], 'Matched W8/W9 reference recipe '+metal)
    # A different summation order checks the balanced expression independently.
    mcontrast = math.fsum([models[CELLS[1]]['energy'], -models[CELLS[0]]['energy'], -anchors['Dy']['energy_kcal_mol'], anchors['La']['energy_kcal_mol']])
    rcontrast = math.fsum([refs[CELLS[1]]['energy'], -refs[CELLS[0]]['energy'], -refs['aquo8_central_Dy']['energy'], refs['aquo8_central_La']['energy']])
    error = mcontrast-rcontrast; old = result['comparison']['energy']
    check(max(abs(mcontrast-old['model']), abs(rcontrast-old['reference']), abs(error-old['error_kcal_mol'])) < 1e-8, 'Independent balanced contrast arithmetic')
    check(old['passed'] == (abs(error) <= 1), 'Unchanged1kcal energy criterion')
    xyz = bound['rows'][CELLS[0]]['geometry']['positions_A']
    check(xyz == bound['rows'][CELLS[1]]['geometry']['positions_A'] and len(xyz) == 28, 'Paired full28 coordinate identity')
    body_results = {}
    for label in ('La', 'Dy', 'paired'):
        if label == 'paired':
            errors = [[(a-b)-(c-d) for a,b,c,d in zip(am,al,rm,rl)] for am,al,rm,rl in zip(models[CELLS[1]]['forces'],models[CELLS[0]]['forces'],refs[CELLS[1]]['forces'],refs[CELLS[0]]['forces'])]
        else:
            gid = 'aquo9_central_'+label
            errors = [[a-b for a,b in zip(aa,bb)] for aa,bb in zip(models[gid]['forces'],refs[gid]['forces'])]
        calc = independent_body_error(xyz, errors); recorded = result['comparison']['rigid'][label]
        check(calc['component_count'] == 57 and max(abs(calc[k]-recorded[k]) for k in ('RMS_kcal_mol_A', 'maximum_component_kcal_mol_A')) < 1e-12 and calc['passed'] == recorded['passed'], 'Independent full rigid-body arithmetic '+label)
        body_results[label] = calc
    check(result['comparison']['declared_comparisons'] == result['comparison']['available_comparisons'] == 4 and result['comparison']['passed'] is False, 'Four complete checks; failure retained')
    blob = (WORK/'run_v1/completion_event.json').read_bytes(); (out/'completion_event.snapshot.json').write_bytes(blob)
    scheduler = json.loads(blob)['records']['17670']
    check(scheduler['state'] == 'COMPLETED' and scheduler['exit_code'] == '0:0' and int(scheduler['allocated_cpus']) == 1, 'Actual scheduler completion/resources')
    cost = {'scheduler_elapsed_seconds': int(scheduler['elapsed_seconds']), 'allocated_CPUs': int(scheduler['allocated_cpus']), 'allocated_CPU_seconds': int(scheduler['elapsed_seconds'])*int(scheduler['allocated_cpus']), 'native_programs': 4, 'summed_native_wall_seconds': native_wall, 'aggregate_native_child_CPU_seconds': result['aggregate_native_child_user_CPU_seconds']+result['aggregate_native_child_system_CPU_seconds'], 'adapter_parent_CPU_seconds': result['parent_CPU_seconds'], 'adapter_wall_seconds': result['elapsed_seconds']}
    review = {'status': 'actual_replay_passed_scientific_screen_failed', 'job_id': '17670', 'utc': dt.datetime.now(dt.timezone.utc).isoformat(), 'implementation': pin(__file__), 'actual_result': result_pin, 'bound': result['bound'], 'release': result['release'], 'source': result['implementation'], 'completion_snapshot': pin(out/'completion_event.snapshot.json'), 'checks': checks, 'check_count': len(checks), 'verified_pins': list(pins.values()), 'new_model_cells': new_evidence, 'saved_W8_anchors': anchors, 'saved_reference_cells': ref_evidence, 'balanced_contrast': {'model_kcal_mol': mcontrast, 'reference_kcal_mol': rcontrast, 'error_kcal_mol': error, 'threshold_kcal_mol': 1, 'passed': abs(error) <= 1}, 'rigid_response': body_results, 'complete_comparisons': 4, 'failed_comparisons': 3, 'native_cost': cost, 'new_native_calls': 0, 'repeated_controls': 0, 'local_threads': 1, 'review_wall_seconds': time.monotonic()-start, 'review_CPU_seconds': time.process_time()-cpu, 'reserved_donor_outcomes_read': False, 'physical_accuracy_qualified': False, 'automatic_refit_authorized': False}
    (out/'REVIEW.json').write_text(json.dumps(review, indent=2, allow_nan=False)+'\n')
    print(json.dumps({'review': pin(out/'REVIEW.json'), 'checks': len(checks), 'contrast': review['balanced_contrast'], 'rigid': {k: {f:v[f] for f in ('RMS_kcal_mol_A','maximum_component_kcal_mol_A','passed')} for k,v in body_results.items()}, 'cost': cost, 'review_wall_seconds':review['review_wall_seconds'], 'review_CPU_seconds': review['review_CPU_seconds']}))


if __name__ == '__main__':
    main()
