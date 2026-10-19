#!/usr/bin/env python3
"""Bounded whole2049 MACE+fixed-FF ALPB derivative check; every stage is Slurm-only."""
import argparse
import copy
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import resource
import socket
import sys
import time
import traceback

PROTOCOL = 'mex_whole2049_before_four_state_directional_check_v1'
PANEL_SHA = 'dd611ca03982e4f2548a88dabf9477e33aba7592fc2b41d832e60d13edd2add7'
COLLECTION_SHA = '4c0c8b7379a2210f64267f0df94b1df72b54902d7233ec4b9614d4e941349d27'
QUALIFICATION_SHA = '15fb2a83b3f57e84ed7fcda0f61094dbad241ec265edd42fc6cce88dc4e55987'
DIRECTIONS_SHA = '8972af45b94e6e1eefdc916030b3fa08fc1fc79a1ae581340a95da2f01c9ee47'
SOLVENT_SHA = 'e236a22dbb57a1105ef74d0542d96827dfbcc661390eba7b785b1150f1e62df8'
EV, HA, BOHR = 23.06054783061903, 627.509474, .529177210903
DIRECTIONS = ('core_backbone', 'boundary_anchor', 'external_water', 'spectator')
MOVING = (['A/84//C'], ['A/83//CA'], ['A/301//O', 'A/301//H1', 'A/301//H2'], ['A/201//ND'])
STATES = ('before_La_aligned', 'before_Dy_aligned', 'before_La_one_Nd_reversed', 'before_Dy_one_Nd_reversed')
STEPS = (.0002, .0001)
LIMITS = {'maximum_MACE_calls': 64, 'maximum_ALPB_calls': 16, 'forward_seconds': 600.,
          'peak_device_fraction': .90, 'absolute_kcal_mol_A': .005, 'relative': .0001}
BIND = ('protocol', 'implementation', 'specification', 'panel_implementation', 'panel_result',
        'old_qualification', 'directions', 'solvent_implementation', 'solvent_config', 'limits')


def require(ok, message):
    if not ok:
        raise ValueError(message)


def pin(path):
    p = Path(path).resolve()
    with p.open('rb') as f:
        sha = hashlib.file_digest(f, 'sha256').hexdigest()
    return {'path': str(p), 'sha256': sha}


def verify(p):
    require(pin(p['path'])['sha256'] == p['sha256'], 'Changed pin: ' + p['path'])
    return Path(p['path'])


def read(p):
    return json.loads(verify(p).read_text())


def write(path, value):
    with Path(path).open('x') as f:
        json.dump(value, f, indent=2, allow_nan=False); f.write('\n')
    return pin(path)


def module(name, p):
    spec = importlib.util.spec_from_file_location(name, verify(p))
    obj = importlib.util.module_from_spec(spec); spec.loader.exec_module(obj)
    return obj


def binding(c):
    return {k: c[k]['sha256'] if k == 'implementation' else c[k] for k in BIND}


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False).encode()).hexdigest()


def fixture_id(direction, h, sign):
    return '%s_h%.4f_%s' % (direction, h, 'minus' if sign == -1 else 'plus')


def prepare(cfg):
    require(os.environ.get('SLURM_JOB_ID', '').isdigit(), 'Slurm required before scientific imports')
    require(cfg['protocol'] == PROTOCOL and cfg['limits'] == LIMITS
            and pin(__file__)['sha256'] == cfg['implementation']['sha256'], 'Frozen protocol differs')
    for k in BIND[1:-1]: verify(cfg[k])
    for k, sha in (('panel_implementation', PANEL_SHA), ('panel_result', COLLECTION_SHA),
                   ('old_qualification', QUALIFICATION_SHA), ('directions', DIRECTIONS_SHA),
                   ('solvent_implementation', SOLVENT_SHA)):
        require(cfg[k]['sha256'] == sha, 'Frozen source differs: ' + k)
    panel = module('whole_force_panel', cfg['panel_implementation'])
    parent = read(cfg['panel_result']); pcfg = read(parent['config'])
    require(parent['status'] == 'complete' and parent['stage'] == 'collect'
            and parent['implementation']['sha256'] == PANEL_SHA
            and parent['total_new_MACE_forwards_across_panel'] == 8, 'Complete whole panel required')
    pm, cells, mc, np = panel.prepare(pcfg)
    panel.preflight_gate(pcfg, pm)
    scout = panel.scout_gate(pcfg, pm, cells, np)
    remaining = read(pcfg['remaining_result'])
    require(parent['binding'] == panel.binding(pcfg) and read(parent['manifest']) == pm
            and parent['parents'] == {'scout': pcfg['scout_result'], 'remaining': pcfg['remaining_result']}
            and remaining['status'] == 'complete' and remaining['rows'] == parent['rows']
            and remaining['native_forward_calls'] == remaining['new_MACE_calls_attempted'] == 7
            and remaining['rows'][0] == dict(scout, reused_from=pcfg['scout_result'])
            and [r['task_id'] for r in parent['rows']] == list(panel.IDS), 'Actual parent panel replay differs')
    centers = {}
    for task, row in zip(pm['tasks'], parent['rows']):
        if task['task_id'] in STATES:
            raw, assembled = panel.admit_saved(row, task, cells['before_' + task['metal']], np)
            centers[task['task_id']] = {'task': task, 'row': row, 'assembled': assembled}
    require(tuple(centers) == STATES, 'Four before centers required')
    require(cells['before_La']['positions_A'] == cells['before_Dy']['positions_A'], 'Central paired coordinates differ')
    ids = pm['source_ids']; x = np.asarray(cells['before_La']['positions_A'])
    static = read(pcfg['solvent_parent']); scfg = read(cfg['solvent_config'])
    require(read(static['build']['parent'])['config'] == cfg['solvent_config'], 'Solvent build parent differs')
    binary = static['build']['binary']; verify(binary); verify(static['build']['fortran'])
    require(static['build']['command'] == scfg['compiler_command'] + [static['build']['fortran']['path'], '-o', binary['path']]
            + scfg['link_arguments'], 'Original solvent build command differs')
    for p in scfg['build_inputs']: verify(p)
    q = np.loadtxt(verify(cells['before_La']['solvent']['charges'])).tolist()
    require(q == np.loadtxt(verify(cells['before_Dy']['solvent']['charges'])).tolist(), 'Paired fixed charges differ')
    old = read(cfg['old_qualification']); oc = read(old['config']); directions = read(cfg['directions'])
    old_preflight = read(old['preflight'])
    require(old['status'] == 'passed' and old['all_declared_passed'] and str(old['slurm_job_id']) == '12952'
            and old_preflight['directions'] == cfg['directions'] and directions['selected_before_fresh_calls']
            and directions['new_molecular_calls'] == 0 and directions['source_ids'] == ids
            and directions['starts'] == oc['starts'] and directions['starts']['La'] == cells['before_La']['point']
            and directions['rule'] == list(DIRECTIONS)
            and oc['settings']['FD_steps_A'] == list(STEPS)
            and oc['settings']['FD_absolute_kcal_mol_A'] == LIMITS['absolute_kcal_mol_A']
            and oc['settings']['FD_relative'] == LIMITS['relative'], 'Original fixture protocol differs')
    candidates = []
    for rp in old['new_points']:
        request = read(rp)
        if request['state']['metal'] == 'La':
            require(request['source_ids'] == ids and request['coordinate_dimensions'] == 6147, 'Old physical order differs')
            if request['all_physical_coordinates_A'] != x.tolist():
                candidates.append((rp, request))
    require(len(candidates) == 16, 'Exactly sixteen archived La displacements required')
    fixtures = []; used = set()
    for direction, moving in zip(DIRECTIONS, MOVING):
        vector = directions['vectors']['La'][direction]; v = np.asarray(vector)
        require(v.shape == (2049, 3) and np.isfinite(v).all()
                and abs(float(np.linalg.norm(v)) - 1.) <= 1e-12
                and [ids[i] for i in np.flatnonzero(np.any(v != 0., axis=1))] == moving, 'Frozen direction identity/norm differs')
        oldcheck = next(c for c in old['checks'] if c['metal'] == 'La' and c['check'] == direction)
        require(oldcheck['direction_A_per_A'] == vector and oldcheck['moving_source_ids'] == moving
                and oldcheck['passed'] and [s['h_A'] for s in oldcheck['steps']] == list(STEPS), 'Archived direction declaration differs')
        for h in STEPS:
            for sign in (-1, 1):
                expected = x-h*v if sign == -1 else x+h*v
                matches = [(rp, r) for rp, r in candidates if np.array_equal(r['all_physical_coordinates_A'], expected)]
                require(len(matches) == 1, 'Archived signed displacement is missing or ambiguous')
                rp, request = matches[0]; require(rp['sha256'] not in used, 'Duplicate displacement reuse'); used.add(rp['sha256'])
                point = read(request['point']); atoms = read(point['physical']); mapping = read(point['mapping'])
                require(point['status'] == 'complete' and point['metal'] == mapping['metal'] == 'La'
                        and mapping['physical_atoms'] == atoms and [a['id'] for a in atoms] == ids
                        and [a['element'] for a in atoms] == cells['before_La']['symbols']
                        and [a['xyz_A'] for a in atoms] == request['all_physical_coordinates_A']
                        and request['native_energies'] == point['native_energies'], 'Archived displacement POINT differs')
                solvent = panel.admit_solvent(point, atoms, mapping, q, binary, np)
                fixtures.append({'fixture_id': fixture_id(direction, h, sign), 'direction': direction, 'h_A': h,
                    'sign': sign, 'old_request': rp, 'point': request['point'], 'physical': point['physical'],
                    'mapping': point['mapping'], 'coordinates_digest': digest(request['all_physical_coordinates_A']),
                    'solvent': solvent})
    require(len(used) == 16, 'Displacement denominator incomplete')
    manifest = {'protocol': PROTOCOL, 'binding': binding(cfg), 'panel_manifest': parent['manifest'],
                'centers': {k: c['row'] for k, c in centers.items()}, 'directions': cfg['directions'],
                'source_ids': ids, 'fixtures': fixtures, 'solvent_binary': binary,
                'solvent_build': static['build'], 'fixed_charges': cells['before_La']['solvent']['charges'],
                'reused_central_MACE': 4, 'reused_central_ALPB': 2, 'reused_displaced_La_ALPB': 16,
                'declared_MACE_displacements': 64, 'declared_new_Dy_ALPB': 16, 'limits': LIMITS}
    return panel, pcfg, mc, np, cells, centers, fixtures, directions['vectors']['La'], q, scfg, manifest


def new_solvent(panel, helper, fixture, cells, q, binary, env, out, np):
    atoms = copy.deepcopy(read(fixture['physical'])); atoms[2045]['element'] = 'Dy'
    require([a['element'] for a in atoms] == cells['before_Dy']['symbols'], 'Displaced Dy inventory differs')
    physical = write(out/'PHYSICAL.json', atoms)
    mapping = {'physical_atoms': atoms, 'source_ids': cells['before_Dy']['source_ids'],
               'metal': 'Dy', 'parent_mapping': cells['before_Dy']['mapping'],
               'archived_La_fixture': fixture, 'scope': 'whole physical only; no core, field or cap update'}
    mp = write(out/'PHYSICAL_MAPPING.json', mapping)
    c = {'physical': atoms, 'x': np.asarray([a['xyz_A'] for a in atoms]), 'q': q, 'mapping_pin': mp}
    sol = helper.solvent_cell(Path(binary['path']), c, out/'native', env, np)
    # A solvent-only gradient view permits the unchanged receipt admission routine.
    # This is not a hybrid POINT or an independent full-energy/gradient evaluation.
    gradient = write(out/'ALPB_GRADIENT.json', {'source_ids': mapping['source_ids'], 'physical_mapping': mp,
        'component_gradients_kcal_mol_A': {'fixed_FF_charge_ALPB':
            (np.asarray(sol['gradient_hartree_per_bohr'])*HA/BOHR).tolist()}})
    view = {'ALPB': pin(out/'native/RESULT.json'), 'mapping': mp, 'assembled_gradient': gradient,
            'native_energies': {'ALPB_Eh': sol['total_hartree']}}
    admitted = panel.admit_solvent(view, atoms, mapping, q, binary, np)
    return {'status': 'complete', 'fixture_id': fixture['fixture_id'], 'physical': physical,
            'mapping': mp, 'gradient_view': gradient, 'admission': admitted}


def run(cfg, world, out, result):
    panel, pcfg, mc, np, cells, centers, fixtures, vectors, q, scfg, manifest = world
    pre = read(cfg['preflight_result'])
    require(pre['status'] == 'preflight_passed' and pre['stage'] == 'preflight'
            and pre['binding'] == binding(cfg) and read(pre['manifest']) == manifest
            and pre['new_MACE_calls_attempted'] == pre['new_ALPB_calls_attempted'] == pre['native_model_loads_attempted'] == 0,
            'Actual identical zero-call preflight required')
    require(int(os.environ['SLURM_CPUS_PER_TASK']) == 12, 'Twelve allocated CPUs required')
    result['preflight_result'] = cfg['preflight_result']
    loader = module('whole_force_native_loader', pcfg['native_helper'])
    solvent_helper = module('whole_force_solvent_helper', cfg['solvent_implementation'])
    env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1',
               BLIS_NUM_THREADS='1', NUMEXPR_NUM_THREADS='1', LD_LIBRARY_PATH=scfg['library_directory'])
    result['native_model_loads_attempted'] = 1; tick = time.monotonic()
    native = panel.native_class(loader)(mc)
    result['model_load_seconds'] = time.monotonic()-tick
    result['device'] = native.device; result['versions'] = native.versions
    try:
        require('H200' in native.device, 'Declared H200 required')
        for direction in DIRECTIONS:
            v = np.asarray(vectors[direction])
            for h in STEPS:
                pair = [next(f for f in fixtures if f['direction'] == direction and f['h_A'] == h and f['sign'] == s) for s in (-1, 1)]
                for tid in STATES:
                    central = centers[tid]; task = central['task']; metal = task['metal']
                    values = []
                    for fixture in pair:
                        fid = fixture['fixture_id']; sol = fixture['solvent']
                        if metal == 'Dy':
                            slot = next(s for s in result['Dy_solvent_slots'] if s['fixture_id'] == fid)
                            if slot['status'] != 'complete':
                                d = out/('Dy_ALPB_' + fid); d.mkdir()
                                slot.update(status='attempted', output_directory=str(d)); slot.pop('reason', None)
                                result['new_ALPB_calls_attempted'] += 1
                                require(result['new_ALPB_calls_attempted'] <= 16, 'ALPB ceiling exceeded')
                                slot.update(new_solvent(panel, solvent_helper, fixture, cells, q, manifest['solvent_binary'], env, d, np))
                            sol = slot['admission']
                        atoms = read(fixture['physical']); positions = [a['xyz_A'] for a in atoms]
                        rid = tid + '__' + fid; slot = next(s for s in result['MACE_slots'] if s['request_id'] == rid)
                        d = out/rid; d.mkdir()
                        template = {k: task[k] for k in ('task_id', 'geometry', 'metal', 'hypothesis', 'charge',
                                                      'multiplicity', 'metal_index', 'atom_count', 'all_electron_count')}
                        identity = {'symbols': cells['before_'+metal]['symbols'], 'positions_A': positions,
                                    'source_ids': manifest['source_ids'], 'state_template': template,
                                    'model_binding': panel.binding(pcfg)}
                        request = {'request_id': rid, 'parent_center': central['row'], 'state_template': template,
                                   'fixture': fixture, 'actual_input': identity, 'input_digest': digest(identity),
                                   'ALPB': sol['result'], 'no_center_cache_key_inherited': True}
                        slot['request'] = write(d/'REQUEST.json', request)
                        slot.update(status='attempted', output_directory=str(d)); slot.pop('reason', None)
                        result['new_MACE_calls_attempted'] += 1
                        require(result['new_MACE_calls_attempted'] <= 64, 'MACE ceiling exceeded')
                        raw = native.evaluate(template, identity['symbols'], positions)
                        slot['raw_MACE'] = write(d/'MACE.json', raw)
                        require(raw['forward_seconds'] <= LIMITS['forward_seconds']
                                and raw['peak_GPU_allocated_bytes'] <= LIMITS['peak_device_fraction']*raw['device_total_memory_bytes'],
                                'Native forward resource gate failed')
                        slot.update(status='complete', ALPB=sol['result'], input_digest=request['input_digest'])
                        values.append((raw['energy_eV'], sol['energy_hartree']))
                    minus, plus = values
                    differences = {'whole_MACE_vacuum': (plus[0]-minus[0])*EV,
                                   'fixed_FF_charge_ALPB': (plus[1]-minus[1])*HA}
                    differences['total'] = math.fsum(differences.values())
                    check = next(c for c in result['checks'] if c['state'] == tid and c['direction'] == direction and c['h_A'] == h)
                    parts = {}; assembled = central['assembled']
                    for component, delta in differences.items():
                        g = assembled['gradient_kcal_mol_A'] if component == 'total' else assembled['component_gradients_kcal_mol_A'][component]
                        analytic = float(np.sum(np.asarray(g)*v)); numeric = delta/(2*h)
                        tolerance = LIMITS['absolute_kcal_mol_A'] + LIMITS['relative']*max(abs(analytic), abs(numeric))
                        residual = numeric-analytic
                        parts[component] = {'analytic': analytic, 'finite_difference': numeric, 'residual': residual,
                                            'tolerance': tolerance, 'passed': math.isfinite(residual) and abs(residual) <= tolerance}
                    check.update(status='passed' if all(p['passed'] for p in parts.values()) else 'failed', components=parts)
                    if check['status'] != 'passed':
                        result['status'] = 'failed_derivative_gate'
                        raise ValueError('Directional derivative gate failed: ' + tid + '/' + direction + '/' + str(h))
        require(result['new_MACE_calls_attempted'] == native.calls == 64 and result['new_ALPB_calls_attempted'] == 16,
                'Complete run call denominator differs')
        result['status'] = 'passed'; result['all_declared_passed'] = True
    finally:
        result['native_forward_calls'] = native.calls
        result['peak_GPU_allocated_bytes'] = native.torch.cuda.max_memory_allocated()
        result['device_total_memory_bytes'] = native.torch.cuda.get_device_properties(0).total_memory


def main():
    p = argparse.ArgumentParser(description=__doc__); p.add_argument('stage', choices=('preflight', 'run'))
    p.add_argument('--config', type=Path, required=True); p.add_argument('--output', type=Path, required=True)
    a = p.parse_args(); require(os.environ.get('SLURM_JOB_ID', '').isdigit(), 'Slurm required')
    out = a.output.resolve(); out.mkdir(parents=True, exist_ok=False); tick = time.monotonic()
    fids = [fixture_id(d, h, s) for d in DIRECTIONS for h in STEPS for s in (-1, 1)]
    result = {'protocol': PROTOCOL, 'stage': a.stage, 'status': 'execution_failure', 'all_declared_passed': False,
              'implementation': pin(__file__), 'config': pin(a.config), 'slurm_job_id': os.environ['SLURM_JOB_ID'],
              'host': socket.gethostname(), 'new_MACE_calls_attempted': 0, 'native_forward_calls': 0,
              'native_model_loads_attempted': 0, 'new_ALPB_calls_attempted': 0,
              'new_DFT_calls': 0, 'new_GFN_calls': 0, 'new_MM_calls': 0, 'new_center_calls': 0, 'optimizer_steps': 0,
              'MACE_slots': [{'request_id': t+'__'+f, 'status': 'unavailable', 'reason': 'not_attempted'} for t in STATES for f in fids],
              'Dy_solvent_slots': [{'fixture_id': f, 'status': 'unavailable', 'reason': 'not_attempted'} for f in fids],
              'checks': [{'state': t, 'direction': d, 'h_A': h, 'status': 'unavailable'} for t in STATES for d in DIRECTIONS for h in STEPS],
              'interpretation': 'Sampled BEFORE-gradient consistency only; neither aggregate hypothesis is selected or physically qualified'}
    try:
        cfg = json.loads(a.config.read_text()); result['binding'] = binding(cfg)
        world = prepare(cfg); result['manifest'] = write(out/'MANIFEST.json', world[-1])
        if a.stage == 'preflight': result['status'] = 'preflight_passed'
        else: run(cfg, world, out, result)
    except Exception as exc:
        if result['status'] != 'failed_derivative_gate': result['status'] = 'execution_failure'
        result['failure'] = f'{type(exc).__name__}: {exc}'; result['traceback'] = traceback.format_exc()
        for slot in result['MACE_slots'] + result['Dy_solvent_slots']:
            if slot['status'] == 'attempted': slot.update(status='unavailable', reason=result['failure'])
    finally:
        result['elapsed_seconds'] = time.monotonic()-tick
        result['process_CPU_seconds'] = resource.getrusage(resource.RUSAGE_SELF).ru_utime + resource.getrusage(resource.RUSAGE_SELF).ru_stime
        result['peak_parent_RSS_KiB'] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        result['allocated_CPUs'] = int(os.environ.get('SLURM_CPUS_PER_TASK', '1'))
        result['native_environment'] = {k: os.environ.get(k) for k in ('CUDA_VISIBLE_DEVICES', 'OMP_NUM_THREADS',
            'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS', 'SLURM_JOB_GPUS')}
        write(out/'RESULT.json', result)
    return 0 if result['status'] in ('preflight_passed', 'passed') else 1


if __name__ == '__main__':
    sys.exit(main())
