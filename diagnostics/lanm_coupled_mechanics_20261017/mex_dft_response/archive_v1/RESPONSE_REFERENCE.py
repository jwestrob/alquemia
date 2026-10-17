#!/usr/bin/env python3
"""Fixed eight-cell Mex embedded reference; saved preparation and one native call."""
import argparse
import copy
from decimal import Decimal
import hashlib
import importlib
import importlib.util
import json
import math
import os
from pathlib import Path
import resource
import sys
import time
import traceback

PROTOCOL = 'mex_saved_terminal_embedded_PBE0_response_eight_v1'
INVENTORY_SHA = 'e4f226f90d0da0fdb9a9aab246536742e26a9aeb72820e7f3b51fa4582fa4a9f'
QUALIFIER_SHA = '74b56bf12de0543251c59dcf7a6bb1bb3b84afdedc078ac1ff3045471dec7f48'
PROPOSAL_SHA = '98c04e3b7118052791cb562531ee8538b616805d7530c72302c548256ca87369'
ALIASES = ('La_terminal', 'Dy_terminal', 'La_cartesian_terminal', 'Dy_cartesian_terminal')
METALS = ('La', 'Dy')
IDS = tuple(a+'__'+m for a in ALIASES for m in METALS)
SCOUT = IDS[0]
COMPONENTS = ('PBE0_embedded', 'MM_ledger', 'fixed_FF_charge_ALPB')
GROUPS = ('retained_bonded', 'MM_Coulomb', 'MM_LJ', 'QM_MM_LJ')
SETTINGS = {'maximum_new_DFT_cells': 8, 'mpi_ranks': 12, 'concurrent_endpoints': 1,
            'first_cost_scout': SCOUT, 'new_other_component_calls': 0,
            'deformation_work_resolution_kcal_mol': 1.,
            'projection_agreement_target_kcal_mol_A': 1.,
            'projection_rule': 'both full physical own-metal restricted-to-broad chords, unit Cartesian L2 norm',
            'native_gradient_coordinate_atol_A': 1e-6, 'native_gradient_energy_atol_Eh': 1e-8}
CAVEATS = {'physical_state_qualified': False, 'global_force_qualified': False,
           'SCF_error_bound_available': False, 'stationary_minimum_claimed': False,
           'affinity_claimed': False, 'spectators': ['A/201//ND', 'A/202//ND', 'A/204//ND'],
           'represented_waters': 167, 'whole_charge_e': 2, 'QM_charge_e': 0,
           'electronic_model': 'spin-free frozen-f; physical La singlet/Dy sextet, effective singlets; no SOC',
           'scaffold_reference_limit': 'The unchanged MM ledger cannot independently validate MM stiffness.'}


def require(ok, message):
    if not ok:
        raise ValueError(message)


def pin(path):
    p = Path(path).resolve()
    return {'path': str(p), 'sha256': hashlib.sha256(p.read_bytes()).hexdigest()}


def verify(record):
    require(record is not None and pin(record['path'])['sha256'] == record['sha256'],
            'Missing/changed pinned artifact: '+str(record))
    return Path(record['path'])


def read(record):
    return json.loads(verify(record).read_text())


def write(path, value):
    with Path(path).open('x') as f:
        json.dump(value, f, indent=2, allow_nan=False)
        f.write('\n')
    return pin(path)


def canonical(value):
    return json.loads(json.dumps(value, allow_nan=False))


def source_hashes(value):
    """Only for code/specification trees: admit verified byte-identical snapshots."""
    if isinstance(value, dict) and set(value) == {'path','sha256'}:
        verify(value)
        return value['sha256']
    require(isinstance(value, dict), 'Unsupported source-pin tree')
    return {k:source_hashes(v) for k,v in value.items()}


def module(name, record):
    spec = importlib.util.spec_from_file_location(name, verify(record))
    obj = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(obj)
    return obj


def binding(cfg):
    return {k: cfg[k] for k in ('protocol', 'implementation', 'specification', 'plan',
            'inventory', 'proposal', 'qualifier', 'recipe_config', 'gradient_helpers', 'settings')}


def inputs(cfg):
    require(os.environ.get('SLURM_JOB_ID', '').isdigit(), 'Slurm required, including saved analysis')
    require(cfg['protocol'] == PROTOCOL and cfg['settings'] == SETTINGS, 'Frozen protocol differs')
    require(pin(__file__)['sha256'] == cfg['implementation']['sha256'], 'Adapter source differs')
    for p in (cfg['implementation'], cfg['specification'], cfg['plan'], cfg['proposal']):
        verify(p)
    require(cfg['inventory']['sha256'] == INVENTORY_SHA and cfg['qualifier']['sha256'] == QUALIFIER_SHA
            and cfg['proposal']['sha256'] == PROPOSAL_SHA,
            'Frozen inventory/derivative implementation differs')
    inv = read(cfg['inventory']); qc = read(cfg['recipe_config'])
    require(inv['plan'] == cfg['plan'] and inv['DFT_recipe_config'] == cfg['recipe_config']
            and qc['implementation']['sha256'] == QUALIFIER_SHA, 'Original recipe lineage differs')
    require(inv['declared_cells'] == 8 and inv['declared_geometries'] == 4
            and [c['task_id'] for c in inv['cells']] == list(IDS)
            and [g['alias'] for g in inv['geometries']] == list(ALIASES)
            and inv['existing_DFT_reusable_cells'] == 0, 'Fixed source denominator differs')
    require(set(cfg['endpoint_results']) == set(IDS), 'All eight result slots must remain explicit')
    qual = module('mex_response_derivative_qualification', cfg['qualifier'])
    original, ref, world, origins, native = qual.inputs(qc)
    for name, record in cfg['gradient_helpers'].items():
        verify(record)
        obj = importlib.import_module(name)
        require(pin(obj.__file__)['sha256'] == record['sha256'], 'Imported gradient helper differs: '+name)
        native[name] = obj
    require(cfg['gradient_helpers']['metal_environment_force_assembly']['sha256'] == pin(world[12].__file__)['sha256'],
            'Original cap-scatter implementation differs')
    manifest = read(inv['reference_manifest'])
    require(manifest == world[7] and manifest['method'] == native['metal_environment_frozen_f'].METHOD,
            'Original analytic recipe differs')
    for name, record in cfg['gradient_helpers'].items():
        require(record['sha256'] == manifest['implementation'][name+'.py']['sha256'],
                'Original native gradient implementation differs: '+name)
    return inv, qc, qual, world, origins, native


def source_replay(inv, out):
    """Use the reviewed pool's saved admission paths, never its scoring entry point."""
    old = read(inv['pool_result']); pcfg = read(old['config'])
    require(old['status'] == 'complete' and old['slurm_job_id'] == '12981'
            and old['declared_cells'] == 14 and old['declared_geometries'] == 7, 'Actual common7 unavailable')
    pool = module('mex_response_saved_pool', old['implementation'])
    require(old['protocol_id'] == pool.PROTOCOL and pcfg['implementation']['sha256'] == old['implementation']['sha256'],
            'Pool implementation/protocol differs')
    ecfg = read(pcfg['evaluator_config'])
    em = module('mex_response_saved_evaluator', ecfg['implementation'])
    sm = module('mex_response_saved_search', pcfg['search_implementation'])
    opt = module('mex_response_saved_optimizer', pcfg['optimizer'])
    ev = em.MexEvaluator(ecfg, out/'saved_evaluator')
    guard, unused = sm.build_guard(ev, opt, out)
    oldpool = module('mex_response_old_common5', ev.cfg['pool_implementation'])
    matrix, by_coordinate, baseline, baseline_pin = pool.baseline(ev, oldpool, guard, out)
    slots, failures = pool.admit_searches(pcfg, ev, sm, opt, guard, matrix, by_coordinate)
    require(all(s['status'] == 'admitted' for s in slots), 'Both actual broad searches must be admitted')
    for actual in old['matrix']:
        for metal, cell in actual['cells'].items():
            if cell.get('cartesian_request'):
                pool.import_request(ev, metal, cell['cartesian_request'])
    for row in matrix:
        actual = next(a for a in old['matrix'] if a['coordinate_sha256'] == row['coordinate_sha256'])
        require(row['aliases'] == actual['aliases'], 'Common7 candidate identity differs')
        for metal in METALS:
            cell = pool.cached_cell(ev, metal, row['_x'])
            require(cell is not None and cell['point'] == actual['cells'][metal]['point']
                    and cell['native_energies'] == actual['cells'][metal]['native_energies'], 'Common7 raw cell replay differs')
            row['cells'][metal] = cell
        pool.update_work(ev, row)
    replay = pool.summarize(oldpool, ev, matrix, pool.ALIASES, slots)
    require(canonical(replay) == old['common7'], 'Common7 native menu replay differs')
    require(ev.new_evaluations == 0 and all(v == 0 for v in ev.counters.values()) and ev.old.native is None,
            'Saved replay must not initialize a model or evaluate components')
    return ev, {'pool': inv['pool_result'], 'baseline_replay': baseline_pin,
                'terminal_slots': slots, 'failed_exact_requests': failures,
                'identity_sha256': ev.identity_hash, 'counters': ev.counters}


def shared_cell(c, ev, origins, qc, native):
    """Replay current geometry, ledger and raw solvent; archived maps stay authoritative."""
    np = ev.np; point = read(c['approximate_point']); trial = read(c['trial']); mp = read(c['mapping'])
    physical = read(c['physical']); x = np.asarray([a['xyz_A'] for a in physical]); ids = [a['id'] for a in physical]
    require(point['status'] == 'complete' and point['metal'] == c['metal'] and point['trial'] == c['trial']
            and point['physical'] == c['physical'] and trial['mapping'] == c['mapping']
            and trial['source_inputs'] == c['source_inputs'] and mp['inputs'] == c['source_inputs']
            and mp['physical_atoms'] == physical and ids == ev.ids, 'Current POINT/source/mapping differs')
    require(x.shape == (2049, 3) and np.isfinite(x).all() and len(set(ids)) == 2049
            and hashlib.sha256(np.asarray(x, dtype='<f8').tobytes()).hexdigest() == c['coordinate_sha256'],
            'Full physical coordinate identity differs')
    require(mp['metal'] == c['metal'] and len(mp['core']) == 209 and len(mp['field']) == 1840
            and len(mp['caps']) == 4 and mp['xyz'] == trial['xyz'] == c['core_xyz']
            and mp['pointcharges'] == read(trial['field'])['source'] == c['pointcharges'], 'Current core/field inventory differs')
    mapper = native['metal_environment_force_assembly']
    regenerated = mapper.prepare_mapping(verify(c['source_inputs']), 'A', c['metal'])
    # Existing mapper admission already checks cap anchors/coordinates (1e-8 A)
    # and Jacobians (1e-12). Retain archived Jacobians for the native scatter.
    for k in ('inputs', 'configuration', 'metal', 'physical_atoms', 'core', 'field', 'xyz', 'pointcharges'):
        require(regenerated[k] == mp[k], 'Archived current mapping differs: '+k)
    for old, new in zip(mp['caps'], regenerated['caps']):
        for k in old:
            if k.startswith('jacobian_'):
                require(np.allclose(old[k], new[k], rtol=0, atol=1e-12), 'Existing cap-Jacobian gate failed')
            else:
                require(old[k] == new[k], 'Cap identity/anchor differs')
    source = read(c['source_inputs'])
    require(source['configurations']['A']['core_mapping'] == c['core_mapping']
            and source['boundary_mapping'] == c['boundary_mapping'], 'Current boundary provenance differs')
    mm = read(c['MM_result']); mc = read(c['MM_config']); ledger = read(c['MM_ledger'])
    dispatch = (ev.cfg['base_implementation'], ev.cfg['pool_implementation'], ev.pcfg['remaining_implementation'])
    for p in (*dispatch, mc['implementation'], mc['classical_helper'], mc['common_helper']):
        verify(p)
    require(mm['status'] == 'complete' and mm['config'] == c['MM_config'] and len(mm['rows']) == 1
            and mc['trials'] == [c['trial']] and mc['ledgers'] == {m: c['MM_ledger'] for m in METALS}
            and mm['stage'] == 'mm' and mm['implementation'] == mc['implementation']
            and mc['implementation']['sha256'] in {p['sha256'] for p in dispatch}
            and ev.cfg['base_implementation']['sha256'] == qc['helpers']['MM']['sha256']
            and mc['classical_helper'] == ledger['reused_builder'], 'MM execution/ledger lineage differs')
    for p in ledger['artifacts'].values():
        verify(p)
    verify(ledger['native_system'])
    mr = mm['rows'][0]
    require(mr['status'] == 'complete' and mr['trial'] == c['trial'] and mr['metal'] == c['metal']
            and mr['task_id'] == trial['task_id'] and set(mr['components_kcal_mol']) == set(GROUPS)
            and set(mr['component_gradients_kcal_mol_A']) == set(GROUPS), 'MM state/groups differ')
    mi = {s:i for i,s in enumerate(mr['source_ids'])}
    require(len(mi) == len(mr['source_ids']) == 2049 and set(mi) == set(ids), 'MM full source IDs differ')
    order = [mi[s] for s in ids]
    require(all(np.asarray(v).shape == (2049,3) and np.isfinite(v).all()
                for v in mr['component_gradients_kcal_mol_A'].values()), 'Full finite MM component gradients required')
    mg = sum(np.asarray(v)[order] for v in mr['component_gradients_kcal_mol_A'].values())
    require(math.fsum(mr['components_kcal_mol'].values()) == mr['total_kcal_mol'], 'MM native energy closure differs')
    sol = read(c['ALPB']); ex = read(sol['execution'])
    require(ex['returncode'] == 0 and ex['command'] == [qc['ALPB_binary']['path'], sol['xyz']['path'], sol['charges']['path']]
            and sol['charges'] == c['ALPB_charges']
            and sol['charges']['sha256'] == origins[c['metal']]['solvent_charges']['sha256'], 'Native fixed-charge solvent differs')
    for k in ('gradient', 'xyz', 'charges', 'execution'):
        verify(sol[k])
    stdout_pin = pin(verify(sol['execution']).parent/'stdout.txt')
    fields = {}
    for line in verify(stdout_pin).read_text().splitlines():
        if line.split() and line.split()[0] in ('COMPONENTS', 'CHARGE_SUM', 'ROUNDING_BOUNDS', 'CAVITY_GLOBAL'):
            key, *values = line.split(); require(key not in fields, 'Duplicated solvent field')
            fields[key] = list(map(float, values))
    require(dict(zip(('Gelec','Ghb','Gsasa','Gshift'), fields['COMPONENTS'])) == sol['components_hartree']
            and fields['CHARGE_SUM'] == [sol['charge_e']] and fields['ROUNDING_BOUNDS'] == [0., 0.]
            and fields['CAVITY_GLOBAL'] == sol['cavity_global']
            and math.fsum(sol['components_hartree'].values()) == sol['total_hartree'], 'Raw solvent energy/precision differs')
    sx = native['affordable_common'].xyz(verify(sol['xyz']))
    require(len(sx) == 2049 and [a[0] for a in sx] == [a['element'] for a in physical]
            and np.array_equal([a[1:] for a in sx], x), 'Solvent full physical geometry differs')
    q = np.loadtxt(verify(sol['charges'])); rawsg = np.loadtxt(verify(sol['gradient']))
    require(q.shape == (2049,) and np.isfinite(q).all() and abs(math.fsum(q)-2.) <= 1e-5
            and rawsg.shape == (2049, 3) and np.isfinite(rawsg).all()
            and np.array_equal(rawsg, sol['gradient_hartree_per_bohr']), 'Raw solvent charge/gradient differs')
    sg = rawsg * native['affordable_common'].HA_TO_KCAL / native['affordable_common'].BOHR_TO_A
    ag = read(c['approximate_assembled_gradient'])
    require(ag['source_ids'] == ids and ag['physical_mapping'] == c['mapping']
            and point['assembled_gradient'] == c['approximate_assembled_gradient'], 'Approximate gradient provenance differs')
    for name, g in (('MM_ledger', mg), ('fixed_FF_charge_ALPB', sg)):
        require(g.shape == (2049,3) and np.isfinite(g).all()
                and np.array_equal(g, ag['component_gradients_kcal_mol_A'][name]), 'Archived shared gradient differs: '+name)
    require(point['native_energies']['MM_groups_kcal_mol'] == mr['components_kcal_mol']
            and point['native_energies']['ALPB_Eh'] == sol['total_hartree'], 'Current native shared energies differ')
    xyz = native['affordable_common'].xyz(verify(c['core_xyz']))
    numbers = native['metal_environment_electronic_state'].NUMBERS
    state = c['reference_state']; original_task = c['reference_origin_task']
    require(original_task['electronic_state'] == state and original_task['metal'] == c['metal']
            and len(xyz) == 209 and [a[0] for a in xyz].count(c['metal']) == 1
            and sum(numbers[a[0]] for a in xyz) == state['all_electron_count']
            and state['charge'] == 0 and state['explicit_electrons'] == 832
            and state['all_electron_count']-state['ecp_core_electrons'] == 832
            and state['effective_multiplicity'] == 1
            and state['physical_multiplicity'] == (1 if c['metal'] == 'La' else 6)
            and state['ecp_core_electrons'] == (46 if c['metal'] == 'La' else 55), 'Native current electronic inventory/state differs')
    return {'MM_result': c['MM_result'], 'MM_config': c['MM_config'], 'MM_ledger': c['MM_ledger'],
            'ALPB': c['ALPB'], 'ALPB_stdout': stdout_pin, 'source_ids': ids,
            'native_energies': {'MM_kcal_mol': mr['total_kcal_mol'], 'MM_groups_kcal_mol': mr['components_kcal_mol'],
                                'ALPB_hartree': sol['total_hartree']},
            'component_gradients_kcal_mol_A': {'MM_ledger': mg.tolist(), 'fixed_FF_charge_ALPB': sg.tolist()},
            'approximate_point': c['approximate_point'], 'approximate_assembled_gradient': c['approximate_assembled_gradient']}


def freeze_directions(inv, np):
    directions = {}
    for metal in METALS:
        before = next(c for c in inv['cells'] if c['task_id'] == metal+'_terminal__'+metal)
        after = next(c for c in inv['cells'] if c['task_id'] == metal+'_cartesian_terminal__'+metal)
        a, b = read(before['physical']), read(after['physical'])
        require([v['id'] for v in a] == [v['id'] for v in b], 'Physical response chord order differs')
        delta = np.asarray([v['xyz_A'] for v in b])-np.asarray([v['xyz_A'] for v in a])
        length = float(np.linalg.norm(delta)); require(math.isfinite(length) and length > 0, 'Degenerate response chord')
        directions[metal] = {'before': before['physical'], 'after': after['physical'],
            'source_ids': [v['id'] for v in a], 'displacement_A': delta.tolist(),
            'Cartesian_L2_length_A': length, 'unit_direction': (delta/length).tolist()}
    return {'rule': SETTINGS['projection_rule'], 'selected_before_new_reference_outcomes': True, 'directions': directions}


def prepare(cfg, out, result, context):
    inv, qc, qual, world, origins, native = context
    ev, replay = source_replay(inv, out)
    result['source_replay'] = write(out/'SOURCE_REPLAY.json', replay)
    result['directions'] = write(out/'DIRECTIONS.json', freeze_directions(inv, ev.np))
    original_manifest = read(inv['reference_manifest']); tasks = []; taskpins = []; shared = {}
    for c in inv['cells']:
        tid = c['task_id']; d = out/tid; d.mkdir()
        shared[tid] = write(d/'SHARED.json', shared_cell(c, ev, origins, qc, native))
        original = next(t for t in original_manifest['tasks'] if t['metal'] == c['metal'])
        require(original == c['reference_origin_task'], 'Original task lineage differs')
        for p, name in ((c['core_xyz'], 'core.xyz'), (c['pointcharges'], 'environment.pc'), (original['input'], 'endpoint.inp')):
            (d/name).write_bytes(verify(p).read_bytes())
        text = (d/'endpoint.inp').read_text()
        require(text.startswith('! '+original_manifest['method']+'\n') and 'DoEQ false' in text
                and '* xyzfile 0 1 core.xyz' in text, 'Exact original embedded analytic deck required')
        task = {'task_id': tid, 'metal': c['metal'], 'charge': 0, 'multiplicity': 1,
                'physical_multiplicity': c['reference_state']['physical_multiplicity'],
                'electronic_state': c['reference_state'], 'task_type': 'analytic_gradient',
                'input': pin(d/'endpoint.inp'), 'xyz': pin(d/'core.xyz'), 'pointcharges': pin(d/'environment.pc'),
                'output_path': str(d/'endpoint.out'), 'engrad_path': str(d/'endpoint.engrad'),
                'pointcharge_gradient_path': str(d/'endpoint.runtime.pcgrad'),
                'geometry_alias': c['geometry_alias'], 'coordinate_sha256': c['coordinate_sha256'],
                'origin_task': original, 'physical': c['physical'], 'mapping': c['mapping'],
                'source_inputs': c['source_inputs'], 'core_mapping': c['core_mapping'], 'boundary_mapping': c['boundary_mapping'],
                'approximate_point': c['approximate_point'], 'shared_components': shared[tid], 'caveats': CAVEATS}
        task['cache_key'] = hashlib.sha256(json.dumps({'protocol': PROTOCOL, 'binding': binding(cfg),
            'task_state': c['reference_state'], 'physical': c['physical'], 'mapping': c['mapping'],
            'input_sha256': task['input']['sha256'], 'core': c['core_xyz'], 'field': c['pointcharges'],
            'orca': qc['orca'], 'assets': qc['basis_assets'], 'native': qc['native_implementation'],
            'MM': c['MM_result'], 'ALPB': c['ALPB']}, sort_keys=True).encode()).hexdigest()
        taskpins.append(write(d/'TASK.json', task)); tasks.append(task)
    for alias in ALIASES:
        pair = [read(t['physical']) for t in tasks if t['geometry_alias'] == alias]
        require(len(pair) == 2 and [a['id'] for a in pair[0]] == [a['id'] for a in pair[1]]
                and ev.np.array_equal([a['xyz_A'] for a in pair[0]], [a['xyz_A'] for a in pair[1]]),
                'Both-metal full physical coordinate pairing differs')
    header = {'protocol_id': PROTOCOL, 'agreement': cfg['specification'], 'plan': cfg['plan'],
        'method': original_manifest['method'], 'orca': qc['orca'], 'assets': qc['basis_assets'],
        'implementation': qc['native_implementation'], 'gradient_helpers': cfg['gradient_helpers'],
        'execution_policy': {k: qc['native_implementation'][v] for k,v in
            (('task_runner','run_orca_task_manifest.py'), ('runtime_renderer','render_orca_runtime_input.py'))},
        'execution_resources': {'mpi_ranks': 12, 'concurrent_tasks': 1}, 'tasks': tasks}
    result.update(status='prepared', tasks=taskpins, shared_components=shared,
        full_manifest=write(out/'FULL_MANIFEST.json', header), input_admission_passed=True,
        execution_authorized=False, future_gate_presence=gate_presence(cfg))


def prepared(cfg, context):
    inv, qc, qual, world, origins, native = context
    p = read(cfg['preparation_result']); pcfg = read(p['config'])
    require(p['status'] == 'prepared' and p['stage'] == 'prepare' and p['protocol'] == PROTOCOL
            and p['input_admission_passed'] and p['new_DFT_calls_attempted'] == 0
            and binding(pcfg) == binding(cfg), 'Actual exact zero-call preparation required')
    full = read(p['full_manifest']); tasks = {read(t)['task_id']: read(t) for t in p['tasks']}
    require(tuple(tasks) == IDS and list(tasks.values()) == full['tasks'], 'Eight-task prepared manifest differs')
    require(full['method'] == read(inv['reference_manifest'])['method'] and full['orca'] == qc['orca']
            and full['assets'] == qc['basis_assets'] and full['implementation'] == qc['native_implementation'], 'Prepared native recipe differs')
    for c in inv['cells']:
        t = tasks[c['task_id']]
        for k in ('physical','mapping','source_inputs','core_mapping','boundary_mapping','approximate_point'):
            require(t[k] == c[k], 'Prepared source provenance differs: '+k); verify(t[k])
        require(t['task_type'] == 'analytic_gradient' and t['electronic_state'] == c['reference_state']
                and t['shared_components'] == p['shared_components'][t['task_id']], 'Prepared state/components differ')
        for a,b in (('xyz',c['core_xyz']), ('pointcharges',c['pointcharges']), ('input',c['reference_origin_task']['input'])):
            require(verify(t[a]).read_bytes() == verify(b).read_bytes(), 'Prepared native bytes differ: '+a)
        read(t['shared_components'])
    require(read(p['directions']) == freeze_directions(inv, world[-1]), 'Frozen response directions differ')
    return p, full, tasks


def gate_presence(cfg):
    return {k: cfg.get(k) is not None for k in ('qualification_result','first_scout_decision',
                                              'response_scout_result','remaining_decision')}


def qualification_gate(cfg, out, context):
    inv, qc, qual, world, origins, native = context
    actual = read(cfg['qualification_result']); ac = read(actual['config'])
    require(actual['protocol'] == qual.PROTOCOL and actual['stage'] == 'collect'
            and actual['status'] == 'complete' and actual['scope'] == 'all_eight'
            and actual['all_eight_directional_derivatives_passed'] is True
            and actual['selected_tasks'] == list(qual.IDS)
            and actual['implementation']['sha256'] == QUALIFIER_SHA, 'Actual full derivative qualification required')
    for k in ('protocol','settings','implementation','plan','specification','origin_result','python','helpers',
              'native_implementation','orca','basis_assets','MM_template','ALPB_binary','ALPB_library_directory'):
        if k in ('implementation','specification','native_implementation'):
            require(source_hashes(ac[k]) == source_hashes(qc[k]), 'Derivative/reference source bytes differ: '+k)
        else:
            require(ac[k] == qc[k], 'Derivative/reference recipe differs: '+k)
    replay = {}; qual.collect(ac, out, replay)
    for k in replay:
        require(replay[k] == actual[k], 'Saved full derivative replay differs: '+k)
    require(replay['status'] == 'complete' and len(replay['checks']) == 4
            and all(v['passed'] for v in replay['checks']), 'Full derivative replay failed')
    return {'qualification_result': cfg['qualification_result'], 'scope': 'all_eight',
            'checks': replay['checks'], 'molecular_wall_seconds': {t: replay['rows'][t]['DFT']['molecular_wall_seconds'] for t in qual.IDS},
            'derivative_tolerance': qual.SETTINGS, 'new_molecular_calls': 0}


def preflight_gate(cfg):
    f = read(cfg['preflight_result'])
    require(f['protocol'] == PROTOCOL and f['stage'] == 'preflight' and f['status'] == 'complete'
            and f['input_admission_passed'] and f['preparation_result'] == cfg['preparation_result']
            and binding(read(f['config'])) == binding(cfg) and f['new_DFT_calls_attempted'] == 0,
            'Exact saved zero-call preflight required')
    return f


def decision_gate(cfg, proof, scope):
    key = 'first_scout_decision' if scope == 'first_cost_scout' else 'remaining_decision'
    d = read(cfg[key]); expected = [SCOUT] if scope == 'first_cost_scout' else list(IDS[1:])
    require(d['status'] == 'approved' and d['scope'] == scope and d['protocol'] == PROTOCOL
            and d['allowed_tasks'] == expected and d['qualification_result'] == cfg['qualification_result']
            and d['preparation_result'] == cfg['preparation_result'] and d['preflight_result'] == cfg['preflight_result']
            and d['qualification_molecular_wall_seconds'] == proof['molecular_wall_seconds']
            and isinstance(d['measured_cost_review'], str) and d['measured_cost_review'].strip()
            and d['cost_evidence'], 'Pinned root measured-cost decision is missing/incompatible')
    for p in d['cost_evidence']:
        verify(p)
    return d


def admit_gradient(manifest_pin, qc, qual, native):
    manifest = read(manifest_pin); require(len(manifest['tasks']) == 1, 'One native endpoint per execution')
    task = manifest['tasks'][0]
    require(task['task_id'] in IDS and task['task_type'] == 'analytic_gradient' and task['engrad_path'], 'Analytic endpoint required')
    row = qual.admit_sp(manifest_pin, qc, native)
    try:
        if row['status'] != 'available':
            return canonical(row)
        receipt = read(row['receipt'])
        require(receipt['expected_engrad_complete'] and receipt['manifest'] == manifest_pin
                and receipt['task_id'] == task['task_id'] and receipt['orca_version'] == '6.1.1', 'Native analytic receipt differs')
        parsed = native['metal_environment_frozen_f'].parse(task, {'electronic_state': task['electronic_state']})
        require(parsed['energy_hartree'] == row['energy_hartree'] and parsed['engrad'] == receipt['artifacts']['engrad'],
                'Analytic/native energy or gradient artifact differs')
        raw = native['affordable_response'].read_engrad(verify(parsed['engrad']))
        pc = pin(task['pointcharge_gradient_path'])
        pg = native['metal_environment_reference'].read_pcgrad(verify(pc), 1840)
        mp = read(task['mapping'])
        require(task['xyz']['sha256'] == mp['xyz']['sha256'] and task['pointcharges']['sha256'] == mp['pointcharges']['sha256'],
                'Native analytic current core/field differs')
        mapped = native['metal_environment_force_assembly'].scatter(mp, raw['gradient_Ha_per_bohr'], pg)
        tokens = [t for line in verify(parsed['engrad']).read_text().splitlines()
                  if line.strip() and not line.lstrip().startswith('#') for t in line.split()]
        pc_tokens = [t for line in verify(pc).read_text().splitlines()[1:] for t in line.split()]
        precision = {}
        for name, values in (('QM', tokens[2:2+3*209]), ('pointcharges', pc_tokens)):
            quantum = [float(Decimal(1).scaleb(Decimal(t.replace('D','E').replace('d','e')).as_tuple().exponent)) for t in values]
            precision[name] = {'minimum_printed_quantum_Ha_per_bohr': min(quantum),
                               'maximum_printed_quantum_Ha_per_bohr': max(quantum), 'is_error_bound': False}
        row.update(parsed, status='available', pointcharge_gradient=pc, pointcharge_gradient_units='Hartree/bohr',
                   mapped_physical_gradient=mapped, printed_gradient_precision=precision, caveats=CAVEATS)
    except Exception as exc:
        row.update(status='unavailable', admission_failure=repr(exc), traceback=traceback.format_exc())
    return canonical(row)


def force_report(task, dft, shared, directions, np):
    ag = read(shared['approximate_assembled_gradient']); ids = shared['source_ids']
    g = np.asarray(dft['mapped_physical_gradient']['gradient_kcal_mol_A'])
    components = {'PBE0_embedded': g, **{k: np.asarray(v) for k,v in shared['component_gradients_kcal_mol_A'].items()}}
    total = sum(components.values()); approximate = np.asarray(ag['gradient_kcal_mol_A'])
    require(total.shape == approximate.shape == (2049,3) and np.isfinite(total).all()
            and np.isfinite(approximate).all(), 'Full physical force unavailable')
    per_atom = []
    for i, sid in enumerate(ids):
        per_atom.append({'source_id': sid, 'reference_force_kcal_mol_A': (-total[i]).tolist(),
            'approximate_force_kcal_mol_A': (-approximate[i]).tolist(),
            'reference_minus_approximate_force_kcal_mol_A': (approximate[i]-total[i]).tolist(),
            'reference_force_norm_kcal_mol_A': float(np.linalg.norm(total[i])),
            'approximate_force_norm_kcal_mol_A': float(np.linalg.norm(approximate[i])),
            'force_discrepancy_norm_kcal_mol_A': float(np.linalg.norm(approximate[i]-total[i]))})
    projections = {}
    for metal, d in directions['directions'].items():
        require(d['source_ids'] == ids, 'Frozen projection atom order differs')
        u = np.asarray(d['unit_direction']); a = float(np.sum(total*u)); b = float(np.sum(approximate*u))
        projections[metal] = {'reference_dE_ds_kcal_mol_A': a, 'approximate_dE_ds_kcal_mol_A': b,
            'reference_minus_approximate_kcal_mol_A': a-b,
            'within_descriptive_agreement_target': abs(a-b) <= SETTINGS['projection_agreement_target_kcal_mol_A'],
            'reference_components_kcal_mol_A': {k:float(np.sum(v*u)) for k,v in components.items()},
            'approximate_components_kcal_mol_A': {k:float(np.sum(np.asarray(v)*u)) for k,v in ag['component_gradients_kcal_mol_A'].items()}}
    return {'quantity': 'gradient=dE/dR; force=-gradient', 'source_ids': ids,
        'gradient_kcal_mol_A': total.tolist(), 'component_gradients_kcal_mol_A': {k:v.tolist() for k,v in components.items()},
        'per_atom': per_atom, 'projections': projections,
        'maximum_reference_force_kcal_mol_A': max(r['reference_force_norm_kcal_mol_A'] for r in per_atom),
        'maximum_force_discrepancy_kcal_mol_A': max(r['force_discrepancy_norm_kcal_mol_A'] for r in per_atom),
        'atom_vector_RMS_reference_force_kcal_mol_A': float(np.sqrt(np.mean(np.sum(total*total, axis=1)))),
        'atom_vector_RMS_force_discrepancy_kcal_mol_A': float(np.sqrt(np.mean(np.sum((total-approximate)**2, axis=1)))),
        'force_error_qualified': False, 'SCF_error_bound_available': False,
        'precision_policy': 'Raw engrad/pcgrad and SCF printouts retained. Printed digits are not an electronic-force error bound.'}


def assemble(task, dft, shared, directions, np, ha):
    values = dict(shared['native_energies'], PBE0_embedded_hartree=dft['energy_hartree'])
    energies = {'PBE0_embedded': values['PBE0_embedded_hartree']*ha, 'MM_ledger': values['MM_kcal_mol'],
                'fixed_FF_charge_ALPB': values['ALPB_hartree']*ha}
    return {'native_energies': values, 'component_energies_kcal_mol': energies,
            'total_energy_kcal_mol': math.fsum(energies.values()),
            'forces': force_report(task, dft, shared, directions, np), 'caveats': CAVEATS}


def replay_endpoint(cfg, record, context, prep, tasks):
    inv, qc, qual, world, origins, native = context
    r = read(record); tid = r.get('task_id'); require(tid in IDS, 'Unknown response endpoint')
    rc = read(r['config'])
    require(r['protocol'] == PROTOCOL and r['stage'] == 'endpoint' and r['status'] == 'complete'
            and binding(rc) == binding(cfg) and r['preparation_result'] == cfg['preparation_result']
            and r['qualification_result'] == rc['qualification_result'] == cfg['qualification_result']
            and r['new_DFT_calls_attempted'] == r['new_DFT_calls_requested'] == 1
            and all(r[k] == 0 for k in ('new_MM_calls','new_ALPB_calls','new_MACE_calls','new_GFN_calls','optimizer_steps')),
            'Admitted response endpoint required')
    proof = read(r['qualification_replay'])
    require(proof['qualification_result'] == cfg['qualification_result'] and proof['scope'] == 'all_eight'
            and len(proof['checks']) == 4 and all(v['passed'] for v in proof['checks']), 'Endpoint derivative-gate provenance differs')
    preflight_gate(rc); decision_gate(rc, proof, 'first_cost_scout')
    if tid != SCOUT:
        decision = decision_gate(rc, proof, 'remaining_seven')
        require(decision['response_scout_result'] == rc['response_scout_result']
                and decision['first_scout_decision'] == rc['first_scout_decision'], 'Endpoint cost-gate lineage differs')
    expected = tasks[tid]; task = read(r['task']); manifest = read(r['manifest'])
    require(manifest['tasks'] == [task] and task['prepared_task'] == next(p for p in prep['tasks'] if read(p)['task_id'] == tid),
            'Executed task/preparation link differs')
    for k in expected:
        if k in ('xyz','input','pointcharges'):
            require(verify(task[k]).read_bytes() == verify(expected[k]).read_bytes(), 'Executed bytes differ: '+k)
        elif k not in ('output_path','engrad_path','pointcharge_gradient_path'):
            require(task[k] == expected[k], 'Executed state/geometry differs: '+k)
    dft = admit_gradient(r['manifest'], qc, qual, native)
    require(dft['status'] == 'available' and dft == r['DFT'], 'Native analytic replay failed')
    shared = read(task['shared_components'])
    assembled = assemble(task, dft, shared, read(prep['directions']), world[-1], qual.HA)
    require(assembled == read(r['assembled']), 'Saved full-gradient/native-energy assembly differs')
    return r, assembled


def endpoint(cfg, out, result, context, tid):
    require(tid in IDS and cfg['endpoint_results'][tid] is None, 'Unknown/already recorded endpoint; no automatic retry')
    inv, qc, qual, world, origins, native = context
    prep, full, tasks = prepared(cfg, context); preflight_gate(cfg)
    proof = qualification_gate(cfg, out, context)
    first = decision_gate(cfg, proof, 'first_cost_scout')
    if tid != SCOUT:
        scout, assembled = replay_endpoint(cfg, cfg['response_scout_result'], context, prep, tasks)
        require(scout['task_id'] == SCOUT and scout['new_DFT_calls_attempted'] == 1, 'Actual admitted first scout required')
        decision = decision_gate(cfg, proof, 'remaining_seven')
        require(decision['response_scout_result'] == cfg['response_scout_result']
                and decision['first_scout_decision'] == cfg['first_scout_decision']
                and decision['response_scout_molecular_wall_seconds'] == scout['DFT']['molecular_wall_seconds'],
                'Actual response-scout cost decision differs')
    result['qualification_replay'] = write(out/'QUALIFICATION_REPLAY.json', proof)
    task = copy.deepcopy(tasks[tid])
    for k, name in (('xyz','core.xyz'), ('pointcharges','environment.pc'), ('input','endpoint.inp')):
        (out/name).write_bytes(verify(task[k]).read_bytes()); task[k] = pin(out/name)
    task.update(output_path=str(out/'endpoint.out'), engrad_path=str(out/'endpoint.engrad'),
                pointcharge_gradient_path=str(out/'endpoint.runtime.pcgrad'),
                prepared_task=next(p for p in prep['tasks'] if read(p)['task_id'] == tid))
    tp = write(out/'TASK.json', task)
    manifest = dict(full, tasks=[task], full_prepared_manifest=prep['full_manifest'], preparation_result=cfg['preparation_result'])
    mp = write(out/'MANIFEST.json', manifest)
    result.update(task_id=tid, task=tp, manifest=mp, preparation_result=cfg['preparation_result'],
        directions=prep['directions'], qualification_result=cfg['qualification_result'],
        first_scout_decision=cfg['first_scout_decision'], remaining_decision=cfg.get('remaining_decision'))
    runtime = module('mex_response_runtime', qc['helpers']['reference_runtime'])
    result['allocation'] = runtime.allocation(native['render_orca_runtime_input'], out)
    write(out/'START.json', result)
    try:
        result['new_DFT_calls_requested'] = 1
        native['run_orca_task_manifest'].run_manifest(verify(mp), orca_path=verify(qc['orca']), workers=1, nprocs=12,
            selected_tasks={tid}, expected_runner_sha256=qc['native_implementation']['run_orca_task_manifest.py']['sha256'])
    except Exception as exc:
        result['runner_failure'] = repr(exc)
    finally:
        result['new_DFT_calls_attempted'] = int((out/'endpoint.out').exists())
        for k in ('xyz','input','pointcharges'):
            verify(task[k])
        result['DFT'] = admit_gradient(mp, qc, qual, native)
    if result['DFT']['status'] == 'available':
        value = assemble(task, result['DFT'], read(task['shared_components']), read(prep['directions']), world[-1], qual.HA)
        result.update(status='complete', assembled=write(out/'ASSEMBLED.json', value))
        result['rows'][tid].update(status='complete', task=tp, assembled=result['assembled'])
    else:
        result['status'] = 'unavailable'
        result['rows'][tid].update(status='unavailable', task=tp, DFT=result['DFT'])


def work(new, old, ha):
    return {'PBE0_embedded': float((Decimal(new['DFT']['printed_energy'])-Decimal(old['DFT']['printed_energy']))*Decimal(str(ha))),
            'MM_ledger': math.fsum(new['native_energies']['MM_groups_kcal_mol'][k]-old['native_energies']['MM_groups_kcal_mol'][k] for k in GROUPS),
            'fixed_FF_charge_ALPB': (new['native_energies']['ALPB_hartree']-old['native_energies']['ALPB_hartree'])*ha}


def approximate_work(value, origin, units):
    return {'MACE_vacuum': (value['MACE_eV']-origin['MACE_eV'])*units.EV,
            'GFN2_field_difference': math.fsum([value['GFN2_embedded_Eh']-origin['GFN2_embedded_Eh'],
                 -(value['GFN2_isolated_Eh']-origin['GFN2_isolated_Eh'])])*units.HA,
            'MM_ledger': math.fsum(value['MM_groups_kcal_mol'][k]-origin['MM_groups_kcal_mol'][k] for k in origin['MM_groups_kcal_mol']),
            'fixed_FF_charge_ALPB': (value['ALPB_Eh']-origin['ALPB_Eh'])*units.HA}


def collect(cfg, out, result, context):
    inv, qc, qual, world, origins, native = context
    prep, full, tasks = prepared(cfg, context)
    rows = result['rows']; costs = []
    static = module('mex_response_units', qc['helpers']['ALPB'])
    require(static.HA == qual.HA, 'Existing native energy-unit conversion differs')
    if any(cfg['endpoint_results'].values()):
        result['qualification_replay'] = write(out/'QUALIFICATION_REPLAY.json', qualification_gate(cfg, out, context))
    for tid in IDS:
        record = cfg['endpoint_results'][tid]
        if record is None:
            rows[tid].update(status='not_available', reason='No actual endpoint pin supplied'); continue
        rows[tid]['result'] = record
        try:
            raw = read(record)
            costs.append({'task_id': tid, 'result': record, 'slurm_job_id': raw.get('slurm_job_id'),
                'status': raw.get('status'), 'molecular_wall_seconds': raw.get('DFT', {}).get('molecular_wall_seconds'),
                'allocation': raw.get('allocation'), 'requested': raw.get('new_DFT_calls_requested'),
                'attempted': raw.get('new_DFT_calls_attempted'), 'cost_unavailable_is_not_zero': True})
            r, a = replay_endpoint(cfg, record, context, prep, tasks)
            require(r['task_id'] == tid, 'Result slot/task differs')
            rows[tid].update(status='complete', DFT=r['DFT'], assembled=r['assembled'], native_energies=a['native_energies'])
        except Exception as exc:
            rows[tid].update(status='unavailable', failure=repr(exc), traceback=traceback.format_exc())
    result['actual_endpoint_costs'] = costs; result['directions'] = prep['directions']
    comparisons = []
    for before in ALIASES[:2]:
        for after in ALIASES[2:]:
            row = {'before': before, 'after': after, 'status': 'unavailable', 'metals': {}}
            for metal in METALS:
                a, b = rows[after+'__'+metal], rows[before+'__'+metal]
                if a['status'] != 'complete' or b['status'] != 'complete':
                    continue
                ref = work(a, b, qual.HA); ca = next(c for c in inv['cells'] if c['task_id'] == after+'__'+metal)
                cb = next(c for c in inv['cells'] if c['task_id'] == before+'__'+metal)
                na, nb = read(ca['approximate_point'])['native_energies'], read(cb['approximate_point'])['native_energies']
                ap = approximate_work(na, nb, static)
                rt, at = math.fsum(ref.values()), math.fsum(ap.values())
                row['metals'][metal] = {'reference_components_kcal_mol': ref, 'approximate_components_kcal_mol': ap,
                    'reference_work_kcal_mol': rt, 'approximate_work_kcal_mol': at, 'reference_minus_approximate_kcal_mol': rt-at,
                    'reference_minus_replaced_electronic_work_kcal_mol': ref['PBE0_embedded']-math.fsum([ap['MACE_vacuum'],ap['GFN2_field_difference']]),
                    'within_reporting_resolution': abs(rt-at) <= SETTINGS['deformation_work_resolution_kcal_mol']}
            if set(row['metals']) == set(METALS):
                d = row['metals']; rd = d['Dy']['reference_work_kcal_mol']-d['La']['reference_work_kcal_mol']
                ad = d['Dy']['approximate_work_kcal_mol']-d['La']['approximate_work_kcal_mol']
                row.update(status='complete', reference_Dy_minus_La_work_kcal_mol=rd,
                    approximate_Dy_minus_La_work_kcal_mol=ad, differential_discrepancy_kcal_mol=rd-ad,
                    differential_within_reporting_resolution=abs(rd-ad) <= SETTINGS['deformation_work_resolution_kcal_mol'])
            comparisons.append(row)
    result['deformation_comparisons'] = comparisons
    own = {m: next(r for r in comparisons if r['before'] == m+'_terminal' and r['after'] == m+'_cartesian_terminal')['metals'].get(m) for m in METALS}
    result['own_metal_path_differential'] = {'status': 'unavailable'}
    if all(own.values()):
        rd = own['Dy']['reference_work_kcal_mol']-own['La']['reference_work_kcal_mol']
        ad = own['Dy']['approximate_work_kcal_mol']-own['La']['approximate_work_kcal_mol']
        result['own_metal_path_differential'] = {'status': 'complete', 'reference_kcal_mol': rd,
            'approximate_kcal_mol': ad, 'discrepancy_kcal_mol': rd-ad,
            'within_reporting_resolution': abs(rd-ad) <= SETTINGS['deformation_work_resolution_kcal_mol']}
    complete = all(r['status'] == 'complete' for r in rows.values())
    result['paired_four_geometry_pool'] = {'status': 'unavailable', 'required_tasks': list(IDS)}
    if complete:
        methods = {}
        for method in ('embedded_PBE0_plus_shared_MM_ALPB','approximate_original_model'):
            menus = {}
            for metal in METALS:
                base = rows['La_terminal__'+metal]
                approximate_base = read(next(c for c in inv['cells'] if c['task_id'] == 'La_terminal__'+metal)['approximate_point'])['native_energies']
                candidates = {}
                for alias in ALIASES:
                    if method == 'embedded_PBE0_plus_shared_MM_ALPB':
                        parts = work(rows[alias+'__'+metal], base, qual.HA)
                    else:
                        point = read(next(c for c in inv['cells'] if c['task_id'] == alias+'__'+metal)['approximate_point'])
                        parts = approximate_work(point['native_energies'], approximate_base, static)
                    candidates[alias] = math.fsum(parts.values())
                minimum = min(candidates.values()); oldmin = min(candidates[a] for a in ALIASES[:2])
                menus[metal] = {'work_from_La_terminal_kcal_mol': candidates, 'minimum_work_kcal_mol': minimum,
                    'old_two_geometry_minimum_kcal_mol': oldmin, 'additional_work_from_old_two_kcal_mol': minimum-oldmin,
                    'exact_ties': [a for a in ALIASES if candidates[a] == minimum],
                    'within_1kcal_resolution': [a for a in ALIASES if candidates[a]-minimum <= 1.],
                    'selection_policy': 'mathematical minimum; exact ties retained; resolution band descriptive'}
            methods[method] = {'metals': menus, 'Dy_minus_La_additional_accommodation_kcal_mol':
                menus['Dy']['additional_work_from_old_two_kcal_mol']-menus['La']['additional_work_from_old_two_kcal_mol']}
        result['paired_four_geometry_pool'] = {'status': 'complete', 'methods': methods, 'caveats': CAVEATS}
    result.update(status='complete' if complete else 'partial', preparation_result=cfg['preparation_result'],
        comparison_targets_are_native_gates=False, new_molecular_calls=0)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--stage', choices=('prepare','preflight','endpoint','collect'), required=True)
    p.add_argument('--task', choices=IDS); p.add_argument('--config', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True); args = p.parse_args()
    require(os.environ.get('SLURM_JOB_ID', '').isdigit(), 'Slurm required')
    out = args.output.resolve(); out.mkdir(parents=True, exist_ok=False); start = time.monotonic()
    cfg = json.loads(args.config.read_text())
    result = {'protocol': PROTOCOL, 'stage': args.stage, 'status': 'execution_failure', 'implementation': pin(__file__),
        'config': pin(args.config), 'settings': SETTINGS, 'slurm_job_id': os.environ['SLURM_JOB_ID'],
        'declared_cells': 8, 'rows': {t: {'task_id':t, 'status':'unavailable'} for t in IDS},
        'new_DFT_calls_requested': 0, 'new_DFT_calls_attempted': 0, 'new_MM_calls': 0, 'new_ALPB_calls': 0,
        'new_MACE_calls': 0, 'new_GFN_calls': 0, 'optimizer_steps': 0, 'caveats': CAVEATS}
    try:
        context = inputs(cfg)
        if args.stage == 'prepare':
            prepare(cfg, out, result, context)
        elif args.stage == 'preflight':
            prep, full, tasks = prepared(cfg, context)
            ev, replay = source_replay(context[0], out)
            for c in context[0]['cells']:
                require(shared_cell(c, ev, context[4], context[1], context[5]) == read(prep['shared_components'][c['task_id']]),
                        'Prepared saved shared components differ')
            result.update(status='complete', input_admission_passed=True, preparation_result=cfg['preparation_result'],
                          source_replay=write(out/'SOURCE_REPLAY.json', replay), future_gate_presence=gate_presence(cfg), execution_authorized=False)
            if cfg.get('qualification_result') is not None:
                result['qualification_replay'] = write(out/'QUALIFICATION_REPLAY.json', qualification_gate(cfg, out, context))
        elif args.stage == 'endpoint':
            endpoint(cfg, out, result, context, args.task)
        else:
            collect(cfg, out, result, context)
        verify(result['config'])
    except Exception as exc:
        result.update(status_before_failure=result['status'], status='execution_failure',
                      failure=repr(exc), traceback=traceback.format_exc())
    result['elapsed_seconds'] = time.monotonic()-start
    result['application_allocated_CPU_seconds'] = result['elapsed_seconds']*int(os.environ.get('SLURM_CPUS_ON_NODE','1'))
    result['parent_peak_RSS_KiB'] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    result['parent_RSS_is_not_MPI_peak'] = True
    write(out/'RESULT.json', result)
    print(json.dumps({'status': result['status'], 'result': pin(out/'RESULT.json')}))
    return 0 if result['status'] in ('prepared','complete','partial') else 1


if __name__ == '__main__':
    raise SystemExit(main())
