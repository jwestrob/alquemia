#!/usr/bin/env python3
"""Eight fixed whole2049 native MACE cells with four saved fixed-FF ALPB cells.

Every stage requires Slurm. No native solvent/MM/GFN call or optimizer exists here.
The pinned native loader is reused unchanged; only size/state checks are generalized.
"""
import argparse
import hashlib
import importlib.metadata
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

PROTOCOL = 'mex_whole2049_fixed_motion_two_spin_hypotheses_v1'
ROOT = Path('/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs')
EV, HA, BOHR = 23.06054783061903, 627.509474, .529177210903
SPINS = {'aligned': {'La': 10, 'Dy': 15},
         'one_Nd_reversed': {'La': 4, 'Dy': 9}}
IDS = tuple(f'{g}_{m}_{h}' for h in SPINS for g in ('before', 'after') for m in ('La', 'Dy'))
SCOUT = 'before_La_aligned'
ALIASES = {'before': 'La_terminal', 'after': 'La_cartesian_terminal'}
NATIVE_SHA = '0e48b8c6c2f271512f2c445f826036d8a68f7dca0aabc5923697ab3ecb7d9896'
POOL_SHA = '9d5439a42f3cd14bb4fdd63f883cf6e807ccb1752735de1d59cdd13bf5938010'
SEARCH_SHA = 'c494798a7f2a8bf17d51b642a2345cb9673f165e52f607a9ee3044bfd6cdda62'
SOURCE_SHA = '4b3b98983af0cb5fcdbf7b8f52e35063308201ce9df777a30484acf5146b7930'
BIND = ('protocol', 'implementation', 'specification', 'MACE_config', 'native_helper',
        'model_inspection', 'source_inputs', 'pool_result', 'search_result', 'solvent_parent', 'runtime')
LIMITS = {'scout_forward_seconds': 600., 'scout_peak_device_fraction': .90,
          'maximum_new_forwards': 8, 'reporting_resolution_kcal_mol': 1.}
ION_IDS = ['A/203//ND', 'A/201//ND', 'A/202//ND', 'A/204//ND']


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
    s = importlib.util.spec_from_file_location(name, verify(p))
    obj = importlib.util.module_from_spec(s); s.loader.exec_module(obj)
    return obj


def binding(c):
    return {k: c[k]['sha256'] if k == 'implementation' else c[k] for k in BIND}


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False).encode()).hexdigest()


def vectors(x, n=2049):
    require(len(x) == n and all(len(r) == 3 and all(math.isfinite(v) for v in r) for r in x),
            'Full finite physical vectors required')
    return x


def state_check(task, symbols, positions, numbers, elements, embeddings):
    require(task['task_id'] in IDS and task['task_id'] == f"{task['geometry']}_{task['metal']}_{task['hypothesis']}"
            and task['geometry'] in ALIASES and task['hypothesis'] in SPINS, 'Task outside declared eight states')
    require(task['charge'] == 2 and task['multiplicity'] == SPINS[task['hypothesis']][task['metal']]
            and task['metal_index'] == 2045 and task['atom_count'] == 2049, 'Declared aggregate charge/spin differs')
    require(len(symbols) == 2049 and symbols[2045:] == [task['metal'], 'Nd', 'Nd', 'Nd']
            and set(symbols[:2045]) <= {'H', 'C', 'N', 'O'} and symbols.count(task['metal']) == 1,
            'Whole physical atom/ion inventory differs')
    vectors(positions)
    require({numbers[s] for s in symbols} <= set(elements), 'Unsupported checkpoint species')
    electrons = sum(numbers[s] for s in symbols) - task['charge']
    require(electrons >= task['multiplicity']-1 and (electrons-task['multiplicity']+1) % 2 == 0,
            'All-electron spin parity differs')
    for name, value in (('total_charge', task['charge']), ('total_spin', task['multiplicity'])):
        spec = embeddings[name]
        require(spec['type'] == 'categorical' and value == int(value)
                and 0 <= value+spec.get('offset', 0) < spec['num_classes'], 'Checkpoint category differs')
    return electrons


def batch_check(state, task):
    require(state['charge'] == 2 and state['spin_multiplicity'] == task['multiplicity']
            and state['head'] == 'omol' and state['selected_metal_index'] == 2045
            and state['atoms'] == 2049 and state['graph_count'] == 1
            and state['metal_neighbor_edge_count'] > 0 and state['charge_prediction_claimed'] is False,
            'Actual native graph state differs')


def admit_solvent(point, atoms, mapping, q, binary, np):
    sol = read(point['ALPB']); ex = read(sol['execution'])
    for key in ('gradient', 'xyz', 'charges'):
        verify(sol[key])
    start_pin = pin(verify(sol['execution']).parent/'START.json'); start = read(start_pin)
    command = [binary['path'], sol['xyz']['path'], sol['charges']['path']]
    require(ex['returncode'] == 0 and ex['command'] == command and start['command'] == command
            and start['xyz'] == sol['xyz'] and start['charges'] == sol['charges']
            and start['physical_mapping'] == point['mapping'], 'Saved native solvent execution differs')
    rows = [line.split() for line in verify(sol['xyz']).read_text().splitlines()[2:] if line.strip()]
    require(int(verify(sol['xyz']).read_text().splitlines()[0]) == len(rows) == 2049
            and all(len(r) == 4 for r in rows)
            and [r[0] for r in rows] == [a['element'] for a in atoms]
            and [[float(v) for v in r[1:]] for r in rows] == [a['xyz_A'] for a in atoms], 'ALPB physical geometry differs')
    rawq = np.loadtxt(verify(sol['charges'])); rawg = np.loadtxt(verify(sol['gradient']))
    require(rawq.shape == (2049,) and np.isfinite(rawq).all() and rawq.tolist() == q
            and abs(math.fsum(rawq)-2.) <= 1e-5 and rawg.shape == (2049,3) and np.isfinite(rawg).all()
            and rawg.tolist() == sol['gradient_hartree_per_bohr'], 'Saved raw solvent charge/gradient differs')
    stdout_pin = pin(verify(sol['execution']).parent/'stdout.txt'); fields = {}
    for line in verify(stdout_pin).read_text().splitlines():
        parts = line.split()
        if parts and parts[0] in ('COMPONENTS', 'CHARGE_SUM', 'ROUNDING_BOUNDS', 'CAVITY_GLOBAL'):
            require(parts[0] not in fields, 'Duplicate solvent native field')
            fields[parts[0]] = list(map(float, parts[1:]))
    require(len(fields['COMPONENTS']) == 4
            and dict(zip(('Gelec','Ghb','Gsasa','Gshift'), fields['COMPONENTS'])) == sol['components_hartree']
            and fields['CHARGE_SUM'] == [sol['charge_e']] and fields['ROUNDING_BOUNDS'] == [0.,0.]
            and fields['CAVITY_GLOBAL'] == sol['cavity_global']
            and math.fsum(sol['components_hartree'].values()) == sol['total_hartree']
            and sol['total_hartree'] == point['native_energies']['ALPB_Eh'], 'Raw solvent energy/precision replay differs')
    require(sol['fixed_charge_rule'] == 'unmodified ff19SB/TIP3P physical charges plus four formal3 ions; no cap population projection',
            'Fixed solvent recipe differs')
    g = rawg*HA/BOHR  # Exact original Mex producer operation order.
    ag = read(point['assembled_gradient'])
    require(ag['source_ids'] == [a['id'] for a in atoms] and ag['physical_mapping'] == point['mapping']
            and np.array_equal(g, ag['component_gradients_kcal_mol_A']['fixed_FF_charge_ALPB']),
            'Archived solvent full-gradient replay differs')
    return {'result': point['ALPB'], 'execution': sol['execution'], 'start': start_pin, 'stdout': stdout_pin,
            'xyz': sol['xyz'], 'charges': sol['charges'], 'gradient': sol['gradient'],
            'energy_hartree': sol['total_hartree'], 'gradient_kcal_mol_A': g.tolist()}


def prepare(cfg):
    require(os.environ.get('SLURM_JOB_ID', '').isdigit(), 'Every stage requires Slurm')
    require(cfg['protocol'] == PROTOCOL and cfg['limits'] == LIMITS and cfg['spin_hypotheses'] == SPINS
            and pin(__file__)['sha256'] == cfg['implementation']['sha256'], 'Frozen protocol/settings differ')
    for k in BIND[1:-1]: verify(cfg[k])
    require(sys.prefix == cfg['runtime']['sys_prefix'] and sys.executable == cfg['runtime']['python_executable'],
            'Qualified native virtual environment differs')
    require(cfg['native_helper']['sha256'] == NATIVE_SHA and cfg['pool_result']['sha256'] == POOL_SHA
            and cfg['search_result']['sha256'] == SEARCH_SHA and cfg['source_inputs']['sha256'] == SOURCE_SHA,
            'Fixed parent/source selection differs')
    import numpy as np
    from ase.data import atomic_numbers
    mc = read(cfg['MACE_config']); inspect = read(cfg['model_inspection'])
    require(pin(sys.executable)['sha256'] == mc['python']['sha256']
            and {k: importlib.metadata.version(k) for k in mc['expected_versions']} == mc['expected_versions'],
            'Qualified Python/software differs')
    for record in mc['helper_pins'].values(): verify(record)
    verify(mc['model']['checkpoint'])
    require(mc['model']['checkpoint'] == inspect['checkpoint'] and inspect['class'] == 'ScaleShiftMACE'
            and mc['model']['embedding_specs'] == inspect['embedding_specs'] and inspect['parameter_dtypes'] == ['torch.float64']
            and mc['model']['dtype'] == 'float64' and mc['model']['execution_adapter'] is None
            and mc['model']['solvent'] is None and not mc['model']['extra_dispersion']
            and mc['model']['spin_input_convention'] == 'multiplicity' and not mc['model']['periodic'],
            'Native model recipe/inspection differs')
    source = read(cfg['source_inputs']); export = read(source['export']); baseatoms = read(export['atoms'])
    require(source['source_id'] == 'Mex_8FNS' and source['waters']['total_source'] == 167
            and source['charge_accounting']['total_e'] == 2 and export['atom_count'] == len(baseatoms) == 2045,
            'Original physical inventory differs')
    verify(source['waters']['inventory'])
    for p in export['forcefields']: verify(p)
    expected_ids = [a['id'] for a in baseatoms]+ION_IDS
    q = [a['charge_e'] for a in baseatoms]+[3.,3.,3.,3.]
    static = read(cfg['solvent_parent']); binary = static['build']['binary']; verify(binary)
    require(static['status'] == 'complete' and static['build']['returncode'] == 0, 'Original ALPB build unavailable')
    pool = read(cfg['pool_result']); run = read(cfg['search_result']); search = read(run['search'])
    require(pool['status'] == 'complete' and str(pool['slurm_job_id']) == '12981'
            and run['status'] == 'complete' and str(run['slurm_job_id']) == '12970' and run['metal'] == 'La'
            and cfg['search_result'] in pool['search_results'] and pool['physical_source_ids'] == expected_ids
            and search['physical_state_admitted'] and search['final_geometry']['pass']
            and search['accepted_steps'] == 80 and not search['converged'], 'Admitted fixed search/pool differs')
    cells = {}; coords = {}
    for geometry, alias in ALIASES.items():
        matches = [r for r in pool['matrix'] if alias in r['aliases']]
        require(len(matches) == 1 and set(matches[0]['cells']) == {'La','Dy'}, 'Fixed paired geometry unavailable')
        row = matches[0]
        for metal, c in row['cells'].items():
            point = read(c['point']); atoms = read(point['physical']); mapping = read(point['mapping'])
            require(point['status'] == 'complete' and point['metal'] == metal and c['physical'] == point['physical']
                    and c['native_energies'] == point['native_energies'] and c['assembled_gradient'] == point['assembled_gradient']
                    and mapping['metal'] == metal and mapping['physical_atoms'] == atoms
                    and [a['id'] for a in atoms] == expected_ids and len(set(expected_ids)) == 2049
                    and not any(a['id'].startswith('cap/') for a in atoms), 'Exact whole physical POINT/map differs')
            require([a['element'] for a in atoms[:2045]] == [a['element'] for a in baseatoms]
                    and [a['element'] for a in atoms[2045:]] == [metal,'Nd','Nd','Nd'], 'Spectators/elements changed')
            if metal == 'La':
                saved = search['start_point' if geometry == 'before' else 'terminal_point']['point']
                require(c['point'] == saved, 'Search start/terminal selection differs')
            key = f'{geometry}_{metal}'; symbols = [a['element'] for a in atoms]; x = [a['xyz_A'] for a in atoms]
            cells[key] = {'point': c['point'], 'physical': point['physical'], 'mapping': point['mapping'],
                          'source_ids': expected_ids, 'symbols': symbols, 'positions_A': x,
                          'old_native_energies': point['native_energies'],
                          'solvent': admit_solvent(point, atoms, mapping, q, binary, np)}
            coords[key] = x
        require(coords[f'{geometry}_La'] == coords[f'{geometry}_Dy'], 'Paired metal coordinates differ')
    tasks = []
    for hypothesis in SPINS:
        for geometry in ALIASES:
            for metal in ('La','Dy'):
                cell = cells[f'{geometry}_{metal}']
                task = {'task_id': f'{geometry}_{metal}_{hypothesis}', 'geometry': geometry, 'metal': metal,
                        'hypothesis': hypothesis, 'charge': 2, 'multiplicity': SPINS[hypothesis][metal],
                        'metal_index': 2045, 'atom_count': 2049, 'physical': cell['physical'],
                        'old_point': cell['point'], 'ALPB': cell['solvent']['result']}
                task['all_electron_count'] = state_check(task, cell['symbols'], cell['positions_A'], atomic_numbers,
                                                       inspect['elements'], inspect['embedding_specs'])
                task['cache_key'] = digest({'binding': binding(cfg), 'task': task, 'symbols': cell['symbols'],
                                            'positions_A': cell['positions_A']})
                tasks.append(task)
    require(tuple(t['task_id'] for t in tasks) == IDS, 'Eight-state denominator differs')
    controls = []
    for name in ('wrong_atom', 'wrong_charge', 'wrong_spin'):
        t = dict(tasks[0]); cell = cells['before_La']; symbols = list(cell['symbols'])
        if name == 'wrong_atom': symbols[2046] = 'La'
        elif name == 'wrong_charge': t['charge'] = 0
        else: t['multiplicity'] = 1
        try:
            state_check(t, symbols, cell['positions_A'], atomic_numbers, inspect['elements'], inspect['embedding_specs'])
        except ValueError as exc: controls.append({'control': name, 'rejected': True, 'reason': str(exc)})
        else: raise ValueError('Negative state control was accepted: '+name)
    manifest = {'protocol': PROTOCOL, 'binding': binding(cfg), 'tasks': tasks, 'source_ids': expected_ids,
                'reused_ALPB': {k:c['solvent'] for k,c in cells.items()}, 'negative_controls': controls,
                'no_caps': True, 'water_count': 167, 'field_MM_GFN_terms': 'absent', 'limits': LIMITS}
    return manifest, cells, mc, np


def preflight_gate(cfg, manifest):
    r = read(cfg['preflight_result'])
    require(r['status'] == 'preflight_passed' and r['stage'] == 'preflight' and r['binding'] == binding(cfg)
            and r['new_MACE_calls_attempted'] == 0 and r['native_model_loads_attempted'] == 0
            and read(r['manifest']) == manifest, 'Actual zero-call preflight differs')


def native_class(helper):
    class WholeNative(helper.NativeMACE):
        # Initialization is inherited byte-for-byte. The model invocation is the same native forward.
        def evaluate(self, task, symbols, positions):
            tick = time.monotonic(); torch = self.torch
            state_check(task, symbols, positions, self.Z, self.model.atomic_numbers.tolist(), self.model.embedding_specs)
            atoms = self.Atoms(symbols, positions=positions, pbc=False)
            atoms.info.update(charge=task['charge'], spin=task['multiplicity'])
            batch = self.calc._atoms_to_batch(atoms)
            state = self.input_batch(self.calc, atoms, task['charge'], task['multiplicity'], batch, 2045)
            batch_check(state, task)
            forward_tick = time.monotonic()
            with torch.enable_grad():
                self.calls += 1
                value = self.model(batch.to_dict(), training=False, compute_force=True,
                                   compute_virials=False, compute_stress=False, compute_displacement=False,
                                   compute_hessian=False, compute_edge_forces=False, compute_atomic_stresses=False)
            energy = float(value['energy'].detach().cpu().item()); force = value['forces'].detach().cpu()
            forward_seconds = time.monotonic()-forward_tick
            require(math.isfinite(energy) and tuple(force.shape) == (2049,3) and bool(torch.isfinite(force).all()),
                    'Nonfinite native full energy/force')
            require([p._version for p in self.model.parameters()] == self.parameter_versions, 'Model parameters changed')
            return {'energy_eV': energy, 'forces_eV_A': force.tolist(), 'input_state_check': state,
                    'force_convention': 'negative derivative of the same native vacuum energy',
                    'forward_seconds': forward_seconds, 'wall_seconds': time.monotonic()-tick,
                    'peak_GPU_allocated_bytes': torch.cuda.max_memory_allocated(),
                    'device_total_memory_bytes': torch.cuda.get_device_properties(0).total_memory}
    return WholeNative


def assemble(task, raw, cell, np):
    batch_check(raw['input_state_check'], task); vectors(raw['forces_eV_A'])
    require(math.isfinite(raw['energy_eV']) and raw['force_convention'] == 'negative derivative of the same native vacuum energy',
            'Raw MACE energy/force convention differs')
    mg = -np.asarray(raw['forces_eV_A'])*EV; sg = np.asarray(cell['solvent']['gradient_kcal_mol_A'])
    parts = {'whole_MACE_vacuum': raw['energy_eV']*EV, 'fixed_FF_charge_ALPB': cell['solvent']['energy_hartree']*HA}
    return {'source_ids': cell['source_ids'], 'native_energies': {'MACE_eV': raw['energy_eV'],
            'ALPB_Eh': cell['solvent']['energy_hartree']}, 'component_energies_kcal_mol': parts,
            'total_energy_kcal_mol': math.fsum(parts.values()),
            'component_gradients_kcal_mol_A': {'whole_MACE_vacuum': mg.tolist(), 'fixed_FF_charge_ALPB': sg.tolist()},
            'gradient_kcal_mol_A': (mg+sg).tolist(), 'forces_kcal_mol_A': (-(mg+sg)).tolist()}


def admit_saved(row, task, cell, np):
    require(row['status'] == 'complete' and row['task'] == task and row['ALPB'] == cell['solvent']['result'],
            'Saved task/solvent provenance differs')
    raw = read(row['raw_MACE']); assembled = read(row['assembled'])
    require(assembled == assemble(task, raw, cell, np), 'Saved native/component/total replay differs')
    return raw, assembled


def scout_gate(cfg, manifest, cells, np):
    r = read(cfg['scout_result'])
    require(r['status'] == 'scout_passed' and r['stage'] == 'scout' and r['binding'] == binding(cfg)
            and r['preflight_result'] == cfg['preflight_result'] and read(r['manifest']) == manifest
            and r['new_MACE_calls_attempted'] == r['native_forward_calls'] == 1
            and r['native_model_loads_attempted'] == 1 and r['scout_continuation_gate']['passed'],
            'Actual admitted one-call scout required')
    rows = {x['task_id']: x for x in r['rows']}
    require(list(rows) == list(IDS) and all(rows[t]['status'] == 'unavailable' for t in IDS if t != SCOUT),
            'Scout touched an undeclared cell')
    raw, _ = admit_saved(rows[SCOUT], manifest['tasks'][0], cells['before_La'], np)
    require(raw['forward_seconds'] <= LIMITS['scout_forward_seconds']
            and raw['peak_GPU_allocated_bytes'] <= LIMITS['scout_peak_device_fraction']*raw['device_total_memory_bytes'],
            'Measured scout continuation budget exceeded')
    return rows[SCOUT]


def run(stage, cfg, manifest, cells, mc, np, out, result):
    preflight_gate(cfg, manifest)
    require(1 <= int(os.environ['SLURM_CPUS_PER_TASK']) <= 12, 'At most twelve allocated CPUs required')
    result['preflight_result'] = cfg['preflight_result']
    if stage == 'remaining':
        result['rows'][0] = dict(scout_gate(cfg, manifest, cells, np), reused_from=cfg['scout_result'])
        result['scout_result'] = cfg['scout_result']
    helper = module('mex_whole_native_loader', cfg['native_helper'])
    result['native_model_loads_attempted'] = 1; tick = time.monotonic()
    native = native_class(helper)(mc)
    result['model_load_seconds'] = time.monotonic()-tick
    result['device'] = native.device; result['versions'] = native.versions; result['model'] = mc['model']
    require('H200' in native.device, 'Declared H200 GPU required')
    try:
        tasks = manifest['tasks'][:1] if stage == 'scout' else manifest['tasks'][1:]
        for task in tasks:
            tid = task['task_id']; idx = IDS.index(tid); cell = cells[f"{task['geometry']}_{task['metal']}"]
            d = out/tid; d.mkdir()
            write(d/'START.json', {'task': task, 'binding': binding(cfg), 'ALPB': cell['solvent']['result'],
                                  'new_native_calls_ceiling_for_cell': 1})
            result['new_MACE_calls_attempted'] += 1
            result['rows'][idx] = {'task_id': tid, 'task': task, 'status': 'unavailable', 'reason': 'native_request_failed'}
            raw = native.evaluate(task, cell['symbols'], cell['positions_A'])
            rp = write(d/'MACE.json', raw); ap = write(d/'ASSEMBLED.json', assemble(task, raw, cell, np))
            result['rows'][idx] = {'task_id': tid, 'task': task, 'status': 'complete', 'raw_MACE': rp,
                                   'ALPB': cell['solvent']['result'], 'assembled': ap}
            if stage == 'scout':
                gate = {'forward_seconds': raw['forward_seconds'], 'peak_GPU_allocated_bytes': raw['peak_GPU_allocated_bytes'],
                        'device_total_memory_bytes': raw['device_total_memory_bytes'],
                        'passed': raw['forward_seconds'] <= LIMITS['scout_forward_seconds']
                        and raw['peak_GPU_allocated_bytes'] <= LIMITS['scout_peak_device_fraction']*raw['device_total_memory_bytes']}
                result['scout_continuation_gate'] = gate
                result['status'] = 'scout_passed' if gate['passed'] else 'scout_budget_blocked'
        if stage == 'remaining': result['status'] = 'complete'
    finally:
        result['native_forward_calls'] = native.calls
        result['peak_GPU_allocated_bytes'] = native.torch.cuda.max_memory_allocated()
        result['device_total_memory_bytes'] = native.torch.cuda.get_device_properties(0).total_memory


def collect(cfg, manifest, cells, np, result):
    preflight_gate(cfg, manifest); scout = scout_gate(cfg, manifest, cells, np)
    parent = read(cfg['remaining_result'])
    require(parent['status'] == 'complete' and parent['stage'] == 'remaining' and parent['binding'] == binding(cfg)
            and parent['preflight_result'] == cfg['preflight_result'] and parent['scout_result'] == cfg['scout_result']
            and read(parent['manifest']) == manifest and parent['new_MACE_calls_attempted'] == parent['native_forward_calls'] == 7
            and parent['rows'][0] == dict(scout, reused_from=cfg['scout_result']), 'Full finite remaining panel unavailable')
    require([r['task_id'] for r in parent['rows']] == list(IDS), 'Saved eight-cell denominator differs')
    data = {}
    for task, row in zip(manifest['tasks'], parent['rows']):
        raw, assembled = admit_saved(row, task, cells[f"{task['geometry']}_{task['metal']}"], np)
        data[task['task_id']] = assembled
    result['rows'] = parent['rows']; result['parents'] = {'scout': cfg['scout_result'], 'remaining': cfg['remaining_result']}
    works = {}; force_changes = {}
    for h in SPINS:
        works[h] = {}; force_changes[h] = {}
        for metal in ('La','Dy'):
            a, b = (data[f'{g}_{metal}_{h}'] for g in ('after','before'))
            parts = {'whole_MACE_vacuum': (a['native_energies']['MACE_eV']-b['native_energies']['MACE_eV'])*EV,
                     'fixed_FF_charge_ALPB': (a['native_energies']['ALPB_Eh']-b['native_energies']['ALPB_Eh'])*HA}
            works[h][metal] = dict(parts, total_kcal_mol=math.fsum(parts.values()))
            force_changes[h][metal] = {}
            for comp in ('whole_MACE_vacuum', 'fixed_FF_charge_ALPB', 'total'):
                ga = a['gradient_kcal_mol_A'] if comp == 'total' else a['component_gradients_kcal_mol_A'][comp]
                gb = b['gradient_kcal_mol_A'] if comp == 'total' else b['component_gradients_kcal_mol_A'][comp]
                df = -np.asarray(ga)+np.asarray(gb)
                force_changes[h][metal][comp] = {'full_force_change_kcal_mol_A': df.tolist(),
                    'RMS_atom_vector_kcal_mol_A': float(np.sqrt(np.mean(np.sum(df*df, axis=1)))),
                    'maximum_atom_vector_kcal_mol_A': float(np.max(np.linalg.norm(df, axis=1)))}
        works[h]['Dy_minus_La'] = {k: works[h]['Dy'][k]-works[h]['La'][k] for k in works[h]['La']}
    sensitivity = {k: works['one_Nd_reversed'][k]['total_kcal_mol']-works['aligned'][k]['total_kcal_mol']
                   for k in ('La','Dy','Dy_minus_La')}
    old_work = {}
    for metal in ('La','Dy'):
        a, b = (cells[f'{g}_{metal}']['old_native_energies'] for g in ('after','before'))
        parts = {'local_MACE': (a['MACE_eV']-b['MACE_eV'])*EV,
                 'GFN_field': math.fsum((a['GFN2_embedded_Eh']-b['GFN2_embedded_Eh'],
                                         -a['GFN2_isolated_Eh']+b['GFN2_isolated_Eh']))*HA,
                 'MM': math.fsum(a['MM_groups_kcal_mol'][k]-b['MM_groups_kcal_mol'][k] for k in a['MM_groups_kcal_mol']),
                 'ALPB': (a['ALPB_Eh']-b['ALPB_Eh'])*HA}
        old_work[metal] = dict(parts, total_kcal_mol=math.fsum(parts.values()))
    result.update(status='complete', work_kcal_mol=works, force_changes=force_changes,
                  spin_hypothesis_work_changes_kcal_mol=sensitivity,
                  spin_dependence_above_reporting_resolution={k:abs(v)>1. for k,v in sensitivity.items()},
                  old_hybrid_work_comparator_kcal_mol=old_work,
                  total_new_MACE_forwards_across_panel=8, reused_ALPB_count=4,
                  interpretation='Aggregate conditioning sensitivity; no local-spin, accuracy, affinity or production claim')


def main():
    p = argparse.ArgumentParser(); p.add_argument('stage', choices=('preflight','scout','remaining','collect'))
    p.add_argument('--config', type=Path, required=True); p.add_argument('--output', type=Path, required=True)
    a = p.parse_args(); require(os.environ.get('SLURM_JOB_ID','').isdigit(), 'Slurm required')
    out = a.output.resolve(); out.mkdir(parents=True, exist_ok=False); tick = time.monotonic()
    result = {'protocol': PROTOCOL, 'status': 'execution_failure', 'stage': a.stage, 'config': pin(a.config),
              'implementation': pin(__file__), 'slurm_job_id': os.environ['SLURM_JOB_ID'], 'host': socket.gethostname(),
              'new_MACE_calls_attempted': 0, 'native_forward_calls': 0, 'native_model_loads_attempted': 0,
              'new_ALPB_calls': 0, 'new_MM_calls': 0, 'new_GFN_calls': 0, 'new_DFT_calls': 0, 'optimizer_steps': 0,
              'rows': [{'task_id': t, 'status': 'unavailable', 'reason': 'not_attempted'} for t in IDS]}
    try:
        cfg = json.loads(a.config.read_text()); result['binding'] = binding(cfg)
        manifest, cells, mc, np = prepare(cfg); result['manifest'] = write(out/'MANIFEST.json', manifest)
        if a.stage == 'preflight': result['status'] = 'preflight_passed'
        elif a.stage == 'collect': collect(cfg, manifest, cells, np, result)
        else: run(a.stage, cfg, manifest, cells, mc, np, out, result)
    except Exception as exc:
        result['status'] = 'execution_failure'; result['failure'] = f'{type(exc).__name__}: {exc}'
        result['traceback'] = traceback.format_exc()
    finally:
        torch = sys.modules.get('torch')
        if result['native_model_loads_attempted'] and torch is not None and torch.cuda.is_initialized():
            result['peak_GPU_allocated_bytes'] = torch.cuda.max_memory_allocated()
            result['device_total_memory_bytes'] = torch.cuda.get_device_properties(0).total_memory
        result['elapsed_seconds'] = time.monotonic()-tick
        result['process_CPU_seconds'] = resource.getrusage(resource.RUSAGE_SELF).ru_utime+resource.getrusage(resource.RUSAGE_SELF).ru_stime
        result['peak_parent_RSS_KiB'] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        result['allocated_CPUs'] = int(os.environ.get('SLURM_CPUS_PER_TASK','1'))
        result['native_environment'] = {k:os.environ.get(k) for k in ('CUDA_VISIBLE_DEVICES','OMP_NUM_THREADS',
            'OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','SLURM_JOB_GPUS','SLURM_STEP_GPUS')}
        write(out/'RESULT.json', result)
    return 0 if result['status'] in ('preflight_passed','scout_passed','complete') else 1


if __name__ == '__main__':
    sys.exit(main())
