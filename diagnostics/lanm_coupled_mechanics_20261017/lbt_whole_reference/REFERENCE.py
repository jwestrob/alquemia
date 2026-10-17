#!/usr/bin/env python3
"""Whole-LBT saved-input preflight and one fixed before-La PBE0 force scout."""
import argparse
import datetime as dt
import hashlib
import importlib
import importlib.util
import json
import math
import os
from pathlib import Path
import re
import resource
import socket
import sys
import time
import traceback

PROTOCOL = 'lbt_whole256_fixed_deformation_reference_scout_v1'
INVENTORY_SHA = '48d9f40f1bc7f5f927131d5f8e1b2ad3d1fa12187f7b8c321b7cb5375efe4baa'
PROPOSAL_SHA = '02b45065533355f06f53316e201e165971e0736f77130dde97dcc718bef2af75'
COMPACT_SHA = '18ed7c0f4141b95cba20191bf36867bc839ae8722945b81dfb74fe3eeeda7f25'
METHOD = 'PBE0 D4 def2-TZVP def2/J RIJCOSX NoAutostart DefGrid3 VeryTightSCF EnGrad'
IDS = ('before_La', 'after_La', 'before_Dy', 'after_Dy')
SCOUT = 'before_La'
PARTS = ('Gelec', 'Ghb', 'Gsasa', 'Gshift')
CAVEATS = {
    'terminal_chemistry': 'Conditional free NH3+/COO- termini and five deprotonated acidic side chains; experimental terminal specification unresolved.',
    'vacuum_dianion': 'Post-added fixed-FF ALPB does not polarize or stabilize the vacuum SCF density self-consistently; finite-basis convergence does not prove electron binding or a stable ground state.',
    'frozen_f': 'Physical La multiplicity1/Dy6; effective singlets with La ECP46/Dy ECP55, spin-free frozen-f representation without SOC.',
    'geometry': 'One fixed 7CCO La-proposed before/after deformation; both endpoints nonstationary, no geometry repair or optimization.',
    'reference_scope': 'Whole vacuum PBE0 replaces all internal electronic/MM terms; only unchanged fixed-FF ALPB is added. No CPCM or external point charges.',
    'force_qualification': 'Native analytic-gradient admission only; no numerical derivative test or statistical force-accuracy threshold in this four-cell proposal.',
    'interpretation': 'Independent conditional force/work diagnostic, not an affinity prediction, an exact MACE training-functional replay, or a unique MM-parameter attribution.',
}


def require(ok, message):
    if not ok:
        raise ValueError(message)


def pin(path):
    p = Path(path).resolve()
    with p.open('rb') as handle:
        digest = hashlib.file_digest(handle, 'sha256').hexdigest()
    return {'path': str(p), 'sha256': digest}


def verify(record):
    require(pin(record['path'])['sha256'] == record['sha256'], 'Pin changed: ' + record['path'])
    return Path(record['path'])


def read(record):
    return json.loads(verify(record).read_text())


def write(path, value):
    with Path(path).open('x') as handle:
        json.dump(value, handle, indent=2, allow_nan=False)
        handle.write('\n')
    return pin(path)


def canonical(value):
    return json.loads(json.dumps(value, allow_nan=False))


def load_module(name, record):
    spec = importlib.util.spec_from_file_location(name, verify(record))
    obj = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(obj)
    return obj


def scientific_binding(cfg):
    return {k: cfg[k] for k in ('protocol', 'plan', 'inventory', 'proposal', 'specification',
                               'python', 'permitted_execution_tasks', 'runtime_reference_helper', 'ALPB_binary')}


def inputs(cfg):
    require(os.environ.get('SLURM_JOB_ID', '').isdigit(), 'Slurm required, including saved preflight/collection')
    require(cfg['protocol'] == PROTOCOL and pin(__file__)['sha256'] == cfg['implementation']['sha256'], 'Adapter differs')
    verify(cfg['implementation']); verify(cfg['plan']); verify(cfg['proposal']); verify(cfg['specification'])
    require(cfg['inventory']['sha256'] == INVENTORY_SHA and cfg['proposal']['sha256'] == PROPOSAL_SHA,
            'Frozen selection/proposal changed')
    require(cfg['permitted_execution_tasks'] == [SCOUT], 'Only the declared before-La scout is executable')
    require(pin(sys.executable)['sha256'] == cfg['python']['sha256'], 'Python executable differs')
    inv = read(cfg['inventory'])
    require(inv['declared_cells'] == 4 and inv['declared_geometries'] == 2
            and [c['cell_id'] for c in inv['cells']] == list(IDS)
            and inv['primary_scout'] == SCOUT and inv['plan'] == cfg['plan'], 'Inventory denominator/plan differs')
    native = inv['native_reference']; compact = read(native['compact_assets_manifest'])
    require(native['compact_assets_manifest']['sha256'] == COMPACT_SHA and native['method'] == METHOD
            and native['metal_assets'] == compact['assets'] and native['orca'] == compact['orca'], 'Native recipe changed')
    require(cfg['runtime_reference_helper'] == native['runtime_reference_helper'], 'Runtime helper differs')
    for record in compact['implementation'].values():
        verify(record)
    verify(compact['orca']); verify(cfg['ALPB_binary'])
    for assets in compact['assets'].values():
        for record in assets.values():
            verify(record)
    directory = verify(compact['implementation']['metal_environment_frozen_f.py']).parent
    sys.path.insert(0, str(directory))
    names = ('metal_environment_frozen_f', 'run_orca_task_manifest', 'ggr_sensitivity',
             'mace_omol_vacuum', 'render_orca_runtime_input', 'affordable_response',
             'affordable_common', 'metal_environment_electronic_state')
    modules = {name: importlib.import_module(name) for name in names}
    for name, obj in modules.items():
        require(pin(obj.__file__)['sha256'] == compact['implementation'][name + '.py']['sha256'],
                'Imported native helper differs: ' + name)
    model = modules['metal_environment_frozen_f']
    require(model.METHOD == METHOD, 'Native analytic parser method differs')
    runtime = load_module('whole_lbt_original_runtime', cfg['runtime_reference_helper'])
    import numpy as np
    return inv, compact, modules, runtime, np


def state(symbols, metal, numbers, charge=-2):
    require(charge == -2 and len(symbols) == 256 and symbols[-1] == metal
            and metal in ('La', 'Dy') and all(s in {'C', 'H', 'N', 'O'} for s in symbols[:-1]),
            'Whole256 atom order/elements/charge differ')
    core = {'La': 46, 'Dy': 55}[metal]
    total = sum(numbers[s] for s in symbols) - charge
    explicit = total - core
    require(explicit > 0 and explicit % 2 == 0, 'Whole effective closed-shell electron parity differs')
    return dict(metal=metal, oxidation_state_hypothesis=3, charge=charge,
                physical_multiplicity=1 if metal == 'La' else 6,
                physical_f_occupation=0 if metal == 'La' else 9,
                effective_multiplicity=1, ecp_core_electrons=core,
                all_electron_count=total, explicit_electrons=explicit,
                alpha_electrons=explicit // 2, beta_electrons=explicit // 2,
                representation='spin_free_4f0_valence_model' if metal == 'La' else 'spin_free_frozen_4f9_valence_model',
                spin_orbit_included=False, physical_state_qualified=False,
                whole_physical_atom_count=256, target_index=255, external_charge_count=0, synthetic_cap_count=0)


def input_text(basis, aux):
    return f'! {METHOD}\n%scf\n Guess PModel\nend\n%basis\n{basis}\n{aux}\nend\n* xyzfile -2 1 physical.xyz\n'


def asset_checks(compact, metal):
    assets = compact['assets'][metal]
    basis, aux = verify(assets['basis']).read_text(), verify(assets['aux']).read_text()
    require(re.findall(r'(?mi)^\s*NewECP\s+(\w+)', basis) == [metal]
            and re.findall(r'(?mi)^\s*N_core\s+(\d+)', basis) == [str(46 if metal == 'La' else 55)]
            and re.findall(r'(?mi)^\s*NewGTO\s+(\w+)', basis) == [metal]
            and re.findall(r'(?mi)^\s*NewAuxJGTO\s+(\w+)', aux) == [metal], 'Declared metal basis/ECP/auxiliary differs')
    template = input_text(basis, aux)
    require(not re.search(r'%pal|%maxcore|pointcharges|CPCM|SMD|COSMO|\bOpt\b|MaxIter|DoEQ', template, re.I),
            'Unexpected whole-vacuum input control')
    return template


def physical_inputs(inv, np, model):
    response = read(inv['hybrid_response_result']); search = read(inv['hybrid_search'])
    require(response['status'] == search['status'] == 'complete'
            and response['source_id'] == '7CCO' and response['metal'] == 'La' and response['region'] == 'baseline'
            and response['search'] == inv['hybrid_search'] and search['physical_state_admitted']
            and not search['converged'] and search['accepted_steps'] == 80
            and search['termination'] == 'accepted_step_limit', 'Selected actual response differs')
    rows = {c['cell_id']: c for c in inv['cells']}
    require(search['start_point']['point'] == rows['before_La']['hybrid_point']
            and search['terminal_point']['point'] == rows['after_La']['hybrid_point'], 'Selected before/after points differ')
    pool = read(inv['hybrid_parent_pool']); whole = read(inv['whole_MACE_pool'])
    require(pool['status'] == whole['status'] == 'complete', 'Parent pools unavailable')
    cells = {}
    for tid in IDS:
        row = rows[tid]; parent = row if row['hybrid_point'] else rows['after_La']
        atoms = read(parent['physical']); positions = read(parent['positions'])
        ids = [a['id'] for a in atoms]
        require(len(atoms) == 256 and len(set(ids)) == 256 and ids[-1] == 'metal'
                and atoms[-1]['element'] == 'M' and positions['source_ids'] == ids
                and all(a['kind'] in ('protein_source', 'terminal_completion') for a in atoms[:-1])
                and [(a['id'], a['element']) for a in atoms if a['kind'] == 'terminal_completion']
                    == [('A/1//H2', 'H'), ('A/1//H3', 'H'), ('A/17//OXT', 'O')]
                and atoms[-1]['kind'] == 'selected_metal', 'Physical inventory/order differs')
        x = np.asarray(positions['positions_A'], dtype=float)
        require(x.shape == (256, 3) and np.isfinite(x).all()
                and np.array_equal(x, [a['xyz_A'] for a in atoms]), 'Physical positions differ')
        symbols = [a['element'] for a in atoms[:-1]] + [row['metal']]
        st = state(symbols, row['metal'], model.NUMBERS)
        charges = read(parent['global_solvent_charges'])
        require(charges['physical_ids'] == ids and len(charges['charges_e']) == 256
                and all(math.isfinite(q) for q in charges['charges_e'])
                and charges['charges_e'][-1] == 3.0 and abs(math.fsum(charges['charges_e']) + 2) <= 1e-5,
                'Whole fixed charges differ')
        if row['hybrid_point']:
            p = read(row['hybrid_point']); task = read(p['task'])
            require(p['status'] == 'complete' and p['positions'] == row['positions']
                    and task['physical_source'] == row['physical'] and p['mapping'] == row['hybrid_mapping']
                    and p['component_receipts']['ALPB'] == row['ALPB']
                    and p['state']['metal'] == row['metal'] and p['state']['whole_charge'] == -2
                    and p['point']['positions_A'] == x.tolist(), 'Saved hybrid/source join differs')
            if tid.startswith('before'):
                selected = [c for c in pool['cells'] if c['source_id'] == '7CCO' and c['region'] == 'baseline'
                            and c['candidate_id'] == 'cartesian_01' and c['metal'] == row['metal']]
                require(len(selected) == 1 and selected[0]['point'] == row['hybrid_point'], 'Before pooled source differs')
        cells[tid] = dict(row=row, atoms=atoms, ids=ids, x=x, symbols=symbols, state=st,
                          charges=charges, physical_parent=parent['physical'])
    require(np.array_equal(cells['before_La']['x'], cells['before_Dy']['x'])
            and np.array_equal(cells['after_La']['x'], cells['after_Dy']['x'])
            and all(c['ids'] == cells[SCOUT]['ids'] for c in cells.values()), 'Paired deformation differs')
    require(len({c['state']['explicit_electrons'] for c in cells.values()}) == 1, 'La/Dy frozen-f explicit counts differ')
    return cells


def admit_alpb(cell, cfg, np, model):
    record = cell['row']['ALPB']
    if record is None:
        return {'status': 'unavailable', 'reason': 'Selected after-La geometry as Dy has no bound solvent receipt; no call here.'}
    raw = read(record)
    for p in raw['artifacts'].values():
        verify(p)
    start = read(raw['artifacts']['START.json']); execution = read(raw['artifacts']['EXECUTION.json'])
    command = [str(verify(cfg['ALPB_binary'])), str(verify(start['xyz'])), str(verify(start['charges']))]
    require(start['command'] == execution['command'] == command and execution['returncode'] == 0,
            'Saved ALPB execution differs')
    physical = read(start['physical_source'])
    require(physical == read(cell['physical_parent']), 'Saved ALPB physical provenance differs')
    atoms = model.xyz(verify(start['xyz']))
    q = [float(v) for v in verify(start['charges']).read_text().split()]
    require([a[0] for a in atoms] == cell['symbols'] and np.array_equal([a[1:] for a in atoms], cell['x'])
            and q == cell['charges']['charges_e'], 'Saved ALPB coordinates/charge/metal differ')
    text = verify(raw['artifacts']['stdout.txt']).read_text()
    def values(label, size):
        found = re.findall(r'^' + label + r'\s+(.+)$', text, re.M)
        require(len(found) == 1, 'Missing native solvent field: ' + label)
        v = [float(x) for x in found[0].split()]
        require(len(v) == size and all(math.isfinite(x) for x in v), 'Invalid native solvent field')
        return v
    parts = dict(zip(PARTS, values('COMPONENTS', 4)))
    require(parts == raw['components_hartree'] and math.fsum(parts.values()) == raw['total_hartree']
            and values('CHARGE_SUM', 1)[0] == raw['charge_e'] and abs(raw['charge_e'] + 2) <= 1e-5
            and values('ROUNDING_BOUNDS', 2) == raw['rounding_bounds_hartree'] == [0., 0.]
            and values('CAVITY_GLOBAL', 4) == raw['cavity_global']
            and raw['gradient_units'] == 'Hartree/Bohr; archived only', 'Saved native solvent scalar replay differs')
    g = np.loadtxt(verify(raw['artifacts']['gradient.dat']))
    require(g.shape == (256, 3) and np.isfinite(g).all(), 'Saved full solvent gradient unavailable')
    return dict(status='available', receipt=record, energy_hartree=raw['total_hartree'],
                components_hartree=parts, gradient_hartree_per_bohr=g.tolist(),
                gradient_kcal_mol_A=(g * model.HA_TO_KCAL / model.BOHR_TO_A).tolist(), new_calls=0)


def directions_and_controls(cells, inv, model, np):
    x0, x1 = cells['before_La']['x'], cells['after_La']['x']
    chord = x1 - x0; length = float(np.linalg.norm(chord))
    require(length > 0 and math.isfinite(length), 'Declared deformation is zero/nonfinite')
    audit = read(inv['force_audit'])
    row = next(r for r in audit['rows'] if r['source_id'] == '7CCO' and r['metal'] == 'La')
    require(audit['status'] == 'complete' and row['status'] == 'available'
            and row['hybrid_point']['point'] == cells[SCOUT]['row']['hybrid_point'], 'Original force audit differs')
    ca = cells[SCOUT]['ids'].index('A/1//CA')
    atom = next(a for a in row['top20_atoms'] if a['id'] == 'A/1//CA')
    vector = np.asarray(atom['forces']['hybrid'], dtype=float)
    require(atom['index'] == ca and not atom['electronic_real'] and np.linalg.norm(vector) > 0,
            'Phe1 CA force target differs')
    local = np.zeros_like(chord); local[ca] = vector / np.linalg.norm(vector)
    controls = []
    for label, symbols, metal, charge in (
        ('wrong_total_charge', cells[SCOUT]['symbols'], 'La', -1),
        ('missing_physical_atom', cells[SCOUT]['symbols'][:-1], 'La', -2),
        ('wrong_target_species', cells[SCOUT]['symbols'], 'Dy', -2),
    ):
        try:
            state(symbols, metal, model.NUMBERS, charge)
        except ValueError as exc:
            controls.append({'name': label, 'rejected': True, 'reason': str(exc)})
        else:
            raise ValueError('Negative state control was admitted: ' + label)
    return dict(physical_ids=cells[SCOUT]['ids'], frozen_before_new_reference_calls=True,
                normalization='Euclidean norm over full physical Cartesian array',
                directions={'full_deformation_chord': (chord / length).tolist(),
                            'initial_hybrid_Phe1_CA_force': local.tolist()},
                selection_source=inv['force_audit'], negative_controls=controls,
                numerical_derivative_qualification=False)


def prepare(cfg, out, inv, compact, cells, model, np):
    tasks = []
    for tid in IDS:
        cell = cells[tid]; directory = out / tid; directory.mkdir()
        xyz = directory / 'physical.xyz'
        xyz.write_text('256\nFixed whole-LBT physical coordinates; no caps or external field\n' + ''.join(
            '%s %.17g %.17g %.17g\n' % (s, *v) for s, v in zip(cell['symbols'], cell['x'])))
        replay = model.xyz(xyz)
        require([a[0] for a in replay] == cell['symbols'] and np.array_equal([a[1:] for a in replay], cell['x']),
                'XYZ serialization changed exact physical input')
        deck = directory / 'endpoint.inp'; deck.write_text(asset_checks(compact, cell['row']['metal']))
        task = dict(task_id=tid, source_id='7CCO', geometry=cell['row']['geometry'], metal=cell['row']['metal'],
                    xyz=pin(xyz), input=pin(deck), output_path=str(directory / 'endpoint.out'),
                    engrad_path=str(directory / 'endpoint.engrad'), task_type='analytic_gradient',
                    charge=-2, multiplicity=1, physical_multiplicity=cell['state']['physical_multiplicity'],
                    electronic_state=cell['state'], physical_parent=cell['physical_parent'],
                    positions=cell['row']['positions'], physical_ids=cell['ids'], atom_count=256,
                    target_index=255, ALPB=cell['row']['ALPB'], caveats=CAVEATS)
        write(directory / 'TASK.json', task); tasks.append(task)
    full = dict(protocol_id=PROTOCOL, agreement=cfg['specification'], plan=cfg['plan'], method=METHOD,
                adapter=pin(__file__), inventory=cfg['inventory'], inputs_config=cfg['_config_pin'],
                compact_manifest=inv['native_reference']['compact_assets_manifest'], assets=compact['assets'],
                orca=compact['orca'], implementation=compact['implementation'], execution_policy=compact['execution_policy'],
                execution_resources={'mpi_ranks': 12, 'concurrent_tasks': 1}, tasks=tasks, caveats=CAVEATS,
                declared_cells=4, execution_allowlist=[SCOUT], no_MM_or_GFN_or_external_field=True)
    full_pin = write(out / 'FULL_MANIFEST.json', full)
    manifest = dict(full, full_manifest=full_pin, tasks=[t for t in tasks if t['task_id'] == SCOUT])
    return full_pin, write(out / 'MANIFEST.json', manifest)


def signatures(manifest):
    return [{k: t[k] for k in ('task_id', 'source_id', 'geometry', 'metal', 'charge', 'multiplicity',
                               'physical_multiplicity', 'electronic_state', 'physical_parent', 'positions',
                               'physical_ids', 'atom_count', 'target_index', 'ALPB', 'caveats')} |
            {'xyz_sha256': t['xyz']['sha256'], 'input_sha256': t['input']['sha256']} for t in manifest['tasks']]


def validate_prepared(full, cells, compact, cfg, model, np):
    require(full['protocol_id'] == PROTOCOL and full['method'] == METHOD and full['inventory'] == cfg['inventory']
            and full['orca'] == compact['orca'] and full['assets'] == compact['assets']
            and full['implementation'] == compact['implementation'] and full['execution_policy'] == compact['execution_policy']
            and full['adapter']['sha256'] == cfg['implementation']['sha256']
            and full['execution_allowlist'] == [SCOUT] and full['caveats'] == CAVEATS
            and [t['task_id'] for t in full['tasks']] == list(IDS), 'Prepared whole-reference manifest differs')
    verify(full['adapter'])
    for task in full['tasks']:
        c = cells[task['task_id']]; atoms = model.xyz(verify(task['xyz']))
        require(task['electronic_state'] == c['state'] and task['charge'] == -2 and task['multiplicity'] == 1
                and task['physical_multiplicity'] == c['state']['physical_multiplicity']
                and task['physical_parent'] == c['physical_parent'] and task['positions'] == c['row']['positions']
                and task['physical_ids'] == c['ids'] and task['ALPB'] == c['row']['ALPB']
                and task['atom_count'] == 256 and task['target_index'] == 255 and task['task_type'] == 'analytic_gradient'
                and task['source_id'] == '7CCO' and task['metal'] == c['row']['metal']
                and [a[0] for a in atoms] == c['symbols'] and np.array_equal([a[1:] for a in atoms], c['x'])
                and verify(task['input']).read_text() == asset_checks(compact, c['row']['metal']),
                'Prepared whole-reference exact state/deck/coordinates differ: ' + task['task_id'])


def preflight_gate(cfg, full, directions, solvent):
    prior = read(cfg['preflight_result']); pcfg = read(prior['config'])
    require(prior['protocol'] == PROTOCOL and prior['stage'] == 'preflight' and prior['status'] == 'preflight_passed'
            and prior['new_DFT_calls_attempted'] == 0 and prior['new_ALPB_calls'] == 0
            and prior['implementation']['sha256'] == cfg['implementation']['sha256']
            and scientific_binding(pcfg) == scientific_binding(cfg), 'Actual unchanged zero-call preflight required')
    require(signatures(read(prior['full_manifest'])) == signatures(full)
            and read(prior['directions']) == directions and prior['solvent_admission'] == solvent,
            'Zero-call preflight inputs/state/solvent do not replay')
    return cfg['preflight_result']


def admit(manifest_pin, compact, modules, runtime):
    model = modules['metal_environment_frozen_f']; collector = modules['ggr_sensitivity']
    residual = modules['mace_omol_vacuum']
    manifest, rows = collector.executed(verify(manifest_pin))
    require([t['task_id'] for t in manifest['tasks']] == [SCOUT], 'Scout execution denominator differs')
    task = manifest['tasks'][0]; row = rows[SCOUT]
    row.update(task=task, task_id=SCOUT, physical_multiplicity=task['physical_multiplicity'])
    try:
        if Path(task['output_path']).exists():
            row['SCF_evidence'] = runtime.observations(task['output_path'], residual)
        if row['status'] != 'complete':
            return canonical(row)
        parsed = model.parse(task, {'electronic_state': task['electronic_state']})
        require(math.isfinite(parsed['energy_hartree']) and len(parsed['gradient_kcal_mol_per_A']) == 256,
                'Full finite native gradient unavailable')
        text = Path(task['output_path']).read_text()
        require('INITIAL GUESS: MODEL POTENTIAL' in text and not re.search(r'Reading point charge file|^\s*(?:CPCM|SMD|COSMO) SOLVATION MODEL', text, re.M),
                'Vacuum/guess scope changed')
        require(row['SCF_evidence']['ConvCheckMode'] == 'Total+1el-Energy' and row['SCF_evidence']['ConvForced'] == 0,
                'Native SCF convergence semantics changed')
        receipt = read(row['receipt']); sidecar = read(receipt['artifacts']['runtime_sidecar'])
        require(receipt['orca_executable'] == compact['orca'] and receipt['returncode'] == 0
                and receipt['normal_termination'] and receipt['scf_converged'], 'Native execution failed')
        policy = sidecar['memory_policy']
        require(policy['workers'] == 1 and policy['ranks_per_worker'] == 12
                and policy['memory_source'] == 'SLURM_MEM_PER_NODE', 'Native allocation-derived memory policy differs')
        row.update(parsed, status='available', runtime_memory_policy=policy,
                   molecular_wall_seconds=(dt.datetime.fromisoformat(receipt['finished_at_utc']) -
                                           dt.datetime.fromisoformat(receipt['started_at_utc'])).total_seconds(), caveats=CAVEATS)
    except Exception as exc:
        row.update(status='unavailable', admission_failure=repr(exc), energy_hartree=None)
    return canonical(row)


def force_comparison(cell, native, solvent, directions, np, model):
    require(native['status'] == solvent['status'] == 'available', 'Native/solvent scout unavailable')
    g = np.asarray(native['gradient_kcal_mol_per_A']); gs = np.asarray(solvent['gradient_kcal_mol_A'])
    require(g.shape == gs.shape == (256, 3) and np.isfinite(g).all() and np.isfinite(gs).all(), 'Nonfinite full reference gradient')
    ref = g + gs
    hp = read(cell['row']['hybrid_point'])['point']
    wp = read(cell['row']['existing_whole_MACE_cell']['point'])
    require(wp['status'] == 'complete' and wp['source_id'] == '7CCO' and wp['metal'] == 'La'
            and wp['ALPB'] == cell['row']['ALPB'], 'Whole-MACE same-coordinate provenance differs')
    wg = np.load(verify(wp['gradient_kcal_mol_A']), allow_pickle=False)
    whole_positions = np.load(verify(wp['positions']), allow_pickle=False)
    require(np.array_equal(whole_positions, cell['x']), 'Whole-MACE coordinate join differs')
    gradient = {'reference_PBE0_vacuum': g, 'reference_PBE0_plus_fixed_ALPB': ref,
                'hybrid_composite': np.asarray(hp['gradient_kcal_mol_A']), 'whole_MACE_composite': wg}
    require(all(v.shape == (256, 3) and np.isfinite(v).all() for v in gradient.values()), 'Saved comparator gradient unavailable')
    def metrics(v):
        norm = np.linalg.norm(v, axis=1)
        return {'maximum_atom_kcal_mol_A': float(np.max(norm)),
                'RMS_atom_vector_kcal_mol_A': float(np.sqrt(np.mean(norm ** 2)))}
    ca = cell['ids'].index('A/1//CA')
    table = [{'index': i, 'physical_id': name, 'resname': cell['atoms'][i].get('resname'),
              'force_kcal_mol_A': {k: (-v[i]).tolist() for k, v in gradient.items()}} for i, name in enumerate(cell['ids'])]
    return dict(quantity='force=-gradient', units='kcal/mol/Angstrom',
                fixed_physical_input=cell['row']['positions'], force_table=table, Phe1_CA=table[ca],
                Phe1=[r for r in table if r['physical_id'].startswith('A/1/')],
                distributions={k: metrics(v) for k, v in gradient.items()},
                reference_minus_model_force_errors={k: metrics(ref - gradient[k]) for k in ('hybrid_composite', 'whole_MACE_composite')},
                force_projections={name: {k: float(np.sum(-v * np.asarray(direction))) for k, v in gradient.items()}
                                   for name, direction in directions['directions'].items()},
                reference_energy_components_kcal_mol={'whole_PBE0_vacuum': native['energy_hartree'] * model.HA_TO_KCAL,
                                                      'fixed_FF_charge_ALPB': solvent['energy_hartree'] * model.HA_TO_KCAL},
                energy_work_comparison='unavailable: only before-La scout is executable',
                numerical_force_accuracy_qualified=False, caveats=CAVEATS)


def run(args):
    out = Path(args.output).resolve(); out.mkdir(parents=True, exist_ok=False)
    start = time.monotonic()
    result = dict(protocol=PROTOCOL, stage=args.stage, status='execution_failure', implementation=pin(__file__),
                  config=pin(args.config), slurm_job_id=os.environ.get('SLURM_JOB_ID'), host=socket.gethostname(),
                  declared_cells=4, selected_tasks=[], new_DFT_calls_attempted=0, new_ALPB_calls=0,
                  native_runner_invocations=0,
                  new_MACE_calls=0, new_GFN_calls=0, new_MM_calls=0, optimizer_steps=0, caveats=CAVEATS,
                  rows={tid: {'task_id': tid, 'status': 'not_requested'} for tid in IDS})
    try:
        cfg = json.loads(Path(args.config).read_text()); cfg['_config_pin'] = pin(args.config)
        inv, compact, modules, runtime, np = inputs(cfg)
        model = modules['metal_environment_frozen_f']
        cells = physical_inputs(inv, np, model)
        solvent = {tid: admit_alpb(cell, cfg, np, model) for tid, cell in cells.items()}
        require(all(solvent[tid]['status'] == 'available' for tid in IDS[:3])
                and solvent['after_Dy']['status'] == 'unavailable', 'Frozen three-available/one-missing solvent inventory changed')
        directions = directions_and_controls(cells, inv, model, np)
        result.update(inventory=cfg['inventory'], plan=cfg['plan'], specification=cfg['specification'],
                      states={tid: c['state'] for tid, c in cells.items()}, solvent_admission=solvent,
                      physical_ids=cells[SCOUT]['ids'], scientific_binding=scientific_binding(cfg))
        result['directions'] = write(out / 'DIRECTIONS.json', directions)
        if args.stage == 'collect':
            prior = read(cfg['scout_result']); pcfg = read(prior['config'])
            require(prior['protocol'] == PROTOCOL and prior['stage'] == 'scout'
                    and prior['implementation']['sha256'] == cfg['implementation']['sha256']
                    and scientific_binding(pcfg) == scientific_binding(cfg)
                    and prior['selected_tasks'] == [SCOUT] and cfg['preflight_result'] == pcfg['preflight_result'],
                    'Declared actual scout required')
            result['scout_result'] = cfg['scout_result']; result['manifest'] = prior['manifest']
            full = read(prior['full_manifest'])
            validate_prepared(full, cells, compact, cfg, model, np)
            result['preflight_result'] = preflight_gate(cfg, full, directions, solvent)
            manifest = read(prior['manifest'])
            require(manifest['full_manifest'] == prior['full_manifest']
                    and manifest['tasks'] == [t for t in full['tasks'] if t['task_id'] == SCOUT]
                    and {k: v for k, v in manifest.items() if k not in ('full_manifest', 'tasks')}
                        == {k: v for k, v in full.items() if k != 'tasks'}, 'Actual one-task manifest differs from frozen whole panel')
            require(read(prior['directions']) == directions and prior['solvent_admission'] == solvent,
                    'Scout direction/solvent evidence differs')
            row = admit(prior['manifest'], compact, modules, runtime)
            require(row == prior['rows'][SCOUT], 'Raw actual scout admission does not replay')
        else:
            result['allocation'] = runtime.allocation(modules['render_orca_runtime_input'], out)
            full_pin, manifest_pin = prepare(cfg, out, inv, compact, cells, model, np)
            validate_prepared(read(full_pin), cells, compact, cfg, model, np)
            result.update(full_manifest=full_pin, manifest=manifest_pin)
            if args.stage == 'preflight':
                result.update(status='preflight_passed', all_four_states_frozen=True, full_four_available=False)
                row = None
            else:
                result['preflight_result'] = preflight_gate(cfg, read(full_pin), directions, solvent)
                result['selected_tasks'] = [SCOUT]
                write(out / 'START.json', result)
                try:
                    result['native_runner_invocations'] += 1
                    modules['run_orca_task_manifest'].run_manifest(verify(manifest_pin), orca_path=verify(compact['orca']),
                        workers=1, nprocs=12, selected_tasks={SCOUT},
                        expected_runner_sha256=compact['execution_policy']['task_runner']['sha256'])
                except Exception as exc:
                    result['runner_error'] = repr(exc)
                finally:
                    task = read(manifest_pin)['tasks'][0]
                    result['new_DFT_calls_attempted'] = int(Path(task['output_path']).exists())
                row = admit(manifest_pin, compact, modules, runtime)
        if row is not None:
            result['rows'][SCOUT] = row
            result['status'] = 'complete' if row['status'] == 'available' else 'unavailable'
            result['full_four_available'] = False
            if row['status'] == 'available':
                comparison = force_comparison(cells[SCOUT], row, solvent[SCOUT], directions, np, model)
                result['force_comparison'] = write(out / 'FORCE_COMPARISON.json', comparison)
                result['Phe1_CA_force_comparison'] = comparison['Phe1_CA']
                result['force_distributions'] = comparison['distributions']
                result['completed_reference_cells'] = 1
    except Exception as exc:
        result.update(status='execution_failure', failure=repr(exc), traceback=traceback.format_exc())
    finally:
        result['elapsed_seconds'] = time.monotonic() - start
        result['application_allocated_CPU_seconds'] = result['elapsed_seconds'] * int(os.environ.get('SLURM_CPUS_ON_NODE', '0'))
        result['peak_parent_RSS_KiB'] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        result['parent_RSS_is_not_MPI_peak'] = True
        result['physical_state_qualified'] = result['force_consistency_qualified'] = False
        write(out / 'RESULT.json', result)
    print(json.dumps({k: result[k] for k in ('status', 'selected_tasks', 'new_DFT_calls_attempted', 'elapsed_seconds')}))
    return 0 if result['status'] in ('complete', 'preflight_passed') else 1


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--stage', choices=('preflight', 'scout', 'collect'), required=True)
    parser.add_argument('--config', required=True)
    parser.add_argument('--output', required=True)
    raise SystemExit(run(parser.parse_args()))
