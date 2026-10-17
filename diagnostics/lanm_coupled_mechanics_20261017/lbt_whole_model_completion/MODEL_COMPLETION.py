#!/usr/bin/env python3
"""Finite whole256 model completion. All stages require Slurm; no geometry search."""
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
from types import SimpleNamespace

PROTOCOL = 'lbt_whole256_fixed_deformation_model_completion_v1'
IDS = ('before_La', 'after_La', 'before_Dy', 'after_Dy')
NEW_MACE = ('after_La', 'after_Dy')
NEW_ALPB = ('after_Dy',)
INVENTORY_SHA = '48d9f40f1bc7f5f927131d5f8e1b2ad3d1fa12187f7b8c321b7cb5375efe4baa'
REFERENCE_SHA = 'c22c989674e280806abace63c67148b0181caee0b7b20bcd85b8e3d75d6c2a86'
NATIVE_SHA = '0e48b8c6c2f271512f2c445f826036d8a68f7dca0aabc5923697ab3ecb7d9896'
SOLVENT_SHA = '563b5e1a32fd3bf94544a42c244c54270c7d5bb497d60747db52d4415d4ff22e'
EV, HA, BOHR = 23.06054783061903, 627.509474, .529177210903
BIND_KEYS = ('protocol', 'implementation', 'specification', 'plan', 'inventory',
             'reference_implementation', 'reference_preflight', 'MACE_config',
             'native_MACE_helper', 'solvent_parent', 'solvent_helper')


def require(ok, message):
    if not ok:
        raise ValueError(message)


def pin(path):
    p = Path(path).resolve()
    with p.open('rb') as stream:
        digest = hashlib.file_digest(stream, 'sha256').hexdigest()
    return {'path': str(p), 'sha256': digest}


def verify(record):
    require(pin(record['path'])['sha256'] == record['sha256'], 'Changed pin: ' + record['path'])
    return Path(record['path'])


def read(record):
    return json.loads(verify(record).read_text())


def write(path, value):
    with Path(path).open('x') as stream:
        json.dump(value, stream, indent=2, allow_nan=False); stream.write('\n')
    return pin(path)


def module(name, record):
    spec = importlib.util.spec_from_file_location(name, verify(record))
    obj = importlib.util.module_from_spec(spec); spec.loader.exec_module(obj)
    return obj


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False).encode()).hexdigest()


def binding(cfg):
    # A byte-identical staged executable is permitted; all scientific inputs retain exact pins.
    return {k: (cfg[k]['sha256'] if k == 'implementation' else cfg[k]) for k in BIND_KEYS}


def input_task(cell, old_task, xyz):
    metal = cell['row']['metal']
    return dict(old_task, task_id=cell['row']['cell_id'], origin_task_id=old_task['task_id'],
                physical_source=cell['physical_parent'], xyz=xyz, metal=metal,
                charge=-2, multiplicity=1 if metal == 'La' else 6)


def admit_mace(raw, cell, h, need_force):
    h.state(raw['input_state_check'], 256, -2, 1 if cell['row']['metal'] == 'La' else 6, 255)
    require(math.isfinite(raw['energy_eV']), 'Nonfinite native MACE energy')
    if need_force:
        require(raw['force_convention'] == 'negative derivative of the same native vacuum energy',
                'MACE force convention differs')
        h.vectors(raw['forces_eV_A'], 256)
    else:
        require(raw['forces_eV_A'] is None and raw['force_convention'] is None,
                'Archived before-Dy is an energy-only request')


def hybrid_context(cells):
    admitted = {}
    for tid in ('before_La', 'after_La', 'before_Dy'):
        r = read(cells[tid]['row']['hybrid_point']); p = r['point']; receipts = r['component_receipts']
        ml, emb, iso, mm, sol = (read(receipts[k]) for k in ('MACE', 'GFN_embedded', 'GFN_isolated', 'MM', 'ALPB'))
        require(emb['status'] == iso['status'] == 'accepted' and mm['status'] == 'complete',
                'Archived hybrid components are unavailable')
        native = {'MACE_eV': ml['energy_eV'], 'GFN_embedded_Eh': emb['energy_hartree'],
                  'GFN_isolated_Eh': iso['energy_hartree'], 'MM_groups_kcal_mol': mm['component_energies_kcal_mol'],
                  'ALPB_Eh': sol['total_hartree']}
        require(set(native['MM_groups_kcal_mol']) == {'retained_bonded', 'MM_LJ', 'MM_Coulomb', 'QM_MM_LJ'}
                and mm['energy_kcal_mol'] == math.fsum(native['MM_groups_kcal_mol'].values()),
                'Archived hybrid four-group MM closure differs')
        components = {'MACE_vacuum': native['MACE_eV']*EV,
                      'GFN2_field_difference': (native['GFN_embedded_Eh']-native['GFN_isolated_Eh'])*HA,
                      'MM_ledger': mm['energy_kcal_mol'], 'fixed_FF_charge_ALPB': native['ALPB_Eh']*HA}
        require(all(math.isfinite(v) for v in components.values()) and native == p['native_energies']
                and components == p['component_energies_kcal_mol']
                and math.fsum(components.values()) == p['total_energy_kcal_mol'],
                'Archived hybrid native/component/total energy closure differs')
        admitted[tid] = {'point': cells[tid]['row']['hybrid_point'], 'native_energies': native,
                         'component_energies_kcal_mol': components, 'receipts': receipts}
    na, nb = (admitted[k]['native_energies'] for k in ('after_La', 'before_La'))
    work = {'MACE_vacuum': (na['MACE_eV']-nb['MACE_eV'])*EV,
            'GFN2_field_difference': math.fsum((na['GFN_embedded_Eh']-nb['GFN_embedded_Eh'],
                                              -na['GFN_isolated_Eh']+nb['GFN_isolated_Eh']))*HA,
            'MM_ledger': math.fsum(na['MM_groups_kcal_mol'][k]-nb['MM_groups_kcal_mol'][k]
                                  for k in na['MM_groups_kcal_mol']),
            'fixed_FF_charge_ALPB': (na['ALPB_Eh']-nb['ALPB_Eh'])*HA}
    work['composite_kcal_mol'] = math.fsum(work.values())
    return {'admitted': admitted, 'La_work': work, 'Dy_work': None}


def prepare(cfg):
    require(os.environ.get('SLURM_JOB_ID', '').isdigit(), 'Slurm required for every stage')
    require(cfg['protocol'] == PROTOCOL and pin(__file__)['sha256'] == cfg['implementation']['sha256'],
            'Model completion implementation/protocol differs')
    for key in BIND_KEYS[1:]:
        verify(cfg[key])
    require(cfg['inventory']['sha256'] == INVENTORY_SHA and cfg['reference_implementation']['sha256'] == REFERENCE_SHA
            and cfg['native_MACE_helper']['sha256'] == NATIVE_SHA and cfg['solvent_helper']['sha256'] == SOLVENT_SHA,
            'Frozen inventory/helpers differ')
    inv = read(cfg['inventory']); ref = module('whole_completion_reference', cfg['reference_implementation'])
    h = module('whole_completion_solvent', cfg['solvent_helper'])
    native = module('whole_completion_native_mace', cfg['native_MACE_helper'])
    common = module('whole_completion_common', inv['native_reference']['native_helpers']['affordable_common.py'])
    import numpy as np
    from ase.data import atomic_numbers
    model = SimpleNamespace(NUMBERS=atomic_numbers, xyz=common.xyz,
                            HA_TO_KCAL=common.HA_TO_KCAL, BOHR_TO_A=common.BOHR_TO_A)
    require(model.HA_TO_KCAL == HA and model.BOHR_TO_A == BOHR, 'Native unit constants differ')
    require(inv['plan'] == cfg['plan'] and inv['declared_cells'] == 4 and inv['declared_geometries'] == 2
            and [c['cell_id'] for c in inv['cells']] == list(IDS), 'Four fixed cells differ')
    cells = ref.physical_inputs(inv, np, model)
    hybrid = hybrid_context(cells)
    rp = read(cfg['reference_preflight']); rpcfg = read(rp['config'])
    require(rp['status'] == 'preflight_passed' and rp['stage'] == 'preflight' and rp['slurm_job_id'] == '13049'
            and rp['implementation']['sha256'] == REFERENCE_SHA and rp['inventory'] == cfg['inventory']
            and rp['plan'] == cfg['plan'] and rp['all_four_states_frozen']
            and rp['states'] == {k: c['state'] for k, c in cells.items()}
            and all(rp[k] == 0 for k in ('new_DFT_calls_attempted', 'new_ALPB_calls', 'new_MACE_calls',
                                         'new_GFN_calls', 'new_MM_calls', 'optimizer_steps')),
            'Actual13049 four-state preparation differs')
    full = read(rp['full_manifest'])
    require([t['task_id'] for t in full['tasks']] == list(IDS), 'Reference full denominator differs')
    for task in full['tasks']:
        c = cells[task['task_id']]; atoms = common.xyz(verify(task['xyz']))
        require(task['physical_parent'] == c['physical_parent'] and task['positions'] == c['row']['positions']
                and task['physical_ids'] == c['ids'] and task['electronic_state'] == c['state']
                and [a[0] for a in atoms] == c['symbols'] and np.array_equal([a[1:] for a in atoms], c['x']),
                'Exact13049 physical geometry/state differs')
    directions = read(rp['directions'])
    require(directions['physical_ids'] == cells['before_La']['ids']
            and directions['selection_source'] == inv['force_audit'] and directions['frozen_before_new_reference_calls']
            and set(directions['directions']) == {'full_deformation_chord', 'initial_hybrid_Phe1_CA_force'},
            'Frozen13049 directions differ')
    for vector in directions['directions'].values():
        a = np.asarray(h.vectors(vector, 256))
        require(abs(float(np.linalg.norm(a)) - 1.) <= 1e-12, 'Frozen full-array direction is not normalized')
    # Consume the exact archived vectors; do not regenerate a direction in a new numerical environment.
    mcfg = read(cfg['MACE_config'])
    require(pin(sys.executable)['sha256'] == mcfg['python']['sha256'], 'MACE Python executable differs')
    for record in mcfg['helper_pins'].values():
        verify(record)
    verify(mcfg['model']['checkpoint'])
    versions = {k: importlib.metadata.version(k) for k in ('torch', 'mace-torch', 'ase')}
    require(versions == mcfg['expected_versions'], 'Pinned native model software differs')
    sp = read(cfg['solvent_parent'])
    require(sp['status'] == 'complete' and sp['slurm_job_id'] == '12715'
            and sp['implementation']['sha256'] == SOLVENT_SHA and sp['build']['returncode'] == 0,
            'Actual native solvent build differs')
    binary = sp['build']['binary']; verify(binary)
    require(binary == rpcfg['ALPB_binary'], 'Reference/native solvent binary differs')
    h.pins_in(sp['build'])
    env = dict(os.environ)
    for k, v in sp['native_environment'].items():
        if v is None: env.pop(k, None)
        else: env[k] = v
    solvent = {k: ref.admit_alpb(c, {'ALPB_binary': binary}, np, model) for k, c in cells.items()}
    require(solvent == rp['solvent_admission'], 'Exact three archived ALPB admissions differ from13049')
    originals, raw_mace, old_tasks = {}, {}, {}
    for tid in ('before_La', 'before_Dy'):
        c = cells[tid]; saved = c['row']['existing_whole_MACE_cell']; p = read(saved['point'])
        require(saved['candidate_id'] == 'cartesian_01' and saved['status'] == p['status'] == 'complete'
                and p['source_id'] == '7CCO' and p['metal'] == c['row']['metal']
                and read(p['physical_source']) == c['atoms'] and p['ALPB'] == c['row']['ALPB']
                and h.xyz(p['xyz']) == (c['symbols'], c['x'].tolist())
                and np.array_equal(np.load(verify(p['positions']), allow_pickle=False), c['x']),
                'Archived whole before point physical lineage differs')
        recipe = p['recipe']; task = read(p['task']); raw = read(p['MACE'])
        require(recipe['model'] == mcfg['model'] and recipe['MACE_versions'] == versions
                and recipe['physical_ids'] == c['ids'] and recipe['symbols'] == c['symbols']
                and recipe['fixed_FF_charges_e'] == c['charges']['charges_e']
                and recipe['physical_charge'] == -2 and recipe['physical_multiplicity'] == c['state']['physical_multiplicity']
                and recipe['target_index'] == 255 and recipe['native_ALPB_binary'] == binary
                and recipe['native_ALPB_build']['linked_libraries'] == sp['build']['runtime_libraries']
                and recipe['protocol_id'] == 'lbt3_matched_fixedFF_response_v1'
                and recipe['field_and_MM_ledger'] == 'absent by whole-peptide partition'
                and recipe['solvent_definition'] == 'native ALPB water fixed ff19SB plus formal3 target; no SCC',
                'Archived whole model recipe differs')
        h.pins_in(recipe['native_ALPB_build'])
        require(p['cache_key'] == digest({'recipe': recipe, 'positions_A': c['x'].tolist()})
                and task['physical_source'] == p['physical_source'] and task['xyz'] == p['xyz']
                and task['charge'] == -2 and task['multiplicity'] == c['state']['physical_multiplicity']
                and task['metal'] == c['row']['metal'] and task['metal_index'] == 255
                and task['atom_count'] == 256 and task['all_electron_count'] == c['state']['all_electron_count']
                and task['representation'] == 'whole' and task['source_cap_map'] is None,
                'Saved whole task/cache state differs')
        need_force = tid == 'before_La'
        require(p['compute_force'] is need_force, 'Archived force-request availability changed')
        admit_mace(raw, c, h, need_force)
        s = solvent[tid]
        native_e = dict(MACE_eV=raw['energy_eV'], field_Eh=0., MM_kcal_mol=0.,
                        ALPB_Eh=s['energy_hartree'], ALPB_parts_Eh=s['components_hartree'])
        require(p['native_energies'] == saved['native_energies'] == native_e
                and p['components'] == {'MACE_eV': raw['energy_eV'], 'fixed_ALPB_hartree': s['energy_hartree']},
                'Archived whole native energy replay differs')
        if need_force:
            g = -np.asarray(raw['forces_eV_A']) * EV + np.asarray(s['gradient_hartree_per_bohr']) * (HA / BOHR)
            require(np.array_equal(g, np.load(verify(p['gradient_kcal_mol_A']), allow_pickle=False)),
                    'Archived complete whole gradient replay differs')
        originals[tid], raw_mace[tid], old_tasks[c['row']['metal']] = p, raw, task
    signatures = {tid: {'physical_parent': c['physical_parent'], 'positions': c['row']['positions'],
                       'physical_ids': c['ids'], 'symbols': c['symbols'], 'positions_A': c['x'].tolist(),
                       'charge': -2, 'multiplicity': c['state']['physical_multiplicity'],
                       'fixed_FF_charges_e': c['charges']['charges_e']} for tid, c in cells.items()}
    return dict(cfg=cfg, inv=inv, ref=ref, h=h, native=native, model=model, np=np, cells=cells,
                reference=rp, directions=directions, hybrid_context=hybrid, mcfg=mcfg, solvent_parent=sp, binary=binary, env=env,
                solvent=solvent, originals=originals, raw_mace=raw_mace, old_tasks=old_tasks,
                signatures=signatures, versions=versions)


def package(w, tid, raw_pin, solvent, output, task_pin=None):
    c, np, h = w['cells'][tid], w['np'], w['h']
    raw = read(raw_pin); has_force = tid != 'before_Dy'
    admit_mace(raw, c, h, has_force)
    require(solvent['status'] == 'available', 'Native solvent missing')
    native = {'MACE_eV': raw['energy_eV'], 'fixed_ALPB_hartree': solvent['energy_hartree']}
    row = dict(cell_id=tid, status='available', source_id='7CCO', metal=c['row']['metal'],
               physical_source=c['physical_parent'], positions=c['row']['positions'], physical_ids=c['ids'],
               charge=-2, physical_multiplicity=c['state']['physical_multiplicity'], target_index=255,
               MACE=raw_pin, ALPB=solvent['receipt'], native_energies=native,
               force_status='available' if has_force else 'unavailable_energy_only_archived_request',
               hybrid_point=c['row']['hybrid_point'], task=task_pin,
               force_arrays=None, force_projections_kcal_mol_A=None)
    if has_force:
        mf = np.asarray(raw['forces_eV_A']) * EV
        sf = -np.asarray(solvent['gradient_hartree_per_bohr']) * (HA / BOHR)
        arrays = {'MACE_vacuum': mf, 'fixed_FF_ALPB': sf, 'composite': mf + sf}
        require(all(a.shape == (256, 3) and np.isfinite(a).all() for a in arrays.values()), 'Invalid full physical forces')
        row['force_arrays'] = write(output / (tid + '_FORCES.json'),
                                   {'physical_ids': c['ids'], 'units': 'kcal/mol/Angstrom',
                                    'quantity': 'force = negative derivative',
                                    'components': {k: a.tolist() for k, a in arrays.items()}})
        row['force_projections_kcal_mol_A'] = {
            name: {k: float(np.sum(a * np.asarray(vector))) for k, a in arrays.items()}
            for name, vector in w['directions']['directions'].items()}
        row['force_residuals_kcal_mol_A'] = {
            k: {'max_atom_vector_norm': float(np.linalg.norm(a, axis=1).max()),
                'atom_vector_RMS': float(np.sqrt(np.mean(np.sum(a*a, axis=1))))}
            for k, a in arrays.items()}
    return row


def compare(w, rows):
    by = {r['cell_id']: r for r in rows}
    if any(by[k]['status'] != 'available' for k in IDS):
        return {'status': 'unavailable', 'reason': 'All four native energy cells are required; no survivor differential.'}
    work = {}
    for metal in ('La', 'Dy'):
        a, b = by['after_' + metal]['native_energies'], by['before_' + metal]['native_energies']
        de, dh = a['MACE_eV'] - b['MACE_eV'], a['fixed_ALPB_hartree'] - b['fixed_ALPB_hartree']
        work[metal] = {'MACE_eV': de, 'fixed_ALPB_hartree': dh,
                       'MACE_kcal_mol': de*EV, 'fixed_ALPB_kcal_mol': dh*HA,
                       'composite_kcal_mol': math.fsum((de*EV, dh*HA))}
    diff = {k: work['Dy'][k] - work['La'][k] for k in work['La']}
    hybrid = w['hybrid_context']['La_work']
    return dict(status='complete', definition='W_M = E_M(after) - E_M(before); delta_W = W_Dy - W_La',
                per_metal_work=work, Dy_minus_La_deformation_work=diff,
                hybrid_La_work=hybrid, hybrid_Dy_work=None,
                hybrid_Dy_unavailable_reason='No existing hybrid after-Dy evaluation; closed branch remains closed.',
                material_reporting_resolution_kcal_mol=1., accuracy_or_preference_claim=False,
                energy_coverage=4, full_force_coverage=3)


def run(w, out, result):
    cfg = w['cfg']; prior = read(cfg['preflight_result'])
    require(prior['status'] == 'preflight_passed' and prior['stage'] == 'preflight'
            and prior['scientific_binding'] == binding(cfg) and prior['signatures'] == w['signatures']
            and prior['directions'] == w['reference']['directions']
            and prior['solvent_admission'] == w['solvent'] and prior['hybrid_context'] == w['hybrid_context']
            and prior['model_loaded'] is False
            and prior['MACE_calls_attempted'] == prior['ALPB_calls_attempted'] == 0,
            'Actual unchanged zero-call preflight required')
    native = None
    for tid in NEW_MACE:
        row = next(r for r in result['rows'] if r['cell_id'] == tid); c = w['cells'][tid]
        directory = out / tid; directory.mkdir()
        xyz = directory / 'whole.xyz'
        xyz.write_text('256\nExact frozen whole-LBT deformation endpoint\n' + ''.join(
            '%s %.17g %.17g %.17g\n' % (s, *v) for s, v in zip(c['symbols'], c['x'])))
        require(w['h'].xyz(pin(xyz)) == (c['symbols'], c['x'].tolist()), 'Exact native XYZ serialization differs')
        task = input_task(c, w['old_tasks'][c['row']['metal']], pin(xyz))
        tp = write(directory / 'TASK.json', task)
        request = {'cell_id': tid, 'status': 'unavailable', 'task': tp, 'MACE_calls_attempted': 0,
                   'ALPB_calls_attempted': 0, 'MACE': None, 'ALPB': c['row']['ALPB']}
        result['requests'].append(request); tick = time.monotonic()
        try:
            if native is None:
                t = time.monotonic(); native = w['native'].NativeMACE(w['mcfg'])
                result['model_load_seconds'] = time.monotonic()-t
                result['model_loaded'] = True; result['GPU_device'] = native.device
            require(result['MACE_calls_attempted'] < 2 and tid in NEW_MACE, 'MACE call ceiling/allowlist')
            before = native.calls; t = time.monotonic()
            try:
                raw = native.evaluate(task, c['symbols'], c['x'], compute_force=True)
            finally:
                count = native.calls-before; request['MACE_calls_attempted'] += count
                result['MACE_calls_attempted'] += count
                result['MACE_call_seconds'] += time.monotonic()-t
            require(count == 1, 'One native MACE request required')
            request['MACE'] = write(directory / 'MACE.json', raw)
            admit_mace(raw, c, w['h'], True)
            solvent = w['solvent'][tid]
            if tid in NEW_ALPB:
                require(result['ALPB_calls_attempted'] == 0, 'One ALPB call ceiling')
                for lib in w['solvent_parent']['build']['runtime_libraries']: verify(lib)
                result['ALPB_calls_attempted'] += 1; request['ALPB_calls_attempted'] += 1
                t = time.monotonic()
                try:
                    w['h'].solvent(w['binary'], {'symbols': c['symbols'], 'x': c['x'].tolist(),
                        'q': c['charges']['charges_e'], 'task': task}, directory/'ALPB', w['env'])
                finally:
                    result['ALPB_call_seconds'] += time.monotonic()-t
                request['ALPB'] = pin(directory/'ALPB'/'RESULT.json')
                candidate = dict(c, row=dict(c['row'], ALPB=request['ALPB']))
                solvent = w['ref'].admit_alpb(candidate, {'ALPB_binary': w['binary']}, w['np'], w['model'])
            completed = package(w, tid, request['MACE'], solvent, out, tp)
            row.clear(); row.update(completed); request['status'] = 'complete'
        except Exception as exc:
            row.update(status='unavailable', reason=str(exc)); request.update(failure=str(exc), traceback=traceback.format_exc())
            raise
        finally:
            request['elapsed_seconds'] = time.monotonic()-tick
            record = write(directory/'REQUEST.json', request); request['record'] = record
    require(result['MACE_calls_attempted'] == 2 and result['ALPB_calls_attempted'] == 1,
            'Actual finite completion denominator differs')
    result['peak_GPU_allocated_bytes'] = native.torch.cuda.max_memory_allocated()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--stage', choices=('preflight', 'run', 'collect'), required=True)
    parser.add_argument('--config', required=True); parser.add_argument('--output', required=True)
    args = parser.parse_args(); out = Path(args.output).resolve(); out.mkdir(parents=True, exist_ok=False)
    cp = pin(args.config); cfg = read(cp); tick = time.monotonic()
    result = dict(protocol=PROTOCOL, stage=args.stage, status='execution_failure', implementation=pin(__file__),
                  config=cp, slurm_job_id=os.environ.get('SLURM_JOB_ID'), host=socket.gethostname(),
                  declared_cells=4, new_MACE_allowlist=list(NEW_MACE), new_ALPB_allowlist=list(NEW_ALPB),
                  MACE_calls_attempted=0, ALPB_calls_attempted=0, new_DFT_calls=0, new_GFN_calls=0,
                  new_MM_calls=0, optimizer_steps=0, new_geometry_count=0, model_loaded=False,
                  model_load_seconds=0., MACE_call_seconds=0., ALPB_call_seconds=0., requests=[],
                  rows=[{'cell_id': k, 'status': 'unavailable', 'reason': 'Not yet admitted.'} for k in IDS],
                  comparison={'status': 'unavailable'}, physical_state_qualified=False, accuracy_claim=False)
    try:
        w = prepare(cfg)
        result.update(scientific_binding=binding(cfg), signatures=w['signatures'],
                      directions=w['reference']['directions'], solvent_admission=w['solvent'],
                      hybrid_context=w['hybrid_context'],
                      physical_ids=w['cells']['before_La']['ids'], versions=w['versions'],
                      native_environment=w['solvent_parent']['native_environment'],
                      model_recipe={k: p['recipe'] for k, p in w['originals'].items()},
                      caveats={'before_Dy_force': 'Archived whole-MACE request is energy-only; no added force call.',
                               'hybrid_after_Dy': 'Unavailable; no new hybrid evaluation authorized.',
                               'scope': 'Same fixed 7CCO La deformation, nonstationary80-step endpoint, q-2 free termini; physical multiplicities1/6.',
                               'spin': 'Physical MACE sextet Dy input is distinct from effective frozen-f DFT singlet.',
                               'interpretation': 'Model sensitivity only, no affinity, accuracy or stationary-state claim.'})
        for tid in ('before_La', 'before_Dy'):
            slot = next(r for r in result['rows'] if r['cell_id'] == tid)
            slot.clear(); slot.update(package(w, tid, w['originals'][tid]['MACE'], w['solvent'][tid], out,
                                              w['originals'][tid]['task']))
        if args.stage == 'preflight':
            result['status'] = 'preflight_passed'
        elif args.stage == 'run':
            run(w, out, result); result['comparison'] = compare(w, result['rows']); result['status'] = 'complete'
        else:
            previous = read(cfg['run_result'])
            require(previous['protocol'] == PROTOCOL and previous['stage'] == 'run'
                    and previous['scientific_binding'] == binding(cfg) and previous['signatures'] == w['signatures'],
                    'Saved run identity differs')
            result['parent_run'] = cfg['run_result']; result['parent_costs'] = {
                k: previous[k] for k in ('MACE_calls_attempted', 'ALPB_calls_attempted', 'elapsed_seconds',
                                        'MACE_call_seconds', 'ALPB_call_seconds', 'model_load_seconds')}
            for tid in NEW_MACE:
                old = next(r for r in previous['rows'] if r['cell_id'] == tid)
                if old['status'] != 'available': continue
                c = w['cells'][tid]; task = read(old['task']); sy, x = w['h'].xyz(task['xyz'])
                require(task == input_task(c, w['old_tasks'][c['row']['metal']], task['xyz'])
                        and sy == c['symbols'] and x == c['x'].tolist(), 'Saved new task differs')
                request = next(r for r in previous['requests'] if r['cell_id'] == tid)
                require(read(request['record']) == {k: v for k, v in request.items() if k != 'record'}
                        and request['status'] == 'complete' and request['task'] == old['task']
                        and request['MACE'] == old['MACE'] and request['ALPB'] == old['ALPB']
                        and request['MACE_calls_attempted'] == 1
                        and request['ALPB_calls_attempted'] == (1 if tid == 'after_Dy' else 0), 'Saved native request differs')
                solvent = w['ref'].admit_alpb(dict(c, row=dict(c['row'], ALPB=old['ALPB'])),
                                             {'ALPB_binary': w['binary']}, w['np'], w['model'])
                slot = next(r for r in result['rows'] if r['cell_id'] == tid)
                slot.clear(); slot.update(package(w, tid, old['MACE'], solvent, out, old['task']))
                require(read(slot['force_arrays']) == read(old['force_arrays'])
                        and slot['native_energies'] == old['native_energies'], 'Saved new energy/full-force replay differs')
            if previous['status'] == 'complete':
                require(previous['MACE_calls_attempted'] == 2 and previous['ALPB_calls_attempted'] == 1,
                        'Completed parent call denominator differs')
                result['comparison'] = compare(w, result['rows']); result['status'] = 'complete'
                require(result['comparison'] == previous['comparison'], 'Saved comparison replay differs')
            else:
                result['status'] = 'unavailable_parent_failure'
                result['comparison'] = {'status': 'unavailable', 'reason': 'Failed parent retained; no survivor differential.'}
    except Exception as exc:
        result.update(status='execution_failure', failure=str(exc), traceback=traceback.format_exc())
    result['elapsed_seconds'] = time.monotonic()-tick
    cpus = int(os.environ.get('SLURM_CPUS_PER_TASK', '0'))
    result['allocation'] = {'cpus_per_task': cpus, 'job_gpus': os.environ.get('SLURM_JOB_GPUS'),
                            'cuda_visible_devices': os.environ.get('CUDA_VISIBLE_DEVICES')}
    result['application_allocated_CPU_seconds'] = cpus*result['elapsed_seconds']
    result['application_allocated_GPU_seconds'] = result['elapsed_seconds'] if args.stage == 'run' and result['model_loaded'] else None
    result['peak_parent_RSS_KiB'] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    write(out/'RESULT.json', result)
    return 0 if result['status'] in ('complete', 'preflight_passed') else 1


if __name__ == '__main__':
    raise SystemExit(main())
