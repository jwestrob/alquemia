#!/usr/bin/env python3
"""Four pinned saved-force comparisons; Slurm only, no molecular executors."""
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

PROTOCOL = 'lbt_saved_partition_force_reference_v1'
STATES = (('baseline', 'La'), ('baseline', 'Dy'), ('expanded', 'La'), ('expanded', 'Dy'))
COMPONENTS = ('MACE_vacuum', 'GFN2_field_difference', 'MM_ledger', 'fixed_FF_charge_ALPB')
REFERENCE_SHA = '9229f7bebb533dd44a024bb80eca10eded2b3a64cc3972b245d65d82b5ecda97'
POOL_SHA = 'd98aafb629185e785a3bd291802e46c374cfeb35b46b52486147c51303f1d592'
POINT_SHAS = ('3fac44a78143ece9cba653e0eaec5aac59852166c32515c400210425d4f2da8f',
              '6995daa6daaed7df8951d0bb6d680ac6bb4480d8831ddf16ba77f91c22e0aad0',
              '2bd337b7ade180e2fb046b32d859eff1acbe1d62acfb0350e5c0138da4b0e338',
              '6da98a9c37c735a8874ddd2e7c5055ab9d9c9f3f69431bbd4d03fef68f5f77c5')
DIRECTION_SHA = '6c306212746400b7bf25635c0548d45f53d17eb5a39cc95782a1fa72e97d52ed'
EV = 23.06054783061903
HA = 627.509474
BOHR = 0.529177210903


def require(ok, message):
    if not ok:
        raise ValueError(message)


def pin(path):
    p = Path(path).resolve()
    with p.open('rb') as f:
        return {'path': str(p), 'sha256': hashlib.file_digest(f, 'sha256').hexdigest()}


def verify(record):
    require(isinstance(record, dict) and set(record) == {'path', 'sha256'}, 'Malformed evidence pin')
    require(pin(record['path'])['sha256'] == record['sha256'], 'Changed input: ' + record['path'])
    return Path(record['path'])


def read(record):
    return json.loads(verify(record).read_text())


def verify_tree(value):
    if isinstance(value, dict):
        if set(value) == {'path', 'sha256'}:
            verify(value)
        else:
            for v in value.values():
                verify_tree(v)
    elif isinstance(value, list):
        for v in value:
            verify_tree(v)


def write(path, value):
    with Path(path).open('x') as f:
        json.dump(value, f, indent=2, allow_nan=False)
        f.write('\n')
    return pin(path)


def array(value, shape=(256, 3)):
    a = np.asarray(value, dtype=float)
    require(a.shape == shape and np.isfinite(a).all(), 'Invalid finite array: ' + str(shape))
    return a


def metrics(value):
    a = array(value)
    n = np.linalg.norm(a, axis=1)
    return {'maximum_atom_kcal_mol_A': float(n.max()),
            'RMS_atom_vector_kcal_mol_A': float(np.sqrt(np.mean(n*n))),
            'mean_atom_norm_kcal_mol_A': float(n.mean()),
            'vector_sum_kcal_mol_A': np.sum(a, axis=0).tolist(),
            'maximum_atom_index': int(np.argmax(n))}


def references(ready, result):
    require(ready['protocol'] == PROTOCOL and ready['implementation']['sha256'] == pin(__file__)['sha256'],
            'Source/protocol differs')
    verify_tree(ready)
    require(ready['reference']['sha256'] == REFERENCE_SHA and ready['pool']['sha256'] == POOL_SHA
            and [x['point']['sha256'] for x in ready['cells']] == list(POINT_SHAS)
            and [(x['region'], x['metal']) for x in ready['cells']] == list(STATES), 'Frozen selection differs')
    ref, pool = read(ready['reference']), read(ready['pool'])
    require(ref['status'] == 'complete' and ref['stage'] == 'collect' and ref['slurm_job_id'] == '13424'
            and ref['full_four_available'] and ref['new_DFT_calls_attempted'] == 0, 'Admitted full reference required')
    require(pool['status'] == 'complete' and pool['stage'] == 'run' and pool['slurm_job_id'] == '13013'
            and pool['declared_hybrid_cells'] == 104 and len(pool['cells']) == 104, 'Complete recovered pool required')
    verify_tree(ref['scientific_binding'])
    require(ref['directions']['sha256'] == DIRECTION_SHA, 'Frozen directions differ')
    directions = read(ref['directions'])
    require(directions['frozen_before_new_reference_calls'] and set(directions['directions']) ==
            {'full_deformation_chord', 'initial_hybrid_Phe1_CA_force'}, 'Frozen direction identity differs')
    dv = {k: array(v) for k, v in directions['directions'].items()}
    spec = importlib.util.spec_from_file_location('saved_partition_scatter', verify(ready['mapping_helper']))
    mapper = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mapper)
    mapper.np = np  # Use only existing scatter; never initialize its molecular dependencies.
    result.update(reference=ready['reference'], pool=ready['pool'], directions=ref['directions'],
                  reference_caveats=ref['caveats'], force_convention='force=-gradient; error=hybrid-reference',
                  units='kcal/mol/Angstrom', closure_policy='Exact original component ordering and array equality; no tolerance relaxation')
    refs = {}
    for metal in ('La', 'Dy'):
        tid = 'before_' + metal
        native = ref['rows'][tid]
        signature = next(s for s in ref['original_signatures'] if s['task_id'] == tid)
        require(native['status'] == 'available' and native['task_id'] == tid
                and native['electronic_state'] == signature['electronic_state']
                and native['task']['physical_parent'] == signature['physical_parent']
                and native['task']['positions'] == signature['positions'], 'Reference native/state linkage differs')
        state = signature['electronic_state']
        require(state['metal'] == metal and state['charge'] == -2 and state['effective_multiplicity'] == 1
                and state['physical_multiplicity'] == (1 if metal == 'La' else 6)
                and state['whole_physical_atom_count'] == 256 and state['target_index'] == 255,
                'Whole reference state differs')
        for k in ('receipt', 'output', 'engrad'):
            verify(native[k])
        verify_tree(native['task'])
        sol = ref['solvent_admission'][tid]
        require(sol['status'] == 'available', 'Reference solvent unavailable')
        solvent = read(sol['receipt'])
        verify_tree(solvent['artifacts'])
        raw_s = array(np.loadtxt(verify(solvent['artifacts']['gradient.dat'])))
        require(np.array_equal(raw_s, sol['gradient_hartree_per_bohr'])
                and np.array_equal(raw_s*HA/BOHR, sol['gradient_kcal_mol_A']), 'Reference solvent replay differs')
        g = array(native['gradient_kcal_mol_per_A'])
        f = -(g + array(sol['gradient_kcal_mol_A']))
        report = read(ref['force_reports'][tid])
        ids = signature['physical_ids']
        require(len(ids) == len(set(ids)) == 256 and ids[3] == 'A/1//CA' and ids[-1] == 'metal'
                and ids == directions['physical_ids']
                and [a['physical_id'] for a in report['force_table']] == ids,
                'Reference physical order differs')
        require(np.array_equal(f, [a['force_kcal_mol_A']['reference_PBE0_plus_fixed_ALPB'] for a in report['force_table']]),
                'Reference force report does not replay')
        for name, d in dv.items():
            require(float(np.sum(f*d)) == report['force_projections'][name]['reference_PBE0_plus_fixed_ALPB'],
                    'Reference frozen projection differs')
        x = read(signature['positions'])['positions_A']
        require([a['xyz_A'] for a in read(signature['physical_parent'])] == x, 'Reference physical coordinates differ')
        refs[metal] = dict(force=f, signature=signature, positions=x, ids=ids, solvent=sol,
                           raw_solvent_gradient=raw_s, report=report, native=native)
    require(refs['La']['positions'] == refs['Dy']['positions'] and refs['La']['ids'] == refs['Dy']['ids'],
            'Reference metals are not paired on the same geometry')
    return pool, refs, dv, mapper


def one_cell(cell, pool, reference, directions, mapper, out):
    region, metal = cell['region'], cell['metal']
    matches = [r for r in pool['cells'] if (r['source_id'], r['candidate_id'], r['region'], r['metal']) ==
               ('7CCO', 'cartesian_01', region, metal)]
    require(len(matches) == 1 and matches[0]['status'] == 'complete' and matches[0]['point'] == cell['point'],
            'Selected pool cell differs')
    row = matches[0]
    raw = read(cell['point']); p = raw['point']; mapping = read(raw['mapping'])
    verify_tree(raw['component_receipts']); verify_tree(raw['source_task'])
    task = read(raw['task']); part = read(task['partition'])
    require(task['source_id'] == '7CCO' and task['region'] == region and task['metal'] == metal,
            'Current task identity differs')
    charges = read(task['global_solvent_charges'])
    require(charges['physical_ids'] == reference['ids'] and charges['total_e'] == -2
            and task['global_solvent_charges']['sha256'] == '2fd973f2798b8c3525c45a195b2580d8e7d08cc32c4f91a0c31b6375b19fa02c',
            'Original full physical fixed charges differ')
    require(raw['status'] == 'complete' and raw['identity_sha256'] == pool['identity_sha256']
            and raw['cache_key'] == row['cache_key'] and (p['source_id'], p['region'], p['metal']) == ('7CCO', region, metal)
            and raw['state'] == p['state'] and p['positions_A'] == read(raw['positions'])['positions_A'] == reference['positions']
            and [a['xyz_A'] for a in mapping['physical_atoms']] == reference['positions']
            and [a['id'] for a in mapping['physical_atoms']] == reference['ids'], 'Physical/model identity differs')
    state = p['state']
    require(state['whole_charge'] == state['QM_charge'] == task['charge'] == -2 and state['metal'] == metal
            and state['MACE_multiplicity'] == task['multiplicity'] == (1 if metal == 'La' else 6)
            and state['GFN_uhf'] == task['gfn2_uhf'] == 0 and state['GFN_effective_multiplicity'] == 1,
            'Hybrid electronic state differs')
    rc = {k: read(v) for k, v in raw['component_receipts'].items() if k != 'MM_execution'}
    ml, mm, emb, iso, sol = (rc[k] for k in ('MACE', 'MM', 'GFN_embedded', 'GFN_isolated', 'ALPB'))
    require(emb['status'] == iso['status'] == 'accepted' and mm['status'] == 'complete', 'Native component unavailable')
    require(ml['task'] == raw['task'] and ml['physical_source'] == task['physical_source'] and ml['xyz'] == task['xyz']
            and ml['input_state_check']['charge'] == -2 and ml['input_state_check']['spin_multiplicity'] == state['MACE_multiplicity']
            and ml['input_state_check']['atoms'] == len(mapping['core'])
            and ml['force_convention'] == 'negative vacuum energy derivative'
            and ml['model']['dtype'] == 'float64' and ml['model']['spin_input_convention'] == 'multiplicity',
            'Raw MACE input/force convention differs')
    for component in (emb, iso):
        verify_tree(component)
    require(mm['positions'] == raw['positions'], 'Raw MM geometry differs')
    require(np.array_equal(sum(np.asarray(mm['component_gradients_kcal_mol_A'][k]) for k in mm['group_order']),
                           mm['gradient_kcal_mol_A']), 'Raw MM group closure differs')
    native = dict(MACE_eV=ml['energy_eV'], GFN_embedded_Eh=emb['energy_hartree'], GFN_isolated_Eh=iso['energy_hartree'],
                  MM_groups_kcal_mol=mm['component_energies_kcal_mol'], ALPB_Eh=sol['total_hartree'])
    require(native == p['native_energies'] == row['native_energies']
            and p['component_energies_kcal_mol'] == row['component_energies_kcal_mol']
            and p['total_energy_kcal_mol'] == row['total_energy_kcal_mol'] == math.fsum(p['component_energies_kcal_mol'].values()),
            'Archived component/native energy closure differs')
    require(raw['component_receipts']['ALPB'] == reference['solvent']['receipt'], 'Common fixed-FF solvent receipt differs')
    verify_tree(sol['artifacts'])
    sg = array(np.loadtxt(verify(sol['artifacts']['gradient.dat'])))
    require(np.array_equal(sg, reference['raw_solvent_gradient']), 'Common native solvent derivative differs')
    z = np.zeros((len(mapping['field']), 3))
    eg = mapper.scatter(mapping, np.asarray(emb['qm_gradient']['hartree_per_bohr']),
                        np.asarray(emb['pointcharge_gradient']['hartree_per_bohr']))*HA/BOHR
    ig = mapper.scatter(mapping, np.asarray(iso['qm_gradient']['hartree_per_bohr']), z)*HA/BOHR
    gradients = {'MACE_vacuum': mapper.scatter(mapping, -np.asarray(ml['forces_eV_A'])*EV, z),
                 'GFN2_field_difference': eg-ig, 'MM_ledger': array(mm['gradient_kcal_mol_A']),
                 'fixed_FF_charge_ALPB': sg*HA/BOHR}
    require(set(p['component_gradients_kcal_mol_A']) == set(COMPONENTS)
            and all(np.array_equal(array(gradients[k]), p['component_gradients_kcal_mol_A'][k]) for k in COMPONENTS)
            and np.array_equal(sum(gradients[k] for k in COMPONENTS), p['gradient_kcal_mol_A']),
            'Exact raw component/full gradient replay differs')
    forces = {k: -gradients[k] for k in COMPONENTS}; force = -array(p['gradient_kcal_mol_A'])
    if region == 'baseline':
        require(np.array_equal(force, [a['force_kcal_mol_A']['hybrid_composite'] for a in reference['report']['force_table']]),
                'Exact original baseline force replay differs')
        for name, d in directions.items():
            require(float(np.sum(force*d)) == reference['report']['force_projections'][name]['hybrid_composite'],
                    'Exact baseline frozen projection differs')
    real = {a['physical_index'] for a in mapping['core'] if 'physical_index' in a}
    require(real == set(part['qm_source_indices']) | {255} and len(real) == (170 if region == 'baseline' else 208)
            and len(mapping['caps']) == 2, 'Partition membership differs')
    roles = {}
    for c in mapping['caps']:
        for key in ('retained', 'omitted'):
            i = c[key+'_index']; require(reference['ids'][i] == c[key+'_source_id'], 'Cap identity differs')
            require((i in real) == (key == 'retained') and i not in roles, 'Cap membership differs')
            roles[i] = key+'_boundary_anchor'
    field = {a['source_index'] for a in mapping['field']}
    require(len(field) == len(mapping['field']) and field == set(range(256))-real-
            {i for i, role in roles.items() if role == 'omitted_boundary_anchor'}, 'Field coverage differs')
    require(roles.get(3) == (None if region == 'baseline' else 'omitted_boundary_anchor'), 'Phe1 CA role differs')
    error = force-reference['force']
    atoms = [dict(index=i, physical_id=pid, electronic_real=i in real, boundary_role=roles.get(i),
                  category='target_metal' if i == 255 else roles.get(i, 'electronic_real_other' if i in real else 'MM_other'),
                  hybrid_force=force[i].tolist(), reference_force=reference['force'][i].tolist(),
                  hybrid_minus_reference=error[i].tolist(), components={k: v[i].tolist() for k, v in forces.items()})
             for i, pid in enumerate(reference['ids'])]
    table = write(out/(region+'_'+metal+'_FORCES.json'), {'units':'kcal/mol/Angstrom', 'quantity':'force=-gradient',
                  'point':cell['point'], 'mapping':raw['mapping'], 'atoms':atoms})
    return dict(status='available', region=region, metal=metal, point=cell['point'], state=state,
                full_table=table, hybrid_metrics=metrics(force), reference_metrics=metrics(reference['force']),
                error_metrics=metrics(error), Phe1_CA=atoms[3], caps=mapping['caps'],
                projections={name: {'hybrid':float(np.sum(force*d)), 'reference':float(np.sum(reference['force']*d)),
                    'hybrid_minus_reference':float(np.sum(error*d)),
                    'Phe1_CA_components':{k:float(np.sum(v[3]*d[3])) for k,v in forces.items()},
                    'Phe1_CA_projection_definition':'Contribution using the original full-direction vector at atom3; no local renormalization'}
                    for name,d in directions.items()},
                checks={'exact_raw_component_replay':True,'exact_full_gradient_closure':True,
                        'exact_baseline_replay':True if region=='baseline' else None,
                        'identical_physical_coordinates':True,'same_raw_ALPB_receipt_and_gradient':True},
                solvent_conversion_order_difference=metrics(gradients['fixed_FF_charge_ALPB']-array(reference['solvent']['gradient_kcal_mol_A'])))


def main():
    a = argparse.ArgumentParser(description=__doc__)
    a.add_argument('--ready', required=True); a.add_argument('--output', required=True)
    args = a.parse_args()
    require(os.environ.get('SLURM_JOB_ID', '').isdigit(), 'Saved-force arithmetic requires Slurm')
    attempts = []
    def tripwire(event, values):
        if event in ('subprocess.Popen', 'os.system', 'os.exec', 'os.posix_spawn', 'os.fork'):
            attempts.append(event); raise RuntimeError('External processes forbidden in saved-only analysis')
    sys.addaudithook(tripwire)
    out = Path(args.output).resolve(); out.mkdir(parents=True, exist_ok=False); start = time.monotonic()
    result = dict(protocol=PROTOCOL, status='execution_failure', implementation=pin(__file__), ready=pin(args.ready),
                  slurm_job_id=os.environ['SLURM_JOB_ID'], declared_cells=4,
                  rows=[dict(region=r, metal=m, status='unavailable_not_admitted') for r,m in STATES],
                  new_molecular_calls={'DFT':0,'MACE':0,'GFN':0,'MM':0,'ALPB':0}, new_geometries=0,
                  optimization_steps=0, external_process_attempts=attempts,
                  expanded_after={'status':'unavailable','reason':'No saved expanded hybrid after geometry; outside this four-cell check.'},
                  accuracy_qualified=False, parameter_change_supported=False, preference_claim=False)
    try:
        ready = read(result['ready'])
        global np
        import numpy as np
        require(np.__version__ == ready['numpy_version'], 'Frozen NumPy replay environment differs')
        pool, refs, directions, mapper = references(ready, result)
        for cell, slot in zip(ready['cells'], result['rows']):
            try:
                slot.update(one_cell(cell, pool, refs[cell['metal']], directions, mapper, out))
            except Exception as exc:
                slot.update(status='unavailable', failure=repr(exc), traceback=traceback.format_exc())
        result['status'] = 'complete' if all(r['status']=='available' for r in result['rows']) else 'unavailable'
        verify(result['ready'])
    except Exception as exc:
        result.update(status='execution_failure', failure=repr(exc), traceback=traceback.format_exc())
    result['elapsed_seconds'] = time.monotonic()-start
    write(out/'RESULT.json', result)
    return 0 if result['status']=='complete' else 1


if __name__ == '__main__':
    raise SystemExit(main())
