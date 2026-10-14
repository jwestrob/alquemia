#!/usr/bin/env python3
"""Metadata freeze or Slurm-only reciprocal metal-position preparation; no energies."""
import argparse
import copy
import hashlib
import json
import os
from pathlib import Path
import resource
import shutil
import socket
import subprocess
import sys
import time
import traceback

PROTOCOL = 'ExaA_reciprocal_PQQ_frame_metal_position_preparation_v1'
INVENTORY_SHA = '757a13523549697fd52cf9e94b3ba2ffdd01c463f3969fd8a5d01c1459d28343'
METALS = ('Ca', 'La')
CANDIDATES = ('origin', 'adaptive_La')
PQQ_NAMES = tuple([f'C{i}' for i in range(1, 15)]+['N1', 'N2']+[f'O{i}' for i in range(1, 9)])
THREADS = ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS', 'BLIS_NUM_THREADS', 'OPENMM_CPU_THREADS')


def pin(path):
    p = Path(path).resolve()
    return {'path': str(p), 'sha256': hashlib.sha256(p.read_bytes()).hexdigest()}


def check(value):
    if pin(value['path']) != value:
        raise ValueError('Changed pinned input: '+value['path'])
    return Path(value['path'])


def read(value):
    return json.loads((check(value) if isinstance(value, dict) else Path(value)).read_text())


def write(path, value):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('x') as f:
        json.dump(value, f, indent=2, allow_nan=False); f.write('\n')


def freeze(inventory, plan, output):
    """Only stdlib path/hash/JSON operations; numerical work belongs to execute."""
    invpin = pin(inventory)
    if invpin['sha256'] != INVENTORY_SHA:
        raise ValueError('Not the frozen two-source inventory')
    inv = read(invpin); sources = []
    for sample in (0, 3):
        cid = f'q9z4j7-pqq-la_model__conditioned_2Ca_v1__seed-1_sample-{sample}'
        row, = [r for r in inv['rows'] if r['case_id'] == cid]
        proposal = read(row['native_proposals']['La'])
        pm = read(proposal['manifest']); batch = read(pm['batch_manifest'])
        task, = [t for t in pm['tasks'] if t['metal'] == 'La']
        sources.append({'sample': sample, 'role': row['role'], 'source_case_id': cid,
                        'record': row, 'proposal_manifest': proposal['manifest'],
                        'worker_manifest': pm['batch_manifest'], 'angular_task': task})
    first = read(sources[0]['worker_manifest'])
    if read(sources[1]['worker_manifest'])['source_modules'] != first['source_modules']:
        raise ValueError('Original physical helpers differ between sources')
    out = Path(output).resolve(); out.mkdir(parents=True, exist_ok=False)
    impl = out/'implementation'; impl.mkdir()
    implementation = {}
    for name in ('PREPARE.py', 'RUN.sbatch'):
        src = Path(__file__).with_name(name); dst = impl/name
        shutil.copyfile(src, dst); implementation[name] = pin(dst)
    m = {'protocol': PROTOCOL, 'inventory': invpin, 'plan': pin(plan), 'sources': sources,
         'implementation': implementation, 'source_snapshot': first['source_snapshot'],
         'source_modules': first['source_modules'], 'python': pin(first['optimizer_python']),
         'original_worker_settings': first['settings'], 'PQQ_names': PQQ_NAMES,
         'alignment': 'unweighted named24heavy proper Kabsch source0_to_source3; exact inverse for reverse',
         'settings': {'kinematic_replay_A': 1e-8, 'proper_rotation_tolerance': 1e-10,
                      'minimum_relative_second_singular_value': 1e-12,
                      'maximum_metal_displacement_A': 0.8},
         'output': str(out/'execution'), 'declared_sources': 2, 'declared_new_cells': 8,
         'reused_original_cells': 8, 'new_energy_calls': 0, 'new_optimizations': 0,
         'scientific_execution': 'unrun', 'modeled_not_observed': True}
    write(out/'manifest.json', m); print(json.dumps(pin(out/'manifest.json')))


def placement():
    job = os.environ.get('SLURM_JOB_ID')
    if not job or not job.isdigit():
        raise RuntimeError('Scientific preparation requires actual Slurm')
    text = subprocess.check_output(['scontrol', 'show', 'job', '-o', job], text=True)
    f = dict(x.split('=', 1) for x in text.split() if '=' in x)
    if f.get('JobState') != 'RUNNING' or f.get('NumNodes') != '1' or f.get('NumCPUs') != '1':
        raise RuntimeError('Expected one allocated CPU on one running node')
    if int(os.environ['SLURM_CPUS_ON_NODE']) != 1 or len(os.sched_getaffinity(0)) < 1:
        raise RuntimeError('Allocation environment/affinity differs')
    if any(os.environ.get(k) != '1' for k in THREADS):
        raise RuntimeError('Expected one numerical thread')
    return {'job_id': job, 'hostname': socket.gethostname(), 'scontrol': text,
            'affinity': sorted(os.sched_getaffinity(0)), 'allocated_CPUs': 1}


def execute(manifest, manifest_sha):
    place = placement()  # Precedes every numerical/scientific import.
    mp = {'path': str(Path(manifest).resolve()), 'sha256': manifest_sha}; m = read(mp)
    if m['protocol'] != PROTOCOL or pin(__file__) != m['implementation']['PREPARE.py']:
        raise ValueError('Wrong preparation snapshot/protocol')
    out = Path(m['output']); out.mkdir(parents=True, exist_ok=False)
    began = time.monotonic(); cpu0 = time.process_time()
    receipt = {'status': 'running', 'manifest': mp, 'placement': place, 'new_energy_calls': 0,
               'new_optimizations': 0}
    rows = [{'case_id': s['source_case_id'], 'source_case_id': s['source_case_id'],
             'sample': s['sample'], 'role': s['role'], 'status': 'unavailable',
             'reason': 'preparation not completed',
             'matrix': {z: {q: {'status': 'unavailable', 'reason': 'preparation not completed'}
                            for q in CANDIDATES} for z in METALS}} for s in m['sources']]
    write(out/'START.json', receipt)
    try:
        check(m['inventory']); check(m['plan']); check(m['python'])
        for p in m['implementation'].values(): check(p)
        for p in m['source_modules'].values(): check(p)
        if m['PQQ_names'] != list(PQQ_NAMES) or m['inventory']['sha256'] != INVENTORY_SHA:
            raise ValueError('Fixed input/alignment names changed')
        sys.path.insert(0, m['source_snapshot'])
        import numpy as np
        import adaptive_metal_proposals as metal_geometry
        import mace_site_kinematics as kinematic_module
        import mace_hybrid as atomic_module
        import affordable_common as common
        for module in (metal_geometry, kinematic_module, atomic_module, common):
            if pin(module.__file__) != m['source_modules'][Path(module.__file__).name]:
                raise ValueError('Wrong physical helper import')
        if metal_geometry.SETTINGS['final_physical_tolerance_A'] != m['original_worker_settings']['final_physical_tolerance_A']:
            raise ValueError('Original physical tolerance changed')
        inv = read(m['inventory']); prepared = {}; names = list(PQQ_NAMES)
        for s in m['sources']:
            r = s['record']; sample = s['sample']
            original, = [v for v in inv['rows'] if v['case_id'] == s['source_case_id']]
            if r != original or r['role'] != ('control' if sample == 0 else 'target'):
                raise ValueError('Source identity differs')
            for key in ('prepared_admission', 'source_preparation', 'state_signature', 'full_source', 'spectator_admission'):
                check(r[key])
            for p in r['maps'].values(): check(p)
            for p in r['native_proposals'].values(): check(p)
            pm = read(s['proposal_manifest']); batch = read(s['worker_manifest'])
            if pm['batch_manifest'] != s['worker_manifest'] or batch['source_modules'] != m['source_modules']:
                raise ValueError('Physical helper lineage differs')
            if s['angular_task'] not in pm['tasks']:
                raise ValueError('Saved four-mode task differs')
            a = read(r['prepared_admission']); state = read(r['state_signature'])
            if a['status'] != 'prepared' or not a['mapping_eligible'] or not a['paired_maps_equal'] or (a['atom_count'], a['cap_count']) != (192, 14):
                raise ValueError('Original source admission unavailable')
            if state != r['state'] or state['endpoint_charges'] != {'Ca': -3, 'La': -2} or state['multiplicities'] != {'Ca': 1, 'La': 1} or state['water_inventory']:
                raise ValueError('Original dry fixed state differs')
            if read(r['maps']['Ca']) != read(r['maps']['La']):
                raise ValueError('Paired physical maps differ')
            mapping = read(r['source_preparation'])['mapping']['source_to_qm']
            if sorted(x['qm_index'] for x in mapping) != list(range(1, 192)):
                raise ValueError('Incomplete source-to-QM membership')
            pqq = [x for x in mapping if x['kind'] == 'opaque_cofactor' and x['element'] != 'H']
            if len(pqq) != 24 or {x['name'] for x in pqq} != set(names) or len({x['fragment'] for x in pqq}) != 1:
                raise ValueError('Exactly24 uniquely named whole-PQQ heavy atoms required')
            byname = {x['name']: x for x in pqq}
            if any(byname[n]['element'] != n[0] for n in names):
                raise ValueError('PQQ name/element mismatch')
            indices = [byname[n]['qm_index'] for n in names]
            kin = kinematic_module.Kinematics(read(r['maps']['Ca'])['context'])
            ids = [x['id'] for x in kin.modes]
            metal_indices = [ids.index('metal_'+axis) for axis in 'xyz']
            for i, k in enumerate(metal_indices):
                mode = kin.modes[k]
                if mode['kind'] != 'metal_translation' or mode['unit'] != 'angstrom' or mode['moving_indices'] != [kin.metal] or not np.array_equal(mode['axis'], np.eye(3)[i]):
                    raise ValueError('Unsupported metal coordinate block')
            task = copy.deepcopy(s['angular_task']); task['active_indices'] = metal_indices+task['active_indices']
            if len(task['active_indices']) != 7 or len(set(task['active_indices'])) != 7:
                raise ValueError('Expected metal3 plus original four donor modes')
            proposal = read(r['native_proposals']['La'])
            if proposal['status'] != 'proposal_available' or not proposal['proposal']['geometry_checks']['pass'] or not proposal['final_geometry']['physical_feasible']:
                raise ValueError('Saved donor proposal unavailable')
            atoms = {}; vectors = {}; original_checks = {}
            for q in CANDIDATES:
                qvec = np.zeros(len(kin.modes)) if q == 'origin' else np.asarray(proposal['proposal']['full_q'])
                if not np.array_equal(qvec[metal_indices], np.zeros(3)):
                    raise ValueError('Original metal was not fixed')
                if not np.array_equal(metal_geometry.full_q(task, qvec[task['active_indices']]), qvec):
                    raise ValueError('Saved motion outside original active subspace')
                for z in METALS:
                    cell = r['matrix'][z][q]; atom = common.xyz(check(cell['xyz']))
                    check(cell['MACE']); check(cell['MACE_request']); check(cell['native_frozen_ALPB_result'])
                    if len(atom) != 192 or atom[0][0] != z or cell['charge'] != state['endpoint_charges'][z] or cell['multiplicity'] != 1:
                        raise ValueError('Original endpoint identity/state differs')
                    atomic_module.check_atoms(atom, cell['charge'])
                    if float(np.max(np.abs(kin.evaluate(qvec)[1]-np.asarray([v[1:] for v in atom])))) > m['settings']['kinematic_replay_A']:
                        raise ValueError('Saved donor-coordinate replay differs')
                    atoms[z, q] = atom
                if atoms['Ca', q][1:] != atoms['La', q][1:] or atoms['Ca', q][0][1:] != atoms['La', q][0][1:]:
                    raise ValueError('Original paired geometry differs')
                original_checks[q] = metal_geometry.final_geometry(kin, task, qvec[task['active_indices']], [v[0] for v in atoms['Ca', q]])
                vectors[q] = qvec
            ca0 = np.asarray([a[1:] for a in atoms['Ca', 'origin']]); ca1 = np.asarray([a[1:] for a in atoms['Ca', 'adaptive_La']])
            if not np.array_equal(ca0[indices+[0]], ca1[indices+[0]]):
                raise ValueError('PQQ or metal changed along original donor path')
            prepared[sample] = {'record': r, 'kin': kin, 'task': task, 'metal_indices': metal_indices,
                                'atoms': atoms, 'vectors': vectors, 'pqq_indices': indices,
                                'pqq': ca0[indices], 'metal': ca0[0], 'original_checks': original_checks}
        if set(prepared) != {0, 3} or prepared[0]['record']['state'] != prepared[3]['record']['state']:
            raise ValueError('Two source states differ')
        # Row-vector convention: X3 ≈ X0 @ R + t; one fit, exact inverse.
        x0, x3 = prepared[0]['pqq'], prepared[3]['pqq']; c0, c3 = x0.mean(axis=0), x3.mean(axis=0)
        u, singular, vt = np.linalg.svd((x0-c0).T@(x3-c3))
        if singular[0] <= 0 or singular[1] <= singular[0]*m['settings']['minimum_relative_second_singular_value']:
            raise ValueError('PQQ alignment is collinear/undefined; planar rank2 is allowed')
        handed = np.eye(3); handed[-1, -1] = 1.0 if np.linalg.det(u@vt) > 0 else -1.0
        rotation = u@handed@vt; translation = c3-c0@rotation
        tol = m['settings']['proper_rotation_tolerance']
        if abs(float(np.linalg.det(rotation))-1) > tol or float(np.max(np.abs(rotation.T@rotation-np.eye(3)))) > tol:
            raise ValueError('Alignment is not a proper rotation')
        reverse_rotation = rotation.T; reverse_translation = -translation@reverse_rotation
        new_points = {3: prepared[0]['metal']@rotation+translation,
                      0: prepared[3]['metal']@reverse_rotation+reverse_translation}
        roundtrip = float(np.max(np.abs((x0@rotation+translation)@reverse_rotation+reverse_translation-x0)))
        if roundtrip > tol:
            raise ValueError('Reciprocal transform roundtrip failed')
        alignment = {'status': 'prepared', 'modeled_not_observed': True, 'named_PQQ_heavy_atoms': names,
                     'source0_to_source3_rotation_row_vectors': rotation.tolist(), 'translation_A': translation.tolist(),
                     'reverse_rotation': reverse_rotation.tolist(), 'reverse_translation_A': reverse_translation.tolist(),
                     'singular_values_A2': singular.tolist(), 'determinant': float(np.linalg.det(rotation)),
                     'fit_RMSD_A': float(np.sqrt(np.mean(np.sum((x0@rotation+translation-x3)**2, axis=1)))),
                     'inverse_roundtrip_max_A': roundtrip, 'fit_calls': 1,
                     'PQQ_indices_by_sample': {str(k): v['pqq_indices'] for k, v in prepared.items()},
                     'transplanted_metal_A_by_host': {str(k): v.tolist() for k, v in new_points.items()},
                     'scope': 'Only mapped metal point transferred; no protein/cofactor atom is transplanted.'}
        write(out/'ALIGNMENT.json', alignment); ap = pin(out/'ALIGNMENT.json')
        for row in rows:
            h = prepared[row['sample']]; r = h['record']; sample = row['sample']; directory = out/f'sample{sample}'; directory.mkdir()
            row.update({k: r[k] for k in ('selection_id', 'maps', 'state_signature', 'source_preparation', 'prepared_admission', 'full_source', 'spectator_admission')})
            row.update(alignment=ap, donor_source_case_id=prepared[3 if sample == 0 else 0]['record']['case_id'],
                       transformed_metal_A=new_points[sample].tolist(), modeled_not_observed=True)
            admission = {'status': 'unavailable', 'pass': False, 'source_case_id': r['case_id'], 'alignment': ap,
                         'old_geometry_checks': h['original_checks'], 'candidate_checks': {}, 'preserved_original_displacement_origin': True}
            try:
                delta = new_points[sample]-h['metal']; displacement = float(np.linalg.norm(delta))
                admission['metal_displacement_A'] = displacement
                if displacement > m['settings']['maximum_metal_displacement_A']+metal_geometry.SETTINGS['final_physical_tolerance_A']:
                    raise ValueError('Transferred metal exceeds original0.8A physical bound')
                for q in CANDIDATES:
                    qvec = h['vectors'][q].copy(); qvec[h['metal_indices']] = delta
                    checks = metal_geometry.final_geometry(h['kin'], h['task'], qvec[h['task']['active_indices']], [v[0] for v in h['atoms']['Ca', q]])
                    rendered = h['kin'].evaluate(qvec)[1]
                    expected = np.asarray([a[1:] for a in h['atoms']['Ca', q]]); expected[0] = new_points[sample]
                    error = float(np.max(np.abs(rendered-expected)))
                    if error > m['settings']['kinematic_replay_A']:
                        raise ValueError('Metal-only kinematic replay changes nonmetal geometry')
                    admission['candidate_checks'][q] = {'pass': True, 'physical_geometry': checks, 'full_q': qvec.tolist(),
                                                        'rendered_vs_exact_saved_max_A': error}
                    for z in METALS:
                        old = r['matrix'][z][q]; original = check(old['xyz']); lines = original.read_text().splitlines(keepends=True)
                        lines[1] = 'Modeled reciprocal source metal position; all nonmetal coordinate lines unchanged\n'
                        lines[2] = z+' '+' '.join(format(float(v), '.17g') for v in new_points[sample])+'\n'
                        xp = directory/f'{z}_{q}.xyz'
                        with xp.open('x') as f: f.write(''.join(lines))
                        actual = common.xyz(xp); atomic = atomic_module.check_atoms(actual, old['charge'])
                        if actual[1:] != h['atoms'][z, q][1:] or xp.read_text().splitlines(keepends=True)[3:] != original.read_text().splitlines(keepends=True)[3:]:
                            raise ValueError('Nonmetal coordinates/lines changed')
                        if not np.array_equal(actual[0][1:], new_points[sample]):
                            raise ValueError('Metal point serialization changed')
                        row['matrix'][z][q] = {'status': 'prepared', 'xyz': pin(xp), 'charge': old['charge'],
                            'multiplicity': old['multiplicity'], 'atom_count': old['atom_count'],
                            'valence_electrons': old['GFN2_valence_electrons'], 'atomic_check': atomic,
                            'original_xyz': old['xyz'], 'original_MACE': old['MACE'], 'original_MACE_eV': old['MACE_eV'],
                            'original_frozen_ALPB': old['native_frozen_ALPB_result'], 'geometry_checks': admission['candidate_checks'][q]}
                for q in CANDIDATES:
                    ca = common.xyz(check(row['matrix']['Ca'][q]['xyz'])); la = common.xyz(check(row['matrix']['La'][q]['xyz']))
                    if ca[1:] != la[1:] or ca[0][1:] != la[0][1:]: raise ValueError('New paired coordinates differ')
                admission.update(status='prepared', **{'pass': True}, nonmetal_lines_unchanged=True, paired_geometry_identical=True)
                row.update(status='prepared', reason=None)
            except Exception as exc:
                admission.update(reason=f'{type(exc).__name__}: {exc}', traceback=traceback.format_exc())
                row['reason'] = admission['reason']
                for z in METALS:
                    for q in CANDIDATES: row['matrix'][z][q].update(status='unavailable', reason=admission['reason'])
            write(directory/'GEOMETRY_ADMISSION.json', admission); row['geometry_admission'] = pin(directory/'GEOMETRY_ADMISSION.json')
            for z in METALS:
                for q in CANDIDATES: row['matrix'][z][q]['geometry_admission'] = row['geometry_admission']
        receipt['status'] = 'complete'
    except BaseException as exc:
        receipt.update(status='failed', error=f'{type(exc).__name__}: {exc}', traceback=traceback.format_exc())
        raise
    finally:
        receipt.update(wall_seconds=time.monotonic()-began, process_CPU_seconds=time.process_time()-cpu0,
                       peak_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
        write(out/'EXECUTION.json', receipt)
        result = {'protocol': PROTOCOL, 'manifest': mp, 'inventory': m['inventory'], 'terminal': True,
                  'status': 'prepared' if all(r['status'] == 'prepared' for r in rows) else 'unavailable',
                  'rows': rows, 'declared_sources': 2, 'declared_new_cells': 8, 'reused_original_cells': 8,
                  'prepared_sources': sum(r['status'] == 'prepared' for r in rows), 'execution': pin(out/'EXECUTION.json'),
                  'new_energy_calls': 0, 'new_optimizations': 0, 'modeled_not_observed': True, 'production_changed': False}
        write(out/'PREPARATION.json', result)
    print(json.dumps({'preparation': pin(out/'PREPARATION.json'), 'prepared_sources': result['prepared_sources'], 'declared_sources': 2}))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__); sub = p.add_subparsers(dest='action', required=True)
    f = sub.add_parser('freeze'); f.add_argument('--inventory', required=True); f.add_argument('--plan', required=True); f.add_argument('--output', required=True)
    e = sub.add_parser('execute'); e.add_argument('--manifest', required=True); e.add_argument('--manifest-sha', required=True)
    a = p.parse_args()
    if a.action == 'freeze': freeze(a.inventory, a.plan, a.output)
    else: execute(a.manifest, a.manifest_sha)
