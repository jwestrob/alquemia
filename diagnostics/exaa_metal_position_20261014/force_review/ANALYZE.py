#!/usr/bin/env python3
"""Freeze metadata or analyze twenty real saved force arrays inside Slurm."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import socket
import subprocess
import sys
import time
import traceback

PROTOCOL = 'ExaA_saved_metal_donor_force_diagnostic_v1'
INVENTORY_SHA = '757a13523549697fd52cf9e94b3ba2ffdd01c463f3969fd8a5d01c1459d28343'
THREADS = ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','BLIS_NUM_THREADS','OPENMM_CPU_THREADS')


def pin(path):
    p = Path(path).resolve()
    return {'path': str(p), 'sha256': hashlib.sha256(p.read_bytes()).hexdigest()}


def check(p):
    if pin(p['path']) != p:
        raise ValueError('Pinned artifact changed: '+p['path'])
    return Path(p['path'])


def read(p):
    return json.loads((check(p) if isinstance(p, dict) else Path(p)).read_text())


def write(p, data):
    p = Path(p); p.parent.mkdir(parents=True, exist_ok=True)
    with p.open('x') as f:
        json.dump(data, f, indent=2, allow_nan=False); f.write('\n')


def freeze(root, output):
    root = Path(root).resolve(); out = Path(output).resolve()
    inv = pin(root/'workspaces/exaa_donor_response_20261014/inventory_v1/INVENTORY.json')
    if inv['sha256'] != INVENTORY_SHA:
        raise ValueError('Not the exact two-source inventory')
    prep = pin(root/'workspaces/exaa_metal_position_20261014/preparation_v1/execution/PREPARATION.json')
    mace = pin(root/'workspaces/exaa_metal_position_20261014/scoring_v1/execution/FINAL_COLLECTION.json')
    pm = read(read(prep)['manifest'])
    sources = []
    for s in read(inv)['rows']:
        proposals = {z: read(p) for z, p in s['native_proposals'].items()}
        proposal_manifest = proposals['Ca']['manifest']
        if proposals['La']['manifest'] != proposal_manifest:
            raise ValueError('Two proposal manifests differ')
        sources.append({'source_case_id': s['case_id'], 'role': s['role'],
                        'proposal_manifest': proposal_manifest,
                        'original_proposals': s['native_proposals']})
    out.mkdir(parents=True, exist_ok=False); impl = out/'implementation'; impl.mkdir()
    pins = {}
    for name in ('ANALYZE.py','RUN.sbatch'):
        dst = impl/name; shutil.copyfile(Path(__file__).with_name(name), dst); pins[name] = pin(dst)
    m = {'protocol': PROTOCOL, 'inventory': inv, 'preparation': prep, 'MACE_collection': mace,
         'plan': pin(Path(__file__).with_name('PLAN.md')), 'implementation': pins,
         'source_snapshot': pm['source_snapshot'], 'source_modules': pm['source_modules'],
         'python': pm['python'], 'sources': sources, 'output': str(out/'execution'),
         'declared_cells': 20, 'coordinate_replay_tolerance_A': 1e-8,
         'saved_angular_gradient_replay_tolerance_kcal_mol_rad': 1e-8,
         'metal_projection_replay_tolerance_eV_A': 1e-12,
         'force_hamiltonian': 'MACE OMOL total vacuum energy only',
         'new_model_calls': 0, 'new_optimizations': 0}
    write(out/'manifest.json', m); print(json.dumps(pin(out/'manifest.json')))


def placement():
    job = os.environ.get('SLURM_JOB_ID')
    if not job or not job.isdigit():
        raise RuntimeError('Numerical saved-force analysis requires Slurm')
    text = subprocess.check_output(['scontrol','show','job','-o',job], text=True)
    d = dict(v.split('=',1) for v in text.split() if '=' in v)
    if (d.get('JobState'), d.get('NumNodes'), d.get('NumCPUs')) != ('RUNNING','1','1'):
        raise RuntimeError('Expected one CPU on one running node')
    if int(os.environ.get('SLURM_CPUS_ON_NODE','0')) != 1 or any(os.environ.get(k) != '1' for k in THREADS):
        raise RuntimeError('Expected one allocated numerical thread')
    return {'job_id': job, 'hostname': socket.gethostname(), 'scontrol': text,
            'affinity': sorted(os.sched_getaffinity(0)), 'allocated_CPUs': 1}


def execute(manifest, sha):
    place = placement()  # Before NumPy or physical-helper imports.
    mp = {'path': str(Path(manifest).resolve()), 'sha256': sha}; m = read(mp)
    if m['protocol'] != PROTOCOL or pin(__file__) != m['implementation']['ANALYZE.py']:
        raise ValueError('Wrong analysis implementation')
    out = Path(m['output']); out.mkdir(parents=True, exist_ok=False)
    began = time.monotonic(); cpu0 = time.process_time(); rows = []
    for s in m['sources']:
        for z in ('Ca','La'):
            for placement_name, names in [('original',('origin','adaptive_Ca','adaptive_La')),
                                           ('transplanted',('origin','adaptive_La'))]:
                for name in names:
                    rows.append({'source_case_id': s['source_case_id'], 'role': s['role'], 'metal': z,
                                 'candidate': name, 'placement': placement_name,
                                 'status': 'unavailable', 'reason': 'analysis not completed',
                                 'donor_projection_status': 'unavailable'})
    write(out/'START.json', {'manifest': mp, 'placement': place, 'declared_cells': len(rows)})
    failure = None; sources_out = []
    try:
        for p in m['implementation'].values(): check(p)
        for p in m['source_modules'].values(): check(p)
        check(m['python']); check(m['plan'])
        if pin(sys.executable) != m['python']:
            raise ValueError('Wrong numerical Python interpreter')
        sys.path.insert(0, m['source_snapshot'])
        import numpy as np
        import adaptive_force_diagnostic as projection
        import mace_site_kinematics as kinematics
        import affordable_common as common
        for module in (projection, kinematics, common):
            if pin(module.__file__) != m['source_modules'][Path(module.__file__).name]:
                raise ValueError('Wrong archived physical helper')
        ev_to_kcal = projection.EV_TO_KCAL
        inv = {r['case_id']: r for r in read(m['inventory'])['rows']}
        prep = {r['source_case_id']: r for r in read(m['preparation'])['rows']}
        scored = {r['source_case_id']: r for r in read(m['MACE_collection'])['rows']}
        for spec in m['sources']:
            cid = spec['source_case_id']; src = inv[cid]; pre = prep[cid]; sc = scored[cid]
            if pre['status'] != 'prepared' or pre['maps'] != src['maps'] or pre['state_signature'] != src['state_signature']:
                raise ValueError('Original/transplanted source provenance differs')
            state = read(src['state_signature']); maps = {z: read(p) for z,p in src['maps'].items()}
            if maps['Ca'] != maps['La'] or state != src['state']:
                raise ValueError('Paired maps/state differs')
            kin = kinematics.Kinematics(maps['Ca']['context'])
            modes = kin.modes; metal_indices = [next(i for i,x in enumerate(modes) if x['id']=='metal_'+a) for a in 'xyz']
            for k,i in enumerate(metal_indices):
                mode = modes[i]
                if mode['unit'] != 'angstrom' or mode['moving_indices'] != [kin.metal] or not np.array_equal(mode['axis'],np.eye(3)[k]):
                    raise ValueError('Metal mode is not pure single-ion translation')
            original0 = np.asarray([a[1:] for a in common.xyz(check(src['matrix']['Ca']['origin']['xyz']))])
            moved0 = np.asarray([a[1:] for a in common.xyz(check(pre['matrix']['Ca']['origin']['xyz']))])
            delta = moved0[0]-original0[0]; length = float(np.linalg.norm(delta))
            if length <= 1e-12 or not np.array_equal(original0[1:],moved0[1:]):
                raise ValueError('Undefined or nonmetal transplant displacement')
            unit = delta/length
            pm = read(spec['proposal_manifest']); tasks = {t['metal']: t for t in pm['tasks']}
            props = {z: read(p) for z,p in spec['original_proposals'].items()}
            q_saved = {'origin': np.zeros(len(modes))}
            for z in ('Ca','La'):
                if props[z]['status'] != 'proposal_available' or props[z]['manifest'] != spec['proposal_manifest']:
                    raise ValueError('Original admitted proposal missing')
                q_saved['adaptive_'+z] = np.asarray(props[z]['proposal']['full_q'])
            sources_out.append({'source_case_id': cid, 'role': src['role'], 'metal_displacement_A': delta.tolist(),
                                'metal_displacement_norm_A': length, 'transplant_unit_vector': unit.tolist(),
                                'maps': src['maps'], 'state_signature': src['state_signature'],
                                'original_proposals': spec['original_proposals'], 'proposal_manifest': spec['proposal_manifest'],
                                'original_boundary_flags': {z: props[z]['boundary_flag'] for z in ('Ca','La')},
                                'old_active_mode_ids': {z: tasks[z]['active_mode_ids'] for z in ('Ca','La')}})
            for row in [r for r in rows if r['source_case_id'] == cid]:
                try:
                    z = row['metal']; name = row['candidate']; translated = row['placement'] == 'transplanted'
                    cell = sc['matrix'][z][name]['transplanted'] if translated else src['matrix'][z][name]
                    receipt_pin = cell['receipt'] if translated else cell['MACE']; native = read(receipt_pin)
                    req = read(native['request']); atoms = common.xyz(check(cell['xyz']))
                    coords = np.asarray([a[1:] for a in atoms]); force = np.load(check(native['forces']), allow_pickle=False)
                    if (native['status'] != 'complete' or req['xyz'] != cell['xyz'] or
                        req['charge'] != state['endpoint_charges'][z] or req['multiplicity'] != 1 or
                        native['model'] != src['matrix'][z]['origin']['MACE_model'] or
                        native['state']['atoms'] != 192 or native['state']['selected_metal_index'] != 0 or
                        native['state']['charge'] != req['charge'] or native['state']['spin_multiplicity'] != 1 or
                        native['force_definition'] != 'negative Cartesian total vacuum OMOL energy gradient, eV/angstrom' or
                        atoms[0][0] != z or coords.shape != (192,3) or force.shape != (192,3) or not np.isfinite(force).all()):
                        raise ValueError('Saved force/coordinate/model/state mismatch')
                    base_atoms = common.xyz(check(src['matrix'][z][name]['xyz']))
                    if atoms[1:] != base_atoms[1:] or not np.array_equal(coords[0],moved0[0] if translated else original0[0]):
                        raise ValueError('Unexpected physical coordinate change')
                    fm = force[0]; norm = float(np.linalg.norm(fm)); parallel = float(fm@unit)
                    row.update(status='available', reason=None, xyz=cell['xyz'], MACE=receipt_pin,
                               request=native['request'], forces=native['forces'],
                               metal_force_eV_A=fm.tolist(), metal_force_norm_eV_A=norm,
                               metal_force_projection_toward_transplant_eV_A=parallel,
                               metal_force_kcal_mol_A=(fm*ev_to_kcal).tolist(),
                               metal_force_norm_kcal_mol_A=norm*ev_to_kcal,
                               metal_force_projection_toward_transplant_kcal_mol_A=parallel*ev_to_kcal,
                               local_energy_derivative_toward_transplant_kcal_mol_A=-parallel*ev_to_kcal)
                    try:
                        q = q_saved[name].copy()
                        if not np.array_equal(q[metal_indices],np.zeros(3)):
                            raise ValueError('Original donor proposal moved metal')
                        if translated: q[metal_indices] = delta
                        rendered, vectors, raw, normalized, lengths = projection.project(kin.data,q,force)
                        residual = float(np.max(np.abs(rendered-coords)))
                        if residual > m['coordinate_replay_tolerance_A']:
                            raise ValueError('Saved coordinate kinematic replay failed')
                        metal_residual = float(np.max(np.abs(-raw[metal_indices]/ev_to_kcal-fm)))
                        if metal_residual > m['metal_projection_replay_tolerance_eV_A']:
                            raise ValueError('Direct/projected metal force differs')
                        own_replay = None
                        if not translated and name in ('origin','adaptive_'+z):
                            point = props[z]['origin' if name=='origin' else 'proposal']
                            if point['forces'] != native['forces'] or point['MACE'] != receipt_pin:
                                raise ValueError('Own saved-gradient force receipt differs')
                            own_replay = float(np.max(np.abs(raw[tasks[z]['active_indices']]-point['gradient_kcal_mol_rad'])))
                            if own_replay > m['saved_angular_gradient_replay_tolerance_kcal_mol_rad']:
                                raise ValueError('Own archived angular gradient does not replay')
                        donor = [{'id': mode['id'], 'kind': mode['kind'], 'unit': mode['unit'],
                                  'active_in_old_search': mode['id'] in tasks[z]['active_mode_ids'],
                                  'energy_derivative_kcal_per_unit': float(raw[i]),
                                  'physical_heavy_norm_A_per_unit': float(lengths[i]),
                                  'normalized_energy_derivative_kcal_mol_A': float(normalized[i])}
                                 for i,mode in enumerate(modes) if mode['unit']=='radian']
                        row.update(donor_projection_status='available', donor_projection_reason=None,
                                   full_q=q.tolist(), donor_modes=donor, coordinate_replay_max_A=residual,
                                   metal_projection_replay_max_eV_A=metal_residual,
                                   own_saved_gradient_replay_max_kcal_mol_rad=own_replay)
                    except Exception as exc:
                        row.update(donor_projection_status='unavailable', donor_projection_reason=str(exc))
                except Exception as exc:
                    row.update(status='unavailable', reason=str(exc))
        if len(rows) != 20:
            raise ValueError('Declared cell denominator changed')
    except Exception as exc:
        failure = {'reason': str(exc), 'traceback': traceback.format_exc()}
    elapsed = time.monotonic()-began
    result = {'protocol': PROTOCOL, 'manifest': mp, 'sources': sources_out, 'rows': rows,
              'terminal': True, 'declared_cells': 20, 'available_metal_force_cells': sum(r['status']=='available' for r in rows),
              'available_donor_projection_cells': sum(r['donor_projection_status']=='available' for r in rows),
              'failure': failure, 'force_hamiltonian': m['force_hamiltonian'],
              'new_model_calls': 0, 'new_optimizations': 0, 'composite_gradient': None,
              'interpretation': 'Saved local forces only; no finite-work, uncertainty or optimized-state claim.'}
    write(out/'RESULT.json',result)
    lines = ['# Saved metal forces', '', 'MACE vacuum forces only; positive projection points toward the actual transplanted ion position.', '',
             '| Host | Placement | Geometry | Metal | Norm (kcal/mol/Å) | Projection (kcal/mol/Å) | Donor replay |',
             '|---|---|---|---|---:|---:|---|']
    for r in rows:
        values = (f"{r['metal_force_norm_kcal_mol_A']:.6f}",f"{r['metal_force_projection_toward_transplant_kcal_mol_A']:.6f}") if r['status']=='available' else ('unavailable','unavailable')
        lines.append(f"| {r['role']} | {r['placement']} | {r['candidate']} | {r['metal']} | {values[0]} | {values[1]} | {r['donor_projection_status']} |")
    (out/'REPORT.md').write_text('\n'.join(lines)+'\n')
    write(out/'EXECUTION.json', {'status':'complete' if failure is None else 'failed', 'manifest':mp,
          'placement':place,'result':pin(out/'RESULT.json'),'wall_seconds':elapsed,
          'process_CPU_seconds':time.process_time()-cpu0,'allocated_CPU_seconds':elapsed,
          'new_model_calls':0,'new_optimizations':0,'failure':failure})
    if failure: raise RuntimeError(failure['reason'])


if __name__ == '__main__':
    p = argparse.ArgumentParser(); p.add_argument('operation',choices=('freeze','execute'))
    p.add_argument('--root'); p.add_argument('--output'); p.add_argument('--manifest'); p.add_argument('--manifest-sha')
    a = p.parse_args()
    if a.operation=='freeze': freeze(a.root,a.output)
    else: execute(a.manifest,a.manifest_sha)
