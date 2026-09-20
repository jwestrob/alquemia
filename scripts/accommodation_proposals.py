"""Native MACE proposes one physical pose; two actual composite energies select it."""
from __future__ import annotations
import argparse
from collections import Counter
import fcntl
import json
import os
from pathlib import Path
import shutil
import time
import traceback
import numpy as np

from affordable_common import HA_TO_KCAL, InvalidArtifact, cache_key, read_json, record, verify, write_new, xyz
from mace_hybrid import EV_TO_KCAL, check_atoms, write_xyz
from mace_site_kinematics import Kinematics
from accommodation_nonlinear import WarmGPU, relative_components, contrast_components
from accommodation_torsion_profiles import role_mode, geometry_check
from compact_solvation import input_text, completed, diagnostics
from compact_solvation_compare import mix_pair
from compact_solvation_scanner import decision

PROTOCOL = 'native_MACE_proposal_primary_composite_selection_v1'
SETTINGS = {'method': 'L-BFGS-B', 'bounds_radian': [-.8, .8], 'maxiter': 24,
            'maxfun': 80, 'maxls': 10, 'gtol': .2, 'ftol': 1e-12,
            'boundary_flag_margin_radian': .02, 'selection_decrease_kcal_mol': .1,
            'origin_MACE_replay_tolerance_kcal_mol': .01,
            'new_overlap_H_A': .55, 'new_overlap_heavy_A': 1.0}
ROLE_ORDER = ('anchor_glutamate', 'extra_acidic_ligand_homolog')
PLM = ('PQQSEQ_83440678cbbd658047c9', 'PQQSEQ_07ab500e3df76b30d71c')


def check_low_source(manifest, task, origin, expected_orca):
    """Reuse only actual primary-recipe outputs on identical parsed coordinates."""
    m = read_json(manifest)
    if m['orca'] != expected_orca:
        raise InvalidArtifact('q0 ORCA executable differs')
    if (task['charge'], task['multiplicity']) != (origin['charge'], origin['multiplicity']):
        raise InvalidArtifact('q0 electronic state differs')
    if xyz(verify(task['xyz'])) != xyz(verify(origin['xyz'])):
        raise InvalidArtifact('q0 archived coordinates differ')
    if verify(task['input']).read_text() != input_text(origin['charge'], 1, task['medium'], 'native'):
        raise InvalidArtifact('q0 is not the exact primary native GFN2 recipe')
    pin = m.get('reused', {}).get(task['task_id']) or completed(manifest, task['task_id'])
    if pin is None or completed(verify(pin['manifest']), pin['task_id']) != pin:
        raise InvalidArtifact('q0 execution receipt unavailable or changed')
    actual_manifest = read_json(verify(pin['manifest']))
    actual = next(t for t in actual_manifest['all_tasks'] if t['task_id'] == pin['task_id'])
    if (actual_manifest['orca'] != expected_orca or
            xyz(verify(actual['xyz'])) != xyz(verify(origin['xyz'])) or
            verify(actual['input']).read_text() != verify(task['input']).read_text()):
        raise InvalidArtifact('q0 reused receipt describes another recipe or geometry')
    audit = diagnostics(pin, actual)
    if audit['charge_sanity_status'] != 'pass':
        raise InvalidArtifact('q0 native state check failed')
    return {'receipt': pin, 'task': actual, 'audit': audit,
            'energy_hartree': pin['energy_hartree']}


def prepare(references, torsion, primary_result, compact, agreement, output, cpu_python, gpu_python):
    refs, old, primary, cm = map(read_json, (references, torsion, primary_result, compact))
    inv = read_json(verify(cm['inventory']))
    if refs['model'] != old['model'] or refs['model'] != inv['model']:
        raise InvalidArtifact('source checkpoints differ')
    tm = read_json(verify(primary['manifest']))
    if tm['design'] != record(torsion):
        raise InvalidArtifact('PLM primary results have another source')
    lm_path = verify(tm['low_manifest']); lm = read_json(lm_path)
    mr = read_json(verify(primary['MACE_result']))
    if mr['model'] != old['model']:
        raise InvalidArtifact('PLM MACE source checkpoint differs')
    cases = refs['cases'] + [c for c in old['cases'] if c['case_id'] in PLM]
    if len(cases) != 30 or len({c['case_id'] for c in cases}) != 30 or len(refs['cases']) != 28:
        raise InvalidArtifact('declared 30-context population differs')
    out = Path(output).resolve(); out.mkdir(parents=True, exist_ok=False)
    tasks, audits = [], []
    for c in cases:
        parent = read_json(verify(c['source']['parent']))
        selector = parent['fixed_core']['requested_roles']['extra_acidic_ligand_homolog']
        rn = selector['resname'] if isinstance(selector, dict) else selector.split(':')[1][:3]
        roles = ['anchor_glutamate'] + (['extra_acidic_ligand_homolog'] if rn == 'ASP' else [])
        expected = [role_mode(parent, role) for role in roles]
        if c['modes'] != dict(zip(roles, expected)):
            raise InvalidArtifact('active roles differ from original chemistry')
        for metal in ('Ca', 'La'):
            ep = c['origins'][metal]
            design_origin = ep['xyz']; origin_difference = 0.0
            if c['case_id'] in PLM:
                point = next(p for p in old['points'] if p['case_id'] == c['case_id'] and
                             p['metal'] == metal and p['point'] == 'origin')
                a, b = xyz(verify(ep['xyz'])), xyz(verify(point['xyz']))
                if [r[0] for r in a] != [r[0] for r in b]:
                    raise InvalidArtifact('executed q0 composition differs from design')
                origin_difference = float(np.max(np.abs(np.array([r[1:] for r in a]) - np.array([r[1:] for r in b]))))
                if origin_difference > 1e-12:
                    raise InvalidArtifact('executed q0 differs from existing mapping tolerance')
                ep = {**ep, 'xyz': point['xyz']}
            kin = Kinematics(read_json(verify(c['maps'][metal]))['context'])
            atoms = xyz(verify(ep['xyz'])); names = [m['id'] for m in kin.modes]
            active = [names.index(name) for name in expected]; zero = np.zeros(len(names))
            if not np.allclose(kin.evaluate(zero)[1], [a[1:] for a in atoms], atol=1e-12, rtol=0):
                raise InvalidArtifact('origin map does not reproduce archived coordinates')
            check_atoms(atoms, ep['charge'])
            origin_checks = geometry_check(kin, zero, [a[0] for a in atoms])
            if not origin_checks['pass']:
                raise InvalidArtifact('origin physical mapping unavailable')
            task = {'task_id': c['case_id'] + '__' + metal, 'case_id': c['case_id'], 'metal': metal,
                    'xyz': ep['xyz'], 'charge': ep['charge'], 'multiplicity': ep['multiplicity'],
                    'mapping': c['maps'][metal], 'source_preparation': c['preparation'],
                    'active_roles': roles, 'active_mode_ids': expected, 'active_indices': active,
                    'mode_count': len(names), 'origin_geometry_checks': origin_checks,
                    'design_origin_xyz': design_origin, 'executed_q0_difference_A': origin_difference,
                    'label_scope': c['label_scope'], 'q0_status': 'unavailable', 'q0': None}
            try:
                lows = {}
                if c['case_id'] in PLM:
                    point_id = c['case_id'] + '__origin__' + metal
                    native = next(r for r in mr['rows'] if r['task_id'] == point_id)
                    if native['status'] != 'complete' or xyz(verify(native['xyz'])) != atoms:
                        raise InvalidArtifact('PLM native MACE origin differs')
                    native_energy = native['energy_eV']; native_pin = primary['MACE_result']
                    for medium in ('vacuum', 'alpb'):
                        lt = next(t for t in lm['all_tasks'] if t['point_id'] == point_id and t['medium'] == medium)
                        lows[medium] = check_low_source(lm_path, lt, ep, refs['orca'])
                else:
                    item = next(x for x in inv['cases'] if x['case_id'] == c['case_id'])
                    nep = item['representations']['context']['endpoints'][metal]
                    native_pin = nep['native_MACE_receipt']; native = read_json(verify(native_pin))
                    if native['energy_eV'] != nep['native_MACE_energy_eV'] or xyz(verify(nep['xyz'])) != atoms:
                        raise InvalidArtifact('reference native MACE origin differs')
                    native_energy = native['energy_eV']
                    for medium in ('vacuum', 'alpb'):
                        lt = next(t for t in cm['all_tasks'] if
                                  (t['case_id'], t['representation'], t['metal'], t['medium']) ==
                                  (c['case_id'], 'context', metal, medium))
                        lows[medium] = check_low_source(compact, lt, ep, refs['orca'])
                if read_json(verify(lows['vacuum']['audit']['parameter_export'])) != read_json(verify(lows['alpb']['audit']['parameter_export'])):
                    raise InvalidArtifact('q0 GFN2 media parameter sets differ')
                task.update(q0_status='available', q0={'native_MACE_receipt': native_pin, 'low': lows,
                    'components': {'MACE_eV': native_energy,
                        'GFN2_vacuum_hartree': lows['vacuum']['energy_hartree'],
                        'GFN2_ALPB_hartree': lows['alpb']['energy_hartree']}})
            except Exception as exc:
                task['q0_reason'] = str(exc)
            tasks.append(task); audits.append({k: task[k] for k in ('task_id', 'q0_status')} | {'reason': task.get('q0_reason')})
    impl = out / 'implementation'; impl.mkdir(); pins = {}
    for p in Path(__file__).parent.glob('*.py'):
        dest = impl / p.name; shutil.copyfile(p, dest); pins[p.name] = record(dest)
    # Released composite bands are read from the completed reference comparison.
    comparison_path = Path(compact).parent / 'comparison_v1.json'
    comparison = read_json(comparison_path)
    m = {'protocol_id': PROTOCOL, 'settings': SETTINGS, 'tasks': tasks, 'cases': cases,
         'references': record(references), 'torsion': record(torsion), 'primary_result': record(primary_result),
         'compact': record(compact), 'agreement': record(agreement), 'frozen_comparison': record(comparison_path),
         'frozen_context_bands': comparison['calibration']['context']['bands'],
         'model': refs['model'], 'software': refs['software'], 'orca': refs['orca'], 'implementation': pins,
         'cpu_python': str(Path(cpu_python).absolute()), 'gpu_python': str(Path(gpu_python).absolute()),
         'cpu_executable': record(cpu_python), 'gpu_executable': record(gpu_python),
         'algorithm_starts': 60, 'maximum_new_GFN2_singlepoints': 120, 'new_DFT_calls': 0,
         'threshold_refitted': False, 'baseline_changed': False,
         'resources': {'GPU_cpus': 32, 'GPUs': 1, 'GPU_host_mem_MiB': 200000,
                       'GFN2_cpus': 64, 'mpi_ranks': 8, 'concurrent_tasks': 8}}
    write_new(out / 'manifest.json', m)
    audit = {'manifest': record(out / 'manifest.json'), 'rows': audits, 'endpoint_denominator': 60,
             'q0_reuse_available': sum(t['q0_status'] == 'available' for t in tasks),
             'new_molecular_calls': 0}
    write_new(out / 'Q0_AUDIT.json', audit)
    return validate(out / 'manifest.json')


def validate(manifest):
    m = read_json(manifest)
    fold_scope = m.get('evaluation_scope') == 'primary225_fold_transfer_v1'
    if fold_scope:
        from accommodation_fold_proposals import validate_scope
        validate_scope(m)
    count = 416 if fold_scope else 60
    cases = 208 if fold_scope else 30
    if m['protocol_id'] != PROTOCOL or m['settings'] != SETTINGS or len(m['tasks']) != count:
        raise InvalidArtifact('frozen proposal policy differs')
    if len(m['cases']) != cases or {(t['case_id'], t['metal']) for t in m['tasks']} != {(c['case_id'], z) for c in m['cases'] for z in ('Ca', 'La')}:
        raise InvalidArtifact('complete declared context/endpoint scope differs')
    sources = ('fold_design', 'fold_comparison', 'proposal_reference', 'sources') if fold_scope else ('references', 'torsion', 'primary_result', 'compact')
    for key in (*sources, 'agreement', 'software', 'orca', 'cpu_executable', 'gpu_executable', 'frozen_comparison'):
        verify(m[key])
    for pin in m['implementation'].values():
        verify(pin)
    from mace_omol import model
    if model(verify(m['software'])) != m['model']:
        raise InvalidArtifact('native model differs')
    for t in m['tasks']:
        for key in ('xyz', 'mapping', 'source_preparation'):
            verify(t[key])
        check_atoms(xyz(verify(t['xyz'])), t['charge'])
        if t['multiplicity'] != 1:
            raise InvalidArtifact('unsupported spin')
    from affordable_common import paired
    for c in m['cases']:
        ca, la = [next(t for t in m['tasks'] if t['case_id'] == c['case_id'] and t['metal'] == z) for z in ('Ca', 'La')]
        paired(verify(la['xyz']), verify(ca['xyz']), la['charge'], ca['charge'])
    return {'status': 'validated', 'manifest': record(manifest), 'starts': count,
            'q0_available': sum(t['q0_status'] == 'available' for t in m['tasks']),
            'new_molecular_calls': 0}


class MACEProposal:
    def __init__(self, manifest, task, gpu):
        self.manifest = record(manifest); self.task = task; self.gpu = gpu
        self.kin = Kinematics(read_json(verify(task['mapping']))['context'])
        self.atoms = xyz(verify(task['xyz']))
        self.directory = Path(manifest).parent / 'proposals' / task['task_id']
        self.directory.mkdir(parents=True, exist_ok=True)
        self.requests = []

    def evaluate(self, active_q, purpose):
        active_q = np.asarray(active_q, dtype=float)
        if active_q.shape != (len(self.task['active_indices']),) or not np.isfinite(active_q).all() or np.max(np.abs(active_q)) > .8 + 1e-12:
            raise InvalidArtifact('trial outside declared domain')
        q = np.zeros(self.task['mode_count']); q[self.task['active_indices']] = active_q
        _, coords, _, jac = self.kin.evaluate(q)
        key = cache_key({'manifest': self.manifest, 'task': self.task['task_id'], 'q': active_q.tolist()})
        directory = self.directory / 'evaluations' / key; path = directory / 'result.json'
        if path.exists():
            r = read_json(path)
            if r['cache_key'] != key or r['manifest'] != self.manifest:
                raise InvalidArtifact('MACE cache differs')
            self.requests.append({'purpose': purpose, 'result': record(path), 'reused': True})
            if r['status'] != 'complete':
                raise InvalidArtifact('previous failed proposal evaluation retained')
            return r
        directory.mkdir(parents=True, exist_ok=False)
        r = {'cache_key': key, 'manifest': self.manifest, 'task_id': self.task['task_id'],
             'active_q_radian': active_q.tolist(), 'full_q': q.tolist(), 'status': 'failed', 'MACE': None}
        try:
            checks = geometry_check(self.kin, q, [a[0] for a in self.atoms]); r['geometry_checks'] = checks
            if not checks['pass']:
                raise InvalidArtifact('unsupported physical proposal geometry')
            xp = directory / 'context.xyz'
            if np.array_equal(active_q, np.zeros_like(active_q)):
                shutil.copyfile(verify(self.task['xyz']), xp)
            else:
                write_xyz(xp, [(a[0], *p) for a, p in zip(self.atoms, coords)])
            request = {'manifest': self.manifest, 'task_id': self.task['task_id'], 'xyz': record(xp),
                       'charge': self.task['charge'], 'multiplicity': self.task['multiplicity']}
            request_path = directory / 'mace_request.json'; write_new(request_path, request)
            native, pin = self.gpu.evaluate(record(request_path)); r['MACE'] = pin
            if native['status'] != 'complete':
                raise InvalidArtifact('native MACE unavailable: ' + native.get('reason', 'unknown'))
            forces = np.load(verify(native['forces']), allow_pickle=False)
            gradient = -np.einsum('mij,ij->m', jac[self.task['active_indices']], forces) * EV_TO_KCAL
            r.update(status='complete', coordinate=record(xp), MACE_eV=native['energy_eV'],
                     gradient_kcal_mol_rad=gradient.tolist(), forces=native['forces'])
        except Exception as exc:
            r.update(reason=str(exc), traceback=traceback.format_exc())
        write_new(path, r); self.requests.append({'purpose': purpose, 'result': record(path), 'reused': False})
        if r['status'] != 'complete':
            raise InvalidArtifact(r['reason'])
        return r


def optimize(manifest, task, gpu):
    from scipy.optimize import minimize
    ev = MACEProposal(manifest, task, gpu); path = ev.directory / 'result.json'
    if path.exists():
        r = read_json(path)
        if r['manifest'] != record(manifest):
            raise InvalidArtifact('proposal belongs to another manifest')
        return record(path)
    start = time.monotonic(); accepted = []
    r = {'task_id': task['task_id'], 'manifest': record(manifest), 'status': 'unavailable',
         'origin': None, 'proposal': None, 'optimizer': None}
    try:
        if task['q0_status'] != 'available':
            raise InvalidArtifact('exact composite origin unavailable: ' + task.get('q0_reason', 'unknown'))
        origin = ev.evaluate(np.zeros(len(task['active_indices'])), 'origin'); r['origin'] = origin
        replay = (origin['MACE_eV'] - task['q0']['components']['MACE_eV']) * EV_TO_KCAL
        r['origin_MACE_replay_delta_kcal_mol'] = replay
        if abs(replay) > SETTINGS['origin_MACE_replay_tolerance_kcal_mol']:
            raise InvalidArtifact('native origin replay failed')
        def objective(q):
            value = ev.evaluate(q, 'optimizer')
            return ((value['MACE_eV'] - origin['MACE_eV']) * EV_TO_KCAL,
                    np.asarray(value['gradient_kcal_mol_rad']))
        opt = minimize(objective, np.zeros(len(task['active_indices'])), method=SETTINGS['method'], jac=True,
                       bounds=[tuple(SETTINGS['bounds_radian'])] * len(task['active_indices']),
                       callback=lambda q: accepted.append(q.tolist()),
                       options={k: SETTINGS[k] for k in ('maxiter', 'maxfun', 'maxls', 'gtol', 'ftol')})
        candidate = ev.evaluate(opt.x, 'final_proposal')
        distance = float(np.min(.8 - np.abs(opt.x)))
        r.update(optimizer={'success': bool(opt.success), 'message': str(opt.message), 'iterations': int(opt.nit),
                           'function_evaluations': int(opt.nfev), 'gradient_evaluations': int(opt.njev)},
                 proposal=candidate, boundary_distance_radian=distance,
                 boundary_flag=distance < SETTINGS['boundary_flag_margin_radian'],
                 max_raw_active_gradient_kcal_mol_rad=float(np.max(np.abs(candidate['gradient_kcal_mol_rad']))),
                 MACE_proposal_work_kcal_mol=(candidate['MACE_eV'] - origin['MACE_eV']) * EV_TO_KCAL)
        if not opt.success:
            raise InvalidArtifact('optimizer unsuccessful: ' + str(opt.message))
        if not candidate['geometry_checks']['pass']:
            raise InvalidArtifact('final proposal geometry unsupported')
        r['status'] = 'proposal_available'
    except Exception as exc:
        r.update(reason=str(exc), traceback=traceback.format_exc())
    r.update(requests=ev.requests, accepted_iterations=accepted, wall_seconds=time.monotonic() - start,
             job_id=os.environ.get('SLURM_JOB_ID'), unconstrained_minimum_claimed=False)
    write_new(path, r); return record(path)


def propose(manifest):
    validate(manifest)
    if not os.environ.get('SLURM_JOB_ID') or int(os.environ['SLURM_CPUS_ON_NODE']) != 32 or int(os.environ['SLURM_MEM_PER_NODE']) != 200000:
        raise InvalidArtifact('declared GPU allocation required')
    root = Path(manifest).parent; lock = (root / 'proposal.lock').open('a+')
    fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    start = time.monotonic(); gpu = None; error = None; results = []
    try:
        gpu = WarmGPU(manifest)
        for t in read_json(manifest)['tasks']:
            pin = optimize(manifest, t, gpu); results.append(pin); r = read_json(verify(pin))
            print(json.dumps({'task': t['task_id'], 'status': r['status'], 'work': r.get('MACE_proposal_work_kcal_mol')}), flush=True)
    except Exception as exc:
        error = str(exc)
    finally:
        try:
            if gpu is not None: gpu.close()
        except Exception as exc:
            error = (error + '; ' if error else '') + str(exc)
        elapsed = time.monotonic() - start
        write_new(root / ('proposal_execution_' + os.environ['SLURM_JOB_ID'] + '.json'),
                  {'manifest': record(manifest), 'results': results, 'error': error, 'wall_seconds': elapsed,
                   'allocated_core_seconds': elapsed * 32, 'allocated_GPU_seconds': elapsed,
                   'job_id': os.environ['SLURM_JOB_ID'], 'GPU_command': gpu.command if gpu else None})
    if error: raise InvalidArtifact(error)
    return low_prepare(manifest)


def write_low_shards(dest, low, shard_count):
    """Partition existing explicit tasks; each executor owns its own contained paths."""
    from affordable_workflow import dry_run
    tasks, reused = low['all_tasks'], low['reused']
    shard_pins, task_manifests = [], {}
    for shard in range(shard_count):
        subset = tasks[shard::shard_count]
        selected = {t['task_id'] for t in subset}
        sm = {**low, 'all_tasks': subset, 'tasks': [t for t in subset if t['task_id'] not in reused],
              'reused': {k: v for k, v in reused.items() if k in selected}, 'shard': shard}
        sp = Path(dest) / ('shard_' + str(shard)) / 'manifest.json'; write_new(sp, sm)
        dry = dry_run(sp) if sm['tasks'] else {'status': 'no_new_tasks'}
        write_new(sp.parent / 'PREFLIGHT.json', dry); shard_pins.append(record(sp))
        task_manifests.update({tid: record(sp) for tid in selected})
    return {**low, 'tasks': [], 'shards': shard_pins, 'task_manifests': task_manifests,
            'execute_master_forbidden': True}


def low_completed(low, master, task_id):
    actual = verify(low['task_manifests'][task_id]) if 'task_manifests' in low else master
    return low['reused'].get(task_id) or completed(actual, task_id)


def low_prepare(manifest):
    from affordable_workflow import dry_run
    m = read_json(manifest); root = Path(manifest).parent
    dest = root / 'proposal_GFN2'; dest.mkdir(parents=True, exist_ok=False)
    tasks, unavailable, reused = [], [], {}
    shard_count = 4 if m.get('evaluation_scope') == 'primary225_fold_transfer_v1' else 1
    for t in m['tasks']:
        rp = root / 'proposals' / t['task_id'] / 'result.json'
        if not rp.exists():
            unavailable.append({'task_id': t['task_id'], 'reason': 'proposal missing'}); continue
        r = read_json(rp)
        if r['manifest'] != record(manifest): raise InvalidArtifact('proposal source differs')
        if r['status'] != 'proposal_available':
            unavailable.append({'task_id': t['task_id'], 'reason': r.get('reason', r['status'])}); continue
        p = r['proposal']
        for medium in ('vacuum', 'alpb'):
            tid = t['task_id'] + '__' + medium
            shard = len(tasks) % shard_count
            task_root = dest / ('shard_' + str(shard)) if shard_count > 1 else dest
            directory = task_root / 'tasks' / tid; directory.mkdir(parents=True)
            xp = directory / 'core.xyz'; shutil.copyfile(verify(p['coordinate']), xp)
            ip = directory / 'endpoint.inp'; ip.write_text(input_text(t['charge'], t['multiplicity'], medium, 'native'))
            task = {'task_id': tid, 'case': t['case_id'], 'case_id': t['case_id'], 'metal': t['metal'],
                    'proposal_task_id': t['task_id'], 'medium': medium, 'charge': t['charge'],
                    'multiplicity': t['multiplicity'], 'xyz': record(xp), 'input': record(ip),
                    'output_path': str(directory / 'endpoint.out'), 'proposal_receipt': record(rp)}
            tasks.append(task)
            if xyz(xp) == xyz(verify(t['xyz'])):
                reused[tid] = t['q0']['low'][medium]['receipt']
    low = {'protocol_id': PROTOCOL, 'stage': 'primary_native_GFN2_proposals', 'proposal_manifest': record(manifest),
           'agreement': m['agreement'], 'orca': m['orca'], 'implementation': m['implementation'],
           'execution_policy': {k: m['implementation'][v] for k, v in
                                [('task_runner', 'run_orca_task_manifest.py'), ('runtime_renderer', 'render_orca_runtime_input.py')]},
           'execution_resources': {'mpi_ranks': 8, 'concurrent_tasks': 8},
           'all_tasks': tasks, 'tasks': [t for t in tasks if t['task_id'] not in reused],
           'reused': reused, 'unavailable_proposals': unavailable, 'new_DFT_calls': 0}
    if len(tasks) > m['maximum_new_GFN2_singlepoints']: raise InvalidArtifact('proposal GFN2 scope exceeded')
    if shard_count > 1:
        low = write_low_shards(dest, low, shard_count)
    path = dest / 'manifest.json'; write_new(path, low)
    dry = ({'status': 'shards_dry_run_pass', 'shards': low['shards']} if shard_count > 1 else
           dry_run(path) if low['tasks'] else {'status': 'no_new_tasks'})
    write_new(dest / 'PREFLIGHT.json', dry)
    return {'manifest': record(path), 'new_GFN2_tasks': len(tasks) - len(reused), 'reused': len(reused),
            'unavailable_proposals': len(unavailable), 'preflight': dry}


def select_endpoint(origin, proposal):
    if origin is None or proposal is None:
        return {'status': 'unavailable', 'selected': None, 'work': None}
    work = relative_components(proposal, origin)
    selected = 'proposal' if work['composite_kcal_mol'] < -SETTINGS['selection_decrease_kcal_mol'] else 'origin'
    return {'status': 'available', 'selected': selected, 'work': work,
            'selected_components': proposal if selected == 'proposal' else origin}


def score(ca, la):
    """Use the released composite score arithmetic and conversion policy."""
    return mix_pair({'Ca': ca['MACE_eV'], 'La': la['MACE_eV']},
                    {'Ca': ca['GFN2_vacuum_hartree'], 'La': la['GFN2_vacuum_hartree']},
                    {'Ca': ca['GFN2_ALPB_hartree'], 'La': la['GFN2_ALPB_hartree']})


def evidence_counts(cases, stage):
    counts = Counter()
    for row in cases:
        value = row['old_band_' + stage + '_decision']; expected = row['expected_class']
        outcome = ('unavailable' if value is None else 'unlabeled' if expected is None else
                   'correct' if value == expected + '-supported' else
                   'wrong' if value in ('Ca-supported', 'La-supported') else 'inconclusive')
        counts[outcome] += 1
    return {'denominator': len(cases), **{k: counts[k] for k in ('correct', 'wrong', 'inconclusive', 'unavailable', 'unlabeled')}}


def collect(manifest, output):
    m = read_json(manifest); root = Path(manifest).parent; lp = root / 'proposal_GFN2/manifest.json'
    low = read_json(lp) if lp.exists() else None; endpoints = []; index = {}; low_results = []
    for t in m['tasks']:
        rp = root / 'proposals' / t['task_id'] / 'result.json'; r = read_json(rp) if rp.exists() else None
        row = {'task_id': t['task_id'], 'case_id': t['case_id'], 'metal': t['metal'],
               'status': 'unavailable', 'selection': None, 'origin_components': None, 'proposal_components': None,
               'proposal_receipt': record(rp) if rp.exists() else None,
               'reason': r.get('reason') if r else 'proposal not run', 'low': {},
               'boundary_flag': r.get('boundary_flag') if r else None,
               'boundary_distance_radian': r.get('boundary_distance_radian') if r else None,
               'optimizer': r.get('optimizer') if r else None,
               'origin_xyz': t['xyz'],
               'proposal_xyz': r['proposal']['coordinate'] if r and r.get('proposal') else None}
        if r and r['manifest'] != record(manifest): raise InvalidArtifact('proposal receipt manifest differs')
        if r and r.get('origin') and t['q0_status'] == 'available':
            row['origin_components'] = {**t['q0']['components'], 'MACE_eV': r['origin']['MACE_eV']}
        from accommodation_torsion_analysis import donor_distances
        row['donor_distances_A'] = {name: donor_distances({'mapping': t['mapping'], 'q': r[name]['full_q']})
                                  if r and r.get(name) else None for name in ('origin', 'proposal')}
        if r and r['status'] == 'proposal_available' and low:
            for medium in ('vacuum', 'alpb'):
                lt = next(x for x in low['all_tasks'] if x['proposal_task_id'] == t['task_id'] and x['medium'] == medium)
                rr = {'status': 'unavailable', 'task_id': lt['task_id'], 'energy_hartree': None,
                      'reused_q0': lt['task_id'] in low['reused']}
                try:
                    pin = low_completed(low, lp, lt['task_id'])
                    if pin is None: raise InvalidArtifact('proposal solver receipt unavailable or failed')
                    original_m = read_json(verify(pin['manifest']))
                    original_t = next(x for x in original_m['all_tasks'] if x['task_id'] == pin['task_id'])
                    if xyz(verify(original_t['xyz'])) != xyz(verify(lt['xyz'])) or verify(original_t['input']).read_text() != verify(lt['input']).read_text():
                        raise InvalidArtifact('proposal result input differs')
                    audit = diagnostics(pin, original_t)
                    if audit['charge_sanity_status'] != 'pass': raise InvalidArtifact('proposal charge audit failed')
                    rr.update(status='complete', **pin, audit=audit)
                except Exception as exc:
                    rr['reason'] = str(exc)
                    op = Path(lt['output_path']); receipt = Path(str(op) + '.execution.json')
                    rr['actual_artifacts'] = [record(p) for p in (op, receipt) if p.exists()]
                row['low'][medium] = rr; low_results.append(rr)
            if all(e['status'] == 'complete' for e in row['low'].values()):
                a, v = row['low']['alpb'], row['low']['vacuum']
                if read_json(verify(a['audit']['parameter_export'])) != read_json(verify(v['audit']['parameter_export'])):
                    row['reason'] = 'proposal medium parameter mismatch'
                else:
                    row['proposal_components'] = {'MACE_eV': r['proposal']['MACE_eV'],
                        'GFN2_ALPB_hartree': a['energy_hartree'], 'GFN2_vacuum_hartree': v['energy_hartree']}
                    row['selection'] = select_endpoint(row['origin_components'], row['proposal_components'])
                    row['status'] = row['selection']['status']; row['reason'] = None
        endpoints.append(row); index[t['case_id'], t['metal']] = row
    frozen = read_json(verify(m.get('fold_comparison', m['frozen_comparison']))); cases = []
    for c in m['cases']:
        pair = {z: index[c['case_id'], z] for z in ('Ca', 'La')}
        row = {'case_id': c['case_id'], 'label_scope': c['label_scope'],
               'expected_class': c.get('expected_class_for_later_report_only'),
               'status': 'unavailable', 'R0': None, 'R_proposal': None, 'R_selected': None,
               'old_band_origin_decision': None, 'old_band_selected_decision': None}
        for key in ('root_case_id', 'biological_group', 'source_conditioning_metal', 'canonical_coordinate_match'):
            if key in c: row[key] = c[key]
        for name, field in (('R0', 'origin_components'), ('R_proposal', 'proposal_components')):
            if all(e[field] is not None for e in pair.values()): row[name] = score(pair['Ca'][field], pair['La'][field])
        if all(e['status'] == 'available' for e in pair.values()):
            row.update(status='available', R_selected=score(pair['Ca']['selection']['selected_components'], pair['La']['selection']['selected_components']))
        for stage, field in (('origin', 'R0'), ('selected', 'R_selected')):
            row['old_band_' + stage + '_decision'] = decision(row[field]['composite_R_model_kcal_mol'] if row[field] else None, m['frozen_context_bands'])
        row['endpoint_selection'] = {z: pair[z]['selection'] for z in pair}
        old = next((r for r in frozen['rows'] if r['case_id'] == c['case_id'] and r['representation'] == 'context'), None)
        row['archived_origin_R_model_kcal_mol'] = old['composite_R_model_kcal_mol'] if old else None
        row['origin_replay_contrast_difference_model_kcal_mol'] = (row['R0']['composite_R_model_kcal_mol'] - old['composite_R_model_kcal_mol']) if row['R0'] and old else None
        row['delta_R_model_kcal_mol'] = row['R_selected']['composite_R_model_kcal_mol'] - row['R0']['composite_R_model_kcal_mol'] if row['R_selected'] and row['R0'] else None
        cases.append(row)
    for c in m.get('unavailable_cases', []):
        cases.append({**c, 'expected_class': c['expected_class_for_later_report_only'],
                      'status': 'unavailable', 'R0': None, 'R_proposal': None, 'R_selected': None,
                      'old_band_origin_decision': None, 'old_band_selected_decision': None,
                      'endpoint_selection': {}, 'delta_R_model_kcal_mol': None})
    evaluations = [read_json(p) for p in root.glob('proposals/*/evaluations/*/result.json')]
    natives = [read_json(verify(e['MACE'])) for e in evaluations if e.get('MACE')]
    groups = {'canonical': [r for r in cases if r['case_id'].endswith('-pqq-la_model')],
              'consumed_crystals': [r for r in cases if r['case_id'] in ('1H4I', '4MAE', '1KB0')],
              'unlabeled_PLM': [r for r in cases if r['case_id'] in PLM]}
    if m.get('evaluation_scope') == 'primary225_fold_transfer_v1':
        groups = {'primary225': cases, **{z + '_conditioned': [r for r in cases if r['source_conditioning_metal'] == z] for z in ('La', 'Ca')}}
    costs = [read_json(p) for p in root.glob('proposal_execution_*.json')]
    result = {'protocol_id': PROTOCOL, 'manifest': record(manifest), 'endpoints': endpoints, 'cases': cases,
              'endpoint_denominator': 2 * len(cases), 'case_denominator': len(cases),
              'prepared_endpoint_denominator': len(m['tasks']),
              'available_endpoints': sum(r['status'] == 'available' for r in endpoints),
              'available_cases': sum(r['status'] == 'available' for r in cases),
              'proposal_status_counts': dict(Counter(read_json(p)['status'] for p in root.glob('proposals/*/result.json'))),
              'MACE_calls_started': sum(r['model_call_started'] for r in natives),
              'MACE_calls_complete': sum(r['status'] == 'complete' for r in natives),
              'GFN2_proposal_complete': sum(r['status'] == 'complete' and not r['reused_q0'] for r in low_results),
              'GFN2_proposal_reused': sum(r['status'] == 'complete' and r['reused_q0'] for r in low_results),
              'counts': {g: {s: evidence_counts(rows, s) for s in ('origin', 'selected')} for g, rows in groups.items()},
              'selection_counts': dict(Counter(e['selection']['selected'] if e['selection'] else 'unavailable' for e in endpoints)),
              'frozen_context_bands': m['frozen_context_bands'], 'proposal_execution_receipts': costs,
              'new_DFT_calls': 0, 'threshold_refitted': False, 'baseline_changed': False,
              'interpretation': 'restricted two-geometry electronic descriptor; no unconstrained minimum or PLM accuracy claim'}
    write_new(output, result)
    report(result, Path(output).with_suffix('.md'))
    return {k: v for k, v in result.items() if k not in ('endpoints', 'cases')}


def report(result, output):
    lines = ['# MACE proposal / composite selection', '',
             f"Available: {result['available_cases']}/{result['case_denominator']} cases, {result['available_endpoints']}/{result['endpoint_denominator']} endpoints.", '',
             '| Stratum | Stage | Correct | Wrong | Inconclusive | Unavailable | Unlabeled | Total |',
             '|---|---|---:|---:|---:|---:|---:|---:|']
    for group, stages in result['counts'].items():
        for stage, c in stages.items():
            lines.append('| ' + group + ' | ' + stage + ' | ' + ' | '.join(str(c[k]) for k in ('correct', 'wrong', 'inconclusive', 'unavailable', 'unlabeled', 'denominator')) + ' |')
    lines += ['', '| Case | Original R | Selected R | Change | Selected old-band call |',
              '|---|---:|---:|---:|---|']
    for c in result['cases']:
        values = [c['R0']['composite_R_model_kcal_mol'] if c['R0'] else None,
                  c['R_selected']['composite_R_model_kcal_mol'] if c['R_selected'] else None, c['delta_R_model_kcal_mol']]
        lines.append('| ' + c['case_id'] + ' | ' + ' | '.join('unavailable' if v is None else f'{v:.6f}' for v in values) + ' | ' + str(c['old_band_selected_decision']) + ' |')
    lines += ['', 'R uses model kcal/mol. Both actually evaluated members are required for energy selection; failures are unavailable. Boundary proposals are flagged, never called unconstrained minima. Old bands are an exploratory transfer check, without refitting. Known references are consumed; PLM truth labels are unknown.', '',
              f"MACE calls started/complete: {result['MACE_calls_started']}/{result['MACE_calls_complete']}. New GFN2 proposal endpoints complete: {result['GFN2_proposal_complete']}; exact q0 proposal reuses: {result['GFN2_proposal_reused']}. New DFT: 0.", '']
    with Path(output).open('x') as handle: handle.write('\n'.join(lines))


def main():
    p = argparse.ArgumentParser(description=__doc__); sub = p.add_subparsers(dest='operation', required=True)
    q = sub.add_parser('prepare')
    for key in ('references', 'torsion', 'primary-result', 'compact', 'agreement', 'output', 'cpu-python', 'gpu-python'):
        q.add_argument('--' + key, required=True)
    for name in ('validate', 'propose', 'low-prepare', 'collect'):
        q = sub.add_parser(name); q.add_argument('--manifest', required=True)
        if name == 'collect': q.add_argument('--output', required=True)
    args = vars(p.parse_args()); op = args.pop('operation').replace('-', '_')
    print(json.dumps(globals()[op](**args), default=str))


if __name__ == '__main__': main()
