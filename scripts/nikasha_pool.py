"""Cross-score a common, source-specific pool using archived physical proposals."""
from __future__ import annotations
import argparse
from collections import Counter
from functools import lru_cache
import fcntl
import json
import os
from pathlib import Path
import shutil
import time
import numpy as np

from affordable_common import InvalidArtifact, energy, read_json, record, verify, write_new, xyz
from accommodation_nonlinear import WarmGPU, relative_components
from accommodation_proposals import PROTOCOL as SOURCE_PROTOCOL, PLM, score
from compact_solvation import input_text, completed, diagnostics
from mace_hybrid import check_atoms, write_xyz
from mace_site_kinematics import Kinematics

PROTOCOL = 'nikasha_common_geometry_native_OMOL_GFN2_ALPB_v1'
SETTINGS = {'coordinate_copy_tolerance_A': 1e-12, 'origin_selection_tolerance_kcal_mol': .1,
            'equivalent_cell_tolerance_kcal_mol': .1,
            'candidate_order': ['origin', 'proposal_Ca', 'proposal_La']}
PILOT = ('1H4I', '4MAE', *PLM)


@lru_cache(maxsize=None)
def checked_json(path, digest):
    return read_json(verify({'path': path, 'sha256': digest}))


def pinned(pin):
    return checked_json(pin['path'], pin['sha256'])


def same_geometry(a, b):
    if len(a) != len(b) or [r[0] for r in a[1:]] != [r[0] for r in b[1:]]:
        return False
    if a[0][0] not in ('Ca', 'La') or b[0][0] not in ('Ca', 'La'):
        return False
    return bool(np.max(np.abs(np.array([r[1:] for r in a]) -
                              np.array([r[1:] for r in b]))) <= SETTINGS['coordinate_copy_tolerance_A'])


def reused_cell(task, endpoint, proposal, stage, parent):
    """Check real immutable outputs once, reusing the existing protocol receipts."""
    position = proposal[stage]; atoms = xyz(verify(position['coordinate']))
    target = endpoint['origin_xyz' if stage == 'origin' else 'proposal_xyz']
    if atoms != xyz(verify(target)):
        raise InvalidArtifact('archived candidate coordinates differ')
    native = pinned(position['MACE'])
    request = pinned(native['request'])
    if (native['status'] != 'complete' or native['model'] != parent['model'] or
            request['charge'] != task['charge'] or request['multiplicity'] != task['multiplicity'] or
            xyz(verify(request['xyz'])) != atoms):
        raise InvalidArtifact('archived MACE state or coordinates differ')
    comp = endpoint['origin_components' if stage == 'origin' else 'proposal_components']
    if comp is None or comp['MACE_eV'] != native['energy_eV']:
        raise InvalidArtifact('archived MACE energy unavailable or differs')
    lows = {}; parameters = []
    for medium, field in (('vacuum', 'GFN2_vacuum_hartree'), ('alpb', 'GFN2_ALPB_hartree')):
        if stage == 'origin':
            source = task['q0']['low'][medium]; receipt = source['receipt']; audit = source['audit']
        else:
            source = endpoint['low'][medium]; receipt = source; audit = source['audit']
            if source['status'] != 'complete': raise InvalidArtifact('archived proposal solvent cell missing')
        lm = pinned(receipt['manifest'])
        lt = next(t for t in lm['all_tasks'] if t['task_id'] == receipt['task_id'])
        execution = pinned(receipt['receipt'])
        if (lm['orca'] != parent['orca'] or execution['returncode'] != 0 or
                not execution['normal_termination'] or not execution['scf_converged'] or
                execution['manifest'] != receipt['manifest'] or execution['task_id'] != lt['task_id'] or
                lt['charge'] != task['charge'] or lt['multiplicity'] != task['multiplicity'] or
                xyz(verify(lt['xyz'])) != atoms or
                verify(lt['input']).read_text() != input_text(task['charge'], task['multiplicity'], medium, 'native') or
                audit['charge_sanity_status'] != 'pass'):
            raise InvalidArtifact('archived primary solvent recipe/state/receipt differs')
        if energy(verify(receipt['output'])) != comp[field]:
            raise InvalidArtifact('archived solvent energy differs from raw output')
        parameters.append(pinned(audit['parameter_export']))
        lows[medium] = {k: receipt[k] for k in ('manifest', 'task_id', 'output', 'receipt', 'energy_hartree')}
    if parameters[0] != parameters[1]: raise InvalidArtifact('archived medium parameter mismatch')
    return {'status': 'complete', 'components': comp, 'xyz': target, 'MACE': position['MACE'],
            'low': lows, 'reused': True, 'stage': stage}


def choose_rows(matrix, candidates):
    """Common-pool row minima and the separate original numerical-origin policy."""
    if any(matrix[z].get(q, {}).get('status') != 'complete' for z in ('Ca', 'La') for q in candidates):
        return {'status': 'unavailable', 'reason': 'required_matrix_cell_unavailable',
                'mathematical': None, 'operational': None, 'rows': None}
    rows = {}; chosen = {'mathematical': {}, 'operational': {}}
    for z in ('Ca', 'La'):
        origin = matrix[z]['origin']['components']
        works = {q: relative_components(matrix[z][q]['components'], origin) for q in candidates}
        minimum = min(candidates, key=lambda q: works[q]['composite_kcal_mol'])
        operational = minimum if works[minimum]['composite_kcal_mol'] < -SETTINGS['origin_selection_tolerance_kcal_mol'] else 'origin'
        rows[z] = {'work_from_origin_kcal_mol': works, 'mathematical_candidate': minimum,
                   'operational_candidate': operational, 'candidate_count_is_not_population': True}
        for mode, candidate in (('mathematical', minimum), ('operational', operational)):
            chosen[mode][z] = matrix[z][candidate]['components']
    return {'status': 'available', 'rows': rows, **{mode: score(v['Ca'], v['La']) for mode, v in chosen.items()}}


def prepare(source, population, agreement, output, shards):
    begin = time.monotonic(); collection = read_json(source); parent = pinned(collection['manifest'])
    if parent['protocol_id'] != SOURCE_PROTOCOL:
        raise InvalidArtifact('unsupported proposal source protocol')
    if population not in ('pilot4', 'remaining26', 'primary225') or shards not in (1, 4):
        raise InvalidArtifact('undeclared source/shard scope')
    old = {c['case_id']: c for c in collection['cases']}
    if len(old) != collection['case_denominator']: raise InvalidArtifact('duplicate source IDs')
    if population == 'primary225':
        if parent.get('evaluation_scope') != 'primary225_fold_transfer_v1' or len(old) != 225:
            raise InvalidArtifact('not the declared primary225 source population')
        ids = sorted(old)
    else:
        if len(old) != 30 or not set(PILOT) <= old.keys(): raise InvalidArtifact('original30 scope differs')
        ids = list(PILOT) if population == 'pilot4' else sorted(set(old) - set(PILOT))
    out = Path(output).resolve(); out.mkdir(parents=True, exist_ok=False)
    endpoints = {(e['case_id'], e['metal']): e for e in collection['endpoints']}
    original_tasks = {(t['case_id'], t['metal']): t for t in parent['tasks']}
    cases, tasks = [], []
    for cid in ids:
        c = {'case_id': cid, 'status': 'unavailable', 'reason': None, 'candidates': [],
             'matrix': {'Ca': {}, 'La': {}}, 'aliases': {}, 'old_result': old[cid]}
        pending = []
        try:
            if old[cid]['status'] != 'available':
                raise InvalidArtifact('required archived proposal/origin cell unavailable')
            ts = {z: original_tasks[cid, z] for z in ('Ca', 'La')}
            eps = {z: endpoints[cid, z] for z in ts}
            proposals = {z: pinned(eps[z]['proposal_receipt']) for z in ts}
            origins = {z: xyz(verify(ts[z]['xyz'])) for z in ts}
            if (origins['Ca'][1:] != origins['La'][1:] or origins['Ca'][0][1:] != origins['La'][0][1:] or
                    ts['La']['charge'] - ts['Ca']['charge'] != 1 or
                    ts['Ca']['multiplicity'] != ts['La']['multiplicity'] or
                    ts['Ca']['source_preparation'] != ts['La']['source_preparation'] or
                    pinned(ts['Ca']['mapping']) != pinned(ts['La']['mapping'])):
                raise InvalidArtifact('incompatible chemical state or unequal origins; wet-origin expansion unsupported')
            for z in ts:
                verify(ts[z]['source_preparation'])
                if proposals[z]['manifest'] != collection['manifest'] or proposals[z]['status'] != 'proposal_available':
                    raise InvalidArtifact('proposal source differs or failed')
                k = Kinematics(pinned(ts[z]['mapping'])['context'])
                coords = k.evaluate(proposals[z]['proposal']['full_q'])[1]
                atoms = xyz(verify(proposals[z]['proposal']['coordinate']))
                if not np.allclose(coords, [a[1:] for a in atoms], atol=1e-12, rtol=0):
                    raise InvalidArtifact('proposal physical cap/atom mapping differs')
                if [a[0] for a in atoms] != [a[0] for a in origins[z]]:
                    raise InvalidArtifact('proposal atom inventory changed')
            raw = [('origin', ts['Ca']['xyz']), *[('proposal_' + z, eps[z]['proposal_xyz']) for z in ('Ca', 'La')]]
            candidate_atoms = {}
            for name, pin in raw:
                atoms = xyz(verify(pin)); match = next((q for q, a in candidate_atoms.items() if same_geometry(atoms, a)), None)
                if match is None:
                    candidate_atoms[name] = atoms; c['candidates'].append({'id': name, 'xyz': pin})
                    match = name
                c['aliases'][name] = {'representative': match, 'source_xyz': pin}
            for z in ts:
                for stage, name in (('origin', 'origin'), ('proposal', 'proposal_' + z)):
                    cell = reused_cell(ts[z], eps[z], proposals[z], stage, parent)
                    q = c['aliases'][name]['representative']
                    if q in c['matrix'][z]:
                        previous = c['matrix'][z][q]
                        delta = relative_components(cell['components'], previous['components'])['composite_kcal_mol']
                        if abs(delta) > SETTINGS['equivalent_cell_tolerance_kcal_mol']:
                            raise InvalidArtifact('numerical-copy reused energies inconsistent')
                        previous.setdefault('equivalent_reused_cells', []).append({'cell': cell, 'delta_kcal_mol': delta})
                    else: c['matrix'][z][q] = cell
                for candidate in c['candidates']:
                    q = candidate['id']
                    if q in c['matrix'][z]: continue
                    atoms = candidate_atoms[q]; paired = [(z, *atoms[0][1:]), *atoms[1:]]
                    check_atoms(paired, ts[z]['charge'])
                    tid = cid + '__' + z + '__at_' + q
                    td = out / 'cells' / tid; td.mkdir(parents=True)
                    xp = td / 'context.xyz'; write_xyz(xp, paired)
                    if xyz(xp) != paired: raise InvalidArtifact('cross-cell XYZ serialization changed coordinates')
                    task = {'task_id': tid, 'case_id': cid, 'metal': z, 'candidate': q,
                            'xyz': record(xp), 'charge': ts[z]['charge'], 'multiplicity': ts[z]['multiplicity'],
                            'source_mapping': ts[z]['mapping'], 'source_preparation': ts[z]['source_preparation']}
                    pending.append(task)
                    c['matrix'][z][q] = {'status': 'pending', 'task_id': tid, 'xyz': task['xyz'], 'reused': False}
            c['status'] = 'prepared'; tasks.extend(pending)
        except (InvalidArtifact, KeyError, OSError, StopIteration) as exc:
            c['reason'] = str(exc)
        cases.append(c)
    impl = out / 'implementation'; impl.mkdir(); pins = {}
    for p in Path(__file__).parent.glob('*.py'):
        target = impl / p.name; shutil.copyfile(p, target); pins[p.name] = record(target)
    m = {'protocol_id': PROTOCOL, 'settings': SETTINGS, 'source': record(source),
         'source_manifest': collection['manifest'], 'agreement': record(agreement),
         'population': population, 'declared_case_ids': ids, 'cases': cases, 'tasks': tasks,
         'model': parent['model'], 'software': parent['software'], 'orca': parent['orca'],
         'cpu_python': parent['cpu_python'], 'gpu_python': parent['gpu_python'],
         'implementation': pins, 'shard_count': shards,
         'resources': {'GPU_cpus': 32, 'GPUs': 1, 'host_mem_MiB': 200000,
                       'GFN2_cpus_per_shard': 64, 'mpi_ranks': 8, 'concurrent_tasks': 8},
         'new_MACE_cells': len(tasks), 'maximum_new_GFN2_calls': 2 * len(tasks), 'new_DFT_calls': 0,
         'source_preparation_wall_seconds': time.monotonic() - begin,
         'new_optimizations': 0, 'production_changed': False}
    path = out / 'manifest.json'; write_new(path, m)
    low_prepare(path)
    result = validate(path); write_new(out / 'PREFLIGHT.json', result); return result


def validate(manifest):
    m = read_json(manifest)
    if m['protocol_id'] != PROTOCOL or m['settings'] != SETTINGS: raise InvalidArtifact('pool protocol changed')
    for key in ('source', 'source_manifest', 'agreement', 'software', 'orca'): verify(m[key])
    for pin in m['implementation'].values(): verify(pin)
    n = {'pilot4': 4, 'remaining26': 26, 'primary225': 225}[m['population']]
    if (len(m['cases']) != n or len(set(m['declared_case_ids'])) != n or
            {c['case_id'] for c in m['cases']} != set(m['declared_case_ids'])):
        raise InvalidArtifact('pool denominator differs')
    if len(m['tasks']) > 2 * n or m['maximum_new_GFN2_calls'] != 2 * len(m['tasks']):
        raise InvalidArtifact('unexpected additional cells')
    for t in m['tasks']:
        atoms = xyz(verify(t['xyz'])); check_atoms(atoms, t['charge'])
        if atoms[0][0] != t['metal'] or t['multiplicity'] != 1: raise InvalidArtifact('unsupported cross state')
    return {'status': 'validated', 'manifest': record(manifest), 'denominator': n,
            'prepared': sum(c['status'] == 'prepared' for c in m['cases']),
            'new_MACE_cells': len(m['tasks']), 'new_GFN2_calls': 2 * len(m['tasks']), 'new_DFT_calls': 0}


def low_prepare(manifest):
    from affordable_workflow import dry_run
    m = read_json(manifest); root = Path(manifest).parent / 'solvent'; shards = []
    for shard in range(m['shard_count']):
        sd = root / ('shard_' + str(shard)); tasks = []
        for cell in m['tasks'][shard::m['shard_count']]:
            for medium in ('vacuum', 'alpb'):
                tid = cell['task_id'] + '__' + medium; td = sd / 'tasks' / tid; td.mkdir(parents=True)
                xp = td / 'core.xyz'; shutil.copyfile(verify(cell['xyz']), xp)
                ip = td / 'endpoint.inp'; ip.write_text(input_text(cell['charge'], cell['multiplicity'], medium, 'native'))
                tasks.append({**cell, 'task_id': tid, 'cell_id': cell['task_id'], 'case': cell['case_id'],
                              'medium': medium, 'xyz': record(xp), 'input': record(ip), 'output_path': str(td / 'endpoint.out')})
        low = {'protocol_id': PROTOCOL, 'pool_manifest': record(manifest), 'agreement': m['agreement'],
               'orca': m['orca'], 'tasks': tasks, 'all_tasks': tasks, 'shard': shard,
               'execution_policy': {'task_runner': m['implementation']['run_orca_task_manifest.py'],
                                    'runtime_renderer': m['implementation']['render_orca_runtime_input.py']},
               'execution_resources': {'mpi_ranks': 8, 'concurrent_tasks': 8}}
        lp = sd / 'manifest.json'; write_new(lp, low)
        check = dry_run(lp) if tasks else {'status': 'no_new_tasks'}
        write_new(sd / 'PREFLIGHT.json', check); shards.append(record(lp))
    write_new(root / 'INDEX.json', {'pool_manifest': record(manifest), 'shards': shards})


def execute_mace(manifest):
    validate(manifest); m = read_json(manifest); root = Path(manifest).parent
    if (not os.environ.get('SLURM_JOB_ID') or int(os.environ['SLURM_CPUS_ON_NODE']) != 32 or
            int(os.environ['SLURM_MEM_PER_NODE']) != 200000): raise InvalidArtifact('declared GPU allocation required')
    lock = (root / 'execute.lock').open('a+'); fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    start = time.monotonic(); worker = None; receipts = []; error = None
    try:
        worker = WarmGPU(manifest)
        for t in m['tasks']:
            td = root / 'cells' / t['task_id']; request = td / 'request.json'; result = td / 'mace_result.json'
            expected = {k: t[k] for k in ('task_id', 'xyz', 'charge', 'multiplicity')} | {'manifest': record(manifest)}
            if request.exists():
                if read_json(request) != expected: raise InvalidArtifact('partial request cache mismatch')
            else:
                write_new(request, expected)
            if result.exists():
                existing = read_json(result)
                if pinned(existing['request']) != expected or existing['model'] != m['model']:
                    raise InvalidArtifact('cell cache mismatch')
                receipts.append(record(result)); continue
            evaluated, pin = worker.evaluate(record(request)); receipts.append(pin)
            print(json.dumps({'task_id': t['task_id'], 'status': evaluated['status']}), flush=True)
    except Exception as exc: error = str(exc)
    finally:
        try:
            if worker is not None: worker.close()
        except Exception as exc:
            error = (error + '; ' if error else '') + str(exc)
        elapsed = time.monotonic() - start
        write_new(root / ('MACE_EXECUTION_' + os.environ['SLURM_JOB_ID'] + '.json'),
                  {'manifest': record(manifest), 'receipts': receipts, 'error': error,
                   'wall_seconds': elapsed, 'allocated_core_seconds': 32 * elapsed,
                   'allocated_GPU_seconds': elapsed, 'job_id': os.environ['SLURM_JOB_ID']})
    if error: raise InvalidArtifact(error)
    return {'cells': len(receipts), 'wall_seconds': elapsed}


def collect(manifest, output):
    m = read_json(manifest); root = Path(manifest).parent; new = {}
    for pin in read_json(root / 'solvent/INDEX.json')['shards']:
        lp = verify(pin); lm = read_json(lp)
        for t in lm['tasks']:
            row = {'status': 'unavailable'}
            try:
                result = completed(lp, t['task_id'])
                if result:
                    audit = diagnostics(result, t)
                    if audit['charge_sanity_status'] != 'pass': raise InvalidArtifact('cross-cell charge audit failed')
                    row = {'status': 'complete', **result, 'audit': audit}
            except Exception as exc:
                row['reason'] = str(exc)
            if row['status'] != 'complete':
                row['artifacts'] = [record(p) for p in (Path(t['output_path']), Path(t['output_path'] + '.execution.json')) if p.exists()]
            new.setdefault(t['cell_id'], {})[t['medium']] = row
    cases = []; tasks = {t['task_id']: t for t in m['tasks']}
    for c in m['cases']:
        matrix = json.loads(json.dumps(c['matrix']))
        for z, cells in matrix.items():
            for q, cell in cells.items():
                if cell['reused']: continue
                rp = root / 'cells' / cell['task_id'] / 'mace_result.json'
                native = read_json(rp) if rp.exists() else None; lows = new.get(cell['task_id'], {})
                cell.update(status='unavailable', low=lows, MACE=record(rp) if native else None)
                if native:
                    try:
                        request = pinned(native['request']); task = tasks[cell['task_id']]
                        if (request['manifest'] != record(manifest) or request['task_id'] != task['task_id'] or
                                request['xyz'] != task['xyz'] or request['charge'] != task['charge'] or
                                request['multiplicity'] != task['multiplicity'] or native['model'] != m['model']):
                            raise InvalidArtifact('fresh cross-cell MACE state/geometry/method differs')
                        verify(request['xyz'])
                    except Exception as exc:
                        cell['reason'] = str(exc); continue
                if native and native['status'] == 'complete' and len(lows) == 2 and all(r['status'] == 'complete' for r in lows.values()):
                    params = [pinned(lows[k]['audit']['parameter_export']) for k in ('vacuum', 'alpb')]
                    if params[0] != params[1]: cell['reason'] = 'medium parameter mismatch'; continue
                    cell.update(status='complete', components={'MACE_eV': native['energy_eV'],
                        'GFN2_vacuum_hartree': lows['vacuum']['energy_hartree'],
                        'GFN2_ALPB_hartree': lows['alpb']['energy_hartree']})
        pool = choose_rows(matrix, [q['id'] for q in c['candidates']]) if c['status'] == 'prepared' else {
            'status': 'unavailable', 'reason': c['reason'], 'mathematical': None, 'operational': None, 'rows': None}
        cases.append({**c, 'matrix': matrix, 'pool': pool})
    native_receipts = [read_json(p) for p in root.glob('cells/*/mace_result.json')]
    result = {'protocol_id': PROTOCOL, 'manifest': record(manifest), 'cases': cases,
              'denominator': len(cases), 'available': sum(c['pool']['status'] == 'available' for c in cases),
              'required_new_MACE_cells': len(m['tasks']), 'required_new_GFN2_calls': 2 * len(m['tasks']),
              'MACE_calls_started': sum(r['model_call_started'] for r in native_receipts),
              'MACE_calls_complete': sum(r['status'] == 'complete' for r in native_receipts),
              'GFN2_complete': sum(r['status'] == 'complete' for v in new.values() for r in v.values()),
              'new_DFT_calls': 0, 'new_optimizations': 0, 'production_changed': False,
              'interpretation': 'within-source finite electronic geometry pool, no equilibrium population or affinity claim'}
    write_new(output, result)
    return {k: v for k, v in result.items() if k != 'cases'}


def main():
    parser = argparse.ArgumentParser(description=__doc__); sub = parser.add_subparsers(dest='operation', required=True)
    p = sub.add_parser('prepare')
    for key in ('source', 'population', 'agreement', 'output'): p.add_argument('--' + key, required=True)
    p.add_argument('--shards', type=int, default=1)
    for name in ('validate', 'execute-mace', 'collect'):
        p = sub.add_parser(name); p.add_argument('--manifest', required=True)
        if name == 'collect': p.add_argument('--output', required=True)
    a = vars(parser.parse_args()); op = a.pop('operation').replace('-', '_')
    print(json.dumps(globals()[op](**a)))


if __name__ == '__main__': main()
