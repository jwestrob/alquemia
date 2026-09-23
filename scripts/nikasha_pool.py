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
ADAPTIVE_PROTOCOL = 'nikasha_adaptive_angular_common_geometry_native_OMOL_GFN2_ALPB_v1'
JOINT_PROTOCOL = 'nikasha_joint_metal_angular_common_geometry_native_OMOL_GFN2_ALPB_v1'
SCALED_PROTOCOL = 'nikasha_scaled_angular_common_geometry_native_OMOL_GFN2_ALPB_v1'
SETTINGS = {'coordinate_copy_tolerance_A': 1e-12, 'origin_selection_tolerance_kcal_mol': .1,
            'equivalent_cell_tolerance_kcal_mol': .1,
            'candidate_order': ['origin', 'proposal_Ca', 'proposal_La']}
PILOT = ('1H4I', '4MAE', *PLM)
ADAPTIVE_SETTINGS = {**SETTINGS, 'candidate_order': SETTINGS['candidate_order'] + ['adaptive_Ca', 'adaptive_La']}
JOINT_SETTINGS = {**SETTINGS, 'candidate_order': ADAPTIVE_SETTINGS['candidate_order'] + ['joint_Ca', 'joint_La']}


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
    if m['protocol_id'] == 'nikasha_union_adaptive_minimal_common_geometry_native_OMOL_GFN2_ALPB_v1':
        from union_adaptive import validate_pool
        return validate_pool(manifest)
    if m['protocol_id'] in ('nikasha_declared_finite_geometry_native_OMOL_GFN2_ALPB_v1',
                            'nikasha_collective_scaffold_geometry_native_OMOL_GFN2_ALPB_v1'):
        from nikasha_finite_candidates import validate as validate_finite
        return validate_finite(manifest)
    expected = {PROTOCOL: SETTINGS, ADAPTIVE_PROTOCOL: ADAPTIVE_SETTINGS, JOINT_PROTOCOL: JOINT_SETTINGS,
                SCALED_PROTOCOL: ADAPTIVE_SETTINGS}
    if m['protocol_id'] not in expected or m['settings'] != expected[m['protocol_id']]:
        raise InvalidArtifact('pool protocol changed')
    for key in ('source', 'source_manifest', 'agreement', 'software', 'orca'): verify(m[key])
    for pin in m['implementation'].values(): verify(pin)
    n = {'pilot4': 4, 'remaining26': 26, 'primary225': 225, 'adaptive4': 4, 'adaptive26': 26,
         'joint4': 4, 'scaled30': 30, 'scaled225': 225}[m['population']]
    allowed = {PROTOCOL: ('pilot4', 'remaining26', 'primary225'),
               ADAPTIVE_PROTOCOL: ('adaptive4', 'adaptive26'), JOINT_PROTOCOL: ('joint4',),
               SCALED_PROTOCOL: ('scaled30', 'scaled225')}
    if m['population'] not in allowed[m['protocol_id']]: raise InvalidArtifact('pool population/protocol mismatch')
    adaptive = m['protocol_id'] in (ADAPTIVE_PROTOCOL, JOINT_PROTOCOL, SCALED_PROTOCOL)
    if adaptive:
        for pin in m.get('base_pools', [m.get('base_pool')]): verify(pin)
        verify(m['adaptive_proposals'])
    if m['protocol_id'] == SCALED_PROTOCOL:
        source = pinned(pinned(m['adaptive_proposals'])['manifest'])
        if m['proposal_protocol_id'] != source['protocol_id'] or m['proposal_settings'] != source['settings']:
            raise InvalidArtifact('scaled proposal policy differs from its actual source')
    if m.get('GFN2_maxiter', 125) not in (125, 500): raise InvalidArtifact('unsupported iteration ceiling')
    if m.get('GFN2_maxiter') == 500:
        qualification = pinned(m['numerical_qualification'])
        if not qualification['both_controls_pass'] or not qualification['formerly_failed_cell_complete']:
            raise InvalidArtifact('iteration-ceiling qualification failed')
    if (len(m['cases']) != n or len(set(m['declared_case_ids'])) != n or
            {c['case_id'] for c in m['cases']} != set(m['declared_case_ids'])):
        raise InvalidArtifact('pool denominator differs')
    if len(m['tasks']) > (4 if adaptive else 2) * n or m['maximum_new_GFN2_calls'] != 2 * len(m['tasks']):
        raise InvalidArtifact('unexpected additional cells')
    for t in m['tasks']:
        atoms = xyz(verify(t['xyz'])); check_atoms(atoms, t['charge'])
        if atoms[0][0] != t['metal'] or t['multiplicity'] != 1: raise InvalidArtifact('unsupported cross state')
        if t.get('native_reuse'): native_reuse(t, m)
    return {'status': 'validated', 'manifest': record(manifest), 'denominator': n,
            'prepared': sum(c['status'] == 'prepared' for c in m['cases']),
            'new_MACE_cells': m['new_MACE_cells'], 'new_GFN2_calls': 2 * len(m['tasks']), 'new_DFT_calls': 0}


def native_reuse(task, manifest):
    native = pinned(task['native_reuse']); request = pinned(native['request'])
    if (native['status'] != 'complete' or native['model'] != manifest['model'] or
            request['charge'] != task['charge'] or request['multiplicity'] != task['multiplicity'] or
            xyz(verify(request['xyz'])) != xyz(verify(task['xyz']))):
        raise InvalidArtifact('reused native candidate state/geometry/model differs')
    return native


def expansion_bases(paths):
    """Join actual disjoint original30 pool collections without inventing a result."""
    paths = list(paths) if isinstance(paths, (list, tuple)) else [paths]
    if not paths: raise InvalidArtifact('missing base collection')
    data = [read_json(p) for p in paths]; manifests = [pinned(d['manifest']) for d in data]
    origins = {}; cases = []
    for path, collection, manifest in zip(paths, data, manifests):
        if (collection['protocol_id'] != manifest['protocol_id'] or
                len(collection['cases']) != collection['denominator'] or
                {c['case_id'] for c in collection['cases']} != set(manifest['declared_case_ids'])):
            raise InvalidArtifact('base collection denominator or protocol differs')
        for c in collection['cases']:
            if c['case_id'] in origins: raise InvalidArtifact('duplicate source in base collections')
            origins[c['case_id']] = record(path); cases.append(c)
    if len(paths) == 1: return data[0], manifests[0], origins
    if (len(paths) != 2 or {m['population'] for m in manifests} != {'pilot4', 'remaining26'} or
            any(m['protocol_id'] != PROTOCOL or m['settings'] != SETTINGS for m in manifests)):
        raise InvalidArtifact('only the declared original30 base split can be joined')
    first = manifests[0]
    keys = ('source', 'source_manifest', 'model', 'software', 'orca', 'cpu_python', 'gpu_python', 'resources')
    if any(any(m[k] != first[k] for k in keys) for m in manifests):
        raise InvalidArtifact('base collections have different physical sources or methods')
    source = pinned(first['source_manifest'])
    ids = [c['case_id'] for c in source['cases']]
    if len(ids) != 30 or set(ids) != set(origins): raise InvalidArtifact('original30 source union differs')
    lookup = {c['case_id']: c for c in cases}
    return ({'protocol_id': PROTOCOL, 'cases': [lookup[cid] for cid in ids]},
            {**first, 'population': 'original30', 'declared_case_ids': ids}, origins)


def prepare_expansion(base_pool, proposals, agreement, output, shards=1, maxiter=125, numerical_qualification=None, proposal_family='angular'):
    """Admit actual mapped adaptive candidates into the existing common pool."""
    if proposal_family == 'angular':
        from adaptive_angular_proposals import PROTOCOL as proposal_protocol, final_geometry
        base_protocol, populations, protocol, settings, prefix = PROTOCOL, ('pilot4', 'remaining26'), ADAPTIVE_PROTOCOL, ADAPTIVE_SETTINGS, 'adaptive'
    elif proposal_family == 'joint':
        from adaptive_metal_proposals import PROTOCOL as proposal_protocol, final_geometry
        base_protocol, populations, protocol, settings, prefix = ADAPTIVE_PROTOCOL, ('adaptive4',), JOINT_PROTOCOL, JOINT_SETTINGS, 'joint'
    elif proposal_family == 'scaled-angular':
        from adaptive_completion import PROTOCOL as proposal_protocol
        from adaptive_angular_proposals import final_geometry
        base_protocol, populations, protocol, settings, prefix = PROTOCOL, ('original30', 'primary225'), SCALED_PROTOCOL, ADAPTIVE_SETTINGS, 'adaptive'
    else:
        raise InvalidArtifact('unsupported physical proposal family')
    started = time.monotonic(); base, bm, base_origins = expansion_bases(base_pool)
    base_paths = list(base_pool) if isinstance(base_pool, (list, tuple)) else [base_pool]
    adaptive = read_json(proposals); am = pinned(adaptive['manifest'])
    if (base['protocol_id'] != base_protocol or bm['population'] not in populations or
            adaptive['protocol_id'] != proposal_protocol or am['protocol_id'] != proposal_protocol or
            am['source_manifest'] != bm['source_manifest'] or
            any(am[k] != bm[k] for k in ('model', 'software', 'orca', 'cpu_python', 'gpu_python')) or
            {c['case_id'] for c in base['cases']} != set(bm['declared_case_ids']) or
            {c['case_id'] for c in adaptive['cases']} != set(bm['declared_case_ids']) or shards not in (1, 4)):
        raise InvalidArtifact('adaptive/base pool population or physical method differs')
    n = {'pilot4': 4, 'remaining26': 26, 'adaptive4': 4, 'original30': 30, 'primary225': 225}[bm['population']]
    if len(base['cases']) != n or maxiter not in (125, 500): raise InvalidArtifact('expansion scope differs')
    qualification = None
    if maxiter == 500:
        if numerical_qualification is None: raise InvalidArtifact('iteration ceiling needs actual qualification')
        qualification = read_json(numerical_qualification); qm = pinned(qualification['manifest'])
        if (qualification['protocol_id'] != 'adaptive_pool_native_GFN2_MaxIter500_diagnostic_v1' or
                not qualification['both_controls_pass'] or not qualification['formerly_failed_cell_complete'] or
                qm['orca'] != bm['orca'] or qm['control_tolerance_kcal_mol'] != .05):
            raise InvalidArtifact('iteration ceiling qualification incompatible')
    if proposal_family != 'scaled-angular' and any(c['pool']['status'] != 'available' for c in base['cases']):
        raise InvalidArtifact('completed original common pool required')
    tasks_by_id = {(t['case_id'], t['metal']): t for t in am['tasks']}
    endpoints = {(e['case_id'], e['metal']): e for e in adaptive['endpoints']}
    if len(endpoints) != 2*n: raise InvalidArtifact('adaptive endpoint denominator differs')
    out = Path(output).resolve(); out.mkdir(parents=True, exist_ok=False)
    cases, tasks = [], []
    for prior in base['cases']:
        cid = prior['case_id']; names = [q['id'] for q in prior['candidates']]
        if prior['status'] == 'prepared' and choose_rows(prior['matrix'], names) != prior['pool']:
            raise InvalidArtifact('original pool algebra changed')
        c = json.loads(json.dumps(prior)); c.pop('pool')
        c.update(status='unavailable', reason=None, prior_pool=prior['pool'], base_collection=base_origins[cid])
        for cells in c['matrix'].values():
            for cell in cells.values(): cell['reused'] = True
        pending = []
        try:
            if prior['pool']['status'] != 'available':
                raise InvalidArtifact('inherited base pool unavailable: ' + str(prior.get('reason') or prior['pool'].get('reason')))
            atoms_by_name = {q['id']: xyz(verify(q['xyz'])) for q in c['candidates']}
            points = {}
            for z in ('Ca', 'La'):
                t = tasks_by_id[cid, z]; e = endpoints[cid, z]
                if e['status'] != 'candidate_available': raise InvalidArtifact('adaptive candidate unavailable: ' + z)
                receipt = pinned(e['proposal_receipt']); point = e['candidate']
                if (receipt['status'] != 'proposal_available' or receipt['manifest'] != adaptive['manifest'] or
                        point != receipt['proposal'] or e['origin_components'] != c['matrix'][z]['origin']['components'] or
                        xyz(verify(t['xyz'])) != xyz(verify(c['matrix'][z]['origin']['xyz'])) or
                        e['mapping'] != t['mapping'] or e['source_preparation'] != t['source_preparation']):
                    raise InvalidArtifact('adaptive origin, state, mapping or receipt changed')
                kin = Kinematics(pinned(t['mapping'])['context']); q = np.asarray(point['full_q'])
                if any(q[i] != 0 for i in set(range(len(q))) - set(t['active_indices'])):
                    raise InvalidArtifact('unselected adaptive coordinate changed')
                atoms = xyz(verify(point['coordinate']))
                final_geometry(kin, t, q[t['active_indices']], [a[0] for a in atoms])
                if (not np.allclose(kin.evaluate(q)[1], [a[1:] for a in atoms], atol=1e-12, rtol=0) or
                        [a[0] for a in atoms] != [a[0] for a in xyz(verify(t['xyz']))]):
                    raise InvalidArtifact('adaptive coordinate or inventory differs from physical source')
                points[z] = point
                name = prefix + '_' + z
                match = next((k for k, a in atoms_by_name.items() if same_geometry(a, atoms)), None)
                if match is None:
                    match = name; atoms_by_name[name] = atoms
                    c['candidates'].append({'id': name, 'xyz': point['coordinate']})
                c['aliases'][name] = {'representative': match, 'source_xyz': point['coordinate']}
            for z in ('Ca', 'La'):
                t = tasks_by_id[cid, z]
                for candidate in c['candidates']:
                    name = candidate['id']
                    if name in c['matrix'][z]: continue
                    atoms = atoms_by_name[name]; paired = [(z, *atoms[0][1:]), *atoms[1:]]
                    check_atoms(paired, t['charge']); tid = cid + '__' + z + '__at_' + name
                    td = out / 'cells' / tid; td.mkdir(parents=True)
                    xp = td / 'context.xyz'; write_xyz(xp, paired)
                    if xyz(xp) != paired: raise InvalidArtifact('adaptive cross-coordinate serialization changed')
                    task = {'task_id': tid, 'case_id': cid, 'metal': z, 'candidate': name,
                            'xyz': record(xp), 'charge': t['charge'], 'multiplicity': t['multiplicity'],
                            'source_mapping': t['mapping'], 'source_preparation': t['source_preparation']}
                    if name == c['aliases'][prefix + '_' + z]['representative']:
                        if xyz(verify(points[z]['coordinate'])) == paired:
                            task['native_reuse'] = points[z]['MACE']; native_reuse(task, bm)
                        else:
                            task['native_reuse_not_used'] = 'numerical-copy alias; evaluate the actual representative coordinates'
                    pending.append(task)
                    c['matrix'][z][name] = {'status': 'pending', 'task_id': tid, 'xyz': task['xyz'], 'reused': False}
            c['status'] = 'prepared'; tasks.extend(pending)
        except (InvalidArtifact, KeyError, OSError, StopIteration) as exc:
            c['reason'] = str(exc)
        cases.append(c)
    impl = out / 'implementation'; impl.mkdir(); pins = {}
    for p in Path(__file__).parent.glob('*.py'):
        target = impl / p.name; shutil.copyfile(p, target); pins[p.name] = record(target)
    m = {k: bm[k] for k in ('source', 'source_manifest', 'model', 'software', 'orca', 'cpu_python', 'gpu_python', 'resources')}
    m.update(protocol_id=protocol, settings=settings, agreement=record(agreement),
             adaptive_proposals=record(proposals),
             proposal_family=proposal_family, population=('scaled' if proposal_family == 'scaled-angular' else prefix)+str(n),
             declared_case_ids=bm['declared_case_ids'], cases=cases, tasks=tasks,
             implementation=pins, shard_count=shards,
             new_MACE_cells=sum(not t.get('native_reuse') for t in tasks), maximum_new_GFN2_calls=2*len(tasks),
             new_DFT_calls=0, new_optimizations=0, source_preparation_wall_seconds=time.monotonic()-started,
             production_changed=False, reference=None)
    if len(base_paths) == 1: m['base_pool'] = record(base_paths[0])
    else: m['base_pools'] = [record(p) for p in base_paths]
    if proposal_family == 'scaled-angular':
        m.update(proposal_protocol_id=am['protocol_id'], proposal_settings=am['settings'])
    m.update(GFN2_maxiter=maxiter, numerical_policy_id='native_GFN2_MaxIter'+str(maxiter)+'_unchanged_convergence_v1',
             numerical_qualification=record(numerical_qualification) if qualification else None)
    path = out / 'manifest.json'; write_new(path, m); low_prepare(path)
    result = validate(path); write_new(out / 'PREFLIGHT.json', result); return result


def low_prepare(manifest):
    from affordable_workflow import dry_run
    m = read_json(manifest); root = Path(manifest).parent / 'solvent'; shards = []
    for shard in range(m['shard_count']):
        sd = root / ('shard_' + str(shard)); tasks = []
        for cell in m['tasks'][shard::m['shard_count']]:
            for medium in ('vacuum', 'alpb'):
                tid = cell['task_id'] + '__' + medium; td = sd / 'tasks' / tid; td.mkdir(parents=True)
                xp = td / 'core.xyz'; shutil.copyfile(verify(cell['xyz']), xp)
                ip = td / 'endpoint.inp'; template = input_text(cell['charge'], cell['multiplicity'], medium, 'native')
                if m.get('GFN2_maxiter') == 500: template = template.replace('%scf\n', '%scf\n MaxIter 500\n')
                ip.write_text(template)
                tasks.append({**cell, 'task_id': tid, 'cell_id': cell['task_id'], 'case': cell['case_id'],
                              'medium': medium, 'xyz': record(xp), 'input': record(ip), 'output_path': str(td / 'endpoint.out')})
        low = {'protocol_id': m['protocol_id'], 'pool_manifest': record(manifest), 'agreement': m['agreement'],
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
        for t in m['tasks']:
            if t.get('native_reuse'):
                native_reuse(t, m); continue
            if worker is None: worker = WarmGPU(manifest)
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
                op = Path(t['output_path']); receipt_path = Path(t['output_path'] + '.execution.json')
                row['artifacts'] = [record(p) for p in (op, receipt_path) if p.exists()]
                if receipt_path.exists():
                    receipt = read_json(receipt_path)
                    row['execution_status'] = {k: receipt.get(k) for k in ('returncode', 'normal_termination', 'scf_converged', 'launch_error')}
                    if op.exists() and 'the SCF has not converged' in op.read_text(errors='replace'):
                        row.setdefault('reason', 'SCF_not_converged')
                    else: row.setdefault('reason', 'execution_or_output_validation_failed')
                else: row.setdefault('reason', 'execution_not_complete')
            new.setdefault(t['cell_id'], {})[t['medium']] = row
    cases = []; tasks = {t['task_id']: t for t in m['tasks']}
    for c in m['cases']:
        matrix = json.loads(json.dumps(c['matrix']))
        if c['status'] != 'prepared':
            cases.append({**c, 'matrix': matrix, 'pool': {'status': 'unavailable', 'reason': c['reason'],
                          'mathematical': None, 'operational': None, 'rows': None}})
            continue
        for z, cells in matrix.items():
            for q, cell in cells.items():
                if cell['reused']: continue
                rp = root / 'cells' / cell['task_id'] / 'mace_result.json'
                task = tasks[cell['task_id']]
                native_pin = task.get('native_reuse') or (record(rp) if rp.exists() else None)
                native = pinned(native_pin) if native_pin else None; lows = new.get(cell['task_id'], {})
                cell.update(status='unavailable', low=lows, MACE=native_pin)
                if native:
                    try:
                        request = pinned(native['request'])
                        if task.get('native_reuse'):
                            native_reuse(task, m)
                        elif (request['manifest'] != record(manifest) or request['task_id'] != task['task_id'] or
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
        pool = choose_rows(matrix, [q['id'] for q in c['candidates']])
        cases.append({**c, 'matrix': matrix, 'pool': pool})
    native_receipts = [read_json(p) for p in root.glob('cells/*/mace_result.json')]
    result = {'protocol_id': m['protocol_id'], 'manifest': record(manifest), 'cases': cases,
              'numerical_policy_id': m.get('numerical_policy_id', 'primary_native_GFN2_MaxIter125_v1'),
              'GFN2_maxiter': m.get('GFN2_maxiter', 125),
              'denominator': len(cases), 'available': sum(c['pool']['status'] == 'available' for c in cases),
              'required_new_MACE_cells': m['new_MACE_cells'], 'required_new_GFN2_calls': 2 * len(m['tasks']),
              'native_candidate_reuses': sum(bool(t.get('native_reuse')) for t in m['tasks']),
              'MACE_calls_started': sum(r['model_call_started'] for r in native_receipts),
              'MACE_calls_complete': sum(r['status'] == 'complete' for r in native_receipts),
              'GFN2_complete': sum(r['status'] == 'complete' for v in new.values() for r in v.values()),
              'new_DFT_calls': 0, 'new_optimizations': 0, 'production_changed': False,
              'interpretation': 'within-source finite electronic geometry pool, no equilibrium population or affinity claim'}
    write_new(output, result)
    return {k: v for k, v in result.items() if k != 'cases'}


def recover(collection, diagnostic, output):
    """Explicit numerical-completion result; preserve the original failed matrix."""
    original = read_json(collection); d = read_json(diagnostic); dm = pinned(d['manifest'])
    if (original['protocol_id'] != ADAPTIVE_PROTOCOL or
            d['protocol_id'] != 'adaptive_pool_native_GFN2_MaxIter500_diagnostic_v1' or
            not d['both_controls_pass'] or not d['formerly_failed_cell_complete'] or
            d['complete'] != 3 or d['denominator'] != 3 or dm['control_tolerance_kcal_mol'] != .05):
        raise InvalidArtifact('declared numerical-completion checks did not pass')
    from mace_hybrid import HA_TO_KCAL
    for task, row in zip(dm['tasks'], d['rows']):
        if (row['task_id'] != task['task_id'] or row['status'] != 'complete' or
                not row['parameters_identical'] or row['charge_sanity_status'] != 'pass' or
                pinned(row['parameter_export']) != pinned(task['source_parameters']) or
                energy(verify(row['output'])) != row['energy_hartree'] or
                verify(task['input']).read_text().replace(' MaxIter 500\n', '') != verify(task['source_task']['input']).read_text() or
                xyz(verify(task['xyz'])) != xyz(verify(task['source_task']['xyz']))):
            raise InvalidArtifact('numerical-completion state/recipe/result differs')
        receipt = pinned(row['receipt'])
        if not receipt['scf_converged'] or not receipt['normal_termination'] or receipt['returncode'] != 0:
            raise InvalidArtifact('numerical-completion receipt failed')
        if task['source_endpoint']:
            if abs((row['energy_hartree'] - task['source_endpoint']['energy_hartree']) * HA_TO_KCAL) > .05:
                raise InvalidArtifact('original control energy was not reproduced')
    targets = [(t, r) for t, r in zip(dm['tasks'], d['rows']) if t['source_endpoint'] is None]
    if len(targets) != 1: raise InvalidArtifact('only declared missing endpoint may be completed')
    task, row = targets[0]; source = task['source_task']
    result = json.loads(json.dumps(original))
    case = next(c for c in result['cases'] if c['case_id'] == source['case_id'])
    cell = case['matrix'][source['metal']][source['candidate']]
    previous = cell['low'][source['medium']]
    if (source['medium'] != 'vacuum' or previous['status'] != 'unavailable' or
            task['source_output'] not in previous['artifacts'] or cell['low']['alpb']['status'] != 'complete'):
        raise InvalidArtifact('numerical replacement does not identify the failed required cell')
    replacement = {k: row[k] for k in ('manifest', 'task_id', 'output', 'receipt', 'energy_hartree')}
    cell['low']['vacuum'] = {**replacement, 'status': 'complete', 'diagnostic': record(diagnostic),
                             'audit': row, 'previous_failed_attempt': previous}
    cell.update(status='complete', components={'MACE_eV': pinned(cell['MACE'])['energy_eV'],
                'GFN2_vacuum_hartree': row['energy_hartree'],
                'GFN2_ALPB_hartree': cell['low']['alpb']['energy_hartree']})
    for c in result['cases']:
        c['primary_pool'] = c['pool']
        if c['status'] == 'prepared': c['pool'] = choose_rows(c['matrix'], [q['id'] for q in c['candidates']])
    result.update(numerical_policy_id='primary_native_GFN2_with_explicit_MaxIter500_completion_v1',
                  primary_GFN2_maxiter=original.get('GFN2_maxiter', 125), GFN2_maxiter=500,
                  recovered_from=record(collection), numerical_diagnostic=record(diagnostic),
                  recovery_overlay_applied=True, additional_GFN2_attempts=3,
                  primary_GFN2_complete=original['GFN2_complete'], GFN2_complete=original['GFN2_complete']+1,
                  available=sum(c['pool']['status'] == 'available' for c in result['cases']),
                  implementation=record(__file__))
    write_new(output, result)
    return {'available': result['available'], 'numerical_policy_id': result['numerical_policy_id'],
            'original_result_overwritten': False, 'additional_GFN2_attempts': 3}


def main():
    parser = argparse.ArgumentParser(description=__doc__); sub = parser.add_subparsers(dest='operation', required=True)
    p = sub.add_parser('prepare')
    for key in ('source', 'population', 'agreement', 'output'): p.add_argument('--' + key, required=True)
    p.add_argument('--shards', type=int, default=1)
    p = sub.add_parser('recover')
    for key in ('collection', 'diagnostic', 'output'): p.add_argument('--' + key, required=True)
    p = sub.add_parser('prepare-expansion')
    p.add_argument('--base-pool', nargs='+', required=True)
    for key in ('proposals', 'agreement', 'output'): p.add_argument('--' + key, required=True)
    p.add_argument('--shards', type=int, default=1)
    p.add_argument('--maxiter', type=int, choices=(125, 500), default=125)
    p.add_argument('--numerical-qualification')
    p.add_argument('--proposal-family', choices=('angular', 'joint', 'scaled-angular'), default='angular')
    for name in ('validate', 'execute-mace', 'collect'):
        p = sub.add_parser(name); p.add_argument('--manifest', required=True)
        if name == 'collect': p.add_argument('--output', required=True)
    a = vars(parser.parse_args()); op = a.pop('operation').replace('-', '_')
    print(json.dumps(globals()[op](**a)))


if __name__ == '__main__': main()
