"""One-rank scalar replay of the fixed 28-source union/adaptive reference pools."""
from __future__ import annotations

import argparse
import copy
import datetime as dt
import json
import os
from pathlib import Path
import re
import shutil

from affordable_common import InvalidArtifact, HA_TO_KCAL, read_json, record, verify, write_new, xyz
from affordable_workflow import dry_run, execute as execute_existing
from adaptive_origin_recovery import snapshot
from compact_solvation import completed, diagnostics
from nikasha_pool import choose_rows, score

ROOT = Path(__file__).resolve().parents[1]
PROTOCOL = 'Nikasha_union_reference28_native_GFN2_rank1_scalar_qualification_v1'
REFERENCE_HASH = '2f142883dcfd8b452d6c75e714b49ee1cac7a29ae317a4f289dcbc740498c5c7'
CANDIDATES = ('origin', 'adaptive_Ca', 'adaptive_La')
METALS = ('Ca', 'La')
MEDIA = ('vacuum', 'alpb')
HOSTS = ('node-224-2t-8gpu-1', 'node-128-512g-8gpu-1')
SETTINGS = {'ranks': 1, 'workers': 32, 'allocated_cpus': 32, 'memory_MiB': 65536,
            'reference_ranks': 8, 'case_count': 28, 'candidate_ids': list(CANDIDATES),
            'new_calls': 336, 'cell_tolerance_kcal_mol': .1, 'pool_R_tolerance_kcal_mol': .2,
            'coordinate_compatibility_A': 1e-12, 'numerical_scope': 'scalar_only'}


def same_coordinates(a, b):
    if len(a) != len(b) or any(x[0] != y[0] for x, y in zip(a, b)):
        raise InvalidArtifact('source elements differ')
    delta = max(abs(x[i]-y[i]) for x, y in zip(a, b) for i in (1, 2, 3))
    if delta > SETTINGS['coordinate_compatibility_A']:
        raise InvalidArtifact('source coordinates differ')
    return delta


def actual_task(low):
    receipt = read_json(verify(low['receipt']))
    mp = verify(receipt['manifest']); manifest = read_json(mp)
    task = next(t for t in manifest.get('all_tasks', manifest['tasks']) if t['task_id'] == receipt['task_id'])
    if receipt['parallelism']['nprocs'] != 8:
        raise InvalidArtifact('comparator did not use eight ranks')
    return mp, manifest, task


def native_reuse(cell, task, model):
    result = read_json(verify(cell['MACE']))
    if result['status'] not in ('complete', 'computed') or result['energy_eV'] != cell['components']['MACE_eV']:
        raise InvalidArtifact('native energy unavailable or different')
    if 'request' in result:
        req = read_json(verify(result['request'])); actual_model = result['model']
        native_xyz, charge, mult = req['xyz'], req['charge'], req['multiplicity']
    else:
        parent = read_json(verify(result['manifest'])); actual_model = parent['model']
        nt = next(t for t in parent['tasks'] if t['task_id'] == result['task_id'])
        native_xyz, charge, mult = nt['xyz'], nt['charge'], nt['spin_multiplicity']
    if actual_model != model or (charge, mult) != (task['charge'], task['multiplicity']):
        raise InvalidArtifact('native model/state mismatch')
    delta = same_coordinates(xyz(verify(native_xyz)), xyz(verify(task['xyz'])))
    return {'result': cell['MACE'], 'xyz': native_xyz, 'energy_eV': result['energy_eV'],
            'maximum_coordinate_difference_A': delta}


def exact_key(task):
    return task['input']['sha256'], task['xyz']['sha256'], task['charge'], task['multiplicity']


def reference_rows(reference):
    if record(reference)['sha256'] != REFERENCE_HASH:
        raise InvalidArtifact('not the declared frozen reference')
    r = read_json(reference)
    if len(r['rows']) != 28 or len({x['case_id'] for x in r['rows']}) != 28 or r['calibration_denominator'] != 25:
        raise InvalidArtifact('reference membership differs')
    return r


def prepare(reference, pilot, agreement, output):
    ref = reference_rows(reference); prior = read_json(pilot)
    if prior['ranks'] != 1 or prior['complete'] != 8:
        raise InvalidArtifact('actual completed one-rank pilot required')
    pm = read_json(verify(prior['manifest'])); pilot_keys = {exact_key(t) for t in pm['tasks']}
    root = Path(output).resolve(); root.mkdir(parents=True, exist_ok=False)
    impl = snapshot(Path(__file__).parent, root/'implementation')
    collections = {}; tasks = []; cases = []; source_manifests = {}; backend = None
    maxiters = {}; atom_counts = []; exact_matches = []
    for rr in ref['rows']:
        cp = verify(rr['collection']); pool = collections.setdefault(str(cp), read_json(cp))
        source = next(c for c in pool['cases'] if c['case_id'] == rr['case_id'])
        if source['pool']['status'] != 'available' or set(source['matrix']['Ca']) != set(CANDIDATES):
            raise InvalidArtifact('required original pool missing')
        replay = choose_rows(source['matrix'], CANDIDATES)
        if replay != source['pool']:
            raise InvalidArtifact('original pool replay differs')
        case = {k: rr[k] for k in ('case_id', 'actual_union_case_id', 'role', 'biological_group', 'expected_class', 'collection')}
        case.update(archived_pool=source['pool'], archived_matrix=source['matrix'], source_case_id=source['case_id'])
        cases.append(case)
        for candidate in CANDIDATES:
            paired = {}
            for metal in METALS:
                cell = source['matrix'][metal][candidate]
                for medium in MEDIA:
                    low = cell['low'][medium]; smp, sm, original = actual_task(low)
                    pin = completed(smp, original['task_id'])
                    if pin is None or pin['energy_hartree'] != low['energy_hartree']:
                        raise InvalidArtifact('archived executed energy differs')
                    audit = diagnostics(pin, original)
                    if backend is None: backend = sm['orca']
                    if backend != sm['orca']: raise InvalidArtifact('mixed ORCA executables')
                    source_manifests[str(smp)] = record(smp)
                    body = verify(original['input']).read_text()
                    for required in ('Native-GFN2-xTB NoAutostart', '%maxcore 2000', 'SmearTemp 300', 'UseXTBMixer true'):
                        if required not in body: raise InvalidArtifact('source native recipe differs')
                    if re.search(r'%pal\b|\bEnGrad\b', body, re.I):
                        raise InvalidArtifact('source is not scalar rank-renderable input')
                    count = re.search(r'\bMaxIter\s+(\d+)', body, re.I)
                    limit = 'explicit_'+count[1] if count else 'default_125'
                    maxiters[limit] = maxiters.get(limit, 0)+1
                    native = native_reuse(cell, original, ref['model'])
                    same_coordinates(xyz(verify(cell['xyz'])), xyz(verify(original['xyz'])))
                    atoms = xyz(verify(original['xyz'])); atom_counts.append(len(atoms))
                    if medium == 'vacuum': paired[metal] = (atoms, original['charge'])
                    else: same_coordinates(paired[metal][0], atoms)
                    task_id = f"{rr['case_id']}__{candidate}__{metal}__{medium}"
                    if exact_key(original) in pilot_keys: exact_matches.append(task_id)
                    directory = root/'tasks'/task_id; directory.mkdir(parents=True)
                    t = {k: original[k] for k in ('charge', 'multiplicity', 'metal', 'medium')}
                    for key, filename in (('input', 'endpoint.inp'), ('xyz', 'core.xyz')):
                        shutil.copyfile(verify(original[key]), directory/filename); t[key] = record(directory/filename)
                    t.update(task_id=task_id, case_id=rr['case_id'], candidate=candidate,
                             output_path=str(directory/'endpoint.out'),
                             archived=pin, archived_audit=audit, original_task=original,
                             original_manifest=record(smp), native=native, SCF_limit_policy=limit)
                    tasks.append(t)
            ca, la = paired['Ca'], paired['La']
            if ca[0][0][0] != 'Ca' or la[0][0][0] != 'La' or la[1]-ca[1] != 1:
                raise InvalidArtifact('metal pair/state mismatch')
            same_coordinates([('La', *ca[0][0][1:]), *ca[0][1:]], la[0])
    if exact_matches: raise InvalidArtifact('exact pilot reuse found: adjust finite declaration before running')
    m = {'protocol_id': PROTOCOL, 'settings': SETTINGS, 'reference': record(reference),
         'pilot': record(pilot), 'agreement': record(agreement), 'implementation': impl,
         'orca': backend, 'source_manifests': list(source_manifests.values()), 'cases': cases,
         'tasks': tasks, 'all_tasks': tasks, 'execution_resources': {'mpi_ranks': 1, 'concurrent_tasks': 32},
         'execution_policy': {k: impl[v] for k, v in [('task_runner','run_orca_task_manifest.py'), ('runtime_renderer','render_orca_runtime_input.py')]},
         'inventory': {'cases': 28, 'cells': len(tasks), 'fresh_calls': len(tasks), 'exact_pilot_reuse': exact_matches,
                       'SCF_limits': maxiters, 'minimum_atoms': min(atom_counts), 'maximum_atoms': max(atom_counts)},
         'production_changed': False}
    mp = root/'manifest.json'; write_new(mp, m); v = validate(mp)
    write_new(root/'PREFLIGHT.json', {'manifest': record(mp), 'inventory': m['inventory'], 'dry_run': v, 'new_molecular_calls': 0})
    return {'manifest': record(mp), 'inventory': m['inventory'], 'new_molecular_calls': 0}


def validate(manifest):
    m = read_json(manifest)
    if m['protocol_id'] != PROTOCOL or m['settings'] != SETTINGS:
        raise InvalidArtifact('frozen scalar panel policy differs')
    ref = reference_rows(verify(m['reference'])); verify(m['agreement']); verify(m['pilot'])
    for pin in m['implementation'].values(): verify(pin)
    if m['execution_resources'] != {'mpi_ranks': 1, 'concurrent_tasks': 32}:
        raise InvalidArtifact('rank/worker settings changed')
    expected = {f"{r['case_id']}__{q}__{z}__{s}" for r in ref['rows'] for q in CANDIDATES for z in METALS for s in MEDIA}
    if len(m['tasks']) != 336 or m['tasks'] != m['all_tasks'] or {t['task_id'] for t in m['tasks']} != expected:
        raise InvalidArtifact('fixed 336-cell membership differs')
    if [c['case_id'] for c in m['cases']] != [r['case_id'] for r in ref['rows']]:
        raise InvalidArtifact('case order/membership differs')
    sources = {str(verify(p)): read_json(verify(p)) for p in m['source_manifests']}
    for t in m['tasks']:
        source = sources[str(verify(t['original_manifest']))]
        old = next(x for x in source.get('all_tasks', source['tasks']) if x['task_id'] == t['original_task']['task_id'])
        if old != t['original_task'] or source['orca'] != m['orca']:
            raise InvalidArtifact('original source task/backend differs')
        if any(t[k] != old[k] for k in ('charge', 'multiplicity', 'metal', 'medium')):
            raise InvalidArtifact('source electronic state changed')
        if any(verify(t[k]).read_bytes() != verify(old[k]).read_bytes() for k in ('input', 'xyz')):
            raise InvalidArtifact('source coordinates/input changed')
        out = Path(t['output_path']).resolve()
        if out.parent != verify(t['input']).parent or not out.is_relative_to(Path(manifest).resolve().parent):
            raise InvalidArtifact('task escaped panel directory')
        nr = read_json(verify(t['native']['result']))
        if nr['status'] not in ('complete', 'computed') or nr['energy_eV'] != t['native']['energy_eV']:
            raise InvalidArtifact('native reuse changed')
    return dry_run(manifest)


def execute(manifest):
    validate(manifest)
    if (os.environ.get('SLURM_JOB_NODELIST') not in HOSTS or
        int(os.environ.get('SLURM_CPUS_ON_NODE','0')) != 32 or
        int(os.environ.get('SLURM_MEM_PER_NODE','0')) != 65536):
        raise InvalidArtifact('declared CPU-only32CPU/64GiB host allocation required')
    return execute_existing(manifest)


def collect(manifest, output):
    validate(manifest); m = read_json(manifest); rows = []
    for t in m['tasks']:
        row = {k: t[k] for k in ('task_id', 'case_id', 'candidate', 'metal', 'medium')}
        op = Path(t['output_path']); rp = Path(str(op)+'.execution.json')
        row.update(status='unavailable', reason=None, energy_hartree=None,
                   available_artifacts=[record(p) for p in (op,rp) if p.exists()])
        try:
            pin = completed(manifest,t['task_id'])
            if pin is None: raise InvalidArtifact('native endpoint incomplete or failed')
            rec = read_json(verify(pin['receipt'])); audit = diagnostics(pin,t)
            if rec['parallelism'] != {'nprocs':1,'omp_threads_per_rank':1,'binding':'none'}:
                raise InvalidArtifact('actual rank/thread/binding differs')
            if read_json(verify(audit['parameter_export'])) != read_json(verify(t['archived_audit']['parameter_export'])):
                raise InvalidArtifact('exported native parameters changed')
            text = verify(pin['output']).read_text(); hits = re.findall(r'SCF CONVERGED AFTER\s+(\d+)\s+CYCLES',text)
            wall = (dt.datetime.fromisoformat(rec['finished_at_utc'].replace('Z','+00:00'))-
                    dt.datetime.fromisoformat(rec['started_at_utc'].replace('Z','+00:00'))).total_seconds()
            delta = (pin['energy_hartree']-t['archived']['energy_hartree'])*HA_TO_KCAL
            row.update(status='complete',**pin,audit=audit,wall_seconds=wall,
                       iterations=int(hits[-1]) if hits else None,delta_vs_archived8_kcal_mol=delta,
                       passes_component_gate=abs(delta)<=SETTINGS['cell_tolerance_kcal_mol'])
        except Exception as exc: row['reason'] = str(exc)
        rows.append(row)
    result = {'protocol_id':PROTOCOL,'manifest':record(manifest),'denominator':336,
              'complete':sum(r['status']=='complete' for r in rows),'rows':rows,'new_calls_in_collection':0}
    write_new(output,result); return {'collection':record(output),'complete':result['complete'],'denominator':336}


def decision(value, bands):
    if value is None: return 'unavailable'
    if value <= bands['Ca_max']: return 'Ca-supported'
    if value >= bands['La_min']: return 'La-supported'
    return 'inconclusive'


def compare(collection, output):
    coll = read_json(collection); mp = verify(coll['manifest']); m = read_json(mp)
    if coll['protocol_id'] != PROTOCOL or m['settings'] != SETTINGS: raise InvalidArtifact('collection policy mismatch')
    ref = reference_rows(verify(m['reference'])); by = {r['task_id']:r for r in coll['rows']}; cases = []
    if len(by) != 336 or set(by) != {t['task_id'] for t in m['tasks']}:
        raise InvalidArtifact('collection membership mismatch')
    tasks = {t['task_id']:t for t in m['tasks']}
    for source in m['cases']:
        matrix = copy.deepcopy(source['archived_matrix']); cellrows = []; geometries = []
        for z in METALS:
            for q in CANDIDATES:
                lows = [by[f"{source['case_id']}__{q}__{z}__{s}"] for s in MEDIA]
                cellrows.extend(lows)
                if any(x['status']!='complete' for x in lows):
                    matrix[z][q] = {'status':'unavailable','reason':'required_rank1_component_unavailable'}; continue
                matrix[z][q]['components'] = {'MACE_eV':tasks[lows[0]['task_id']]['native']['energy_eV'],
                    'GFN2_vacuum_hartree':lows[0]['energy_hartree'],'GFN2_ALPB_hartree':lows[1]['energy_hartree']}
                matrix[z][q]['low'] = dict(zip(MEDIA,lows)); matrix[z][q]['reused'] = False
        pool = choose_rows(matrix,CANDIDATES); old = source['archived_pool']; variants = {}
        for mode in ('mathematical','operational'):
            new_R = pool[mode]['composite_R_model_kcal_mol'] if pool['status']=='available' else None
            old_R = old[mode]['composite_R_model_kcal_mol']; delta = new_R-old_R if new_R is not None else None
            variants[mode] = {'archived8_R':old_R,'rank1_R':new_R,'delta_R_kcal_mol':delta,
                'passes_pool_R_gate':abs(delta)<=SETTINGS['pool_R_tolerance_kcal_mol'] if delta is not None else None,
                'archived8_decision':decision(old_R,ref['variants'][mode]['bands']),
                'rank1_decision':decision(new_R,ref['variants'][mode]['bands'])}
        for q in CANDIDATES:
            old_score = score(source['archived_matrix']['Ca'][q]['components'], source['archived_matrix']['La'][q]['components'])
            new_score = score(matrix['Ca'][q]['components'],matrix['La'][q]['components']) if all(matrix[z][q]['status']=='complete' for z in METALS) else None
            geometries.append({'candidate':q,'archived8':old_score,'rank1':new_score,
                'delta_R_kcal_mol':new_score['composite_R_model_kcal_mol']-old_score['composite_R_model_kcal_mol'] if new_score else None})
        complete = all(x['status']=='complete' for x in cellrows)
        cases.append({**{k:source[k] for k in ('case_id','role','biological_group','expected_class')},
            'status':pool['status'],'archived8_pool':old,'rank1_pool':pool,'variants':variants,'geometries':geometries,
            'maximum_component_delta_kcal_mol':max(abs(x['delta_vs_archived8_kcal_mol']) for x in cellrows) if complete else None,
            'all_numerical_gates_pass':complete and all(x['passes_component_gate'] for x in cellrows) and all(v['passes_pool_R_gate'] for v in variants.values())})
    result = {'protocol_id':PROTOCOL,'collection':record(collection),'reference':m['reference'],
              'case_denominator':28,'complete':sum(c['status']=='available' for c in cases),'cell_denominator':336,
              'complete_cells':coll['complete'],'cases':cases,'all_numerical_gates_pass':all(c['all_numerical_gates_pass'] for c in cases),
              'new_molecular_calls':0,'reference_refitted':False,'production_changed':False,'native_gradient_qualified':False}
    write_new(output,result); return {'comparison':record(output),'complete':result['complete'],
                                    'all_numerical_gates_pass':result['all_numerical_gates_pass']}


def main():
    p=argparse.ArgumentParser(description=__doc__); sub=p.add_subparsers(dest='op',required=True)
    for op,fields in {'prepare':('reference','pilot','agreement','output'),'validate':('manifest',),
                      'execute':('manifest',),'collect':('manifest','output'),'compare':('collection','output')}.items():
        q=sub.add_parser(op)
        for f in fields:q.add_argument('--'+f,required=True,type=Path)
    a=vars(p.parse_args());print(json.dumps(globals()[a.pop('op')](**a),indent=2))

if __name__=='__main__':main()
