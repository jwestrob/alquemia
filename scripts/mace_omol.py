"""Versioned MACE-OMOL vacuum descriptor, using the established task executor."""
from __future__ import annotations
import argparse
import copy
import importlib.metadata
import json
import math
import os
from pathlib import Path
import resource
import time
import traceback
import numpy as np
from affordable_common import InvalidArtifact, cache_key, read_json, record, verify, write_new, xyz
from mace_hybrid import EV_TO_KCAL, Z, accepted_attempt, check_atoms, rotation, write_xyz
from mace_canonical_run import inventory, all_tasks
from mace_global_benchmark import snapshot

SCHEMA = 'alquemia.mace_omol.v1'
PROTOCOL = 'mace_omol_0_100m_vacuum_descriptor_v1'
COMPONENT = 'MACE_OMOL_total_vacuum_energy'
CHECKPOINT_SHA = '9b64b4fd5153ca578c694abc57806d8111050de6ff652e695c9b525bc4d36469'
TOL = {'energy_kcal_mol': .01, 'force_max_eV_A': .001}
UNAVAILABLE = {'density_coefficients': None, 'predicted_atomic_charge_status': 'not_produced_by_this_model',
               'solvent_status': 'not_part_of_vacuum_descriptor', 'S_kcal_mol': None,
               'reference': None, 'calibrated_class': None, 'relaxation_correction_kcal_mol': None,
               'entropy_correction_kcal_mol': None}


def model(software):
    s = read_json(software); inspect = read_json(verify(s['inspection']))
    if s['checkpoint']['sha256'] != CHECKPOINT_SHA or inspect['checkpoint'] != s['checkpoint']:
        raise InvalidArtifact('only declared 100M OMOL checkpoint permitted')
    if not set(Z.values()) <= set(inspect['elements']) or inspect['heads'] != ['omol']:
        raise InvalidArtifact('required elements/head unsupported')
    return {'family': 'MACE_OMOL', 'checkpoint': s['checkpoint'], 'head': 'omol', 'dtype': 'float64',
            'periodic': False, 'solvent': None, 'energy_component': COMPONENT,
            'spin_input_convention': 'multiplicity', 'embedding_specs': inspect['embedding_specs'],
            'r_max_A': inspect['r_max_A'], 'execution_adapter': None, 'extra_dispersion': False}


def task_base(source):
    t = copy.deepcopy(source); t.update(kind='core', variant='primary', energy_component=COMPONENT)
    return t


def qualification_tasks(data):
    tasks = []
    for t in all_tasks(data):
        if t['case_id'] not in ('1H4I', '4MAE'):
            continue
        for variant in (('primary', 'repeat', 'rotate', 'translate') if t['case_id']=='1H4I' else ('primary',)):
            new = task_base(t); new.update(task_id=t['task_id']+'_'+variant, variant=variant,
                                           source_xyz=t['xyz'], rotation_matrix=(rotation() if variant=='rotate' else np.eye(3)).tolist(),
                                           translation_A=[10., -7., 3.] if variant=='translate' else [0., 0., 0.])
            tasks.append(new)
    return tasks


def common(source_inventory, software, agreement, output, stage):
    data = inventory(source_inventory); out = Path(output).resolve(); out.mkdir(parents=True, exist_ok=False)
    # Reuse the canonical package's audited dependency closure, plus this worker.
    names = ('mace_omol.py', 'mace_canonical_run.py', 'mace_canonical.py', 'mace_curvature.py',
             'affordable_response.py', 'mace_mechanics_run.py', 'mace_short_engine.py',
             'mace_omol_readout.py', 'mace_omol_coordination.py')
    pins = snapshot(out, names)
    m = {'schema_version': SCHEMA, 'protocol_id': PROTOCOL, 'stage': stage,
         'inventory': record(source_inventory), 'software': record(software), 'agreement': record(agreement),
         'implementation': pins, 'model': model(software), 'tolerances': TOL, 'reused': {},
         'compute_budget': None, 'evidence_use': 'consumed_retrospective_development', **UNAVAILABLE}
    return data, out, m


def seal(out, m):
    for t in m['tasks']:
        t['cache_key'] = cache_key({'task': t, 'model': m['model'], 'software': m['software'], 'implementation': m['implementation']})
    write_new(out/'manifest.json', m)
    return validate(out/'manifest.json')


def prepare_qualification(source_inventory, software, agreement, output):
    data, out, m = common(source_inventory, software, agreement, output, 'qualification')
    tasks = qualification_tasks(data)
    for t in tasks:
        if t['variant'] in ('rotate', 'translate'):
            rows = xyz(verify(t['source_xyz'])); coords = np.array([a[1:] for a in rows])
            center = coords[next(i for i,a in enumerate(rows) if a[0]==t['metal'])]
            moved = (coords-center)@np.array(t['rotation_matrix']).T+center+np.array(t['translation_A'])
            p = out/(t['task_id']+'.xyz'); write_xyz(p, [(a[0], *v) for a,v in zip(rows, moved)])
            t['xyz'] = record(p)
    m['tasks'] = tasks
    return seal(out, m)


def qualified(path):
    saved = read_json(path); mp = verify(saved['manifest']); validate(mp)
    current = collect(mp)
    if saved != current or not current['numerical_gate_pass']:
        raise InvalidArtifact('complete actual OMOL numerical qualification required')
    return current, read_json(mp)


def benchmark_tasks(data, mechanics):
    from mace_mechanics_run import preparation
    p = preparation(mechanics)
    tasks = [task_base(t) for t in all_tasks(data)]
    for name in ('GGR_extended', 'GGR_connected', 'ALPHA_1F6S', 'ALPHA_6IP9'):
        ref = p['cases'][name]; state = read_json(verify(ref))
        for metal in ('La', 'Ca'):
            endpoint = state['grids']['center']['endpoints'][metal]
            tasks.append({'task_id':name+'_'+metal, 'case_id':name, 'kind':'core', 'variant':'primary',
                          'metal':metal, 'expected_class':state['evidence']['direction'],
                          'evaluation_role':'retrospective_nonPQQ_direction', 'evidence':state['evidence'],
                          'preparation':ref, 'explicit_waters':state['explicit_waters'],
                          'energy_component':COMPONENT, **endpoint})
    return tasks


def reuse_primary(data, qualification):
    c, m = qualified(qualification); tasks = {t['task_id']:t for t in m['tasks']}; reused = {}
    for t in all_tasks(data):
        if t['case_id'] not in ('1H4I','4MAE'):
            continue
        key = t['task_id']+'_primary'; source = tasks[key]
        if any(source[k] != t[k] for k in ('xyz','charge','spin_multiplicity','case_id','metal')):
            raise InvalidArtifact('qualification cache does not match frozen benchmark input')
        reused[t['task_id']] = {'source_collection':record(qualification), 'source_task_id':key}
    return reused


def prepare_benchmark(qualification, mechanics, agreement, output):
    qc, qm = qualified(qualification)
    data, out, m = common(verify(qm['inventory']), verify(qm['software']), agreement, output, 'benchmark')
    reused = reuse_primary(data, qualification)
    m.update(qualification=record(qualification), mechanics=record(mechanics), reused=reused,
             tasks=[t for t in benchmark_tasks(data, mechanics) if t['task_id'] not in reused])
    return seal(out, m)


def validate(manifest):
    m = read_json(manifest)
    if m.get('stage','').startswith('coordination_'):
        from mace_omol_coordination import validate as validate_coordination
        return validate_coordination(manifest)
    if m['schema_version'] != SCHEMA or m['protocol_id'] != PROTOCOL or m['tolerances'] != TOL:
        raise InvalidArtifact('OMOL protocol/acceptance changed')
    software = verify(m['software']); s = read_json(software)
    if m['model'] != model(software):
        raise InvalidArtifact('OMOL Hamiltonian/checkpoint changed')
    for ref in [m['agreement'], *m['implementation'].values(), s['python'], s['requirements'], s['checkpoint'],
                s['inspection'], s['download_receipt'], *read_json(verify(s['backend_source_inventory']))['files']]:
        verify(ref)
    data = inventory(verify(m['inventory']))
    if m['stage']=='qualification':
        expected = qualification_tasks(data)
        if len(m['tasks']) != 10 or m['reused']:
            raise InvalidArtifact('finite qualification inventory changed')
    elif m['stage']=='benchmark':
        qc, qm = qualified(verify(m['qualification']))
        reused = reuse_primary(data, verify(m['qualification']))
        if m['reused'] != reused or m['model'] != qm['model'] or m['software'] != qm['software']:
            raise InvalidArtifact('qualified OMOL method/cache changed')
        expected = [t for t in benchmark_tasks(data, verify(m['mechanics'])) if t['task_id'] not in reused]
        if len(m['tasks']) != 60 or len(reused) != 4:
            raise InvalidArtifact('finite60-task benchmark inventory changed')
    elif m['stage']=='readout':
        from mace_omol_readout import source, tasks
        _, parent = source(verify(m['source_collection']))
        expected = tasks(parent)
        if (m['reused'] or m['model'] != parent['model'] or m['software'] != parent['software']
                or m['inventory'] != parent['inventory']):
            raise InvalidArtifact('readout replay changed the source method or inventory')
    else:
        raise InvalidArtifact('unsupported OMOL stage')
    if len(expected) != len(m['tasks']):
        raise InvalidArtifact('expected task count differs')
    for t, e in zip(m['tasks'], expected):
        if any(t.get(k) != v for k,v in e.items() if k!='xyz'):
            raise InvalidArtifact('OMOL input differs from frozen source task')
        if t['variant'] in ('primary', 'repeat'):
            if t['xyz'] != e['xyz']:
                raise InvalidArtifact('original core bytes changed')
        else:
            source = xyz(verify(t['source_xyz'])); actual = xyz(verify(t['xyz']))
            coords = np.array([a[1:] for a in source]); center = coords[next(i for i,a in enumerate(source) if a[0]==t['metal'])]
            wanted = (coords-center)@np.array(e['rotation_matrix']).T+center+np.array(e['translation_A'])
            if [a[0] for a in source] != [a[0] for a in actual] or np.max(np.abs(wanted-np.array([a[1:] for a in actual]))) > 1e-12:
                raise InvalidArtifact('rigid geometry differs from declared transform')
        check_atoms(xyz(verify(t['xyz'])), t['charge'])
        payload = {k:v for k,v in t.items() if k!='cache_key'}
        if t['cache_key'] != cache_key({'task':payload, 'model':m['model'], 'software':m['software'], 'implementation':m['implementation']}):
            raise InvalidArtifact('OMOL scientific cache identity changed')
    return {'status':'pass', 'tasks':len(m['tasks']), 'manifest':record(manifest)}


def input_batch(calc, atoms, charge, multiplicity):
    batch = calc._atoms_to_batch(atoms)
    q = float(batch['total_charge'].item()); spin = float(batch['total_spin'].item())
    if q != charge or spin != multiplicity or calc.head != 'omol':
        raise InvalidArtifact('actual OMOL batch charge/multiplicity/head differs')
    for name, value in (('total_charge', q), ('total_spin', spin)):
        spec = calc.models[0].embedding_specs[name]; index = value+spec.get('offset', 0)
        if spec['type'] != 'categorical' or value != int(value) or not 0 <= index < spec['num_classes']:
            raise InvalidArtifact('state outside checkpoint categorical embedding')
    edges=batch['edge_index']
    return {'charge':q, 'spin_multiplicity':spin, 'head':calc.head,
            'metal_neighbor_edge_count':int(((edges[0]==0)|(edges[1]==0)).sum().item()),
            'atoms':len(atoms), 'graph_count':int(batch['ptr'].numel()-1), 'charge_prediction_claimed':False}


def worker(manifest, task_id, output, memory_mode):
    import torch
    from ase import Atoms
    from mace.calculators import mace_omol
    m = read_json(manifest); t = next(t for t in m['tasks'] if t['task_id']==task_id)
    out = Path(output); start = time.monotonic()
    result = {'task_id':task_id, 'cache_key':t['cache_key'], 'manifest':record(manifest),
              'status':'unavailable', 'memory_mode':memory_mode, 'energy_component':COMPONENT, **UNAVAILABLE}
    try:
        if not os.environ.get('SLURM_JOB_ID') or not torch.cuda.is_available() or memory_mode!='native':
            raise InvalidArtifact('native OMOL requires allocated CUDA')
        torch.set_num_threads(int(os.environ['SLURM_CPUS_PER_TASK'])); torch.set_default_dtype(torch.float64)
        torch.cuda.reset_peak_memory_stats(); props = torch.cuda.get_device_properties(0)
        result['device'] = {'name':props.name, 'total_memory_bytes':props.total_memory,
                            'CUDA_VISIBLE_DEVICES':os.environ.get('CUDA_VISIBLE_DEVICES')}
        calc = mace_omol(model=str(verify(m['model']['checkpoint'])), device='cuda', default_dtype='float64')
        model_obj = calc.models[0]
        if type(model_obj).__name__ != 'ScaleShiftMACE' or dict(model_obj.embedding_specs) != m['model']['embedding_specs']:
            raise InvalidArtifact('loaded model type/state embedding mismatch')
        rows = xyz(verify(t['xyz']))
        if not {Z[a[0]] for a in rows} <= set(model_obj.atomic_numbers.tolist()):
            raise InvalidArtifact('checkpoint lacks a required element')
        atoms = Atoms([a[0] for a in rows], positions=[a[1:] for a in rows], pbc=False)
        atoms.info.update(charge=t['charge'], spin=t['spin_multiplicity']); atoms.calc = calc
        result['input_state_check'] = input_batch(calc, atoms, t['charge'], t['spin_multiplicity'])
        if result['input_state_check']['graph_count'] != 1:
            raise InvalidArtifact('unexpected graph batching')
        if t.get('require_isolated_metal') and result['input_state_check']['metal_neighbor_edge_count'] != 0:
            raise InvalidArtifact('detached metal still has graph neighbor edges')
        versions = {name:p._version for name,p in model_obj.named_parameters()}
        capture = None
        if t.get('capture_native_readout'):
            from mace_omol_readout import NativeCapture
            capture = NativeCapture(model_obj)
        result['model_load_seconds'] = time.monotonic()-start
        torch.cuda.synchronize(); evaluation = time.monotonic()
        value = float(atoms.get_potential_energy()); forces = np.asarray(atoms.get_forces(), dtype=np.float64)
        torch.cuda.synchronize(); result['evaluation_seconds'] = time.monotonic()-evaluation
        if not math.isfinite(value) or forces.shape != (len(atoms),3) or not np.isfinite(forces).all():
            raise InvalidArtifact('invalid OMOL energy/analytic forces')
        result['parameter_versions_unchanged'] = versions == {name:p._version for name,p in model_obj.named_parameters()}
        if not result['parameter_versions_unchanged']:
            raise InvalidArtifact('model parameters mutated during inference')
        fp = out/'forces_eV_A.npy'; np.save(fp, forces)
        result.update(status='computed', energy_eV=value, forces=record(fp), charge_check=None,
                      force_definition='negative_Cartesian_gradient_of_total_vacuum_OMOL_energy')
        if capture is not None:
            result['native_readout'] = capture.save(calc, value, out)
    except Exception as exc:
        result.update(status='failed', reason=str(exc), exception_type=type(exc).__name__); traceback.print_exc()
    finally:
        result.update(wall_seconds=time.monotonic()-start, peak_host_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                      peak_cuda_allocated_bytes=torch.cuda.max_memory_allocated() if torch.cuda.is_available() else None,
                      peak_cuda_reserved_bytes=torch.cuda.max_memory_reserved() if torch.cuda.is_available() else None,
                      slurm_job_id=os.environ.get('SLURM_JOB_ID'), allocated_cpus=os.environ.get('SLURM_CPUS_PER_TASK'),
                      allocated_host_mem_MiB=os.environ.get('SLURM_MEM_PER_NODE'),
                      versions={n:importlib.metadata.version(n) for n in ('torch','mace-torch','ase','e3nn')})
        write_new(out/'result.json', result)
    print(json.dumps({k:result.get(k) for k in ('task_id','status','reason','evaluation_seconds')}), flush=True)
    return result


def accepted_state(result, task):
    state = result.get('input_state_check', {})
    if task.get('require_isolated_metal') and state.get('metal_neighbor_edge_count') != 0:
        return False
    if task.get('capture_native_readout'):
        from mace_omol_readout import checked_arrays
        checked_arrays(result, task)
    return (result.get('energy_component')==COMPONENT and result.get('density_coefficients') is None
            and result.get('parameter_versions_unchanged') is True and state.get('charge')==task['charge']
            and state.get('spin_multiplicity')==task['spin_multiplicity'] and state.get('head')=='omol'
            and state.get('graph_count')==1 and state.get('charge_prediction_claimed') is False
            and math.isfinite(result['energy_eV']))


def collect(manifest):
    if read_json(manifest).get('stage','').startswith('coordination_'):
        from mace_omol_coordination import collect as collect_coordination
        return collect_coordination(manifest)
    mp = Path(manifest).resolve(); m = read_json(mp); rows = {}; attempts = []
    for t in m['tasks']:
        found = []
        for a in sorted((mp.parent/'execution'/t['task_id']).glob('attempt_*')):
            r = accepted_attempt(a,t,mp)
            if r is not None: found.append(r)
            attempts.append({'task_id':t['task_id'], 'path':str(a), 'accepted':r is not None,
                             'receipt':record(a/'receipt.json') if (a/'receipt.json').exists() else None})
        rows[t['task_id']] = found[-1] if found else {'status':'unavailable', 'energy_eV':None}
    complete = all(r['status']=='computed' for r in rows.values()); checks = []
    if complete and m['stage']=='qualification':
        for variant in ('repeat','rotate','translate'):
            deltas = {}
            for metal in ('La','Ca'):
                p = rows[f'1H4I_{metal}_primary']; r = rows[f'1H4I_{metal}_{variant}']
                de = (r['energy_eV']-p['energy_eV'])*EV_TO_KCAL; deltas[metal]=de
                matrix = rotation() if variant=='rotate' else np.eye(3)
                df = float(np.max(np.abs(np.load(verify(r['forces']))@matrix-np.load(verify(p['forces'])))))
                checks.append({'name':metal+'_'+variant, 'energy_error_kcal_mol':de, 'force_error_eV_A':df,
                               'pass':abs(de)<=TOL['energy_kcal_mol'] and df<=TOL['force_max_eV_A']})
            checks.append({'name':'R_'+variant, 'error_kcal_mol':deltas['Ca']-deltas['La'],
                           'pass':abs(deltas['Ca']-deltas['La'])<=TOL['energy_kcal_mol']})
    reused = {}
    if m['stage']=='benchmark':
        source, _ = qualified(verify(m['qualification']))
        for key, entry in m['reused'].items():
            reused[key] = {'source_collection':entry['source_collection'], 'source_task_id':entry['source_task_id'],
                           'result':source['rows'][entry['source_task_id']]}
    result = {'status':'complete' if complete else 'incomplete', 'manifest':record(mp), 'rows':rows,
            'attempts':attempts, 'checks':checks,
            'numerical_gate_pass':complete and bool(checks) and all(c['pass'] for c in checks), **UNAVAILABLE}
    if m['stage']=='benchmark':
        result['reused_rows'] = reused
        result['numerical_gate_pass'] = None
        result['qualified_source'] = m['qualification']
    return result


if __name__=='__main__':
    p = argparse.ArgumentParser(description=__doc__); sub = p.add_subparsers(dest='command',required=True)
    a = sub.add_parser('prepare-qualification')
    for k in ('source-inventory','software','agreement','output'):a.add_argument('--'+k,required=True)
    a = sub.add_parser('prepare-benchmark')
    for k in ('qualification','mechanics','agreement','output'):a.add_argument('--'+k,required=True)
    a = sub.add_parser('prepare-readout')
    for k in ('collection','agreement','output'):a.add_argument('--'+k,required=True)
    args = vars(p.parse_args()); command = args.pop('command')
    from mace_omol_readout import prepare as prepare_readout
    print(json.dumps({'prepare-qualification':prepare_qualification, 'prepare-benchmark':prepare_benchmark,
                      'prepare-readout':prepare_readout}[command](**args),indent=2))
