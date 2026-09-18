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
             'mace_omol_readout.py', 'mace_omol_coordination.py', 'mace_omol_intact.py',
             'mace_omol_edges.py', 'mace_omol_edge_run.py', 'mace_omol_edge_report.py', 'mace_file_checks.py',
             'mace_omol_products.py', 'mace_omol_panel.py', 'mace_omol_panel_prepare.py',
             'mace_omol_intact_inventory.py', 'mace_omol_backbone_audit.py',
             'mace_omol_locality.py', 'mace_omol_spectator.py', 'mace_omol_ablation.py',
             'mace_omol_ablation_run.py', 'mace_omol_ablation_panel.py', 'mace_omol_panel_report.py',
             'mace_omol_prepared.py', 'mace_omol_factorization.py', 'mace_omol_mask_calibration.py',
             'mace_omol_source_prepare.py', 'mace_omol_multisite.py',
             'mace_omol_gradients.py', 'mace_omol_gradient_run.py', 'mace_omol_gradient_worker.py', 'mace_omol_response.py',
             'mace_omol_neutral.py', 'mace_omol_neutral_worker.py', 'mace_omol_neutral_run.py',
             'mace_omol_vacuum_hybrid.py', 'mace_omol_vacuum.py', 'mace_omol_matched_h.py',
             'mace_omol_hybrid_transfer.py','mace_omol_hybrid_transfer_minimum.py','mace_omol_hybrid_response.py','mace_omol_hybrid_minimum.py','mace_metal_response.py','mace_metal_minimum.py','mace_bounded_response.py')
    if stage == 'coupled_donor_path':
        names += ('mace_site_path.py', 'mace_site_kinematics.py')
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
    if read_json(manifest).get('stage')=='coupled_donor_path':
        from mace_site_path import validate as path_validate
        return path_validate(manifest)
    if read_json(manifest).get('stage')=='matched_hybrid_GGR_transfer':
        from mace_omol_hybrid_transfer import validate as validate_transfer
        return validate_transfer(manifest)
    if read_json(manifest).get('stage')=='matched_hybrid_response_validation':
        from mace_omol_hybrid_minimum import validate as validate_response
        return validate_response(manifest)
    if read_json(manifest).get('stage')=='matched_hybrid_response':
        from mace_omol_hybrid_response import validate as validate_response
        return validate_response(manifest)
    if read_json(manifest).get('stage')=='matched_H_core':
        from mace_omol_matched_h import validate_mace
        return validate_mace(manifest)
    if read_json(manifest).get('stage')=='vacuum_context':
        from mace_omol_vacuum_hybrid import validate as validate_vacuum_context
        return validate_vacuum_context(manifest)
    if read_json(manifest).get('stage')=='shared_neutral_core':
        from mace_omol_neutral_run import validate as neutral_validate
        return neutral_validate(manifest)
    if read_json(manifest).get('stage')=='masked_core_response':
        from mace_omol_response import validate as validate_response
        return validate_response(manifest)
    if read_json(manifest).get('stage') in ('masked_gradient_core','masked_gradient_full'):
        from mace_omol_gradient_run import validate as validate_gradients
        return validate_gradients(manifest)
    m = read_json(manifest)
    if m.get('stage')=='ablation_prepared':
        from mace_omol_prepared import validate as validate_prepared
        return validate_prepared(manifest)
    if m.get('stage')=='ablation_canonical':
        from mace_omol_ablation_panel import validate as validate_ablation_panel
        return validate_ablation_panel(manifest)
    if m.get('stage')=='ablation_development':
        from mace_omol_ablation_run import validate as validate_ablation
        return validate_ablation(manifest)
    if m.get('stage')=='spectator_intact':
        from mace_omol_spectator import validate as validate_spectator
        return validate_spectator(manifest)
    if m.get('stage','').startswith('panel_'):
        from mace_omol_panel import validate as validate_panel
        return validate_panel(manifest)
    if m.get('stage','').startswith(('edge_','cpu_','product_')):
        from mace_omol_edge_run import validate as validate_edge
        return validate_edge(manifest)
    if m.get('stage','').startswith('intact_'):
        from mace_omol_intact import validate as validate_intact
        return validate_intact(manifest)
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


def input_batch(calc, atoms, charge, multiplicity, batch=None, metal_index=0):
    batch = calc._atoms_to_batch(atoms) if batch is None else batch
    q = float(batch['total_charge'].item()); spin = float(batch['total_spin'].item())
    if q != charge or spin != multiplicity or calc.head != 'omol':
        raise InvalidArtifact('actual OMOL batch charge/multiplicity/head differs')
    for name, value in (('total_charge', q), ('total_spin', spin)):
        spec = calc.models[0].embedding_specs[name]; index = value+spec.get('offset', 0)
        if spec['type'] != 'categorical' or value != int(value) or not 0 <= index < spec['num_classes']:
            raise InvalidArtifact('state outside checkpoint categorical embedding')
    edges=batch['edge_index']
    return {'charge':q, 'spin_multiplicity':spin, 'head':calc.head,
            'selected_metal_index':metal_index,
            'metal_neighbor_edge_count':int(((edges[0]==metal_index)|(edges[1]==metal_index)).sum().item()),
            'atoms':len(atoms), 'graph_count':int(batch['ptr'].numel()-1), 'charge_prediction_claimed':False}


def worker(manifest, task_id, output, memory_mode):
    if read_json(manifest).get('stage')=='shared_neutral_core':
        from mace_omol_neutral_worker import worker as neutral_worker
        return neutral_worker(manifest,task_id,output,memory_mode)
    if read_json(manifest).get('stage') in ('coupled_donor_path','masked_gradient_core','masked_gradient_full','masked_core_response','matched_hybrid_GGR_transfer','matched_hybrid_response','matched_hybrid_response_validation'):
        from mace_omol_gradient_worker import worker as gradient_worker
        return gradient_worker(manifest,task_id,output,memory_mode)
    import torch
    from ase import Atoms
    from mace.calculators import mace_omol
    m = read_json(manifest); t = next(t for t in m['tasks'] if t['task_id']==task_id)
    device=t.get('execution_device','cuda')
    out = Path(output); start = time.monotonic()
    result = {'task_id':task_id, 'cache_key':t['cache_key'], 'manifest':record(manifest),
              'status':'unavailable', 'memory_mode':memory_mode, 'energy_component':t['energy_component'], **UNAVAILABLE}
    try:
        if (not os.environ.get('SLURM_JOB_ID') or memory_mode!='native' or device not in ('cuda','cpu')
                or (device=='cuda' and not torch.cuda.is_available())):
            raise InvalidArtifact('OMOL requires its declared allocated device')
        torch.set_num_threads(int(os.environ['SLURM_CPUS_PER_TASK'])); torch.set_default_dtype(torch.float64)
        result['execution_device']=device
        if device=='cuda':
            torch.cuda.reset_peak_memory_stats(); props = torch.cuda.get_device_properties(0)
            result['device'] = {'name':props.name, 'total_memory_bytes':props.total_memory,
                                'CUDA_VISIBLE_DEVICES':os.environ.get('CUDA_VISIBLE_DEVICES')}
        else:
            import platform
            result['device']={'name':platform.processor(),'kind':'CPU','threads':torch.get_num_threads(),
                              'node':platform.node(),'cpuinfo':Path('/proc/cpuinfo').read_text().split('\n\n')[0]}
        calc = mace_omol(model=str(verify(m['model']['checkpoint'])), device=device, default_dtype='float64')
        model_obj = calc.models[0]
        if type(model_obj).__name__ != 'ScaleShiftMACE' or dict(model_obj.embedding_specs) != m['model']['embedding_specs']:
            raise InvalidArtifact('loaded model type/state embedding mismatch')
        rows = xyz(verify(t['xyz']))
        if not {Z[a[0]] for a in rows} <= set(model_obj.atomic_numbers.tolist()):
            raise InvalidArtifact('checkpoint lacks a required element')
        atoms = Atoms([a[0] for a in rows], positions=[a[1:] for a in rows], pbc=False)
        atoms.info.update(charge=t['charge'], spin=t['spin_multiplicity']); atoms.calc = calc
        native_batch=calc._atoms_to_batch(atoms) if t.get('energy_only') else None
        metal_index=t.get('metal_index',0)
        if rows[metal_index][0]!=t['metal']:
            raise InvalidArtifact('selected metal source index has the wrong element')
        result['input_state_check'] = input_batch(calc, atoms, t['charge'], t['spin_multiplicity'],native_batch,metal_index)
        if result['input_state_check']['graph_count'] != 1:
            raise InvalidArtifact('unexpected graph batching')
        if t.get('require_isolated_metal') and result['input_state_check']['metal_neighbor_edge_count'] != 0:
            raise InvalidArtifact('detached metal still has graph neighbor edges')
        if 'spectator_index' in t:
            i=t['spectator_index'];edges=native_batch['edge_index']
            count=int(((edges[0]==i)|(edges[1]==i)).sum().item())
            result['input_state_check'].update(spectator_index=i,spectator_neighbor_edge_count=count)
            if rows[i][0]!='Na' or count:
                raise InvalidArtifact('declared sodium spectator has graph interactions or wrong identity')
        versions = {name:p._version for name,p in model_obj.named_parameters()}
        buffer_versions={name:(id(b),b._version) for name,b in model_obj.named_buffers()}
        if t.get('charge_feature_adapter'):
            from mace_omol_ablation import install as install_ablation, ADAPTER as ABLATION, COMPONENT as DESCRIPTOR
            if (t['charge_feature_adapter']!=ABLATION or t['energy_component']!=DESCRIPTOR
                    or not t.get('energy_only') or m['model'].get('charge_feature_adapter')!=ABLATION):
                raise InvalidArtifact('unsupported descriptor feature adapter')
            install_ablation(model_obj)
        if t.get('edge_adapter'):
            from mace_omol_edges import install, ADAPTER
            from mace_omol_products import install as install_products, ADAPTER as PRODUCT_ADAPTER
            if not t.get('energy_only'):
                raise InvalidArtifact('edge adapter supports energy only')
            if t['edge_adapter']['id']==ADAPTER:
                install(model_obj,t['edge_adapter']['chunk_size'])
            elif t['edge_adapter']['id']==PRODUCT_ADAPTER:
                install_products(model_obj,t['edge_adapter']['chunk_size'],t['edge_adapter']['product_chunk_size'])
            else:raise InvalidArtifact('unsupported exact execution adapter')
        capture = None
        if t.get('capture_native_readout'):
            from mace_omol_readout import NativeCapture
            capture = NativeCapture(model_obj)
        result['model_load_seconds'] = time.monotonic()-start
        if device=='cuda':torch.cuda.synchronize()
        evaluation = time.monotonic()
        if t.get('energy_only'):
            if capture is None or list(model_obj.heads)!=['omol'] or calc.energy_units_to_eV!=1.:
                raise InvalidArtifact('energy-only path requires captured native single-head eV output')
            with torch.no_grad():
                result['gradient_computation_enabled']=torch.is_grad_enabled()
                native=model_obj(native_batch.to_dict(),training=False,compute_force=False,
                                 compute_virials=False,compute_stress=False,compute_displacement=False,
                                 compute_hessian=False,compute_edge_forces=False,compute_atomic_stresses=False)
                value=float(native['energy'].item())
                e0=model_obj.atomic_energies_fn(native_batch['node_attrs'])[:,0].detach().cpu().numpy()
            if native['forces'] is not None:
                raise InvalidArtifact('energy-only native call unexpectedly returned forces')
            forces=None
        else:
            value = float(atoms.get_potential_energy()); forces = np.asarray(atoms.get_forces(), dtype=np.float64)
        if device=='cuda':torch.cuda.synchronize()
        result['evaluation_seconds'] = time.monotonic()-evaluation
        if not math.isfinite(value) or (forces is not None and (forces.shape != (len(atoms),3) or not np.isfinite(forces).all())):
            raise InvalidArtifact('invalid OMOL energy/analytic forces')
        result['parameter_versions_unchanged'] = versions == {name:p._version for name,p in model_obj.named_parameters()}
        if not result['parameter_versions_unchanged']:
            raise InvalidArtifact('model parameters mutated during inference')
        if t.get('charge_feature_adapter'):
            from mace_omol_ablation import receipt as ablation_receipt
            result['charge_feature_adapter']=ablation_receipt(model_obj)
            result['output_semantics']='energy_like_descriptor_not_quantum_endpoint'
            result['output_unit']='eV_equivalent_model_units'
        if t.get('edge_adapter'):
            from mace_omol_edges import receipt, ADAPTER
            if t['edge_adapter']['id']==ADAPTER:
                result['execution_adapter']=receipt(model_obj,t['edge_adapter']['chunk_size'])
            else:
                from mace_omol_products import receipt as product_receipt
                result['execution_adapter']=product_receipt(model_obj,t['edge_adapter']['chunk_size'],t['edge_adapter']['product_chunk_size'])
            result['buffer_versions_unchanged']=buffer_versions=={name:(id(b),b._version) for name,b in model_obj.named_buffers()}
            if not result['buffer_versions_unchanged']:
                raise InvalidArtifact('model buffers mutated during edge inference')
        if t.get('energy_only'):
            result.update(status='computed',energy_eV=value,forces=None,charge_check=None,force_definition=None,
                          energy_only=True,analytic_gradient_status='not_requested',
                          native_readout=capture.save_native(value,out,e0))
        else:
            fp = out/'forces_eV_A.npy'; np.save(fp, forces)
            result.update(status='computed', energy_eV=value, forces=record(fp), charge_check=None,
                          force_definition='negative_Cartesian_gradient_of_total_vacuum_OMOL_energy')
        if capture is not None and not t.get('energy_only'):
            result['native_readout'] = capture.save(calc, value, out)
    except Exception as exc:
        result.update(status='failed', reason=str(exc), exception_type=type(exc).__name__); traceback.print_exc()
    finally:
        result.update(wall_seconds=time.monotonic()-start, peak_host_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                      peak_cuda_allocated_bytes=torch.cuda.max_memory_allocated() if device=='cuda' and torch.cuda.is_available() else None,
                      peak_cuda_reserved_bytes=torch.cuda.max_memory_reserved() if device=='cuda' and torch.cuda.is_available() else None,
                      slurm_job_id=os.environ.get('SLURM_JOB_ID'), allocated_cpus=os.environ.get('SLURM_CPUS_PER_TASK'),
                      allocated_host_mem_MiB=os.environ.get('SLURM_MEM_PER_NODE'),
                      versions={n:importlib.metadata.version(n) for n in ('torch','mace-torch','ase','e3nn')})
        write_new(out/'result.json', result)
    print(json.dumps({k:result.get(k) for k in ('task_id','status','reason','evaluation_seconds')}), flush=True)
    return result


def accepted_state(result, task):
    if task.get('neutral_feature_experiment'):
        from mace_omol_neutral_worker import accepted
        return accepted(result,task)
    if task.get('descriptor_gradient_experiment'):
        from mace_omol_gradient_worker import accepted
        return accepted(result,task)
    state = result.get('input_state_check', {})
    expected_component=COMPONENT
    if task.get('charge_feature_adapter'):
        from mace_omol_ablation import ADAPTER as ABLATION, COMPONENT as DESCRIPTOR
        expected_component=DESCRIPTOR
        if (task['charge_feature_adapter']!=ABLATION or not task.get('energy_only')
                or result.get('charge_feature_adapter')!={'id':ABLATION,'calls':1,'rows':state.get('atoms'),
                    'output_max_abs':0.,'physical_charge_input_preserved':True,'spin_embedding_modified':False,
                    'output_semantics':'energy_like_descriptor_not_quantum_endpoint'}
                or result.get('output_semantics')!='energy_like_descriptor_not_quantum_endpoint'
                or result.get('output_unit')!='eV_equivalent_model_units'):
            return False
    elif result.get('charge_feature_adapter') is not None:return False
    if task.get('energy_component')!=expected_component:return False
    if result.get('execution_device','cuda')!=task.get('execution_device','cuda'):return False
    adapter=result.get('execution_adapter')
    if task.get('edge_adapter'):
        wanted=task['edge_adapter']
        if (not adapter or any(adapter.get(k)!=v for k,v in wanted.items())
                or result.get('buffer_versions_unchanged') is not True or adapter.get('gradient_support') is not False
                or len(adapter.get('layers',[]))!=3):return False
        if any(layer['atoms']!=state.get('atoms') or layer['edges']!=layer['processed_edges']
               or layer['chunk_size']!=wanted['chunk_size'] for layer in adapter['layers']):return False
        if 'product_chunk_size' in wanted:
            from mace_omol_products import ADAPTER as PRODUCT_ADAPTER
            layers=adapter.get('product_layers',[])
            if wanted['id']!=PRODUCT_ADAPTER or len(layers)!=3:return False
            if any(layer['atoms']!=state.get('atoms') or layer['atoms']!=layer['processed_atoms']
                   or layer['chunk_size']!=wanted['product_chunk_size']
                   or layer['batches']!=(layer['atoms']+layer['chunk_size']-1)//layer['chunk_size']
                   for layer in layers):return False
        elif adapter.get('product_layers') is not None:return False
    elif adapter is not None:return False
    if state.get('selected_metal_index',0)!=task.get('metal_index',0):
        return False
    if task.get('energy_only'):
        if (result.get('energy_only') is not True or result.get('forces') is not None
                or result.get('force_definition') is not None or result.get('analytic_gradient_status')!='not_requested'
                or result.get('gradient_computation_enabled') is not False):
            return False
    elif result.get('energy_only'):
        return False
    if task.get('require_isolated_metal') and state.get('metal_neighbor_edge_count') != 0:
        return False
    if 'spectator_index' in task and (state.get('spectator_index')!=task['spectator_index']
                                    or state.get('spectator_neighbor_edge_count')!=0):
        return False
    if task.get('capture_native_readout'):
        from mace_omol_readout import checked_arrays
        checked_arrays(result, task)
    return (result.get('energy_component')==expected_component and result.get('density_coefficients') is None
            and result.get('parameter_versions_unchanged') is True and state.get('charge')==task['charge']
            and state.get('spin_multiplicity')==task['spin_multiplicity'] and state.get('head')=='omol'
            and state.get('graph_count')==1 and state.get('charge_prediction_claimed') is False
            and math.isfinite(result['energy_eV']))


def collect(manifest):
    if read_json(manifest).get('stage')=='coupled_donor_path':
        from mace_site_path import collect as path_collect
        return path_collect(manifest)
    if read_json(manifest).get('stage')=='matched_hybrid_GGR_transfer':
        from mace_omol_hybrid_transfer import collect as collect_transfer
        return collect_transfer(manifest)
    if read_json(manifest).get('stage') in ('matched_hybrid_response','matched_hybrid_response_validation'):
        from mace_omol_hybrid_response import collect as collect_response
        return collect_response(manifest)
    if read_json(manifest).get('stage')=='matched_H_core':
        from mace_omol_matched_h import collect_mace
        return collect_mace(manifest)
    if read_json(manifest).get('stage')=='vacuum_context':
        from mace_omol_vacuum_hybrid import collect as collect_vacuum_context
        return collect_vacuum_context(manifest)
    if read_json(manifest).get('stage')=='shared_neutral_core':
        from mace_omol_neutral_run import collect as neutral_collect
        return neutral_collect(manifest)
    if read_json(manifest).get('stage')=='masked_core_response':
        from mace_omol_response import collect as collect_response
        return collect_response(manifest)
    if read_json(manifest).get('stage') in ('masked_gradient_core','masked_gradient_full'):
        from mace_omol_gradient_run import collect as collect_gradients
        return collect_gradients(manifest)
    if read_json(manifest).get('stage')=='ablation_prepared':
        from mace_omol_prepared import collect as collect_prepared
        return collect_prepared(manifest)
    if read_json(manifest).get('stage')=='ablation_canonical':
        from mace_omol_ablation_panel import collect as collect_ablation_panel
        return collect_ablation_panel(manifest)
    if read_json(manifest).get('stage')=='ablation_development':
        from mace_omol_ablation_run import collect as collect_ablation
        return collect_ablation(manifest)
    if read_json(manifest).get('stage')=='spectator_intact':
        from mace_omol_spectator import collect as collect_spectator
        return collect_spectator(manifest)
    if read_json(manifest).get('stage','').startswith('panel_'):
        from mace_omol_panel import collect as collect_panel
        return collect_panel(manifest)
    if read_json(manifest).get('stage','').startswith(('edge_','cpu_','product_')):
        from mace_omol_edge_run import collect as collect_edge
        return collect_edge(manifest)
    if read_json(manifest).get('stage','').startswith('intact_'):
        from mace_omol_intact import collect as collect_intact
        return collect_intact(manifest)
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
