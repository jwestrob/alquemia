"""Allocated derivatives of the separately versioned shared-neutral descriptor."""
from __future__ import annotations
import importlib.metadata
import json
import math
import os
from pathlib import Path
import resource
import time
import traceback
import numpy as np
from affordable_common import InvalidArtifact, read_json, record, verify, write_new, xyz
from mace_hybrid import check_atoms
from mace_omol_neutral import COMPONENT, receipt as feature_receipt, install as install_feature, accepted_receipt
from mace_omol_gradients import ADAPTER, install, receipt


def worker(manifest,task_id,output,memory_mode):
    import torch
    from torch.utils.checkpoint import set_checkpoint_early_stop
    from ase import Atoms
    from mace.calculators import mace_omol
    from mace_omol import UNAVAILABLE, input_batch
    from mace_omol_readout import NativeCapture
    m=read_json(manifest);t=next(t for t in m['tasks'] if t['task_id']==task_id)
    out=Path(output);start=time.monotonic();grad=not t['energy_only']
    r={'task_id':task_id,'cache_key':t['cache_key'],'manifest':record(manifest),'status':'unavailable',
       'memory_mode':memory_mode,'energy_component':COMPONENT,**UNAVAILABLE}
    try:
        if not os.environ.get('SLURM_JOB_ID') or memory_mode!='native' or not torch.cuda.is_available():
            raise InvalidArtifact('gradient pilot requires its allocated GPU')
        torch.set_num_threads(int(os.environ['SLURM_CPUS_PER_TASK']));torch.set_default_dtype(torch.float64)
        torch.cuda.reset_peak_memory_stats();props=torch.cuda.get_device_properties(0)
        r.update(execution_device='cuda',device={'name':props.name,'total_memory_bytes':props.total_memory,
                                               'CUDA_VISIBLE_DEVICES':os.environ.get('CUDA_VISIBLE_DEVICES')})
        calc=mace_omol(model=str(verify(m['model']['checkpoint'])),device='cuda',default_dtype='float64')
        obj=calc.models[0]
        if type(obj).__name__!='ScaleShiftMACE' or dict(obj.embedding_specs)!=m['model']['embedding_specs']:
            raise InvalidArtifact('native model/state embedding changed')
        for parameter in obj.parameters():parameter.requires_grad_(False)
        coords=xyz(verify(t['xyz']))
        if check_atoms(coords,t['charge'])!=t['state'] or coords[t['metal_index']][0]!=t['metal']:
            raise InvalidArtifact('input state/selected metal differs')
        atoms=Atoms([a[0] for a in coords],positions=[a[1:] for a in coords],pbc=False)
        atoms.info.update(charge=t['charge'],spin=t['spin_multiplicity']);batch=calc._atoms_to_batch(atoms)
        r['input_state_check']=input_batch(calc,atoms,t['charge'],t['spin_multiplicity'],batch,t['metal_index'])
        parameters={n:p._version for n,p in obj.named_parameters()}
        buffers={n:(id(b),b._version) for n,b in obj.named_buffers()}
        install_feature(obj)
        if t['derivative_backend']=='checkpointed':install(obj)
        elif t['derivative_backend']!='native':raise InvalidArtifact('undeclared derivative backend')
        capture=NativeCapture(obj);r['model_load_seconds']=time.monotonic()-start
        if list(obj.heads)!=['omol'] or calc.energy_units_to_eV!=1.:
            raise InvalidArtifact('expected native single-head eV model')
        torch.cuda.synchronize();evaluation=time.monotonic()
        # TorchScript can wrap the private early-stop exception as RuntimeError.
        # Complete each recomputation instead; retain the same mathematical graph.
        with torch.set_grad_enabled(grad), set_checkpoint_early_stop(False):
            r['gradient_computation_enabled']=torch.is_grad_enabled()
            native=obj(batch.to_dict(),training=False,compute_force=grad,compute_virials=False,
                       compute_stress=False,compute_displacement=False,compute_hessian=False,
                       compute_edge_forces=False,compute_atomic_stresses=False)
            value=float(native['energy'].item())
            forces=native['forces'].detach().cpu().numpy().copy() if grad else None
            e0=obj.atomic_energies_fn(batch['node_attrs'])[:,0].detach().cpu().numpy()
        torch.cuda.synchronize();r['evaluation_seconds']=time.monotonic()-evaluation
        if not math.isfinite(value) or (grad and (forces.shape!=(len(coords),3) or not np.isfinite(forces).all())):
            raise InvalidArtifact('nonfinite/malformed derivative output')
        if not grad and native['forces'] is not None:raise InvalidArtifact('energy-only task unexpectedly returned forces')
        r['parameter_versions_unchanged']=parameters=={n:p._version for n,p in obj.named_parameters()}
        r['buffer_versions_unchanged']=buffers=={n:(id(b),b._version) for n,b in obj.named_buffers()}
        if not r['parameter_versions_unchanged'] or not r['buffer_versions_unchanged']:
            raise InvalidArtifact('model parameters or buffers changed')
        r.update(charge_feature_adapter=feature_receipt(obj),output_semantics='energy_like_descriptor_not_quantum_endpoint',
                 output_unit='eV_equivalent_model_units',energy_eV=value,energy_only=not grad,
                 native_readout=capture.save_native(value,out,e0),derivative_backend=t['derivative_backend'],
                 descriptor_gradient_experiment=t['descriptor_gradient_experiment'],
                 execution_adapter=receipt(obj,grad) if t['derivative_backend']=='checkpointed' else None,
                 parameter_gradients_disabled=True,checkpoint_early_stop=False,
                 coordinate_derivative_method='analytic_autograd' if grad else None)
        if grad:
            fp=out/'negative_descriptor_gradient_eV_per_A.npy';np.save(fp,forces)
            gp=out/'descriptor_gradient_eV_per_A.npy';np.save(gp,-forces)
            r.update(forces=record(fp),gradient=record(gp),force_definition='negative_Cartesian_gradient_of_shared_neutral_descriptor',
                     gradient_definition='Cartesian_gradient_of_shared_neutral_descriptor',gradient_unit='eV_equivalent_per_A',
                     analytic_gradient_status='computed_descriptor_only_not_validated_physical_force')
        else:r.update(forces=None,gradient=None,force_definition=None,gradient_definition=None,gradient_unit=None,
                      analytic_gradient_status='not_requested')
        r['status']='computed'
    except Exception as exc:
        r.update(status='failed',reason=str(exc),exception_type=type(exc).__name__);traceback.print_exc()
    finally:
        r.update(wall_seconds=time.monotonic()-start,peak_host_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                 peak_cuda_allocated_bytes=torch.cuda.max_memory_allocated() if torch.cuda.is_available() else None,
                 peak_cuda_reserved_bytes=torch.cuda.max_memory_reserved() if torch.cuda.is_available() else None,
                 versions={k:importlib.metadata.version(k) for k in ('mace-torch','torch','ase','e3nn')},
                 slurm_job_id=os.environ.get('SLURM_JOB_ID'),allocated_cpus=os.environ.get('SLURM_CPUS_PER_TASK'),
                 allocated_host_mem_MiB=os.environ.get('SLURM_MEM_PER_NODE'))
        write_new(out/'result.json',r)
    print(json.dumps({k:r.get(k) for k in ('task_id','status','reason','evaluation_seconds','peak_cuda_allocated_bytes')}),flush=True)
    return r


def accepted(result,task):
    from mace_omol_readout import checked_arrays
    grad=not task['energy_only'];s=result.get('input_state_check',{})
    if (result.get('descriptor_gradient_experiment')!=task['descriptor_gradient_experiment']
            or result.get('energy_component')!=COMPONENT or result.get('derivative_backend')!=task['derivative_backend']
            or result.get('output_unit')!='eV_equivalent_model_units'
            or result.get('output_semantics')!='energy_like_descriptor_not_quantum_endpoint'
            or result.get('gradient_computation_enabled')!=grad or result.get('energy_only')!=task['energy_only']
            or result.get('parameter_versions_unchanged') is not True or result.get('buffer_versions_unchanged') is not True
            or result.get('parameter_gradients_disabled') is not True
            or result.get('checkpoint_early_stop') is not False
            or s.get('charge')!=task['charge'] or s.get('spin_multiplicity')!=1 or s.get('head')!='omol'
            or s.get('selected_metal_index')!=task['metal_index'] or s.get('atoms')!=task['state']['atoms'] or s.get('graph_count')!=1):
        return False
    if not accepted_receipt(result.get('charge_feature_adapter'),task,s['atoms']):return False
    adapter=result.get('execution_adapter')
    if task['derivative_backend']=='checkpointed':
        if (not adapter or adapter.get('id')!=ADAPTER or adapter.get('chunk_size')!=1024
                or adapter.get('analytic_gradients_requested')!=grad or adapter.get('coordinate_graph_preserved') is not True
                or adapter.get('aggregation')!='scatter_add_expanded_receiver_index'
                or adapter.get('early_stop') is not False):return False
        for key in ('interactions','products'):
            if len(adapter.get(key,[]))!=3:return False
            for layer in adapter[key]:
                if layer['atoms']!=s['atoms'] or layer['kernel_invocations_including_recomputations']<layer['chunks']:return False
    elif adapter is not None:return False
    if grad:
        if (result.get('coordinate_derivative_method')!='analytic_autograd'
                or result.get('force_definition')!='negative_Cartesian_gradient_of_shared_neutral_descriptor'
                or result.get('gradient_unit')!='eV_equivalent_per_A'):return False
        f=np.load(verify(result['forces']));g=np.load(verify(result['gradient']))
        if f.shape!=(s['atoms'],3) or not np.isfinite(f).all() or not np.array_equal(-f,g):return False
    elif result.get('forces') is not None or result.get('gradient') is not None:return False
    checked_arrays(result,task)
    return math.isfinite(result['energy_eV'])
