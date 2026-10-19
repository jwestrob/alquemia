#!/usr/bin/env python3
"""V2: exact admitted component authority; one regional forward and collection."""
import argparse
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import resource
import socket
import sys
import time
import traceback
from types import SimpleNamespace

PROTOCOL = 'lbt_whole_GFN_regional_vacuum_delta_v1'
CPU_SHA = '46c8d49c8a0fd25adc52853ce8a4bd218c42a331070a05cf7a72a2c864df7bdf'
MACE_SHA = '686eba84fd3dfcf28724d2d4e7957528c987445d7137b18c88a07dfec055f661'
AUTHORITY_SHA = '00add2344edda82078053da5572349b2ae74f6e2052e7b657499fdcb63093d1b'
V1_SHA = '05414274fc4139d7e274f56f1bf513f76c20bed696397890a6eeee99b652dc2c'
IDS = tuple(s+'_'+t for s in ('7CCO', '7CCN') for t in ('before_La', 'after_La', 'before_Dy', 'after_Dy'))
MISSING = '7CCO_after_Dy'
COMPONENTS = ('GFN_full', 'MACE_region', 'negative_GFN_region', 'fixed_FF_ALPB')
HA = 627.509474
EV = 23.06054783061903
BOHR = .529177210903


def require(ok, message):
    if not ok: raise ValueError(message)


def pin(path):
    p = Path(path).resolve()
    with p.open('rb') as stream: sha = hashlib.file_digest(stream, 'sha256').hexdigest()
    return {'path': str(p), 'sha256': sha}


def verify(record):
    require(pin(record['path'])['sha256'] == record['sha256'], 'Changed pin: '+record['path'])
    return Path(record['path'])


def read(record):
    return json.loads(verify(record).read_text())


def write(path, value):
    with Path(path).open('x') as stream:
        json.dump(value, stream, indent=2, allow_nan=False); stream.write('\n')
    return pin(path)


def module(name, record):
    spec = importlib.util.spec_from_file_location(name, verify(record))
    obj = importlib.util.module_from_spec(spec); spec.loader.exec_module(obj)
    return obj


def binding(cfg):
    return {k: cfg[k]['sha256'] if k == 'implementation' else cfg[k] for k in
            ('protocol', 'implementation', 'specification', 'design', 'inventory',
             'cpu_implementation', 'cpu_config', 'cpu_result', 'hybrid_implementation',
             'preflight_result', 'scout_result', 'authority')}


def signatures(manifest):
    return [{k: v for k, v in r.items() if k not in ('whole_xyz', 'region_xyz', 'mapping_replay')}
            | {'whole_xyz_sha256': r['whole_xyz']['sha256'], 'region_xyz_sha256': r['region_xyz']['sha256']}
            for r in manifest['cells']]


def saved_authority(cfg, cpu, cc, parent, out, result):
    """Admit saved inputs and components without reconstructing archived cap sums."""
    require(cfg['authority']['sha256'] == AUTHORITY_SHA, 'Frozen recovery authority differs')
    authority = read(cfg['authority'])
    base = module('whole_base_preparation_frozen', cc['base_implementation'])
    base.verify_tree(cfg); base.verify_tree(authority)
    require(authority['original_implementation']['sha256'] == V1_SHA
            and authority['CPU_result'] == cfg['cpu_result']
            and authority['preflight_result'] == cfg['preflight_result']
            and authority['CPU_manifest'] == parent['manifest'], 'Recovery parent pins differ')
    failed = read(authority['failed_GPU_result']); fc = read(failed['config'])
    require(failed['implementation']['sha256'] == V1_SHA and failed['status'] == 'execution_failure'
            and failed['new_molecular_calls'] == failed['regional_MACE_calls_attempted'] == 0
            and not failed['model_loaded'] and failed['failure'] == "ValueError('Eight prepared parent signatures differ')"
            and fc == read(authority['original_GPU_config']) and fc['cpu_result'] == cfg['cpu_result'],
            'Original zero-call failure identity differs')
    bc = read(cc['base_config']); pf = read(cfg['preflight_result'])
    scout = read(cfg['scout_result']); sc = read(scout['config'])
    require(bc['implementation']['sha256'] == cpu.BASE_SHA and bc['protocol'] == base.PROTOCOL
            and bc['limits'] == base.LIMITS and bc['inventory']['sha256'] == base.INPUT_SHA
            and bc['inventory'] == cfg['inventory'] and bc['design'] == cfg['design'], 'Original scientific declaration differs')
    base.verify_tree(bc)
    require(pin(sys.executable)['sha256'] == bc['python']['sha256']
            and str(Path(sys.prefix).resolve()) == bc['python_prefix'], 'Bound Python environment differs')
    require(pf['stage'] == 'preflight' and pf['status'] == 'preflight_passed' and not pf['model_loaded']
            and pf['scientific_binding'] == base.binding(bc) and pf['new_molecular_calls'] == 0
            and pf['manifest'] == authority['preflight_manifest'], 'Actual original preflight authority unavailable')
    require(scout['protocol'] == base.PROTOCOL and scout['stage'] == 'run-scout' and scout['status'] == 'scout_passed'
            and scout['scientific_binding'] == base.binding(bc) and sc['preflight_result'] == cfg['preflight_result']
            and scout['scout_native_admitted'] and scout['scout_cost_gate_passed']
            and scout['scout_wall_seconds'] <= base.LIMITS['scout_wall_seconds']
            and scout['new_molecular_calls'] == scout['whole_GFN_calls_attempted'] == 1
            and not scout['model_loaded'], 'Actual one-call scout admission/cost unavailable')
    pm = read(authority['preflight_manifest']); cm = read(authority['CPU_manifest'])
    fm = read(authority['failed_GPU_manifest'])
    require(signatures(pm) == signatures(cm) and [r['cell_id'] for r in cm['cells']] == list(IDS)
            and cm['scientific_binding'] == pm['scientific_binding'] == base.binding(bc), 'Exact admitted manifest authority differs')
    # Validate the demonstrated failure category; this is not a relaxed replay gate.
    def raw_signature(m):
        rows = json.loads(json.dumps(signatures(m)))
        for r in rows:
            if r['local_components'] is not None:
                r['local_components'].pop('gradient_difference_kcal_mol_A')
                r['local_components'].pop('component_gradients_kcal_mol_A')
        return rows
    require(raw_signature(pm) == raw_signature(fm) and signatures(pm) != signatures(fm)
            and {k:v for k,v in pm.items() if k != 'cells'} == {k:v for k,v in fm.items() if k != 'cells'},
            'Failure was not confined to the declared derived component representation')
    inv = read(cfg['inventory']); pins = inv['helpers']; hc = read(inv['parents']['qualified_hybrid_config'])
    base.verify_tree(inv)
    require([r['cell_id'] for r in inv['cells']] == list(IDS) and all(hc[k] == v for k,v in pins.items()),
            'Exact eight-state/helper inventory differs')
    runtime = pins['GFN_runtime']; profile = read(pins['GFN_profile'])
    require(pins['GFN_profile']['sha256'] == base.PROFILE_SHA and runtime['parameters']['sha256'] == base.PARAM_SHA,
            'Corrected full La/Dy parameter profile differs')
    import numpy as np
    sys.path.insert(0, str(verify(pins['common_helper']).parent))
    common = module('authority_common', pins['common_helper'])
    gradient = module('authority_gradient', pins['gradient_helper'])
    mapper = module('authority_mapper', pins['preparation_implementation'])
    mapper.np = np; mapper.cap_jacobians = gradient.cap_jacobians
    admit = module('authority_native_admission', pins['admission'])
    backend = module('authority_native_backend', pins['GFN_backend'])
    solvent = module('authority_solvent_state', pins['ALPB_helper'])
    whole = module('authority_whole_reference', bc['whole_reference_helper'])
    require(common.HA_TO_KCAL == HA and common.BOHR_TO_A == BOHR
            and read(pins['GFN_parent_config'])['runtime'] == runtime
            and read(pins['GFN_parent_config'])['tolerances'] == admit.TOL, 'Original units/native gates differ')
    qual = read(pins['GFN_qualification']); qm = read(qual['manifest']); qc = read(qual['config'])
    require(qual['status'] == 'passed' and qual['slurm_job_id'] == '12488' and qual['accepted'] == qual['attempted'] == 8
            and qm['profile'] == pins['GFN_profile'] and all(r['passed'] for r in qual['parity']+qual['solvent_component_parity'])
            and qc['staging']['standalone_settings'] == {'accuracy': .0001, 'maxiter': 500,
                'electronic_temperature_K': 300, 'alpb_reference': 'gsolv'}, 'Qualified native lineage unavailable')
    for r in qual['rows']:
        if r['task']['backend'] == 'standalone':
            require(r['status'] == 'accepted' and all(r['task'][k] == v for k,v in runtime.items()), 'Qualified runtime differs')
    alpb_parent = read(pins['ALPB_evidence']); mace_parent = read(pins['MACE_evidence']); mcfg = read(pins['MACE_config'])
    require(alpb_parent['status'] == 'complete' and alpb_parent['build']['binary'] == pins['ALPB_binary'], 'Fixed-FF solvent lineage differs')
    base.verify_tree(alpb_parent['build']); base.verify_tree(mcfg['helper_pins'])
    require(mace_parent['status'] == 'complete' and mace_parent['model_parameters_unchanged']
            and mace_parent['model'] == pins['model'] == mcfg['model']
            and mace_parent['versions'] == pins['expected_MACE_versions'], 'Native MACE recipe differs')
    h = dict(np=np, pins=pins, runtime=runtime, profile=profile, common=common, gradient=gradient,
             mapper=mapper, admit=admit, backend=backend, solvent=solvent, whole=whole)
    recovery = module('authority_saved_certificate', bc['recovery_implementation'])
    rcfg = read(bc['recovery_config']); partial = read(rcfg['original_partial_result'])
    cert_row = next(r for r in inv['cells'] if r['cell_id'] == '7CCN_after_La'); cert = read(cert_row['solver_exception_certificate'])
    evidence, _, _ = recovery.certify_scc(rcfg, out, partial, cert['original_failed_slot'])
    require(evidence == cert['saved_scc_qualification'] and cert['same_Hamiltonian'] and not cert['identical_solver_protocol']
            and cert['old_failure_remains_unavailable'] and evidence['primary'] == cert_row['local_isolated_GFN2']
            and cm['saved_solver_certificate'] == cert_row['solver_exception_certificate'], 'Certified local solver exception differs')
    cells = {}; rows = {r['cell_id']:r for r in cm['cells']}
    for row in inv['cells']:
        cid = row['cell_id']; m = rows[cid]; metal = row['metal']; mult = 1 if metal == 'La' else 6
        require(all(m[k] == row[k] for k in ('cell_id','source_id','geometry','metal','physical','positions','mapping')),
                'Admitted coordinate/map/input declaration differs')
        base.verify_tree(m)
        atoms = read(row['physical']); ids = [a['id'] for a in atoms]; x = base.shape([a['xyz_A'] for a in atoms],256,np)
        xp = verify(row['positions']); saved_x = np.load(xp,allow_pickle=False) if xp.suffix == '.npy' else read(row['positions'])['positions_A']
        mapping = read(row['mapping']); symbols = [a['element'] for a in atoms[:-1]]+[metal]
        require(np.array_equal(x,saved_x) and ids == m['physical_ids'] and len(set(ids)) == 256 and ids[-1] == 'metal'
                and set(symbols[:-1]) <= {'H','C','N','O'} and atoms[-1]['element'] == 'M', 'Exact full physical input differs')
        require(mapping['physical_atoms'] == atoms and len(mapping['core']) == 172 and len(mapping['caps']) == 2
                and len(mapping['field']) == 84 and mapping['core'][0]['physical_index'] == 255, 'Exact archived map differs')
        stask = row.get('source_task',row.get('source_task_template_La')); base.verify_tree(stask)
        require(stask['charge'] == -2 and stask['region'] == 'baseline', 'Baseline target/charge differs')
        charges = read(stask['global_solvent_charges'])
        require(charges['physical_ids'] == ids and charges['total_e'] == -2
                and abs(math.fsum(charges['charges_e'])+2) <= admit.TOL['charge_closure_e'], 'Fixed FF charges differ')
        syq = [metal]+[a['element'] for a in mapping['core'][1:]]; xq = base.shape([a['xyz_A'] for a in mapping['core']],172,np)
        for kind,sy,n,xx,ne in (('whole',symbols,256,x,750),('region',syq,172,xq,538)):
            nsy,nx = admit.xyz(m[kind+'_xyz'])
            require(nsy == sy and np.array_equal(nx,xx), 'Exact saved native XYZ differs')
            electrons = sum(sum(profile['element'][s]['refocc']) for s in sy)+2
            all_e = sum(admit.Z[s] for s in sy)+2
            expected = {'atoms':n,'charge':-2,'GFN_uhf':0,'GFN_effective_multiplicity':1,'GFN_valence_electrons':ne,
                'physical_all_electrons':all_e,'MACE_physical_multiplicity':mult,'physical_state_qualified':False}
            require(electrons == ne and (all_e-mult+1)%2 == 0 and m['state'][kind] == expected, 'Saved GFN/physical state differs')
        cell = dict(row=row,atoms=atoms,ids=ids,x=x,symbols=symbols,mapping=mapping,core_symbols=syq,core_x=xq,
                    state=m['state'],charges=charges,physical_parent=row['physical'])
        cell['solvent'] = whole.admit_alpb(cell,{'ALPB_binary':pins['ALPB_binary']},np,common)
        require(cell['solvent']['status'] == 'available' and cell['solvent'] == m['ALPB'], 'Exact original solvent admission differs')
        local = m['local_components']
        if cid != MISSING:
            point = read(row['local_POINT']); ml = read(row['local_MACE']); nt = read(row['task'])
            require(point['status'] == 'complete' and point['mapping'] == row['mapping'] and point['task'] == row['task']
                    and point['component_receipts']['MACE'] == row['local_MACE']
                    and point['component_receipts']['GFN_isolated'] == row['local_isolated_GFN2']
                    and point['point']['positions_A'] == x.tolist() and nt['physical_source']['sha256'] == row['physical']['sha256']
                    and nt['metal'] == metal and nt['charge'] == -2 and nt['multiplicity'] == mult, 'Cached local POINT/task differs')
            nsy,nx = admit.xyz(nt['xyz'])
            require(nsy == syq and np.array_equal(nx,xq) and ml['task'] == row['task'] and ml['xyz'] == nt['xyz']
                    and ml['physical_source'] == nt['physical_source'], 'Regional MACE native input differs')
            solvent.state(ml['input_state_check'],172,-2,mult,0); base.shape(ml['forces_eV_A'],172,np)
            require(ml['model'] == pins['model'] and ml['force_convention'] == 'negative vacuum energy derivative'
                    and math.isfinite(ml['energy_eV']), 'Regional MACE recipe/energy differs')
            native = base.native_replay(row['local_isolated_GFN2'],{'source_id':row['source_id'],'region':'baseline',
                'metal':metal,'kind':'isolated','atom_count':172,'charge':-2,'uhf':0,'multiplicity':1,
                'physical_multiplicity':mult,'valence_electrons':538,'pointcharge_count':0,'field_charge_e':0.},h,cid=='7CCN_after_La')
            nst = read(native['task']); base.shape(native['qm_gradient']['hartree_per_bohr'],172,np)
            require(nst['xyz']['sha256'] == nt['xyz']['sha256'] and nst['physical_source'] == nt['physical_source']
                    and local['MACE'] == row['local_MACE'] and local['GFN'] == row['local_isolated_GFN2']
                    and local['mapping'] == row['mapping'] and local['solver_exception_certificate'] == row['solver_exception_certificate']
                    and ml['energy_eV'] == local['MACE_eV'] == point['point']['native_energies']['MACE_eV']
                    and native['energy_hartree'] == local['GFN_hartree'] == point['point']['native_energies']['GFN_isolated_Eh'],
                    'Admitted local component/raw provenance differs')
            for a in local['component_gradients_kcal_mol_A'].values(): base.shape(a,256,np)
            base.shape(local['gradient_difference_kcal_mol_A'],256,np)
            require(set(local['component_gradients_kcal_mol_A']) == {'MACE_region','negative_GFN_region'}
                    and local['old_hybrid_gradient_not_used'], 'Stored regional component schema differs')
        else:
            parent_point = read(row['map_parent_POINT'])
            require(local is None and row['local_MACE'] is row['local_isolated_GFN2'] is row['local_POINT'] is None
                    and parent_point['status'] == 'complete' and parent_point['mapping'] == row['mapping']
                    and parent_point['point']['positions_A'] == x.tolist() and parent_point['state']['metal'] == 'La',
                    'Missing after-Dy original map parent differs')
        cell['local'] = local; cells[cid] = cell
    for sid in ('7CCO','7CCN'):
        for geo in ('before','after'):
            a,b = cells[sid+'_'+geo+'_La'],cells[sid+'_'+geo+'_Dy']
            require(a['ids'] == b['ids'] and np.array_equal(a['x'],b['x']) and a['charges'] == b['charges']
                    and a['mapping'] == b['mapping'], 'Paired metal geometry/map differs')
    require(parent['scout_reused_exactly'] == scout['scout_native_result'] == parent['rows'][0]['whole_GFN'],
            'Original scout was not reused exactly')
    sn = read(scout['scout_native_result']); se = read(sn['execution'])
    require(se['wall_seconds'] == scout['scout_wall_seconds'] and parent['scout_cost'] == se, 'Saved scout cost differs')
    result.update(manifest=authority['CPU_manifest'], recovery_authority=cfg['authority'],
        preserved_failed_GPU=authority['failed_GPU_result'],
        component_authority='Exact admitted CPU/preflight local arrays; no cross-host cap-scatter replay or new tolerance')
    return base,cells,h,rows


def prepare(cfg, out, result):
    require(cfg['protocol'] == PROTOCOL and cfg['implementation']['sha256'] == pin(__file__)['sha256']
            and cfg['cpu_implementation']['sha256'] == CPU_SHA
            and cfg['hybrid_implementation']['sha256'] == MACE_SHA, 'Completion source declaration differs')
    cpu = module('whole_base_cpu_frozen', cfg['cpu_implementation'])
    cc = read(cfg['cpu_config']); parent = read(cfg['cpu_result'])
    require(cc['implementation']['sha256'] == CPU_SHA and parent['implementation']['sha256'] == CPU_SHA
            and parent['config'] == cfg['cpu_config'], 'CPU source/config binding differs')
    verify(parent['implementation'])
    require(parent['protocol'] == cpu.PROTOCOL and parent['stage'] == 'run-development-cpu'
            and parent['status'] == 'development_cpu_complete', 'Complete actual CPU development panel required')
    for key in ('design', 'inventory', 'preflight_result', 'scout_result'):
        require(cfg[key] == cc[key], 'CPU scientific binding differs: '+key)
    for key in ('inventory', 'preflight_result', 'scout_result'):
        require(parent[key] == cfg[key], 'CPU result lineage differs: '+key)
    require(parent['new_molecular_calls'] == 4 and parent['whole_GFN_calls_attempted'] == 3
            and parent['regional_GFN_calls_attempted'] == 1 and parent['regional_MACE_calls_attempted'] == 0
            and not parent['model_loaded'] and parent['cumulative_new_molecular_calls_including_scout'] == 5,
            'Exact CPU call budget differs')
    require(all(parent[k] == 0 for k in ('ALPB_calls', 'MM_calls', 'embedded_GFN_calls', 'DFT_calls', 'optimizer_steps')),
            'Unexpected CPU molecular component')
    require([r['cell_id'] for r in parent['rows']] == list(IDS)
            and all(r['status'] == 'not_requested' and r['whole_GFN'] is None and r['regional_GFN'] is None
                    for r in parent['rows'][4:]), 'Eight CPU slots or unrequested challenge differ')
    requests = parent['requests']
    require([(r['cell_id'], r['scope']) for r in requests] == list(cpu.CALLS)
            and all(r['status'] == 'accepted' for r in requests), 'Four default-solver requests differ')
    d = out/'saved_authority_admission'; d.mkdir()
    base, cells, h, rows = saved_authority(cfg,cpu,cc,parent,d,result)
    native = {}
    for r in parent['rows'][:4]:
        cid = r['cell_id']; m = rows[cid]
        require(all(r[k] == m[k] for k in ('source_id', 'geometry', 'metal', 'state')),
                'CPU row state differs')
        require(r['status'] == ('whole_and_regional_GFN_available' if cid == MISSING else 'whole_GFN_available')
                and (r['regional_GFN'] is None) == (cid != MISSING), 'CPU component availability differs')
        native[cid] = {}
        for kind in ('whole', 'region') if cid == MISSING else ('whole',):
            record = r['whole_GFN' if kind == 'whole' else 'regional_GFN']; n = 256 if kind == 'whole' else 172
            expected = {'source_id': m['source_id'], 'metal': m['metal'], 'kind': 'isolated',
                'region': 'whole' if kind == 'whole' else 'baseline', 'atom_count': n, 'charge': -2,
                'uhf': 0, 'multiplicity': 1, 'physical_multiplicity': m['state'][kind]['MACE_physical_multiplicity'],
                'valence_electrons': m['state'][kind]['GFN_valence_electrons'], 'pointcharge_count': 0,
                'field_charge_e': 0., 'physical_source': m['physical'], 'physical_ids': m['physical_ids']}
            raw = base.native_replay(record, expected, h)
            task = read(raw['task'])
            require(task['xyz']['sha256'] == m[kind+'_xyz']['sha256'], 'Exact CPU native geometry differs')
            if kind == 'region': require(task['mapping'] == m['mapping'], 'Native regional map differs')
            if cid != IDS[0]:
                request = next(x for x in requests if (x['cell_id'], x['scope']) == (cid, kind))
                require(request['native_result'] == record and request['task'] == raw['task']
                        and request['execution'] == raw['execution']
                        and task['config'] == parent['config'], 'CPU request/native receipt differs')
            native[cid][kind] = {'pin': record, 'raw': raw}
    result['rows'] = [{k: r[k] for k in ('cell_id', 'source_id', 'geometry', 'metal', 'physical', 'positions',
                                        'mapping', 'physical_ids', 'state')} | {'status': 'not_requested'}
                      for r in rows.values()]
    result['parent_results'] = {'preflight': cfg['preflight_result'], 'scout': cfg['scout_result'], 'CPU': cfg['cpu_result']}
    result['parent_costs'] = {'scout_native': parent['scout_cost'], 'CPU_result': cfg['cpu_result'],
        'CPU_native_executions': [read(x['execution']) for x in requests],
        'CPU_elapsed_seconds': parent['elapsed_seconds'], 'CPU_allocation': parent['allocation'],
        'scheduler_utilization_not_inferred': True}
    return base, cells, h, rows, native


def mace_task(row, result):
    return {'task_id': MISSING+'__region_MACE', 'cell_id': MISSING, 'source_id': '7CCO',
        'region': 'baseline', 'metal': 'Dy', 'atom_count': 172, 'charge': -2, 'multiplicity': 6,
        'target_index': 0, 'physical_source': row['physical'], 'physical_ids': row['physical_ids'],
        'mapping': row['mapping'], 'xyz': row['region_xyz'], 'protocol': PROTOCOL,
        'config': result['config'], 'cpu_result': result['parent_results']['CPU'],
        'physical_state_qualified': False, 'derivative_qualified': False}


def run_gpu(cfg, cells, h, rows, out, result):
    require(cfg.get('gpu_result') is None, 'A prior GPU result cannot be rerun')
    cell = cells[MISSING]; row = rows[MISSING]
    require(cell['local'] is None and row['local_components'] is None
            and row['state']['region']['MACE_physical_multiplicity'] == 6, 'Only missing after-Dy regional MACE is allowed')
    d = out/'regional_MACE'; d.mkdir()
    (d/'core.xyz').write_bytes(verify(row['region_xyz']).read_bytes())
    task = mace_task(row, result); task['xyz'] = pin(d/'core.xyz')
    task_pin = write(d/'TASK.json', task)
    request = {'task': task_pin, 'cell_id': MISSING, 'status': 'unavailable', 'attempted': False}
    result['request'] = request
    hybrid = module('hybrid_forward_frozen', cfg['hybrid_implementation'])
    warm = module('native_MACE_loader_frozen', h['pins']['native_MACE_helper'])
    mcfg = read(h['pins']['MACE_config'])
    context = SimpleNamespace(native=None, warm=warm, mcfg=mcfg, alpb_helper=h['solvent'],
        cfg={'model': h['pins']['model']}, cost={'MACE_model_loads': 0, 'MACE_calls_attempted': 0})
    raw = {'status': 'unavailable', 'attempted': False, 'task': task_pin,
           'implementation': result['implementation'], 'native_forward': cfg['hybrid_implementation'],
           'native_loader': h['pins']['native_MACE_helper'], 'native_config': h['pins']['MACE_config']}
    try:
        # Exactly the frozen 172-atom native call; no HybridEvaluator initialization or mixed evaluator.
        value = hybrid.HybridEvaluator.mace(context, task, cell['core_symbols'], cell['core_x'])
        require(context.native.calls == context.cost['MACE_calls_attempted'] == 1, 'One-forward ceiling differs')
        raw.update(value, status='accepted', attempted=True)
        raw['forces_artifact'] = write(d/'FORCES_eV_A.json', raw['forces_eV_A'])
        raw['native_runtime'] = {'versions': context.native.versions, 'device': context.native.device,
            'parameter_versions_unchanged': True, 'native_forward_count': context.native.calls,
            'peak_allocated_bytes': context.native.torch.cuda.max_memory_allocated(),
            'peak_reserved_bytes': context.native.torch.cuda.max_memory_reserved(),
            'total_device_bytes': context.native.torch.cuda.get_device_properties(0).total_memory}
        request.update(status='accepted', attempted=True)
    except Exception as exc:
        raw.update(status='unavailable', failure=repr(exc), traceback=traceback.format_exc())
        raise
    finally:
        count = context.cost['MACE_calls_attempted']
        raw['attempted'] = request['attempted'] = bool(count)
        result['new_molecular_calls'] = result['regional_MACE_calls_attempted'] = count
        result['model_loaded'] = context.native is not None
        result['MACE_cost'] = dict(context.cost)
        raw['cost'] = dict(context.cost)
        raw['artifacts'] = {p.name: pin(p) for p in d.iterdir() if p.is_file()}
        result['regional_MACE'] = request['native_result'] = write(d/'RESULT.json', raw)
        result['rows'][3].update(status=('regional_MACE_available' if raw['status'] == 'accepted' else
            'unavailable_attempted' if count else 'unavailable_before_forward'), regional_MACE=result['regional_MACE'])
    result['rows'][3].update(status='regional_MACE_available', regional_MACE=result['regional_MACE'])
    result.update(status='development_gpu_complete', cumulative_new_molecular_calls=6,
                  cumulative_whole_GFN_calls=4, cumulative_regional_GFN_calls=1, cumulative_regional_MACE_calls=1)


def saved_mace(cfg, rows, h, result):
    parent = read(cfg['gpu_result']); pc = read(parent['config'])
    require(parent['implementation']['sha256'] == cfg['implementation']['sha256']
            and pc['implementation']['sha256'] == cfg['implementation']['sha256']
            and binding(pc) == binding(cfg) == parent['scientific_binding'], 'GPU source/scientific binding differs')
    verify(parent['implementation'])
    require(parent['protocol'] == PROTOCOL and parent['stage'] == 'run-development-gpu'
            and parent['status'] == 'development_gpu_complete' and parent['model_loaded']
            and parent['new_molecular_calls'] == parent['regional_MACE_calls_attempted'] == 1
            and parent['cumulative_new_molecular_calls'] == 6, 'Complete one-forward GPU result required')
    require(all(parent[k] == 0 for k in ('whole_GFN_calls_attempted', 'regional_GFN_calls_attempted',
            'ALPB_calls', 'MM_calls', 'embedded_GFN_calls', 'DFT_calls', 'optimizer_steps')), 'Unexpected GPU component calls')
    require([r['cell_id'] for r in parent['rows']] == list(IDS)
            and all(r['status'] == 'not_requested' for r in parent['rows'][4:]), 'GPU eight-slot challenge boundary differs')
    raw = read(parent['regional_MACE']); hbase = rows[MISSING]; task = read(raw['task'])
    expected = mace_task(hbase, {'config': parent['config'], 'parent_results': {'CPU': cfg['cpu_result']}})
    require({k: v for k, v in task.items() if k != 'xyz'} == {k: v for k, v in expected.items() if k != 'xyz'}
            and task['xyz']['sha256'] == expected['xyz']['sha256'], 'Saved regional MACE input/state differs')
    for p in raw['artifacts'].values(): verify(p)
    require(raw['status'] == 'accepted' and raw['attempted'] and raw['model'] == h['pins']['model']
            and raw['force_convention'] == 'negative vacuum energy derivative'
            and raw['native_forward'] == cfg['hybrid_implementation']
            and raw['native_loader'] == h['pins']['native_MACE_helper']
            and raw['native_config'] == h['pins']['MACE_config']
            and raw['cost']['MACE_calls_attempted'] == raw['native_runtime']['native_forward_count'] == 1,
            'Saved native MACE lineage/forward differs')
    require(raw['native_runtime']['versions'] == h['pins']['expected_MACE_versions']
            and raw['native_runtime']['parameter_versions_unchanged'], 'MACE software/parameter provenance differs')
    h['solvent'].state(raw['input_state_check'], 172, -2, 6, 0)
    f = h['np'].asarray(raw['forces_eV_A'], dtype=float)
    require(f.shape == (172, 3) and h['np'].isfinite(f).all() and math.isfinite(raw['energy_eV'])
            and read(raw['forces_artifact']) == raw['forces_eV_A'], 'Saved raw force/energy replay differs')
    require(parent['request']['native_result'] == parent['regional_MACE']
            and parent['request']['task'] == raw['task'] and parent['request']['status'] == 'accepted',
            'GPU request/native result differs')
    result['parent_results']['GPU'] = cfg['gpu_result']
    result['parent_costs']['GPU_result'] = cfg['gpu_result']
    result['parent_costs']['GPU_native'] = {'cost': raw['cost'], 'runtime': raw['native_runtime'],
        'forward_wall_seconds': raw['wall_seconds'], 'elapsed_seconds': parent['elapsed_seconds'], 'allocation': parent['allocation']}
    return parent['regional_MACE'], raw


def collect(cfg, cells, h, rows, native, out, result):
    new_mace_pin, new_mace = saved_mace(cfg, rows, h, result)
    np = h['np']; completed = []
    result['partial_assembled_artifacts'] = []
    for cid in IDS[:4]:
        cell = cells[cid]; m = rows[cid]; full = native[cid]['whole']; local = cell['local']
        if cid == MISSING:
            mace_pin, mace = new_mace_pin, new_mace
            isolated_pin, isolated = native[cid]['region']['pin'], native[cid]['region']['raw']
            gm = -np.asarray(mace['forces_eV_A'], dtype=float)*EV
            gi = np.asarray(isolated['qm_gradient']['hartree_per_bohr'], dtype=float)*HA/BOHR
            zeros = np.zeros((len(cell['mapping']['field']), 3))
            gm_full = h['mapper'].scatter(cell['mapping'], gm, zeros)
            minus_gi_full = -h['mapper'].scatter(cell['mapping'], gi, zeros)
        else:
            mace_pin, isolated_pin = local['MACE'], local['GFN']
            mace, isolated = read(mace_pin), read(isolated_pin)
            gm_full = np.asarray(local['component_gradients_kcal_mol_A']['MACE_region'], dtype=float)
            minus_gi_full = np.asarray(local['component_gradients_kcal_mol_A']['negative_GFN_region'], dtype=float)
        require(m.get('local_components') is None or local['solver_exception_certificate'] is None,
                'No solver exception belongs to the four 7CCO development cells')
        n = {'GFN_full_Eh': full['raw']['energy_hartree'], 'MACE_region_eV': mace['energy_eV'],
             'GFN_region_Eh': isolated['energy_hartree'], 'ALPB_Eh': cell['solvent']['energy_hartree']}
        require(all(math.isfinite(v) for v in n.values()), 'Native energy ledger nonfinite')
        components = {'GFN_full': np.asarray(full['raw']['qm_gradient']['hartree_per_bohr'], dtype=float)*HA/BOHR,
            'MACE_region': gm_full, 'negative_GFN_region': minus_gi_full,
            'fixed_FF_ALPB': np.asarray(cell['solvent']['gradient_kcal_mol_A'], dtype=float)}
        require(all(a.shape == (256, 3) and np.isfinite(a).all() for a in components.values()), 'Full256 component gradient unavailable')
        candidate = sum(components[k] for k in COMPONENTS)
        base = components['GFN_full']+components['fixed_FF_ALPB']
        energies = (n['GFN_full_Eh']*HA, n['MACE_region_eV']*EV, -n['GFN_region_Eh']*HA, n['ALPB_Eh']*HA)
        row = {k: m[k] for k in ('cell_id', 'source_id', 'geometry', 'metal', 'physical', 'positions',
                                 'mapping', 'physical_ids', 'state')}
        row.update(status='complete', component_receipts={'GFN_full': full['pin'], 'MACE_region': mace_pin,
            'GFN_region': isolated_pin, 'ALPB': cell['solvent']['receipt']}, solver_exception_certificate=None,
            native_energies=n, component_gradients_kcal_mol_A={k: a.tolist() for k, a in components.items()},
            candidate_gradient_kcal_mol_A=candidate.tolist(), base_gradient_kcal_mol_A=base.tolist(),
            candidate_force_kcal_mol_A=(-candidate).tolist(), base_force_kcal_mol_A=(-base).tolist(),
            candidate_energy_kcal_mol=math.fsum(energies), base_energy_kcal_mol=math.fsum((energies[0], energies[3])),
            physical_state_qualified=False, derivative_qualified=False, new_molecular_calls=0)
        row['assembled_artifact'] = write(out/(cid+'.json'), row)
        result['partial_assembled_artifacts'].append({'cell_id': cid, 'artifact': row['assembled_artifact']})
        completed.append(row)
    # No four-cell comparison can escape from an incomplete or failed assembly.
    result['rows'][:4] = completed
    result.update(status='complete_7CCO', new_molecular_calls=0, call_accounting={
        'development_total': 6, 'whole_GFN': 4, 'regional_GFN': 1, 'regional_MACE': 1,
        'new_collect_calls': 0, 'remaining_conditional_7CCN_whole_GFN': 4, 'overall_ceiling': 10,
        'repeated_cached_components': 0, 'ALPB': 0, 'MM': 0, 'embedded_GFN': 0, 'DFT': 0},
        comparison_status='Not computed here; frozen saved comparator owns work/force triage',
        challenge_status='not_requested; requires separate frozen work triage and root continuation',
        missing_comparators={'before_Dy_whole_MACE_force': True, 'after_Dy_old_hybrid': True})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--stage', choices=('run-development-gpu', 'collect_7CCO'), required=True)
    parser.add_argument('--config', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    require(os.environ.get('SLURM_JOB_ID', '').isdigit(), 'Slurm required before scientific imports')
    out = args.output.resolve(); out.mkdir(parents=True, exist_ok=False); tick = time.monotonic()
    result = {'protocol': PROTOCOL, 'stage': args.stage, 'status': 'execution_failure',
        'implementation': pin(__file__), 'config': pin(args.config), 'slurm_job_id': os.environ['SLURM_JOB_ID'],
        'host': socket.gethostname(), 'declared_cells': list(IDS),
        'rows': [{'cell_id': cid, 'status': 'not_requested'} for cid in IDS],
        'new_molecular_calls': 0, 'whole_GFN_calls_attempted': 0, 'regional_GFN_calls_attempted': 0,
        'regional_MACE_calls_attempted': 0, 'model_loaded': False, 'ALPB_calls': 0, 'MM_calls': 0,
        'embedded_GFN_calls': 0, 'DFT_calls': 0, 'optimizer_steps': 0,
        'physical_state_qualified': False, 'derivative_qualified': False,
        'limitations': ['Conditional terminal chemistry and physical spin states remain unqualified.',
            'GFN2 UHF0 at 300K can have fractional occupations; no ground-state or derivative qualification.',
            'Regional vacuum delta is not mutual polarization; no absolute cross-element affinity claim.',
            'The 7CCN challenge, before-Dy whole-MACE force and after-Dy old hybrid remain unavailable here.']}
    try:
        cfg = read(result['config'])
        result.update(scientific_binding=binding(cfg), inventory=cfg['inventory'],
                      preflight_result=cfg['preflight_result'], scout_result=cfg['scout_result'])
        _, cells, h, rows, native = prepare(cfg, out, result)
        if args.stage == 'run-development-gpu': run_gpu(cfg, cells, h, rows, out, result)
        else: collect(cfg, cells, h, rows, native, out, result)
    except Exception as exc:
        result.update(status='execution_failure', failure=repr(exc), traceback=traceback.format_exc())
    result.update(elapsed_seconds=time.monotonic()-tick, peak_parent_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        allocation={k: os.environ.get(k) for k in ('SLURM_CPUS_PER_TASK', 'SLURM_CPUS_ON_NODE', 'SLURM_JOB_CPUS_PER_NODE',
                    'SLURM_JOB_GPUS', 'SLURM_STEP_GPUS', 'CUDA_VISIBLE_DEVICES', 'SLURM_MEM_PER_CPU', 'SLURM_MEM_PER_NODE')})
    write(out/'RESULT.json', result)
    print(json.dumps({'status': result['status'], 'result': pin(out/'RESULT.json'), 'new_molecular_calls': result['new_molecular_calls']}))
    return 0 if result['status'] in ('development_gpu_complete', 'complete_7CCO') else 1


if __name__ == '__main__':
    sys.exit(main())
