#!/usr/bin/env python3
"""UNBOUND source draft; no execution approval or configuration is supplied.

Every stage requires Slurm and actual admitted scout/collector/cost-decision
pins. Config pin keys are BIND_KEYS below, plus preflight_result for run and
endpoint_results (task ID -> singleton RESULT pin) for saved collection.
The root decision schema is enforced in decision_gate(). No stage submits jobs.
"""
import argparse
import datetime as dt
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import re
import resource
import socket
import time
import traceback

PROTOCOL = 'lbt_whole256_reference_continuation_draft_v1'
DECISION_PROTOCOL = 'lbt_whole256_reference_continuation_decision_v1'
IDS = ('before_La', 'after_La', 'before_Dy', 'after_Dy')
REMAINING = ('after_La', 'before_Dy', 'after_Dy')
V1_SHA = 'c22c989674e280806abace63c67148b0181caee0b7b20bcd85b8e3d75d6c2a86'
ORIGINAL_CONFIG_SHA = 'dc6127e4fad9ffef90e30ed5770e00fd764a1ba6c1b57cbbfb94b6bdbf038ad8'
ORIGINAL_PREFLIGHT_SHA = '151d8d5f5fc0d4d4729933b64d98e4962e21518504475facb1d61c37d67c554e'
MODEL_COLLECTION_SHA = '64f27e852d2403600ef7a73616c7e52cf50cfbace75689a5b932fde926a06ef0'
MODEL_SOURCE_SHA = 'f35e2225473b1df10288c532a930bdd49186eaaf2e3bab40c65026c44024d8bc'
AFTER_DY_ALPB_SHA = '26f188d829d2921b585ec53db13ac3e72707b9e7a9e825b6e9c9c4615b34477d'
BIND_KEYS = ('protocol', 'implementation', 'specification', 'original_implementation',
             'original_scout_config', 'original_preflight', 'scout_result',
             'scout_collection_result', 'continuation_decision', 'model_completion_result')


def require(ok, message):
    if not ok:
        raise ValueError(message)


def pin(path):
    p = Path(path).resolve()
    with p.open('rb') as f:
        sha = hashlib.file_digest(f, 'sha256').hexdigest()
    return {'path': str(p), 'sha256': sha}


def verify(record):
    require(isinstance(record, dict) and set(record) == {'path', 'sha256'}, 'Required evidence pin missing')
    require(pin(record['path'])['sha256'] == record['sha256'], 'Changed pin: ' + record['path'])
    return Path(record['path'])


def read(record):
    return json.loads(verify(record).read_text())


def write(path, value):
    with Path(path).open('x') as f:
        json.dump(value, f, indent=2, allow_nan=False); f.write('\n')
    return pin(path)


def binding(cfg):
    return {k: cfg[k]['sha256'] if k == 'implementation' else cfg[k] for k in BIND_KEYS}


def decision_gate(cfg, scout, collection):
    """Root writes this evidence only after actual native/force/cost review."""
    decision = read(cfg['continuation_decision'])
    require(decision['protocol'] == DECISION_PROTOCOL and decision['status'] == 'approved_remaining_three'
            and decision['implementation']['sha256'] == cfg['implementation']['sha256']
            and decision['specification'] == cfg['specification']
            and decision['approved_tasks'] == list(REMAINING) and decision['maximum_new_DFT_calls'] == 3
            and decision['allow_retry'] is False, 'Continuation decision scope/source differs')
    verify(decision['implementation'])
    for key in ('original_implementation', 'original_scout_config', 'original_preflight',
                'scout_result', 'scout_collection_result', 'model_completion_result'):
        require(decision[key] == cfg[key], 'Continuation decision evidence differs: ' + key)
    cost = decision['measured_cost']
    require(cost['native_receipt'] == scout['rows']['before_La']['receipt']
            and cost['force_comparison'] == collection['force_comparison']
            and isinstance(cost['scheduler_receipts'], list) and cost['scheduler_receipts'],
            'Actual scout force/cost evidence required')
    for p in [cost['native_receipt'], cost['force_comparison'], *cost['scheduler_receipts']]:
        verify(p)
    return decision


def inputs(cfg):
    require(os.environ.get('SLURM_JOB_ID', '').isdigit(), 'Every stage requires Slurm')
    require(cfg['protocol'] == PROTOCOL and pin(__file__)['sha256'] == cfg['implementation']['sha256'],
            'Draft implementation/protocol differs')
    # Fail on absent future evidence before importing any scientific helper.
    for key in BIND_KEYS[1:]:
        verify(cfg[key])
    require(cfg['original_implementation']['sha256'] == V1_SHA
            and cfg['original_scout_config']['sha256'] == ORIGINAL_CONFIG_SHA
            and cfg['original_preflight']['sha256'] == ORIGINAL_PREFLIGHT_SHA
            and cfg['model_completion_result']['sha256'] == MODEL_COLLECTION_SHA, 'Frozen parent differs')
    scout, collected = read(cfg['scout_result']), read(cfg['scout_collection_result'])
    require(scout['slurm_job_id'] == '13051' and scout['stage'] == 'scout' and scout['status'] == 'complete'
            and scout['selected_tasks'] == ['before_La'] and scout['new_DFT_calls_attempted'] == 1
            and scout['rows']['before_La']['status'] == 'available', 'Admitted actual13051 required')
    require(collected['stage'] == 'collect' and collected['status'] == 'complete'
            and collected['scout_result'] == cfg['scout_result']
            and collected['rows']['before_La']['status'] == 'available'
            and collected['new_DFT_calls_attempted'] == collected['native_runner_invocations'] == 0,
            'Original saved scout collector required')
    decision = decision_gate(cfg, scout, collected)
    spec = importlib.util.spec_from_file_location('whole_continuation_original', verify(cfg['original_implementation']))
    old = importlib.util.module_from_spec(spec); spec.loader.exec_module(old)
    ocfg = read(cfg['original_scout_config'])
    require(ocfg['implementation']['sha256'] == V1_SHA and ocfg['preflight_result'] == cfg['original_preflight'],
            'Original bound scout configuration differs')
    inv, compact, modules, runtime, np = old.inputs(ocfg)
    model = modules['metal_environment_frozen_f']
    cells = old.physical_inputs(inv, np, model)
    original_solvent = {tid: old.admit_alpb(c, ocfg, np, model) for tid, c in cells.items()}
    require(all(original_solvent[k]['status'] == 'available' for k in IDS[:3])
            and original_solvent['after_Dy']['status'] == 'unavailable', 'Original solvent inventory differs')
    pf = read(cfg['original_preflight']); full = read(pf['full_manifest']); directions = read(pf['directions'])
    require(pf['slurm_job_id'] == '13049' and pf['all_four_states_frozen'], 'Original preflight identity differs')
    old.validate_prepared(full, cells, compact, ocfg, model, np)
    old.preflight_gate(ocfg, full, directions, original_solvent)
    require(directions['physical_ids'] == cells['before_La']['ids']
            and directions['frozen_before_new_reference_calls']
            and set(directions['directions']) == {'full_deformation_chord', 'initial_hybrid_Phe1_CA_force'},
            'Original frozen direction identity differs')
    for vector in directions['directions'].values():
        a = np.asarray(vector)
        require(a.shape == (256, 3) and np.isfinite(a).all(), 'Invalid archived direction')
    for saved, stage in ((scout, 'scout'), (collected, 'collect')):
        saved_cfg = read(saved['config']); verify(saved['implementation'])
        require(saved['protocol'] == old.PROTOCOL and saved['stage'] == stage
                and saved['implementation']['sha256'] == saved_cfg['implementation']['sha256'] == V1_SHA
                and old.scientific_binding(saved_cfg) == old.scientific_binding(ocfg)
                and saved['preflight_result'] == cfg['original_preflight']
                and saved['solvent_admission'] == original_solvent and read(saved['directions']) == directions,
                'Scout/collector scientific binding differs')
        sf = read(scout['full_manifest']); old.validate_prepared(sf, cells, compact, ocfg, model, np)
        require(old.signatures(sf) == old.signatures(full), 'Scout full-panel signatures differ')
        manifest = read(saved['manifest'])
        require(saved['manifest'] == scout['manifest'] and manifest['full_manifest'] == scout['full_manifest']
                and manifest['tasks'] == [t for t in sf['tasks'] if t['task_id'] == 'before_La']
                and {k: v for k, v in manifest.items() if k not in ('full_manifest', 'tasks')}
                    == {k: v for k, v in sf.items() if k != 'tasks'}, 'Original singleton manifest differs')
    origin = old.admit(scout['manifest'], compact, modules, runtime)
    require(origin['status'] == 'available' and origin == scout['rows']['before_La'] == collected['rows']['before_La']
            and read(scout['force_comparison']) == read(collected['force_comparison']),
            'Actual original scout native/collector replay differs')
    completion = read(cfg['model_completion_result']); verify(completion['implementation'])
    require(completion['status'] == 'complete' and completion['stage'] == 'collect'
            and completion['slurm_job_id'] == '13142' and completion['implementation']['sha256'] == MODEL_SOURCE_SHA
            and completion['scientific_binding']['reference_preflight'] == cfg['original_preflight']
            and completion['scientific_binding']['inventory'] == ocfg['inventory']
            and completion['directions'] == pf['directions'], 'Completed model adjunct differs')
    whole = {row['cell_id']: row for row in completion['rows']}
    require(set(whole) == set(IDS) and len(completion['rows']) == 4
            and all(row['status'] == 'available' for row in whole.values()), 'Four whole-model energies required')
    for tid, cell in cells.items():
        expected = dict(physical_parent=cell['physical_parent'], positions=cell['row']['positions'],
                        physical_ids=cell['ids'], symbols=cell['symbols'], positions_A=cell['x'].tolist(),
                        charge=-2, multiplicity=cell['state']['physical_multiplicity'],
                        fixed_FF_charges_e=cell['charges']['charges_e'])
        require(completion['signatures'][tid] == expected
                and whole[tid]['physical_source'] == cell['physical_parent']
                and whole[tid]['positions'] == cell['row']['positions']
                and whole[tid]['metal'] == cell['row']['metal'] and whole[tid]['physical_ids'] == cell['ids'],
                'Same-coordinate whole-model comparison differs')
    supplemental = whole['after_Dy']['ALPB']
    require(supplemental['sha256'] == AFTER_DY_ALPB_SHA, 'Supplemental after-Dy solvent differs')
    # Original cells/full manifest/13049 evidence are never mutated.
    after_dy = dict(cells['after_Dy'], row=dict(cells['after_Dy']['row'], ALPB=supplemental))
    solvent = dict(original_solvent, after_Dy=old.admit_alpb(after_dy, ocfg, np, model))
    require(solvent['after_Dy']['status'] == 'available', 'Supplemental after-Dy raw admission failed')
    for tid in IDS:
        require(whole[tid]['ALPB'] == solvent[tid]['receipt'], 'Model/reference solvent receipt differs')
    return dict(cfg=cfg, old=old, ocfg=ocfg, inv=inv, compact=compact, modules=modules, runtime=runtime,
                np=np, model=model, cells=cells, original_preflight=pf, original_full=full,
                original_solvent=original_solvent, directions=directions, origin=origin,
                solvent=solvent, whole=whole, completion=completion, decision=decision)


def prepare_task(w, tid, directory):
    require(tid in REMAINING, 'Only the three remaining cells are executable; before-La repeat forbidden')
    directory.mkdir()
    original = next(t for t in w['original_full']['tasks'] if t['task_id'] == tid)
    task = dict(original)
    for key, name in (('xyz', 'physical.xyz'), ('input', 'endpoint.inp')):
        target = directory / name; target.write_bytes(verify(original[key]).read_bytes()); task[key] = pin(target)
    task.update(output_path=str(directory/'endpoint.out'), engrad_path=str(directory/'endpoint.engrad'))
    require(w['old'].signatures({'tasks': [task]}) == w['old'].signatures({'tasks': [original]}),
            'Remaining native task differs from original13049')
    write(directory/'TASK.json', task)
    manifest = dict(w['original_full'], protocol_id=PROTOCOL, agreement=w['cfg']['specification'],
                    adapter=pin(__file__), original_full_manifest=w['original_preflight']['full_manifest'],
                    inputs_config=w['cfg']['_config_pin'], original_inputs_config=w['original_full']['inputs_config'],
                    continuation_decision=w['cfg']['continuation_decision'], scientific_binding=binding(w['cfg']),
                    execution_allowlist=list(REMAINING), tasks=[task])
    return write(directory/'MANIFEST.json', manifest)


def validate_manifest(w, record, tid):
    m = read(record)
    require(tid in REMAINING and len(m['tasks']) == 1 and m['tasks'][0]['task_id'] == tid,
            'Exactly one remaining native task required')
    expected = dict(w['original_full'], protocol_id=PROTOCOL, agreement=w['cfg']['specification'],
                    adapter=m['adapter'], original_full_manifest=w['original_preflight']['full_manifest'],
                    inputs_config=m['inputs_config'], original_inputs_config=w['original_full']['inputs_config'],
                    continuation_decision=w['cfg']['continuation_decision'], scientific_binding=binding(w['cfg']),
                    execution_allowlist=list(REMAINING), tasks=m['tasks'])
    require(m == expected and m['adapter']['sha256'] == w['cfg']['implementation']['sha256'],
            'Singleton native manifest/protocol differs')
    verify(m['adapter'])
    require(binding(read(m['inputs_config'])) == binding(w['cfg']), 'Singleton configuration differs')
    task = m['tasks'][0]; original = next(t for t in w['original_full']['tasks'] if t['task_id'] == tid)
    require(w['old'].signatures({'tasks': [task]}) == w['old'].signatures({'tasks': [original]}),
            'Singleton physical/state/deck signature differs')
    for key in ('xyz', 'input'):
        require(verify(task[key]).read_bytes() == verify(original[key]).read_bytes(), 'Native input bytes differ')
    d = verify(task['input']).parent
    require(verify(task['xyz']).parent == d and Path(task['output_path']) == d/'endpoint.out'
            and Path(task['engrad_path']) == d/'endpoint.engrad', 'Native task output scope differs')
    return m


def admit(w, manifest_pin, tid):
    """V1 native admission, generalized only from SCOUT to one named task."""
    validate_manifest(w, manifest_pin, tid)
    modules, compact, runtime = w['modules'], w['compact'], w['runtime']
    manifest, rows = modules['ggr_sensitivity'].executed(verify(manifest_pin))
    require([t['task_id'] for t in manifest['tasks']] == [tid], 'Native execution denominator differs')
    task = manifest['tasks'][0]; row = rows[tid]
    row.update(task=task, task_id=tid, physical_multiplicity=task['physical_multiplicity'])
    try:
        if Path(task['output_path']).exists():
            row['SCF_evidence'] = runtime.observations(task['output_path'], modules['mace_omol_vacuum'])
        if row['status'] != 'complete':
            return w['old'].canonical(row)
        parsed = w['model'].parse(task, {'electronic_state': task['electronic_state']})
        require(math.isfinite(parsed['energy_hartree']) and len(parsed['gradient_kcal_mol_per_A']) == 256,
                'Full finite native gradient unavailable')
        text = Path(task['output_path']).read_text()
        require('INITIAL GUESS: MODEL POTENTIAL' in text
                and not re.search(r'Reading point charge file|^\s*(?:CPCM|SMD|COSMO) SOLVATION MODEL', text, re.M),
                'Vacuum/guess scope changed')
        require(row['SCF_evidence']['ConvCheckMode'] == 'Total+1el-Energy'
                and row['SCF_evidence']['ConvForced'] == 0, 'Native SCF convergence semantics changed')
        receipt = read(row['receipt']); sidecar = read(receipt['artifacts']['runtime_sidecar'])
        require(receipt['orca_executable'] == compact['orca'] and receipt['returncode'] == 0
                and receipt['normal_termination'] and receipt['scf_converged'], 'Native execution failed')
        policy = sidecar['memory_policy']
        require(policy['workers'] == 1 and policy['ranks_per_worker'] == 12
                and policy['memory_source'] == 'SLURM_MEM_PER_NODE', 'Native allocation-derived memory policy differs')
        row.update(parsed, status='available', runtime_memory_policy=policy,
                   molecular_wall_seconds=(dt.datetime.fromisoformat(receipt['finished_at_utc']) -
                                           dt.datetime.fromisoformat(receipt['started_at_utc'])).total_seconds(),
                   caveats=w['old'].CAVEATS)
    except Exception as exc:
        row.update(status='unavailable', admission_failure=repr(exc), energy_hartree=None)
    return w['old'].canonical(row)


def force_report(w, tid, native):
    np, cell, model = w['np'], w['cells'][tid], w['model']
    g = np.asarray(native['gradient_kcal_mol_per_A'])
    gs = np.asarray(w['solvent'][tid]['gradient_kcal_mol_A'])
    forces = {'reference_PBE0_vacuum': -g, 'reference_PBE0_plus_fixed_ALPB': -(g+gs)}
    missing = {}; whole = w['whole'][tid]
    if tid == 'before_Dy':
        require(whole['force_arrays'] is None and whole['force_status'] == 'unavailable_energy_only_archived_request',
                'Before-Dy model force availability changed')
        missing['whole_MACE_composite'] = 'Archived energy-only request; no force call authorized.'
    else:
        a = read(whole['force_arrays'])
        require(a['physical_ids'] == cell['ids'] and a['units'] == 'kcal/mol/Angstrom'
                and a['quantity'] == 'force = negative derivative', 'Whole-model force provenance differs')
        raw = read(whole['MACE']); mf = np.asarray(raw['forces_eV_A'])*23.06054783061903
        sf = -np.asarray(w['solvent'][tid]['gradient_hartree_per_bohr'])*(model.HA_TO_KCAL/model.BOHR_TO_A)
        require(np.array_equal(np.asarray(a['components']['MACE_vacuum']), mf)
                and np.array_equal(np.asarray(a['components']['fixed_FF_ALPB']), sf)
                and np.array_equal(np.asarray(a['components']['composite']), mf+sf), 'Whole raw force replay differs')
        forces['whole_MACE_composite'] = np.asarray(a['components']['composite'])
    if cell['row']['hybrid_point'] is None:
        missing['hybrid_composite'] = 'No after-Dy hybrid cell; closed response remains closed.'
    else:
        p = read(cell['row']['hybrid_point'])['point']
        require(p['positions_A'] == cell['x'].tolist() and p['metal'] == cell['row']['metal'],
                'Hybrid force geometry/state differs')
        forces['hybrid_composite'] = -np.asarray(p['gradient_kcal_mol_A'])
    require(all(v.shape == (256, 3) and np.isfinite(v).all() for v in forces.values()), 'Invalid full force array')
    def metrics(v):
        n = np.linalg.norm(v, axis=1)
        return {'maximum_atom_kcal_mol_A': float(n.max()), 'RMS_atom_vector_kcal_mol_A': float(np.sqrt(np.mean(n*n)))}
    table = [dict(index=i, physical_id=name, resname=cell['atoms'][i].get('resname'),
                  force_kcal_mol_A={k: v[i].tolist() for k, v in forces.items()}) for i, name in enumerate(cell['ids'])]
    return dict(cell_id=tid, quantity='force=-gradient', units='kcal/mol/Angstrom',
                fixed_physical_input=cell['row']['positions'], force_table=table,
                Phe1_CA=table[cell['ids'].index('A/1//CA')],
                Phe1=[r for r in table if r['physical_id'].startswith('A/1/')],
                unavailable_comparators=missing, distributions={k: metrics(v) for k, v in forces.items()},
                reference_minus_model_force_errors={k: metrics(forces['reference_PBE0_plus_fixed_ALPB']-v)
                    for k, v in forces.items() if k in ('whole_MACE_composite', 'hybrid_composite')},
                force_projections={name: {k: float(np.sum(v*np.asarray(direction))) for k, v in forces.items()}
                    for name, direction in w['directions']['directions'].items()},
                numerical_force_accuracy_qualified=False, caveats=w['old'].CAVEATS)


def comparison(w, rows):
    if not all(rows[k]['status'] == 'available' for k in IDS):
        return {'status': 'unavailable', 'reason': 'All four reference cells required; no survivor differential.'}
    works = {}
    for metal in ('La', 'Dy'):
        after, before = 'after_'+metal, 'before_'+metal
        d = rows[after]['energy_hartree']-rows[before]['energy_hartree']
        s = w['solvent'][after]['energy_hartree']-w['solvent'][before]['energy_hartree']
        works[metal] = dict(PBE0_vacuum_Eh=d, fixed_FF_ALPB_Eh=s,
                           PBE0_vacuum_kcal_mol=d*w['model'].HA_TO_KCAL,
                           fixed_FF_ALPB_kcal_mol=s*w['model'].HA_TO_KCAL,
                           composite_kcal_mol=math.fsum((d*w['model'].HA_TO_KCAL, s*w['model'].HA_TO_KCAL)))
    return dict(status='complete', definition='W_M=E_M(after)-E_M(before); delta_W=W_Dy-W_La',
                per_metal_work=works, Dy_minus_La_work={k: works['Dy'][k]-works['La'][k] for k in works['La']},
                whole_model_work=w['completion']['comparison']['per_metal_work'],
                whole_model_Dy_minus_La_work=w['completion']['comparison']['Dy_minus_La_deformation_work'],
                hybrid_La_work=w['completion']['comparison']['hybrid_La_work'], hybrid_Dy_work=None,
                reporting_resolution_kcal_mol=1., resolution_is_not_error_bar=True,
                accuracy_or_preference_claim=False, caveats=w['old'].CAVEATS)


def preflight_gate(w):
    p = read(w['cfg']['preflight_result'])
    require(p['protocol'] == PROTOCOL and p['stage'] == 'preflight' and p['status'] == 'preflight_passed'
            and p['scientific_binding'] == binding(w['cfg'])
            and p['original_signatures'] == w['old'].signatures(w['original_full'])
            and p['directions'] == w['original_preflight']['directions']
            and p['original_solvent_admission'] == w['original_solvent']
            and p['solvent_admission'] == w['solvent']
            and p['new_DFT_calls_attempted'] == p['native_runner_invocations'] == 0,
            'Actual unchanged zero-call continuation preflight required')
    require(set(p['manifests']) == set(REMAINING), 'Preflight singleton inventory differs')
    for tid, record in p['manifests'].items():
        validate_manifest(w, record, tid)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--stage', choices=('preflight', 'run', 'collect'), required=True)
    parser.add_argument('--task', choices=REMAINING)
    parser.add_argument('--config', required=True); parser.add_argument('--output', required=True)
    args = parser.parse_args()
    require(os.environ.get('SLURM_JOB_ID', '').isdigit(), 'Every stage requires Slurm')
    out = Path(args.output).resolve(); out.mkdir(parents=True, exist_ok=False); tick = time.monotonic()
    result = dict(protocol=PROTOCOL, stage=args.stage, status='execution_failure', implementation=pin(__file__),
                  config=pin(args.config), slurm_job_id=os.environ['SLURM_JOB_ID'], host=socket.gethostname(),
                  declared_cells=4, remaining_allowlist=list(REMAINING), requested_task=args.task, selected_tasks=[],
                  rows={k: {'task_id': k, 'status': 'not_requested'} for k in IDS},
                  new_DFT_calls_attempted=0, native_runner_invocations=0, new_ALPB_calls=0,
                  new_MACE_calls=0, new_GFN_calls=0, new_MM_calls=0, optimizer_steps=0, new_geometry_count=0,
                  full_four_available=False,
                  comparison={'status': 'unavailable'}, force_reports={}, parent_costs=[])
    try:
        require((args.stage == 'run') == (args.task is not None), '--task is required only for a singleton run')
        cfg = read(result['config']); cfg['_config_pin'] = result['config']; w = inputs(cfg)
        result.update(scientific_binding=binding(cfg), original_signatures=w['old'].signatures(w['original_full']),
                      original_solvent_admission=w['original_solvent'], solvent_admission=w['solvent'],
                      directions=w['original_preflight']['directions'], caveats=w['old'].CAVEATS,
                      scout_result=cfg['scout_result'], scout_collection_result=cfg['scout_collection_result'],
                      continuation_decision=cfg['continuation_decision'], model_completion_result=cfg['model_completion_result'])
        result['rows']['before_La'] = w['origin']
        if args.stage == 'preflight':
            result['manifests'] = {tid: prepare_task(w, tid, out/tid) for tid in REMAINING}
            for tid, record in result['manifests'].items(): validate_manifest(w, record, tid)
            result['status'] = 'preflight_passed'
        elif args.stage == 'run':
            preflight_gate(w); result['preflight_result'] = cfg['preflight_result']
            result['allocation'] = w['runtime'].allocation(w['modules']['render_orca_runtime_input'], out)
            record = prepare_task(w, args.task, out/args.task); validate_manifest(w, record, args.task)
            result.update(manifest=record, selected_tasks=[args.task]); write(out/'START.json', result)
            try:
                result['native_runner_invocations'] += 1
                w['modules']['run_orca_task_manifest'].run_manifest(verify(record), orca_path=verify(w['compact']['orca']),
                    workers=1, nprocs=12, selected_tasks={args.task},
                    expected_runner_sha256=w['compact']['execution_policy']['task_runner']['sha256'])
            except Exception as exc:
                result['runner_error'] = repr(exc)
            finally:
                result['new_DFT_calls_attempted'] = int(Path(read(record)['tasks'][0]['output_path']).exists())
            result['rows'][args.task] = admit(w, record, args.task)
            result['status'] = 'complete' if result['rows'][args.task]['status'] == 'available' else 'unavailable'
        else:
            endpoints = cfg.get('endpoint_results', {})
            require(isinstance(endpoints, dict) and set(endpoints) <= set(REMAINING), 'Only remaining singleton results accepted')
            seen_jobs = set()
            for tid in REMAINING:
                if tid not in endpoints:
                    result['rows'][tid] = {'task_id': tid, 'status': 'unavailable', 'reason': 'No singleton result supplied'}
                    continue
                parent = read(endpoints[tid])
                pcfg = read(parent['config']); verify(parent['implementation'])
                require(parent['protocol'] == PROTOCOL and parent['stage'] == 'run'
                        and parent['implementation']['sha256'] == cfg['implementation']['sha256']
                        and binding(pcfg) == binding(cfg) and parent['requested_task'] == tid
                        and parent['slurm_job_id'] not in seen_jobs, 'Singleton result identity/repeated allocation differs')
                seen_jobs.add(parent['slurm_job_id'])
                require(parent['new_DFT_calls_attempted'] in (0, 1) and parent['native_runner_invocations'] in (0, 1),
                        'Singleton execution budget differs')
                result['parent_costs'].append(dict(result=endpoints[tid], slurm_job_id=parent['slurm_job_id'],
                    new_DFT_calls_attempted=parent['new_DFT_calls_attempted'],
                    native_runner_invocations=parent['native_runner_invocations'], elapsed_seconds=parent['elapsed_seconds'],
                    application_allocated_CPU_seconds=parent['application_allocated_CPU_seconds']))
                if parent['status'] != 'complete':
                    result['rows'][tid] = {'task_id': tid, 'status': 'unavailable', 'parent_result': endpoints[tid],
                                           'reason': 'Failed singleton retained', 'native_evidence': parent['rows'][tid],
                                           'parent_failure': parent.get('failure'), 'runner_error': parent.get('runner_error')}
                    continue
                require(parent['scientific_binding'] == binding(cfg) and parent['selected_tasks'] == [tid]
                        and pcfg['preflight_result'] == parent['preflight_result'], 'Run scientific/preflight binding differs')
                pw = dict(w, cfg=pcfg); preflight_gate(pw)
                native = admit(w, parent['manifest'], tid)
                require(native == parent['rows'][tid] and native['status'] == 'available'
                        and parent['new_DFT_calls_attempted'] == parent['native_runner_invocations'] == 1,
                        'Raw singleton native admission does not replay')
                result['rows'][tid] = native
            result['comparison'] = comparison(w, result['rows'])
            result['status'] = 'complete' if result['comparison']['status'] == 'complete' else 'unavailable'
        if args.stage != 'preflight':
            for tid, native in result['rows'].items():
                if native['status'] == 'available':
                    result['force_reports'][tid] = write(out/(tid+'_FORCES.json'), force_report(w, tid, native))
        result['full_four_available'] = all(result['rows'][k]['status'] == 'available' for k in IDS)
    except Exception as exc:
        result.update(status='execution_failure', failure=repr(exc), traceback=traceback.format_exc())
        result['comparison'] = {'status': 'unavailable', 'reason': 'Execution/admission failed; no survivor comparison'}
    result.update(elapsed_seconds=time.monotonic()-tick, peak_parent_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                  parent_RSS_is_not_MPI_peak=True, physical_state_qualified=False, force_consistency_qualified=False)
    result['application_allocated_CPU_seconds'] = result['elapsed_seconds']*int(os.environ.get('SLURM_CPUS_ON_NODE', '0'))
    result['scheduler_costs'] = 'Use actual completion/accounting receipts; application time is not scheduler elapsed or utilization.'
    write(out/'RESULT.json', result)
    return 0 if result['status'] in ('complete', 'preflight_passed') else 1


if __name__ == '__main__':
    raise SystemExit(main())
