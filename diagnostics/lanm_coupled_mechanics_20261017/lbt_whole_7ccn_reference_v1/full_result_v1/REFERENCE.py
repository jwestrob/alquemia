#!/usr/bin/env python3
"""Four fixed 7CCN whole-peptide reference cells; all stages require Slurm."""
import argparse
import datetime as dt
import hashlib
import importlib
import importlib.util
import json
import math
import os
from pathlib import Path
import re
import resource
import socket
import sys
import time
import traceback

PROTOCOL = 'lbt_whole_7ccn_reference_v1'
INVENTORY_SHA = '4203d1aabfc4f312b8627469e0b94450927e06b73d38294b3db6a5a6739cb0e2'
REFERENCE_SHA = 'c22c989674e280806abace63c67148b0181caee0b7b20bcd85b8e3d75d6c2a86'
IDS = ('before_La', 'after_La', 'before_Dy', 'after_Dy')
REMAINING = IDS[1:]
EV, HA, BOHR = 23.06054783061903, 627.509474, 0.529177210903
BIND_KEYS = ('protocol', 'implementation', 'specification', 'inventory', 'plan', 'python', 'reference_helper')
CAVEATS = {
    'source': 'Lu-conditioned 7CCN; same peptide as consumed7CCO, a distinct prospective deformation, not a blind protein benchmark.',
    'geometry': 'Exact prior whole-MACE La start and 80-step terminal; terminal bound/admitted but nonconverged. No new geometry or optimization.',
    'terminal_chemistry': 'Conditional free NH3+/COO- and five deprotonated acidic side chains; experimental terminal specification unresolved.',
    'vacuum_dianion': 'Post-added fixed-FF ALPB does not polarize the vacuum density; finite-basis SCF convergence does not prove electron binding or ground-state stability.',
    'frozen_f': 'Physical La1/Dy6; effective singlets, La ECP46/Dy ECP55, no SOC or physical-state qualification.',
    'force': 'Native analytic admission only, no finite-difference or force-accuracy qualification. Both Dy MACE force comparisons unavailable; no new MACE requests.',
    'interpretation': 'Per-metal deformation works and Dy-minus-La work difference; no absolute cross-method energies, affinity, classifier or preference claims.',
}


def require(ok, message):
    if not ok: raise ValueError(message)


def pin(path):
    p = Path(path).resolve()
    with p.open('rb') as h: sha = hashlib.file_digest(h, 'sha256').hexdigest()
    return {'path': str(p), 'sha256': sha}


def verify(record):
    require(pin(record['path'])['sha256'] == record['sha256'], 'Pin differs: '+record['path'])
    return Path(record['path'])


def read(record): return json.loads(verify(record).read_text())


def write(path, value):
    with Path(path).open('x') as h: json.dump(value, h, indent=2, allow_nan=False); h.write('\n')
    return pin(path)


def module(name, record):
    spec = importlib.util.spec_from_file_location(name, verify(record))
    obj = importlib.util.module_from_spec(spec); spec.loader.exec_module(obj)
    return obj


def binding(cfg):
    return {k: cfg[k]['sha256'] if k == 'implementation' else cfg[k] for k in BIND_KEYS}


def pins_in(obj):
    if isinstance(obj, dict):
        if set(obj) == {'path', 'sha256'}: verify(obj)
        else:
            for v in obj.values(): pins_in(v)
    elif isinstance(obj, list):
        for v in obj: pins_in(v)


def inputs(cfg):
    require(os.environ.get('SLURM_JOB_ID', '').isdigit(), 'Slurm required for every stage')
    require(cfg['protocol'] == PROTOCOL and pin(__file__)['sha256'] == cfg['implementation']['sha256'], 'Implementation differs')
    for k in BIND_KEYS[1:]: verify(cfg[k])
    require(cfg['inventory']['sha256'] == INVENTORY_SHA and cfg['reference_helper']['sha256'] == REFERENCE_SHA,
            'Frozen selection/helper differs')
    require(pin(sys.executable)['sha256'] == cfg['python']['sha256'], 'Python executable differs')
    inv = read(cfg['inventory']); old = module('sevenccn_whole_reference_primitives', cfg['reference_helper'])
    require(inv['plan'] == cfg['plan'] and inv['source_id'] == '7CCN' and inv['source_conditioned_metal'] == 'Lu'
            and inv['declared_cells'] == 4 and inv['declared_geometries'] == 2
            and [c['cell_id'] for c in inv['cells']] == list(IDS), 'Four-cell selection differs')
    native = inv['native_reference']; compact = read(native['compact_assets_manifest'])
    require(native['compact_assets_manifest']['sha256'] == old.COMPACT_SHA and native['method'] == old.METHOD
            and native['metal_assets'] == compact['assets'] and native['orca'] == compact['orca'], 'Native recipe differs')
    pins_in(compact['implementation']); pins_in(compact['assets']); verify(compact['orca'])
    sys.path.insert(0, str(verify(compact['implementation']['metal_environment_frozen_f.py']).parent))
    names = ('metal_environment_frozen_f', 'run_orca_task_manifest', 'ggr_sensitivity', 'mace_omol_vacuum',
             'render_orca_runtime_input', 'affordable_response', 'affordable_common', 'metal_environment_electronic_state')
    mods = {n: importlib.import_module(n) for n in names}
    for n, obj in mods.items():
        require(pin(obj.__file__)['sha256'] == compact['implementation'][n+'.py']['sha256'], 'Imported helper differs: '+n)
    runtime = module('sevenccn_native_runtime', native['runtime_reference_helper'])
    model = mods['metal_environment_frozen_f']
    require(model.METHOD == old.METHOD and model.HA_TO_KCAL == HA and model.BOHR_TO_A == BOHR, 'Native method/units differ')
    import numpy as np
    search = read(inv['search']); run = read(inv['search_result']); pool = read(inv['whole_MACE_pool']); prep = read(inv['preparation'])
    require(run['status'] == pool['status'] == search['status'] == 'complete' and run['slurm_job_id'] == '12968'
            and len(run['rows']) == 1 and run['rows'][0]['search'] == inv['search']
            and search['source_id'] == '7CCN' and search['metal'] == 'La' and search['physical_state_admitted']
            and not search['converged'] and search['accepted_steps'] == 80 and search['termination'] == 'accepted_step_limit'
            and search['final_geometry']['pass'] and search['final_geometry']['state_compatible']
            and search['original_source'] == inv['original_source'], 'Actual fixed search selection/state differs')
    require(prep['source_id'] == '7CCN' and prep['source_conditioned_metal'] == 'Lu' and prep['atoms'] == 256
            and prep['charge'] == -2 and prep['protein_charge'] == -5 and prep['sequence'] == 'FIDTNNDGWIEGDELLA'
            and prep['peptide_bond_count'] == 16 and prep['caps'] == prep['waters'] == prep['spectators'] == []
            and prep['target_has_no_covalent_graph_edges'] and prep['covalent_bonds'] == inv['covalent_bonds']
            and prep['native_peptide_system'] == inv['native_peptide_system'], 'Prepared chemistry/topology differs')
    verify(inv['original_source']); verify(inv['covalent_bonds']); verify(inv['native_peptide_system'])
    cells, solvent = {}, {}
    for row in inv['cells']:
        tid, metal = row['cell_id'], row['metal']; p = read(row['point']); atoms = read(row['physical'])
        joined = [r for r in pool['paired_matrix'] if r['source_id'] == '7CCN' and r['candidate_id'] == row['candidate_id']]
        require(len(joined) == 1 and joined[0]['cells'][metal]['point'] == row['point']
                and joined[0]['cells'][metal]['status'] == 'complete', 'Exact common-pool cell differs')
        require(p['status'] == 'complete' and p['source_id'] == '7CCN' and p['metal'] == metal
                and all(p[k] == row[k] for k in ('positions', 'MACE', 'ALPB', 'task', 'xyz'))
                and p['physical_source'] == row['physical'], 'POINT child provenance differs')
        ids = [a['id'] for a in atoms]; x = np.load(verify(row['positions']), allow_pickle=False)
        require(len(atoms) == len(set(ids)) == 256 and ids[-1] == 'metal' and atoms[-1]['element'] == 'M'
                and all(a['kind'] in ('protein_source', 'terminal_completion') for a in atoms[:-1])
                and [(a['id'], a['element']) for a in atoms if a['kind'] == 'terminal_completion']
                    == [('A/1//H2','H'), ('A/1//H3','H'), ('A/17//OXT','O')]
                and x.shape == (256,3) and np.isfinite(x).all() and np.array_equal(x, [a['xyz_A'] for a in atoms]),
                'Whole physical atom/order/coordinate inventory differs')
        symbols = [a['element'] for a in atoms[:-1]]+[metal]; st = old.state(symbols, metal, model.NUMBERS)
        recipe = p['recipe']; task = read(p['task']); raw = read(p['MACE']); q = recipe['fixed_FF_charges_e']
        require(recipe['physical_ids'] == ids and recipe['symbols'] == symbols and recipe['physical_charge'] == -2
                and recipe['physical_multiplicity'] == st['physical_multiplicity'] and recipe['target_index'] == 255
                and recipe['protocol_id'] == 'lbt3_matched_fixedFF_response_v1'
                and recipe['field_and_MM_ledger'] == 'absent by whole-peptide partition'
                and recipe['solvent_definition'] == 'native ALPB water fixed ff19SB plus formal3 target; no SCC'
                and len(q) == 256 and all(math.isfinite(v) for v in q) and q[-1] == 3 and abs(math.fsum(q)+2) <= 1e-5,
                'Whole recipe/charge/state differs')
        require(task['physical_source'] == row['physical'] and task['xyz'] == row['xyz'] and task['metal'] == metal
                and task['source_id'] == '7CCN' and task['charge'] == -2 and task['multiplicity'] == st['physical_multiplicity']
                and task['metal_index'] == 255 and task['atom_count'] == 256 and task['all_electron_count'] == st['all_electron_count']
                and task['representation'] == 'whole' and task['source_cap_map'] is None, 'Whole native task differs')
        xyz = model.xyz(verify(row['xyz']))
        require([a[0] for a in xyz] == symbols and np.array_equal([a[1:] for a in xyz], x), 'Saved native XYZ differs')
        key = hashlib.sha256(json.dumps({'recipe':recipe,'positions_A':x.tolist()}, sort_keys=True,
                                       separators=(',', ':'), allow_nan=False).encode()).hexdigest()
        require(key == p['cache_key'], 'Saved whole cache identity differs')
        sc = raw['input_state_check']
        require(sc['atoms'] == 256 and sc['charge'] == -2 and sc['spin_multiplicity'] == st['physical_multiplicity']
                and sc['selected_metal_index'] == 255 and sc['graph_count'] == 1 and sc['head'] == 'omol'
                and sc['charge_prediction_claimed'] is False
                and sc['metal_neighbor_edge_count'] > 0 and math.isfinite(raw['energy_eV']), 'Saved MACE native state differs')
        need = metal == 'La'
        require(p['compute_force'] is need and row['whole_MACE_force_available'] is need, 'Saved MACE force availability differs')
        c = dict(row=row, atoms=atoms, ids=ids, x=x, symbols=symbols, state=st, charges={'charges_e':q},
                 physical_parent=row['physical'], point=p, raw=raw)
        solvent[tid] = old.admit_alpb(c, {'ALPB_binary':recipe['native_ALPB_binary']}, np, model)
        sol = solvent[tid]
        require(sol['status'] == 'available', 'All four archived solvents required')
        native_e = dict(MACE_eV=raw['energy_eV'], field_Eh=0., MM_kcal_mol=0., ALPB_Eh=sol['energy_hartree'],
                        ALPB_parts_Eh=sol['components_hartree'])
        require(native_e == p['native_energies'] == row['native_energies'] == joined[0]['cells'][metal]['native_energies']
                and p['components'] == {'MACE_eV':raw['energy_eV'], 'fixed_ALPB_hartree':sol['energy_hartree']}, 'Saved native energies differ')
        if need:
            f = np.asarray(raw['forces_eV_A']); g = np.load(verify(p['gradient_kcal_mol_A']), allow_pickle=False)
            require(raw['force_convention'] == 'negative derivative of the same native vacuum energy'
                    and f.shape == g.shape == (256,3) and np.isfinite(f).all() and np.isfinite(g).all()
                    and np.array_equal(-f*EV + np.asarray(sol['gradient_hartree_per_bohr'])*(HA/BOHR), g),
                    'Saved whole raw/composite gradient replay differs')
            c['MACE_force'] = f*EV; c['composite_force'] = -g
        else:
            require(raw['forces_eV_A'] is None and raw['force_convention'] is None
                    and p.get('gradient_kcal_mol_A') is None, 'Energy-only Dy force evidence differs')
        cells[tid] = c
    first = cells['before_La']; recipe0 = first['point']['recipe']
    for c in cells.values():
        require(c['ids'] == first['ids'] and c['charges'] == first['charges']
                and all(c['point']['recipe'][k] == recipe0[k] for k in ('model','MACE_versions','native_ALPB_build','native_ALPB_binary')),
                'Matched whole model/solvent lineage differs')
    verify(recipe0['model']['checkpoint']); pins_in(recipe0['native_ALPB_build']); verify(recipe0['native_ALPB_binary'])
    require(search['start_point'] == first['row']['point'] and search['terminal_point'] == cells['after_La']['row']['point']
            and np.array_equal(first['x'],cells['before_Dy']['x'])
            and np.array_equal(cells['after_La']['x'],cells['after_Dy']['x']), 'Fixed paired search deformation differs')
    chord = cells['after_La']['x']-first['x']; length = float(np.linalg.norm(chord))
    require(math.isfinite(length) and length > 0, 'Deformation chord unavailable')
    directions = {'physical_ids':first['ids'], 'source_search':inv['search'], 'normalization':'Full 256x3 Euclidean norm; no alignment or geometry change',
                  'directions':{'full_deformation_chord':(chord/length).tolist()}, 'chord_length_A':length}
    return dict(cfg=cfg, inv=inv, old=old, compact=compact, modules=mods, runtime=runtime, model=model,
                np=np, cells=cells, solvent=solvent, directions=directions, search=search)


def task_signature(t):
    return {k:v for k,v in t.items() if k not in ('xyz','input','output_path','engrad_path')} | {
        'xyz_sha256':t['xyz']['sha256'], 'input_sha256':t['input']['sha256']}


def prepare(w, out, tids):
    tasks=[]; cfg=w['cfg']; old=w['old']
    for tid in tids:
        c=w['cells'][tid]; d=out/tid; d.mkdir(); xyz=d/'physical.xyz'; deck=d/'endpoint.inp'
        xyz.write_text('256\nFixed 7CCN physical coordinates; no caps or external field\n'+''.join(
            '%s %.17g %.17g %.17g\n'%(s,*v) for s,v in zip(c['symbols'],c['x'])))
        deck.write_text(old.asset_checks(w['compact'], c['row']['metal']))
        t=dict(task_id=tid,source_id='7CCN',geometry=c['row']['geometry'],metal=c['row']['metal'],
               xyz=pin(xyz),input=pin(deck),output_path=str(d/'endpoint.out'),engrad_path=str(d/'endpoint.engrad'),
               task_type='analytic_gradient',charge=-2,multiplicity=1,physical_multiplicity=c['state']['physical_multiplicity'],
               electronic_state=c['state'],physical_parent=c['physical_parent'],positions=c['row']['positions'],
               physical_ids=c['ids'],atom_count=256,target_index=255,ALPB=c['row']['ALPB'],caveats=CAVEATS)
        write(d/'TASK.json',t); tasks.append(t)
    compact=w['compact']
    manifest=dict(protocol_id=PROTOCOL,agreement=cfg['specification'],plan=cfg['plan'],method=old.METHOD,
        adapter=pin(__file__),inventory=cfg['inventory'],inputs_config=cfg['_config_pin'],scientific_binding=binding(cfg),
        compact_manifest=w['inv']['native_reference']['compact_assets_manifest'],assets=compact['assets'],orca=compact['orca'],
        implementation=compact['implementation'],execution_policy=compact['execution_policy'],
        execution_resources={'mpi_ranks':12,'concurrent_tasks':1},tasks=tasks,caveats=CAVEATS,declared_cells=4)
    rec=write(out/'MANIFEST.json',manifest); validate_manifest(w,rec,tids); return rec


def validate_manifest(w, rec, tids):
    m=read(rec); cfg=w['cfg']; compact=w['compact']
    require([t['task_id'] for t in m['tasks']] == list(tids) and m['protocol_id'] == PROTOCOL
            and m['scientific_binding'] == binding(cfg) and binding(read(m['inputs_config'])) == binding(cfg)
            and m['adapter']['sha256'] == cfg['implementation']['sha256'] and m['inventory'] == cfg['inventory']
            and m['agreement'] == cfg['specification'] and m['plan'] == cfg['plan'] and m['method'] == w['old'].METHOD
            and m['compact_manifest'] == w['inv']['native_reference']['compact_assets_manifest']
            and all(m[k] == compact[k] for k in ('orca','assets','implementation','execution_policy'))
            and m['execution_resources'] == {'mpi_ranks':12,'concurrent_tasks':1} and m['declared_cells'] == 4
            and m['caveats'] == CAVEATS, 'Native manifest identity differs')
    verify(m['adapter'])
    for t in m['tasks']:
        c=w['cells'][t['task_id']]; xyz=w['model'].xyz(verify(t['xyz'])); d=verify(t['input']).parent
        require(t['source_id'] == '7CCN' and t['metal'] == c['row']['metal'] and t['geometry'] == c['row']['geometry']
                and t['electronic_state'] == c['state'] and t['charge'] == -2 and t['multiplicity'] == 1
                and t['physical_multiplicity'] == c['state']['physical_multiplicity'] and t['task_type'] == 'analytic_gradient'
                and t['physical_parent'] == c['physical_parent'] and t['positions'] == c['row']['positions']
                and t['physical_ids'] == c['ids'] and t['ALPB'] == c['row']['ALPB'] and t['atom_count'] == 256
                and t['target_index'] == 255 and t['caveats'] == CAVEATS
                and [a[0] for a in xyz] == c['symbols'] and w['np'].array_equal([a[1:] for a in xyz],c['x'])
                and verify(t['input']).read_text() == w['old'].asset_checks(compact,c['row']['metal'])
                and verify(t['xyz']).parent == d and Path(t['output_path']) == d/'endpoint.out'
                and Path(t['engrad_path']) == d/'endpoint.engrad', 'Exact whole native task differs')
    return m


def preflight_gate(w):
    p=read(w['cfg']['preflight_result']); cfg=w['cfg']; pcfg=read(p['config'])
    require(p['status'] == 'preflight_passed' and p['stage'] == 'preflight' and p['protocol'] == PROTOCOL
            and p['implementation']['sha256'] == cfg['implementation']['sha256'] and binding(pcfg) == binding(cfg)
            and p['solvent_admission'] == w['solvent'] and read(p['directions']) == w['directions']
            and all(p[k] == 0 for k in ('new_DFT_calls_attempted','native_runner_invocations','new_ALPB_calls','new_MACE_calls','new_MM_calls','new_GFN_calls')),
            'Actual unchanged zero-call preflight required')
    validate_manifest(w,p['manifest'],IDS)
    return p


def admit(w, manifest_pin, tid):
    validate_manifest(w, manifest_pin, [tid])
    compact, modules, runtime = w['compact'], w['modules'], w['runtime']
    model = modules['metal_environment_frozen_f']; collector = modules['ggr_sensitivity']
    residual = modules['mace_omol_vacuum']
    manifest, rows = collector.executed(verify(manifest_pin))
    require([t['task_id'] for t in manifest['tasks']] == [tid], 'Scout execution denominator differs')
    task = manifest['tasks'][0]; row = rows[tid]
    row.update(task=task, task_id=tid, physical_multiplicity=task['physical_multiplicity'])
    try:
        if Path(task['output_path']).exists():
            row['SCF_evidence'] = runtime.observations(task['output_path'], residual)
        if row['status'] != 'complete':
            return w['old'].canonical(row)
        parsed = model.parse(task, {'electronic_state': task['electronic_state']})
        require(math.isfinite(parsed['energy_hartree']) and w['np'].asarray(parsed['gradient_kcal_mol_per_A']).shape == (256,3)
                and w['np'].isfinite(parsed['gradient_kcal_mol_per_A']).all(),
                'Full finite native gradient unavailable')
        text = Path(task['output_path']).read_text()
        require('INITIAL GUESS: MODEL POTENTIAL' in text and not re.search(r'Reading point charge file|^\s*(?:CPCM|SMD|COSMO) SOLVATION MODEL', text, re.M),
                'Vacuum/guess scope changed')
        require(row['SCF_evidence']['ConvCheckMode'] == 'Total+1el-Energy' and row['SCF_evidence']['ConvForced'] == 0,
                'Native SCF convergence semantics changed')
        receipt = read(row['receipt']); sidecar = read(receipt['artifacts']['runtime_sidecar'])
        require(receipt['orca_executable'] == compact['orca'] and receipt['returncode'] == 0
                and receipt['normal_termination'] and receipt['scf_converged'] and receipt['expected_engrad_complete'], 'Native execution failed')
        policy = sidecar['memory_policy']
        require(policy['workers'] == 1 and policy['ranks_per_worker'] == 12
                and policy['memory_source'] == 'SLURM_MEM_PER_NODE', 'Native allocation-derived memory policy differs')
        row.update(parsed, status='available', runtime_memory_policy=policy,
                   molecular_wall_seconds=(dt.datetime.fromisoformat(receipt['finished_at_utc']) -
                                           dt.datetime.fromisoformat(receipt['started_at_utc'])).total_seconds(), caveats=CAVEATS)
    except Exception as exc:
        row.update(status='unavailable', admission_failure=repr(exc), energy_hartree=None)
    return w['old'].canonical(row)



def saved_run(w, record, tid):
    p=read(record); pcfg=read(p['config']); verify(p['implementation'])
    require(p['protocol'] == PROTOCOL and p['stage'] == 'run' and p['requested_task'] == tid
            and p['implementation']['sha256'] == w['cfg']['implementation']['sha256']
            and binding(pcfg) == binding(w['cfg'])
            and p['new_DFT_calls_attempted'] in (0,1) and p['native_runner_invocations'] in (0,1), 'Singleton parent differs')
    if p['status'] != 'complete':
        return p, {'task_id':tid,'status':'unavailable','parent_result':record,'reason':'Failed singleton retained',
                   'native_evidence':p['rows'][tid],'failure':p.get('failure'),'runner_error':p.get('runner_error')}
    require(p['preflight_result'] == pcfg['preflight_result'] == w['cfg']['preflight_result']
            and tid in pcfg['permitted_execution_tasks'], 'Actual singleton preflight/phase differs')
    preflight_gate(dict(w,cfg=pcfg))
    require(p['selected_tasks'] == [tid] and p['new_DFT_calls_attempted'] == p['native_runner_invocations'] == 1
            and p['scientific_binding'] == binding(w['cfg']), 'Complete singleton call identity differs')
    row=admit(w,p['manifest'],tid)
    require(row['status'] == 'available' and row == p['rows'][tid], 'Saved native admission does not replay')
    return p,row


def continuation_gate(w):
    cfg=w['cfg']; scout,row=saved_run(w,cfg['scout_result'],'before_La'); decision=read(cfg['continuation_decision'])
    require(row['status'] == 'available' and scout['status'] == 'complete', 'Admitted actual before-La scout required')
    require(decision['protocol'] == PROTOCOL and decision['status'] == 'approved_remaining_three'
            and decision['scientific_binding'] == binding(cfg) and decision['scout_result'] == cfg['scout_result']
            and decision['preflight_result'] == cfg['preflight_result'] and decision['approved_tasks'] == list(REMAINING)
            and decision['scout_native_receipt'] == row['receipt'] and decision['measured_cost_reviewed'] is True
            and isinstance(decision['cost_review_reason'],str) and decision['cost_review_reason'].strip(),
            'Actual scout admission and explicit measured-cost decision required')
    costs=read(decision['scout_scheduler_cost_receipt'])
    require(str(costs['job_id']) == str(scout['slurm_job_id']) and costs['state'] == 'COMPLETED'
            and costs['allocated_CPUs'] == 12 and isinstance(costs['elapsed_seconds'],(int,float))
            and math.isfinite(costs['elapsed_seconds']) and costs['elapsed_seconds'] > 0,
            'Actual scout scheduler cost receipt differs')
    return scout,row


def force_report(w,tid,row):
    np=w['np']; c=w['cells'][tid]; sol=w['solvent'][tid]
    g=np.asarray(row['gradient_kcal_mol_per_A']); gs=np.asarray(sol['gradient_kcal_mol_A'])
    require(g.shape == gs.shape == (256,3) and np.isfinite(g).all() and np.isfinite(gs).all(), 'Finite full reference gradient required')
    force={'whole_PBE0_vacuum':-g,'whole_PBE0_plus_fixed_ALPB':-(g+gs)}
    if tid.endswith('_La'):
        force.update(whole_MACE_vacuum=c['MACE_force'],whole_MACE_plus_fixed_ALPB=c['composite_force'])
    def metrics(v):
        norms=np.linalg.norm(v,axis=1)
        return {'maximum_atom_kcal_mol_A':float(np.max(norms)), 'RMS_atom_vector_kcal_mol_A':float(np.sqrt(np.mean(norms**2)))}
    errors={}
    if tid.endswith('_La'):
        errors={kind:metrics(force['whole_MACE_'+kind]-force['whole_PBE0_'+kind]) for kind in ('vacuum','plus_fixed_ALPB')}
    return dict(task_id=tid,quantity='force=-gradient',units='kcal/mol/Angstrom',physical=c['physical_parent'],positions=c['row']['positions'],
        force_table=[dict(index=i,physical_id=name,resname=c['atoms'][i].get('resname'),
                          force_kcal_mol_A={k:v[i].tolist() for k,v in force.items()}) for i,name in enumerate(c['ids'])],
        distributions={k:metrics(v) for k,v in force.items()}, MACE_minus_PBE0_force_errors=errors,
        force_projections={name:{k:float(np.sum(v*np.asarray(d))) for k,v in force.items()}
                           for name,d in w['directions']['directions'].items()},
        MACE_force_comparison_status='available' if tid.endswith('_La') else 'unavailable: archived Dy MACE request energy-only',
        numerical_force_accuracy_qualified=False,caveats=CAVEATS)


def comparison(w,rows):
    if any(rows[k]['status'] != 'available' for k in IDS):
        return {'status':'unavailable','required_cells':list(IDS),'reason':'All four native references required; no survivor work comparison'}
    work={}
    for metal in ('La','Dy'):
        a,b='after_'+metal,'before_'+metal
        dft=(rows[a]['energy_hartree']-rows[b]['energy_hartree'])*HA
        mace=(w['cells'][a]['raw']['energy_eV']-w['cells'][b]['raw']['energy_eV'])*EV
        sol=(w['solvent'][a]['energy_hartree']-w['solvent'][b]['energy_hartree'])*HA
        error=mace-dft
        work[metal]={'PBE0_vacuum':dft,'MACE_vacuum':mace,'fixed_FF_ALPB':sol,
                     'PBE0_plus_fixed_ALPB':math.fsum((dft,sol)),'MACE_plus_fixed_ALPB':math.fsum((mace,sol)),
                     'MACE_minus_PBE0_work_error':error,'within_1_kcal_reporting_resolution':abs(error)<=1.}
    difference={k:work['Dy'][k]-work['La'][k] for k in ('PBE0_vacuum','MACE_vacuum','fixed_FF_ALPB','PBE0_plus_fixed_ALPB','MACE_plus_fixed_ALPB')}
    difference['MACE_minus_PBE0_work_error']=difference['MACE_vacuum']-difference['PBE0_vacuum']
    difference['within_1_kcal_reporting_resolution']=abs(difference['MACE_minus_PBE0_work_error'])<=1.
    return {'status':'complete','work_definition':'after minus before at fixed 7CCN-La selected geometry pair',
            'units':'kcal/mol','per_metal':work,'Dy_minus_La_work_difference':difference,'material_resolution_kcal_mol':1.,
            'resolution_is_not_native_admission_or_force_accuracy_gate':True,'caveats':CAVEATS}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--stage',choices=('preflight','run','collect'),required=True)
    parser.add_argument('--task',choices=IDS); parser.add_argument('--config',required=True); parser.add_argument('--output',required=True)
    args=parser.parse_args(); require(os.environ.get('SLURM_JOB_ID','').isdigit(),'Slurm required')
    require((args.stage == 'run') == (args.task is not None),'Only singleton run takes --task')
    out=Path(args.output).resolve(); out.mkdir(parents=True,exist_ok=False); tick=time.monotonic()
    result=dict(protocol=PROTOCOL,stage=args.stage,status='execution_failure',implementation=pin(__file__),config=pin(args.config),
        slurm_job_id=os.environ['SLURM_JOB_ID'],host=socket.gethostname(),requested_task=args.task,selected_tasks=[],declared_cells=4,
        rows={k:{'task_id':k,'status':'not_requested'} for k in IDS},new_DFT_calls_attempted=0,native_runner_invocations=0,
        new_ALPB_calls=0,new_MACE_calls=0,new_GFN_calls=0,new_MM_calls=0,optimizer_steps=0,new_geometry_count=0,
        full_four_available=False,comparison={'status':'unavailable'},force_reports={},parent_costs=[],caveats=CAVEATS)
    try:
        cfg=read(result['config']); cfg['_config_pin']=result['config']; w=inputs(cfg)
        result.update(scientific_binding=binding(cfg),solvent_admission=w['solvent'],states={k:c['state'] for k,c in w['cells'].items()},
                      selection_evidence=w['inv']['search_evidence'],directions=write(out/'DIRECTIONS.json',w['directions']))
        allowed=cfg['permitted_execution_tasks']
        require(allowed in (['before_La'],list(REMAINING)), 'Explicit scout or remaining-three allowlist required')
        if args.stage == 'preflight':
            result['manifest']=prepare(w,out,IDS); result['status']='preflight_passed'
        elif args.stage == 'run':
            require(args.task in allowed,'Requested task outside declared phase; repeat scout forbidden in continuation')
            pre=preflight_gate(w); result['preflight_result']=cfg['preflight_result']
            if args.task != 'before_La':
                continuation_gate(w); result.update(scout_result=cfg['scout_result'],continuation_decision=cfg['continuation_decision'])
            else:
                require(cfg.get('scout_result') is None and cfg.get('continuation_decision') is None,
                        'Scout already bound; repeat forbidden')
            require(not cfg.get('endpoint_results'), 'Existing endpoint binding forbids execution/retry')
            result['allocation']=w['runtime'].allocation(w['modules']['render_orca_runtime_input'],out)
            record=prepare(w,out,[args.task]); t=read(record)['tasks'][0]
            original=next(x for x in read(pre['manifest'])['tasks'] if x['task_id']==args.task)
            require(task_signature(t)==task_signature(original),'Singleton differs from frozen four-cell preflight')
            result.update(manifest=record,selected_tasks=[args.task]); write(out/'START.json',result)
            try:
                result['native_runner_invocations']+=1
                w['modules']['run_orca_task_manifest'].run_manifest(verify(record),orca_path=verify(w['compact']['orca']),
                    workers=1,nprocs=12,selected_tasks={args.task},expected_runner_sha256=w['compact']['execution_policy']['task_runner']['sha256'])
            except Exception as exc: result['runner_error']=repr(exc)
            finally: result['new_DFT_calls_attempted']=int(Path(t['output_path']).exists())
            result['rows'][args.task]=admit(w,record,args.task)
            result['status']='complete' if result['rows'][args.task]['status']=='available' else 'unavailable'
        else:
            preflight_gate(w); result['preflight_result']=cfg['preflight_result']; endpoints=cfg.get('endpoint_results',{})
            require(isinstance(endpoints,dict) and set(endpoints)<=set(IDS),'Only four declared singleton endpoints accepted')
            if any(k in endpoints for k in REMAINING): continuation_gate(w)
            jobs=set()
            for tid in IDS:
                if tid not in endpoints:
                    result['rows'][tid]={'task_id':tid,'status':'unavailable','reason':'No singleton result supplied'}; continue
                parent,row=saved_run(w,endpoints[tid],tid)
                require(parent['slurm_job_id'] not in jobs,'Duplicate parent allocation'); jobs.add(parent['slurm_job_id'])
                if tid in REMAINING:
                    pcfg=read(parent['config'])
                    require(pcfg['scout_result']==cfg['scout_result'] and pcfg['continuation_decision']==cfg['continuation_decision'],
                            'Remaining run scout/cost decision differs')
                    continuation_gate(dict(w,cfg=pcfg))
                result['rows'][tid]=row
                result['parent_costs'].append(dict(result=endpoints[tid],slurm_job_id=parent['slurm_job_id'],
                    elapsed_seconds=parent['elapsed_seconds'],application_allocated_CPU_seconds=parent['application_allocated_CPU_seconds'],
                    new_DFT_calls_attempted=parent['new_DFT_calls_attempted'],native_runner_invocations=parent['native_runner_invocations']))
            result['comparison']=comparison(w,result['rows'])
            result['status']='complete' if result['comparison']['status']=='complete' else 'unavailable'
        for tid,row in result['rows'].items():
            if row['status']=='available': result['force_reports'][tid]=write(out/(tid+'_FORCES.json'),force_report(w,tid,row))
        result['full_four_available']=all(result['rows'][k]['status']=='available' for k in IDS)
    except Exception as exc:
        result.update(status='execution_failure',failure=repr(exc),traceback=traceback.format_exc())
        result['comparison']={'status':'unavailable','reason':'Execution/admission failure retained'}
    result.update(elapsed_seconds=time.monotonic()-tick,peak_parent_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                  parent_RSS_is_not_MPI_peak=True,physical_state_qualified=False,force_consistency_qualified=False)
    result['application_allocated_CPU_seconds']=result['elapsed_seconds']*int(os.environ.get('SLURM_CPUS_ON_NODE','0'))
    result['scheduler_costs']='Actual completion/accounting receipts remain authoritative; application timing is not utilization.'
    write(out/'RESULT.json',result)
    return 0 if result['status'] in ('complete','preflight_passed') else 1


if __name__=='__main__': raise SystemExit(main())
