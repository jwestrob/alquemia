#!/usr/bin/env python3
"""Saved four-cell development comparison; eight declared slots, Slurm only.

No molecular evaluator, optimizer or subprocess is called. The existing scaffold
reader supplies exact admitted reference forces and their original directions.
"""
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

PROTOCOL = 'whole_base_delta_saved_development_comparison_v1'
LOCAL_IDS = ('before_La', 'after_La', 'before_Dy', 'after_Dy')
IDS = tuple(s+'_'+c for s in ('7CCO', '7CCN') for c in LOCAL_IDS)
COMPONENTS = ('GFN_full', 'MACE_region', 'negative_GFN_region', 'fixed_FF_ALPB')
HA, EV, BOHR = 627.509474, 23.06054783061903, 0.529177210903
THRESHOLD = 1.0
FROZEN = {
    'specification': 'b18386f62d9dd04e92f6baeaba10eb1c3a9d58a1e303886f45cd8a246546db91',
    'inventory': 'd94e1cadcfa12e15013cac8e5fbc1679a10c4c566ab64e47ccae4dec09a4cc6a',
    'reference_result': '9229f7bebb533dd44a024bb80eca10eded2b3a64cc3972b245d65d82b5ecda97',
    'whole_MACE_result': '64f27e852d2403600ef7a73616c7e52cf50cfbace75689a5b932fde926a06ef0',
    'reference_reader': 'd736a58917262f9977c6a3152b5e93feaa2183f9abb55689fef97cc657fef515',
    'reference_reader_config': 'bf88807d5acd6d5e2bd20ee0859b101cc871cbb850f964768dc9a756cbd1fd22',
    'directions': '6c306212746400b7bf25635c0548d45f53d17eb5a39cc95782a1fa72e97d52ed',
}


def require(ok, message):
    if not ok:
        raise ValueError(message)


def pin(path):
    path = Path(path).resolve()
    with path.open('rb') as f:
        return {'path': str(path), 'sha256': hashlib.file_digest(f, 'sha256').hexdigest()}


def verify(record):
    require(isinstance(record, dict) and set(record) == {'path', 'sha256'}, 'Unbound or malformed evidence pin')
    require(pin(record['path'])['sha256'] == record['sha256'], 'Changed evidence: '+record['path'])
    return Path(record['path'])


def read(record):
    return json.loads(verify(record).read_text())


def write(path, value):
    with Path(path).open('x') as f:
        json.dump(value, f, indent=2, allow_nan=False)
        f.write('\n')
    return pin(path)


def module(name, record):
    spec = importlib.util.spec_from_file_location(name, verify(record))
    obj = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(obj)
    return obj


def vector(value, np):
    a = np.asarray(value, dtype=float)
    require(a.shape == (256, 3) and np.isfinite(a).all(), 'Expected finite full256 Cartesian array')
    return a


def xyz(record):
    lines = verify(record).read_text().splitlines()
    n = int(lines[0]); rows = [line.split() for line in lines[2:] if line.strip()]
    require(len(rows) == n and all(len(row) == 4 for row in rows), 'Native XYZ layout differs')
    return [row[0] for row in rows], [[float(v) for v in row[1:]] for row in rows]


def references(cfg, np):
    for key, digest in FROZEN.items():
        require(cfg[key]['sha256'] == digest, 'Frozen evidence differs: '+key)
        verify(cfg[key])
    h = module('saved_scaffold_reference_reader', cfg['reference_reader'])
    require((h.HA, h.EV, h.BOHR) == (HA, EV, BOHR), 'Original units differ')
    rc = read(cfg['reference_reader_config'])
    require(rc['reference_result']['sha256'] == cfg['reference_result']['sha256']
            and rc['model_result']['sha256'] == cfg['whole_MACE_result']['sha256']
            and rc['directions'] == cfg['directions'], 'Reference reader input pins differ')
    # Only this saved-data function is used; do not call main(), mm_inputs(),
    # run_cell(), query(), or any other molecular execution path in the module.
    ref, whole, cells, directions = h.saved_inputs(rc, np)
    require(ref['slurm_job_id'] == '13424' and ref['stage'] == 'collect'
            and ref['new_DFT_calls_attempted'] == 0, 'Admitted full reference required')
    inv = read(cfg['inventory'])
    require([r['cell_id'] for r in inv['cells']] == list(IDS), 'Eight-state inventory order differs')
    by = {r['cell_id']: r for r in inv['cells']}
    for tid, cell in cells.items():
        declared = by['7CCO_'+tid]
        require(declared['physical'] == cell['model_row']['physical_source']
                and declared['positions'] == cell['model_row']['positions']
                and declared['ALPB'] == cell['solvent']['receipt'], 'Candidate/reference physical or ALPB join differs')
    return h, ref, whole, cells, directions, by


def candidate_rows(cfg, declared, cells, result, np):
    parent = read(cfg['candidate_result'])
    verify(cfg['candidate_implementation'])
    require(parent['implementation']['sha256'] == cfg['candidate_implementation']['sha256'],
            'Candidate producer source differs')
    verify(parent['implementation'])
    pc = read(parent['config'])
    verify(pc['specification'])
    require(pc['implementation']['sha256'] == parent['implementation']['sha256']
            and parent['inventory'] == pc['inventory'] == cfg['inventory']
            and pc['design'] == cfg['specification'], 'Candidate scientific input binding differs')
    require(parent['protocol'] == 'lbt_whole_GFN_regional_vacuum_delta_v1'
            and parent['stage'] == 'collect_7CCO', 'Expected saved first-source collection')
    require(parent['new_molecular_calls'] == 0, 'Saved collection must not execute molecular calls')
    rows = parent['rows']
    require(isinstance(rows, list) and len(rows) == 8 and [r['cell_id'] for r in rows] == list(IDS),
            'Eight-cell result denominator/order differs')
    result['parent_evidence'] = {'result': cfg['candidate_result'], 'implementation': cfg['candidate_implementation'],
                               'config': parent['config'], 'preflight': parent['preflight_result']}
    result['parent_cost_evidence'] = {k: parent[k] for k in
        ('parent_results', 'parent_costs', 'actual_cost', 'call_accounting', 'scheduler_costs') if k in parent}
    for r, slot in zip(rows, result['rows']):
        slot.update(producer_status=r['status'], producer_result=cfg['candidate_result'])
    require(all(r['status'] == 'not_requested' for r in rows[4:]), '7CCN must remain unexecuted at development triage')
    if parent['status'] != 'complete_7CCO' or not all(r['status'] == 'complete' for r in rows[:4]):
        result['status'] = 'unavailable_parent'
        result['comparison'] = {'status': 'unavailable', 'reason': 'Complete admitted four-cell candidate/base parent required; no survivor work.'}
        result['challenge_triage'] = {'eligible': None, 'status': 'not_evaluable', 'automatic_submission': False}
        return None
    preflight = read(parent['preflight_result'])
    require(preflight['stage'] == 'preflight' and preflight['status'] == 'preflight_passed'
            and preflight['new_molecular_calls'] == 0 and not preflight['model_loaded'],
            'Actual zero-call eight-state preflight required')
    manifest = read(preflight['manifest'])
    require(len(manifest['cells']) == 8 and [r['cell_id'] for r in manifest['cells']] == list(IDS),
            'Eight-state preflight manifest differs')
    require(manifest['scientific_binding']['inventory'] == cfg['inventory'], 'Preflight inventory binding differs')
    prepared = {r['cell_id']: r for r in manifest['cells']}
    recipe = read(cfg['inventory'])['helpers']
    by = {}
    for r in rows[:4]:
        cid = r['cell_id']; d = declared[cid]; p = prepared[cid]; tid = cid[len('7CCO_'):]
        physical = read(r['physical']); mapping = read(r['mapping'])
        for key in ('source_id', 'geometry', 'metal', 'physical', 'positions', 'mapping'):
            require(r[key] == d[key] == p[key], 'Candidate coordinate/state identity differs: '+key)
        require(r['physical_ids'] == p['physical_ids'] == cells[tid]['ids']
                and physical == mapping['physical_atoms'] == cells[tid]['atoms']
                and r['state'] == p['state'], 'Candidate physical or electronic state differs')
        mult = 1 if r['metal'] == 'La' else 6
        for label, n, ne in (('whole', 256, 750), ('region', 172, 538)):
            state = r['state'][label]
            require(state['atoms'] == n and state['charge'] == -2 and state['GFN_uhf'] == 0
                    and state['GFN_effective_multiplicity'] == 1 and state['GFN_valence_electrons'] == ne
                    and state['MACE_physical_multiplicity'] == mult and not state['physical_state_qualified'],
                    'Candidate state convention differs')
        receipts = r['component_receipts']
        require(set(receipts) == {'GFN_full', 'MACE_region', 'GFN_region', 'ALPB'}
                and receipts['ALPB'] == d['ALPB'] == cells[tid]['solvent']['receipt'], 'Component ledger or ALPB pin differs')
        raw = {k: read(v) for k, v in receipts.items()}
        for key, n in (('GFN_full', 256), ('GFN_region', 172)):
            native = raw[key]
            require(native['status'] == 'accepted' and native['attempted'], 'Unadmitted native GFN component')
            task = read(native['task']); execution = read(native['execution'])
            require(task['kind'] == 'isolated' and task['atom_count'] == n and task['charge'] == -2
                    and task['source_id'] == '7CCO' and task['metal'] == r['metal']
                    and task['uhf'] == 0 and task['multiplicity'] == 1
                    and task['physical_multiplicity'] == mult and task['pointcharge_count'] == 0
                    and task['field_charge_e'] == 0 and execution['returncode'] == 0
                    and execution['owned_process_group_empty'], 'Native GFN task/process differs')
            require(task['runtime'] == recipe['GFN_runtime'] and task['full_corrected_profile'] == recipe['GFN_profile']
                    and task['parameter_loading_source'] == recipe['parameter_loading_source'], 'Native GFN parameter lineage differs')
            if key == 'GFN_full':
                symbols = [a['element'] for a in physical[:-1]]+[r['metal']]
                positions = [a['xyz_A'] for a in physical]
                expected_ne = 750
            else:
                symbols = [r['metal']]+[a['element'] for a in mapping['core'][1:]]
                positions = [a['xyz_A'] for a in mapping['core']]
                expected_ne = 538
            require(xyz(task['xyz']) == (symbols, positions) and task['valence_electrons'] == expected_ne,
                    'Raw native coordinate/element/electron identity differs')
            verify(task['xyz']); verify(native['qm_gradient']['artifact']); verify(native['xtbout'])
            for artifact in native['artifacts'].values():
                verify(artifact)
            a = np.asarray(native['qm_gradient']['hartree_per_bohr'], dtype=float)
            require(a.shape == (n, 3) and np.isfinite(a).all(), 'Native gradient shape differs')
            require('solver_control' not in task, 'No solver exception belongs in the four 7CCO cells')
        if d['local_POINT'] is not None:
            require(receipts['MACE_region'] == d['local_MACE'] and receipts['GFN_region'] == d['local_isolated_GFN2'],
                    'Successful cached local component was replaced')
        require(raw['MACE_region']['model'] == recipe['model']
                and raw['MACE_region']['force_convention'] == 'negative vacuum energy derivative',
                'Regional MACE recipe/force convention differs')
        mf = np.asarray(raw['MACE_region']['forces_eV_A'], dtype=float)
        require(mf.shape == (172, 3) and np.isfinite(mf).all(), 'Regional native MACE force missing')
        require(r['solver_exception_certificate'] is None, '7CCO has no certified solver exception')
        n = {'GFN_full_Eh': raw['GFN_full']['energy_hartree'], 'MACE_region_eV': raw['MACE_region']['energy_eV'],
             'GFN_region_Eh': raw['GFN_region']['energy_hartree'], 'ALPB_Eh': raw['ALPB']['total_hartree']}
        require(n == r['native_energies'] and all(math.isfinite(v) for v in n.values()), 'Raw native energy ledger differs')
        components = r['component_gradients_kcal_mol_A']
        require(set(components) == set(COMPONENTS), 'Gradient component ledger differs')
        g = {k: vector(components[k], np) for k in COMPONENTS}
        # Preserve the producer's declared left-to-right component arithmetic.
        candidate = sum(g[k] for k in COMPONENTS)
        base = g['GFN_full']+g['fixed_FF_ALPB']
        require(np.array_equal(candidate, vector(r['candidate_gradient_kcal_mol_A'], np))
                and np.array_equal(base, vector(r['base_gradient_kcal_mol_A'], np)), 'Candidate/base exact internal gradient closure differs')
        require(np.array_equal(g['GFN_full'], np.asarray(raw['GFN_full']['qm_gradient']['hartree_per_bohr'])*HA/BOHR)
                and np.array_equal(g['fixed_FF_ALPB'], cells[tid]['solvent']['gradient_kcal_mol_A']),
                'Whole native or common solvent gradient replay differs')
        energies = (n['GFN_full_Eh']*HA, n['MACE_region_eV']*EV, -n['GFN_region_Eh']*HA, n['ALPB_Eh']*HA)
        require(r['candidate_energy_kcal_mol'] == math.fsum(energies)
                and r['base_energy_kcal_mol'] == math.fsum((energies[0], energies[3])), 'Candidate/base exact energy closure differs')
        by[cid] = r
        slot = next(s for s in result['rows'] if s['cell_id'] == cid)
        slot.update(status='admitted', physical=r['physical'], positions=r['positions'], mapping=r['mapping'],
                    state=r['state'], component_receipts=receipts, native_energies=n,
                    force_coverage=256, source_state='7CCO La-conditioned; paired target states',
                    physical_state_qualified=False, derivative_qualified=False)
    return by


def energy_work(rows):
    """Subtract native scalars before conversion; do not subtract large totals."""
    work = {}
    for metal in ('La', 'Dy'):
        a = rows['7CCO_after_'+metal]['native_energies']
        b = rows['7CCO_before_'+metal]['native_energies']
        raw = {k: a[k]-b[k] for k in ('GFN_full_Eh', 'MACE_region_eV', 'GFN_region_Eh', 'ALPB_Eh')}
        components = {'GFN_full': raw['GFN_full_Eh']*HA,
                      'MACE_region': raw['MACE_region_eV']*EV,
                      'negative_GFN_region': -raw['GFN_region_Eh']*HA,
                      'fixed_FF_ALPB': raw['ALPB_Eh']*HA}
        work[metal] = {'native_before_after_differences': raw,
                       'component_work_kcal_mol': components,
                       'candidate_kcal_mol': math.fsum(components[k] for k in COMPONENTS),
                       'base_only_kcal_mol': math.fsum(components[k] for k in ('GFN_full', 'fixed_FF_ALPB')),
                       'regional_correction_kcal_mol': math.fsum(components[k] for k in ('MACE_region', 'negative_GFN_region'))}
    return work


def work_comparison(rows, ref, whole):
    work = energy_work(rows)
    reference = {}
    whole_work = {}
    for metal in ('La', 'Dy'):
        a, b = 'after_'+metal, 'before_'+metal
        de = ref['rows'][a]['energy_hartree']-ref['rows'][b]['energy_hartree']
        ds = ref['solvent_admission'][a]['energy_hartree']-ref['solvent_admission'][b]['energy_hartree']
        r = {'PBE0_vacuum_Eh': de, 'fixed_FF_ALPB_Eh': ds,
             'PBE0_vacuum_kcal_mol': de*HA, 'fixed_FF_ALPB_kcal_mol': ds*HA,
             'composite_kcal_mol': math.fsum((de*HA, ds*HA))}
        require(r == ref['comparison']['per_metal_work'][metal], 'Original reference work replay differs')
        reference[metal] = r
        wr = {x['cell_id']: x for x in whole['rows']}
        em = wr[a]['native_energies']['MACE_eV']-wr[b]['native_energies']['MACE_eV']
        es = wr[a]['native_energies']['fixed_ALPB_hartree']-wr[b]['native_energies']['fixed_ALPB_hartree']
        wm = {'MACE_eV': em, 'fixed_ALPB_hartree': es, 'MACE_kcal_mol': em*EV,
              'fixed_ALPB_kcal_mol': es*HA, 'composite_kcal_mol': math.fsum((em*EV, es*HA))}
        require(wm == whole['comparison']['per_metal_work'][metal], 'Original whole-MACE work replay differs')
        whole_work[metal] = wm
        v = work[metal]
        v['discrepancy_from_PBE0_kcal_mol'] = {
            'candidate': v['candidate_kcal_mol']-r['composite_kcal_mol'],
            'base_only': v['base_only_kcal_mol']-r['composite_kcal_mol'],
            'whole_MACE': wm['composite_kcal_mol']-r['composite_kcal_mol']}
        v['candidate_absolute_discrepancy_kcal_mol'] = abs(v['discrepancy_from_PBE0_kcal_mol']['candidate'])
        v['candidate_development_gate_passed'] = v['candidate_absolute_discrepancy_kcal_mol'] <= THRESHOLD
    passed = all(work[m]['candidate_development_gate_passed'] for m in ('La', 'Dy'))
    return {'status': 'complete', 'definition': 'work=after-before at the same fixed source; differential=Dy work-La work',
            'per_metal_work': work, 'reference_work': reference, 'whole_MACE_work': whole_work,
            'Dy_minus_La_work_kcal_mol': {
                name: work['Dy'][name]-work['La'][name]
                for name in ('candidate_kcal_mol', 'base_only_kcal_mol', 'regional_correction_kcal_mol')},
            'reference_Dy_minus_La_work_kcal_mol': reference['Dy']['composite_kcal_mol']-reference['La']['composite_kcal_mol'],
            'whole_MACE_Dy_minus_La_work_kcal_mol': whole_work['Dy']['composite_kcal_mol']-whole_work['La']['composite_kcal_mol'],
            'old_hybrid_La_work': whole['comparison']['hybrid_La_work'], 'old_hybrid_Dy_work': None,
            'challenge_triage': {'eligible': passed, 'threshold_absolute_per_metal_kcal_mol': THRESHOLD,
                'requires_both_metals': True, 'rule': 'abs(candidate work-PBE0 work)<=1 kcal/mol for both La and Dy',
                'next_stage': '7CCN conditionally eligible for root scheduling' if passed else '7CCN candidate challenge closed',
                'consumed_development_fixture': True, 'not_reference_uncertainty_or_force_accuracy': True,
                'automatic_submission': False, 'independent_reference_status': '7CCN pending; not consumed here'}}


def force_comparison(rows, cells, directions, h, cfg, out, np):
    """Same full256, named-terminal and chord definitions as scaffold audit14787."""
    vectors = {}
    summaries = {}
    components = {}
    for metal in ('La', 'Dy'):
        endpoints = {}
        for label in ('before', 'after'):
            tid = label+'_'+metal
            cell, r = cells[tid], rows['7CCO_'+tid]
            endpoints[label] = dict(cell['models'])
            endpoints[label]['candidate'] = -vector(r['candidate_gradient_kcal_mol_A'], np)
            endpoints[label]['base_only'] = -vector(r['base_gradient_kcal_mol_A'], np)
        ids = cells['before_'+metal]['ids']
        summaries[metal] = {}
        for name in ('reference_PBE0_plus_fixed_ALPB', 'reference_PBE0_vacuum',
                     'candidate', 'base_only', 'whole_MACE_composite', 'hybrid_composite'):
            fields = {}
            for label in ('before', 'after'):
                if name in endpoints[label]:
                    f = endpoints[label][name]
                    fields['force_'+label] = f
                    if name not in ('reference_PBE0_plus_fixed_ALPB', 'reference_PBE0_vacuum'):
                        fields['error_'+label] = f-endpoints[label]['reference_PBE0_plus_fixed_ALPB']
            available = 'force_before' in fields and 'force_after' in fields
            if available:
                fields['force_change'] = fields['force_after']-fields['force_before']
                if 'error_before' in fields:
                    fields['error_change'] = fields['error_after']-fields['error_before']
            summaries[metal][name] = {'transition_status': 'available' if available else 'unavailable',
                'before_available': 'force_before' in fields, 'after_available': 'force_after' in fields,
                'summaries': {k: h.describe(v, ids, directions, np) for k, v in fields.items()}}
            vectors[metal+'_'+name] = {k: v.tolist() for k, v in fields.items()}
        components[metal] = {}
        for name in COMPONENTS:
            before = -vector(rows['7CCO_before_'+metal]['component_gradients_kcal_mol_A'][name], np)
            after = -vector(rows['7CCO_after_'+metal]['component_gradients_kcal_mol_A'][name], np)
            fields = {'force_before': before, 'force_after': after, 'force_change': after-before}
            components[metal][name] = {k: h.describe(v, ids, directions, np) for k, v in fields.items()}
            vectors[metal+'_component_'+name] = {k: v.tolist() for k, v in fields.items()}
    arrays = write(out/'FORCE_TRANSITIONS.json', {'physical_ids': cells['before_La']['ids'],
        'units': 'kcal/mol/Angstrom', 'force_definition': 'force=-gradient',
        'error_definition': 'model-reference_PBE0_plus_fixed_ALPB',
        'change_definition': 'after-before at the fixed physical coordinates',
        'directions': cfg['directions'], 'arrays': vectors})
    return {'full256_arrays': arrays, 'per_metal_models': summaries, 'component_forces': components,
            'physical_atom_count': 256, 'reference_reader': cfg['reference_reader'],
            'error_definition': 'model-reference_PBE0_plus_fixed_ALPB',
            'change_definition': 'after-before',
            'missing_comparators': {'before_Dy_whole_MACE_force': 'Archived energy-only request',
                                   'after_Dy_hybrid_force_and_work': 'Closed hybrid branch; no paired cell'},
            'force_accuracy_qualified': False, 'not_a_Hessian_or_unique_stiffness_decomposition': True}


def report(result, out):
    lines = ['# Whole-base regional-delta: consumed 7CCO comparison', '',
             'Eight declared states; four development cells and four conditional 7CCN challenge slots.',
             'Work is after minus before. Energy units are kcal/mol; force units are kcal/mol/Angstrom.', '']
    if result['status'] != 'complete':
        lines += ['**No complete work comparison or challenge triage is admitted.**', '',
                  result.get('failure', result.get('comparison', {}).get('reason', 'Incomplete inputs.')), '']
    else:
        comp = result['comparison']
        lines += ['## Deformation work', '', '| Quantity | La | Dy |', '|---|---:|---:|']
        for key in COMPONENTS:
            lines.append('| '+key+' | '+' | '.join(f"{comp['per_metal_work'][m]['component_work_kcal_mol'][key]:.6f}" for m in ('La', 'Dy'))+' |')
        for label, key in (('Candidate', 'candidate_kcal_mol'), ('Base only', 'base_only_kcal_mol')):
            lines.append('| '+label+' | '+' | '.join(f"{comp['per_metal_work'][m][key]:.6f}" for m in ('La', 'Dy'))+' |')
        for label, key in (('Whole PBE0 + same ALPB', 'reference_work'), ('Whole MACE + same ALPB', 'whole_MACE_work')):
            lines.append('| '+label+' | '+' | '.join(f"{comp[key][m]['composite_kcal_mol']:.6f}" for m in ('La', 'Dy'))+' |')
        lines += ['', '## Work discrepancies from the whole-PBE0 reference', '',
                  '| Model minus reference | La | Dy |', '|---|---:|---:|']
        for model in ('candidate', 'base_only', 'whole_MACE'):
            lines.append('| '+model+' | '+' | '.join(f"{comp['per_metal_work'][m]['discrepancy_from_PBE0_kcal_mol'][model]:.6f}" for m in ('La', 'Dy'))+' |')
        triage = comp['challenge_triage']
        lines += ['', '**Development triage: '+('passed' if triage['eligible'] else 'failed')+'.** '+triage['next_stage']+'.',
                  'The frozen rule requires absolute candidate work error <=1 kcal/mol for both metals. It is not a reference error bar or force-accuracy criterion.',
                  'Base-only results are reported independently and do not control this candidate gate.', '',
                  '## Force offsets and changes', '',
                  'Errors are model minus PBE0 plus the same ALPB. Change is after minus before; all256 atoms remain represented.', '',
                  '| Metal / model / field | Full256 vector RMS | Phe1-force projection | Chord projection |',
                  '|---|---:|---:|---:|']
        for metal, models in result['force_comparison']['per_metal_models'].items():
            for name in ('candidate', 'base_only', 'whole_MACE_composite', 'hybrid_composite'):
                record = models[name]
                for field in ('error_before', 'error_after', 'error_change'):
                    if field not in record['summaries']:
                        lines.append(f'| {metal} / {name} / {field} | unavailable | — | — |')
                        continue
                    s = record['summaries'][field]
                    lines.append(f"| {metal} / {name} / {field} | {s['atom_vector_RMS_kcal_mol_A']:.6f} | "
                                 f"{s['projections_kcal_mol_A']['initial_hybrid_Phe1_CA_force']:.6f} | "
                                 f"{s['projections_kcal_mol_A']['full_deformation_chord']:.6f} |")
        lines += ['', 'Full vectors, original named terminal-atom summaries and every component force change are preserved in RESULT.json and FORCE_TRANSITIONS.json.']
    lines += ['', '## Denominator and limits', '',
              '7CCN remains pending and conditionally gated here; no pending reference outcome is read or used.',
              'Before-Dy whole-MACE forces and after-Dy old-hybrid forces/work remain unavailable. No missing array is replaced with zero.',
              'The candidate uses a mechanically embedded vacuum regional correction; it has no mutual high/low polarization. No MM or embedded-GFN term is added.',
              'Conditional terminal chemistry, vacuum-dianion/frozen-f/SCF limitations and nonstationary endpoints persist. No independent candidate derivative or force-accuracy qualification, optimizer, affinity or preference claim follows.',
              'Two endpoints do not determine a Hessian or uniquely separate equilibrium-position and stiffness errors.', '',
              'This process makes zero molecular calls. Parent allocation/call records remain separate from this saved analysis.', '']
    (out/'REPORT.md').write_text('\n'.join(lines))
    return pin(out/'REPORT.md')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    require(os.environ.get('SLURM_JOB_ID', '').isdigit(), 'Actual Slurm allocation required for all replay/arithmetic')
    out = args.output.resolve(); out.mkdir(parents=True, exist_ok=False)
    tick, cpu = time.monotonic(), time.process_time()
    result = {'protocol': PROTOCOL, 'status': 'execution_failure', 'config': pin(args.config),
              'implementation': pin(__file__), 'slurm_job_id': os.environ['SLURM_JOB_ID'], 'host': socket.gethostname(),
              'declared_cells': list(IDS), 'new_molecular_calls': 0, 'new_MACE_calls': 0,
              'new_GFN_calls': 0, 'new_DFT_calls': 0, 'new_ALPB_calls': 0, 'new_MM_calls': 0,
              'optimizer_steps': 0, 'new_geometries': 0, 'model_loaded': False,
              'rows': [{'cell_id': cid, 'status': 'unavailable' if cid.startswith('7CCO') else 'pending_conditional_challenge',
                        'reference_status': 'admitted_whole_PBE0' if cid.startswith('7CCO') else 'not_read_pending_challenge'} for cid in IDS],
              'comparison': {'status': 'unavailable'},
              'challenge_triage': {'eligible': None, 'status': 'not_evaluated', 'automatic_submission': False},
              'physical_state_qualified': False, 'derivative_qualified': False, 'force_accuracy_qualified': False,
              'accuracy_or_preference_claim': False}
    code = 1
    try:
        cfg = read(result['config'])
        require(cfg['protocol'] == PROTOCOL and cfg['implementation']['sha256'] == result['implementation']['sha256'],
                'Comparator source/protocol differs')
        require(cfg['challenge_threshold_kcal_mol'] == THRESHOLD and cfg['declared_cells'] == list(IDS),
                'Frozen development rule/denominator differs')
        require(pin(sys.executable)['sha256'] == cfg['python']['sha256'], 'Pinned comparison Python differs')
        verify(cfg['python'])
        import numpy as np
        h, ref, whole, cells, directions, declared = references(cfg, np)
        result.update(reference_result=cfg['reference_result'], whole_MACE_result=cfg['whole_MACE_result'],
                      inventory=cfg['inventory'], specification=cfg['specification'], directions=cfg['directions'],
                      reference_caveats=ref['caveats'],
                      reference_SCF_evidence={k: r['SCF_evidence'] for k, r in ref['rows'].items()},
                      native_reference_states={k: r['electronic_state'] for k, r in ref['rows'].items()},
                      units={'energy': 'native Eh/eV retained; converted kcal/mol', 'gradient': 'kcal/mol/Angstrom',
                             'GFN_ALPB_native_gradient': 'Hartree/bohr', 'MACE_native_force': 'eV/Angstrom'},
                      numpy_version=np.__version__)
        rows = candidate_rows(cfg, declared, cells, result, np)
        if rows is not None:
            result['comparison'] = work_comparison(rows, ref, whole)
            result['force_comparison'] = force_comparison(rows, cells, directions, h, cfg, out, np)
            result['challenge_triage'] = result['comparison']['challenge_triage']
            result['status'] = 'complete'
            code = 0
    except Exception as exc:
        result.update(status='execution_failure', failure=repr(exc), traceback=traceback.format_exc())
        result['comparison'] = {'status': 'unavailable', 'reason': 'Admission/replay failure; no work or challenge conclusion.'}
        result['challenge_triage'] = {'eligible': None, 'status': 'not_evaluable', 'automatic_submission': False}
    result.update(elapsed_seconds=time.monotonic()-tick, process_CPU_seconds=time.process_time()-cpu,
                  peak_parent_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                  allocated_CPUs=int(os.environ.get('SLURM_CPUS_PER_TASK', '1')))
    result['report'] = report(result, out)
    write(out/'RESULT.json', result)
    return code


if __name__ == '__main__':
    raise SystemExit(main())
