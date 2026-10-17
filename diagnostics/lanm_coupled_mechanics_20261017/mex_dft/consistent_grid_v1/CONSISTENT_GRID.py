#!/usr/bin/env python3
"""Separate three-cell Mex consistent-COSX test; immutable original helpers.

Common config keys are BIND_KEYS. Run needs preflight_result; +/- additionally
need central_result and continuation_decision. Collect takes endpoint_results
with all three IDS keys, using null for unavailable cells. No stage submits.
"""
import argparse
import copy
from decimal import Decimal
import hashlib
import importlib
import importlib.util
import json
import math
import os
from pathlib import Path
import re
import resource
import sys
import time
import traceback

PROTOCOL = 'mex_consistent_COSX_La_O97_three_v1'
DECISION_PROTOCOL = 'mex_consistent_COSX_displacements_decision_v1'
CENTER = 'MM_O97_La_central'
PAIR = ('MM_O97_La_plus', 'MM_O97_La_minus')
IDS = (CENTER, *PAIR)
BIND_KEYS = ('protocol', 'implementation', 'specification', 'root_decision',
             'diagnosis_result', 'failed_collector', 'qualifier', 'recipe_config')
FROZEN = {
    'specification': 'c18dc071b6622624cedc1dec5fa8e0a0a71a0d11ac0abb2988f8c86a7ff19259',
    'root_decision': '2a2d6917a4e00967b27b8e2a64e860dfdca0bc16a4fb7bff518fed3dff457b5e',
    'diagnosis_result': '08e4d217dbe5a6391afadc09ddc8aa6a84c51e2eb7d316b17e6b9362b2a508a4',
    'failed_collector': '334bb6b4374538a5f293f1e6175340309566de92e42af0588485bca5dfe20953',
    'qualifier': '74b56bf12de0543251c59dcf7a6bb1bb3b84afdedc078ac1ff3045471dec7f48',
    'recipe_config': 'd861ad46247254fb3b050eb527fbbc564130cc24c81ac7cea36e7d14abe1d72a'}
OLD_METHOD = '%method\n DoEQ false\nend\n'
NEW_METHOD = '%method\n DoEQ false\n AngularGridX 2,4,4\n IntAccX -1,-1,-1\nend\n'
SETTINGS = dict(maximum_new_native_calls=3, mpi_ranks=12, workers=1,
                h_A=.001, absolute_kcal_mol_A=.005, relative=.0001,
                no_retries=True, new_MM_calls=0, new_ALPB_calls=0)
COMPONENTS = ('PBE0_embedded', 'MM_ledger', 'fixed_FF_charge_ALPB')
HELPERS = ('metal_environment_frozen_f', 'affordable_response',
           'metal_environment_reference', 'metal_environment_force_assembly')


def require(ok, message):
    if not ok:
        raise ValueError(message)


def pin(path):
    p = Path(path).resolve()
    return {'path': str(p), 'sha256': hashlib.sha256(p.read_bytes()).hexdigest()}


def verify(record):
    require(isinstance(record, dict) and set(record) == {'path', 'sha256'}, 'Required pin missing')
    require(pin(record['path'])['sha256'] == record['sha256'], 'Changed pin: ' + record['path'])
    return Path(record['path'])


def read(record):
    return json.loads(verify(record).read_text())


def write(path, value):
    with Path(path).open('x') as stream:
        json.dump(value, stream, indent=2, allow_nan=False); stream.write('\n')
    return pin(path)


def canonical(value):
    return json.loads(json.dumps(value, allow_nan=False))


def module(name, record):
    spec = importlib.util.spec_from_file_location(name, verify(record))
    obj = importlib.util.module_from_spec(spec); spec.loader.exec_module(obj)
    return obj


def binding(cfg):
    return {k: cfg[k]['sha256'] if k == 'implementation' else cfg[k] for k in BIND_KEYS}


def grid_fingerprints(text):
    """Read physical grid properties, excluding timing/MPI batch partitioning."""
    markers = ('DFT GRID GENERATION', 'COSX GRID GENERATION', 'GRIDX 1', 'GRIDX 2', 'GRIDX 3')
    matches = [list(re.finditer(r'^\s*'+re.escape(s)+r'\s*$', text, re.M)) for s in markers]
    require(all(len(m) == 1 for m in matches), 'Expected one set of XC/COSX grid headers')
    pos = [m[0].start() for m in matches]
    require(pos == sorted(pos), 'Grid header order differs')
    end = re.search(r'^Grids setup in\b', text[pos[-1]:], re.M)
    require(end is not None, 'Complete final grid header unavailable')
    blocks = {'XC': text[pos[0]:pos[1]], 'grid1':text[pos[2]:pos[3]],
              'grid2':text[pos[3]:pos[4]], 'grid3':text[pos[4]:pos[4]+end.start()]}
    patterns = {
        'IntAcc': r'General Integration Accuracy\s+IntAcc\s+\.\.\.\s*(.+)',
        'radial': r'Radial Grid Type\s+RadialGrid\s+\.\.\.\s*(.+)',
        'angular': r'Angular Grid \(max\. ang\.\)\s+AngularGrid\s+\.\.\.\s*(.+)',
        'pruning': r'Angular grid pruning method\s+GridPruning\s+\.\.\.\s*(.+)',
        'weights': r'Weight generation scheme\s+WeightScheme\s*\.\.\.\s*(.+)',
        'BFCut': r'Basis function cutoff\s+BFCut\s+\.\.\.\s*(.+)',
        'WCut': r'Integration weight cutoff\s+WCut\s+\.\.\.\s*(.+)',
        'partial_contraction': r'Partially contracted basis set\s+\.\.\.\s*(.+)',
        'rotational_invariance': r'Rotationally invariant grid construction\s+\.\.\.\s*(.+)',
        'points': r'Total number of grid points\s+\.\.\.\s*(\d+)'}
    result = {}
    for name, block in blocks.items():
        row = {}
        for key, pattern in patterns.items():
            hits = re.findall(r'^\s*'+pattern+r'\s*$', block, re.M)
            require(len(hits) == 1, 'Missing/duplicate grid field: '+name+'/'+key)
            row[key] = hits[0].strip()
        require(block.count('Angular grids for H and He will be reduced by one unit') == 1,
                'H/He angular reduction differs')
        row['H_He_reduction'] = True
        if name != 'XC':
            hits = re.findall(r'^UseSFitting\s+\.\.\.\s*(\S+)\s*$', block, re.M)
            require(hits == ['on'], 'COSX fitting differs'); row['UseSFitting'] = 'on'
        result[name] = row
    return result


def validate_grid(text, original):
    observed = grid_fingerprints(text)
    expected = dict(original, grid2=original['grid3'])
    require(observed == expected, 'Executed grid fingerprints differ from frozen grid2=grid3 candidate')
    require('Recomputing exchange energy using gridx3' in text, 'Final integration disabled/absent')
    validate_iteration_limit(text)
    return observed


def validate_iteration_limit(text):
    require(re.findall(r'Maximum # iterations\s+MaxIter\s+\.\.\.\s*(\d+)', text) == ['125'],
            'Native SCF iteration limit differs')


def inputs(cfg, out):
    require(os.environ.get('SLURM_JOB_ID','').isdigit(), 'Every stage requires Slurm')
    require(cfg['protocol'] == PROTOCOL and cfg['implementation']['sha256'] == pin(__file__)['sha256'], 'Protocol/source differs')
    for key in BIND_KEYS[1:]:
        verify(cfg[key])
    require(all(cfg[k]['sha256'] == v for k,v in FROZEN.items()), 'Frozen parent differs')
    decision = read(cfg['root_decision']); diagnosis = read(cfg['diagnosis_result'])
    require(decision['execution_specification']['sha256'] == cfg['specification']['sha256']
            and decision['saved_diagnosis'] == cfg['diagnosis_result']
            and decision['resources']['maximum_new_native_calls'] == 3
            and decision['resources']['no_retries'] and decision['original_failed_qualification_preserved'],
            'Root finite design differs')
    for key in ('draft','execution_specification','saved_diagnosis','scheduler_policy'):
        verify(decision[key])
    require(diagnosis['status'] == 'complete' and diagnosis['original_failure'] == cfg['failed_collector']
            and diagnosis['qualification_status'] == 'failed_derivative_gate', 'Saved diagnosis identity differs')
    failed = read(cfg['failed_collector']); qc = read(cfg['recipe_config'])
    require(failed['config'] == cfg['recipe_config'] and failed['status'] == 'failed_derivative_gate'
            and failed['implementation']['sha256'] == cfg['qualifier']['sha256'], 'Original failed recipe differs')
    qual = module('consistent_grid_original_qualification', cfg['qualifier'])
    original, ref, world, origins, native = qual.inputs(qc)
    prep, full, oldtasks = qual.prepared(qc)
    # Reuse the unchanged raw collector for both old shared components. Its
    # failed derivative is required to replay; it is not promoted or repaired.
    replay = {}; qual.collect(qc, out, replay)
    require(all(failed[k] == v for k,v in replay.items()), 'Original failed collector raw replay differs')
    gradient_helpers = {}
    for name in HELPERS:
        record = world[7]['implementation'][name+'.py']; verify(record)
        obj = importlib.import_module(name)
        require(pin(obj.__file__)['sha256'] == record['sha256'], 'Gradient helper differs: '+name)
        native[name] = obj; gradient_helpers[name] = record
    require(pin(world[12].__file__)['sha256'] == gradient_helpers['metal_environment_force_assembly']['sha256'],
            'Original scatter implementation differs')
    origin = origins['La']; central = origin['reference_evidence']['reference_task']
    source_tasks = {CENTER: dict(central, mapping=origin['physical_mapping']), **{t:oldtasks[t] for t in PAIR}}
    original_outputs = [origin['reference_evidence']['reference_row']['output'],
                        *[failed['rows'][tid]['DFT']['output'] for tid in PAIR]]
    texts = [verify(p).read_text() for p in original_outputs]
    headers = [grid_fingerprints(t) for t in texts]
    require(all(h == headers[0] for h in headers), 'Original three native grids differ')
    h = headers[0]
    require([(h[k]['IntAcc'],h[k]['angular'],h[k]['points']) for k in ('XC','grid1','grid2','grid3')]
            == [('4.959','6 (Lebedev-590)','2341093'),('4.020','2 (Lebedev-110)','256114'),
                ('4.338','3 (Lebedev-194)','549243'),('4.871','4 (Lebedev-302)','1026059')],
            'Original grid inventory differs')
    parser_controls = []
    for record, text in zip(original_outputs, texts):
        validate_iteration_limit(text)
        rejected = False
        try:
            validate_grid(text, h)
        except ValueError as exc:
            require(str(exc) == 'Executed grid fingerprints differ from frozen grid2=grid3 candidate',
                    'Legacy negative parser control failed unexpectedly')
            rejected = True
        require(rejected, 'Old grid2 unexpectedly admitted as new candidate')
        parser_controls.append(dict(output=record, original_grid_and_MaxIter_parse='passed',
                                    rejected_as_new_candidate=True))
    templates = {}
    for tid, task in source_tasks.items():
        for k in ('input','xyz','pointcharges','mapping'):
            verify(task[k])
        text = verify(task['input']).read_text()
        require(text.count(OLD_METHOD) == 1 and text.count('%method') == 1, 'Original method block differs')
        templates[tid] = text.replace(OLD_METHOD, NEW_METHOD)
        require(task['electronic_state'] == central['electronic_state']
                and task['xyz']['sha256'] == central['xyz']['sha256'], 'Fixed La core/state differs')
    require(templates[PAIR[0]] == templates[PAIR[1]] == templates[CENTER].replace(' VeryTightSCF EnGrad\n',' VeryTightSCF\n',1),
            'Only declared EnGrad distinction permitted')
    return dict(cfg=cfg, qc=qc, qual=qual, native=native, np=world[-1], origins=origins, origin=origin,
                prep=prep, full=full, source_tasks=source_tasks, templates=templates,
                original_grids=h, failed=failed, diagnosis=diagnosis, gradient_helpers=gradient_helpers,
                fixtures=read(prep['directions']), parser_controls=parser_controls)


def task_identity(task):
    ignored = {'output_path','engrad_path','pointcharge_gradient_path','prepared_task'}
    return {k:(v['sha256'] if k in ('input','xyz','pointcharges') else v)
            for k,v in task.items() if k not in ignored}


def task_base(w, tid):
    source = w['source_tasks'][tid]
    task = dict(source, task_id=tid, task_type='analytic_gradient' if tid == CENTER else 'single_point',
                source_task=source, numerical_protocol=PROTOCOL)
    task['input'] = {'path':'unwritten', 'sha256':hashlib.sha256(w['templates'][tid].encode()).hexdigest()}
    task['cache_key'] = hashlib.sha256(json.dumps(dict(binding=binding(w['cfg']), task=task_identity(task)),
                                                 sort_keys=True).encode()).hexdigest()
    return task


def stage_task(w, tid, directory):
    require(tid in IDS, 'Unknown three-cell task'); directory.mkdir()
    task = task_base(w, tid)
    for key, name in (('xyz','core.xyz'),('pointcharges','environment.pc')):
        (directory/name).write_bytes(verify(w['source_tasks'][tid][key]).read_bytes()); task[key] = pin(directory/name)
    (directory/'endpoint.inp').write_text(w['templates'][tid]); task['input'] = pin(directory/'endpoint.inp')
    task['output_path'] = str(directory/'endpoint.out')
    if tid == CENTER:
        task.update(engrad_path=str(directory/'endpoint.engrad'), pointcharge_gradient_path=str(directory/'endpoint.runtime.pcgrad'))
    return task


def manifest(w, tasks):
    return dict(w['full'], protocol_id=PROTOCOL, agreement=w['cfg']['specification'], plan=w['cfg']['specification'],
                original_preparation=w['qc']['preparation_result'], source_manifest=w['prep']['full_manifest'],
                implementation=w['qc']['native_implementation'], adapter=w['cfg']['implementation'],
                scientific_binding=binding(w['cfg']), inputs_config=w['cfg']['_config_pin'],
                execution_allowlist=list(IDS), tasks=tasks)


def validate_task(w, task, tid):
    require(task_identity(task) == task_identity(task_base(w,tid)), 'Prepared numerical/source task differs')
    for key in ('xyz','pointcharges'):
        require(verify(task[key]).read_bytes() == verify(w['source_tasks'][tid][key]).read_bytes(), 'Prepared source bytes differ')
    require(verify(task['input']).read_text() == w['templates'][tid], 'Prepared numerical deck differs')
    directory = verify(task['input']).parent
    require(Path(task['output_path']) == directory/'endpoint.out'
            and verify(task['xyz']).parent == verify(task['pointcharges']).parent == directory, 'Task output scope differs')
    if tid == CENTER:
        require(task['engrad_path'] == str(directory/'endpoint.engrad')
                and task['pointcharge_gradient_path'] == str(directory/'endpoint.runtime.pcgrad'), 'Central gradient paths differ')


def preflight(w, out, result):
    tasks = [stage_task(w, tid, out/tid) for tid in IDS]
    for tid, task in zip(IDS, tasks):
        validate_task(w,task,tid)
    result.update(status='prepared', full_manifest=write(out/'FULL_MANIFEST.json',manifest(w,tasks)),
                  task_identities={t['task_id']:task_identity(t) for t in tasks},
                  original_grid_fingerprints=w['original_grids'], native_grid_syntax_executed=False,
                  legacy_parser_controls=w['parser_controls'],
                  original_failed_collector_replay='exact', shared_component_replay='passed',
                  gradient_helpers=w['gradient_helpers'], directions=w['prep']['directions'])


def preflight_gate(w):
    cfg = w['cfg']; saved = read(cfg['preflight_result']); scfg = read(saved['config'])
    require(saved['status'] == 'prepared' and saved['stage'] == 'preflight'
            and saved['scientific_binding'] == binding(cfg) == binding(scfg)
            and saved['new_DFT_calls_attempted'] == saved['native_runner_invocations'] == 0
            and saved['original_failed_collector_replay'] == 'exact'
            and saved['shared_component_replay'] == 'passed'
            and saved['original_grid_fingerprints'] == w['original_grids']
            and saved['legacy_parser_controls'] == w['parser_controls']
            and saved['gradient_helpers'] == w['gradient_helpers'] and saved['directions'] == w['prep']['directions'],
            'Actual passed unchanged zero-call preflight required')
    verify(saved['implementation']); require(saved['implementation']['sha256'] == cfg['implementation']['sha256'], 'Preflight source differs')
    full = read(saved['full_manifest'])
    require([t['task_id'] for t in full['tasks']] == list(IDS), 'Three prepared slots required')
    expected = manifest(w, full['tasks']); expected['inputs_config'] = full['inputs_config']; expected['adapter'] = full['adapter']
    require(full == expected and full['adapter']['sha256'] == cfg['implementation']['sha256'], 'Prepared manifest differs')
    verify(full['adapter']); require(binding(read(full['inputs_config'])) == binding(cfg), 'Manifest source config differs')
    for t in full['tasks']:
        validate_task(w,t,t['task_id'])
        require(saved['task_identities'][t['task_id']] == task_identity(t), 'Prepared identity receipt differs')
    return saved, full


def admit(w, manifest_pin, tid):
    m = read(manifest_pin)
    require(len(m['tasks']) == 1 and m['tasks'][0]['task_id'] == tid and tid in IDS, 'One declared native task required')
    task = m['tasks'][0]; validate_task(w, task, tid)
    expected = manifest(w, m['tasks']); expected.update(inputs_config=m['inputs_config'],adapter=m['adapter'])
    require(m == expected and m['adapter']['sha256'] == w['cfg']['implementation']['sha256'], 'Executed manifest differs')
    verify(m['adapter']); require(binding(read(m['inputs_config'])) == binding(w['cfg']), 'Execution binding differs')
    row = w['qual'].admit_sp(manifest_pin, w['qc'], w['native'])
    try:
        if row['status'] != 'available':
            return canonical(row)
        text = verify(row['output']).read_text()
        row['grid_fingerprints'] = validate_grid(text, w['original_grids'])
        diag = module('consistent_grid_saved_energy_parser', w['diagnosis']['implementation'])
        energies = diag.extract(row['output'])
        require(energies['fields']['final']['value'] == row['printed_energy'], 'Final energy parser differs')
        energies['printed_closures'] = {
            'post_minus_last_cycle_minus_printed_COSX': diag.scalar_closure(energies['fields'], {'post_grid':1,'last_cycle':-1,'cosx_printed':-1}),
            'final_minus_post_minus_printed_D4': diag.scalar_closure(energies['fields'], {'final':1,'post_grid':-1,'D4_printed':-1})}
        row['final_grid_energy_diagnostics'] = energies
        receipt = read(row['receipt'])
        require(receipt['manifest'] == manifest_pin and receipt['task_id'] == tid
                and receipt['artifacts']['xyz'] == task['xyz']
                and receipt['artifacts']['template_input'] == task['input'], 'Native receipt/task link differs')
        verify(task['pointcharges'])
        if tid == CENTER:
            require(receipt['expected_engrad_complete'] and 'Point charge gradient' in text,
                    'Native analytic/field gradient incomplete')
            parsed = w['native']['metal_environment_frozen_f'].parse(task, {'electronic_state':task['electronic_state']})
            require(parsed['energy_hartree'] == row['energy_hartree']
                    and parsed['engrad'] == receipt['artifacts']['engrad'], 'Native analytic energy/receipt differs')
            raw = w['native']['affordable_response'].read_engrad(verify(parsed['engrad']))
            pc = pin(task['pointcharge_gradient_path'])
            pg = w['native']['metal_environment_reference'].read_pcgrad(verify(pc), 1840)
            mp = read(task['mapping']); np = w['np']
            require(task['xyz']['sha256'] == mp['xyz']['sha256']
                    and task['pointcharges']['sha256'] == mp['pointcharges']['sha256']
                    and [a['id'] for a in mp['physical_atoms']] == w['origin']['source_ids'], 'Fresh central map/core/field differs')
            mapped = w['native']['metal_environment_force_assembly'].scatter(mp, raw['gradient_Ha_per_bohr'], pg)
            g = np.asarray(mapped['gradient_kcal_mol_A'])
            require(g.shape == (2049,3) and np.isfinite(g).all(), 'Full finite physical gradient required')
            f = w['fixtures']['MM_O97']
            require(f['id'] == 'A/97//O' and f['index'] == 1019
                    and [i for i,a in enumerate(mp['field']) if a['id'] == f['id']] == [832]
                    and all(a.get('source_id') != f['id'] for a in mp['core'])
                    and all(1019 not in (a['retained_index'],a['omitted_index']) for a in mp['caps']), 'Direct O97 PC-only map differs')
            direct = pg[832]*(627.509474/.529177210903)
            require(np.array_equal(direct,g[1019]), 'Direct native PC-row/assembled gradient differs')
            tokens = [v for line in verify(pc).read_text().splitlines()[1:] for v in line.split()]
            quanta = [float(Decimal(1).scaleb(Decimal(t.replace('D','E')).as_tuple().exponent)) for t in tokens]
            row.update(parsed, status='available', pointcharge_gradient=pc, pointcharge_gradient_units='Hartree/bohr',
                       mapped_physical_gradient=mapped, direct_O97_gradient_kcal_mol_A=direct.tolist(),
                       direct_O97_projection_kcal_mol_A=float(np.dot(direct,f['direction'])),
                       pcgrad_printed_precision=dict(minimum_quantum_Ha_bohr=min(quanta),
                                                     maximum_quantum_Ha_bohr=max(quanta), is_error_bound=False))
    except Exception as exc:
        row.update(status='unavailable', admission_failure=repr(exc), traceback=traceback.format_exc())
    return canonical(row)


def central_assembly(w, row):
    np = w['np']; old = w['origin']
    require(row['status'] == 'available' and row['task']['task_id'] == CENTER, 'Fresh central gradient unavailable')
    components = {k:np.asarray(old['component_gradients_kcal_mol_A'][k]) for k in COMPONENTS}
    components['PBE0_embedded'] = np.asarray(row['mapped_physical_gradient']['gradient_kcal_mol_A'])
    require(all(a.shape == (2049,3) and np.isfinite(a).all() for a in components.values()), 'Central component gradient inventory differs')
    total = sum(components.values())
    native = dict(old['native_energies'], PBE0_embedded_hartree=row['energy_hartree'])
    return dict(metal='La', source_ids=old['source_ids'], physical_atoms=old['physical_atoms'],
                physical_mapping=old['physical_mapping'], electronic_state=old['electronic_state'],
                native_energies=native, component_gradients_kcal_mol_A={k:v.tolist() for k,v in components.items()},
                gradient_kcal_mol_A=total.tolist(), force_kcal_mol_A=(-total).tolist(),
                SCF=row['SCF'], classical_reused=old['classical_reused'], ALPB_reused=old['ALPB_reused'],
                fresh_central_receipt=row['receipt'], fresh_central_pcgrad=row['pointcharge_gradient'],
                shared_origin_result=old['result'], quantity='gradient=dE/dR; force=-gradient',
                QM_MM_classical_Coulomb_included=False, physical_state_qualified=False)


def saved_result(w, record, tid):
    parent = read(record); pcfg = read(parent['config'])
    verify(parent['implementation'])
    require(parent['protocol'] == PROTOCOL and parent['stage'] == 'run' and parent['task_id'] == tid
            and parent['scientific_binding'] == binding(w['cfg']) == binding(pcfg)
            and parent['implementation']['sha256'] == w['cfg']['implementation']['sha256']
            and parent['preflight_result'] == w['cfg']['preflight_result']
            and parent['native_runner_invocations'] == 1 and parent['new_DFT_calls_attempted'] == 1
            and all(parent[k] == 0 for k in ('new_MM_calls','new_ALPB_calls','new_MACE_calls','new_GFN_calls')),
            'Saved endpoint scope/source/call accounting differs')
    if parent['status'] != 'complete':
        return parent, dict(status='unavailable', preserved_result=record, failure=parent.get('error',parent.get('DFT')))
    task = read(parent['task']); m = read(parent['manifest'])
    require(task == m['tasks'][0] and task['prepared_task'] == dict(full_manifest=read(w['cfg']['preflight_result'])['full_manifest'],task_id=tid),
            'Executed task/preflight link differs')
    if tid in PAIR:
        require(parent['central_result'] == pcfg['central_result'] == w['cfg']['central_result']
                and parent['continuation_decision'] == pcfg['continuation_decision'] == w['cfg']['continuation_decision'],
                'Displacement central/cost binding differs')
    row = admit(w,parent['manifest'],tid)
    require(row == parent['DFT'] and row['status'] == 'available', 'Saved raw native admission differs')
    if tid == CENTER:
        require(central_assembly(w,row) == read(parent['assembled']), 'Saved fresh central assembly differs')
    return parent, row


def continuation_gate(w):
    cfg = w['cfg']; parent, central = saved_result(w,cfg['central_result'],CENTER)
    require(parent['status'] == 'complete' and central['status'] == 'available', 'Actual admitted fresh central required')
    decision = read(cfg['continuation_decision'])
    require(decision['protocol'] == DECISION_PROTOCOL and decision['status'] == 'approved_two_displacements'
            and decision['scientific_binding'] == binding(cfg)
            and decision['preflight_result'] == cfg['preflight_result']
            and decision['central_result'] == cfg['central_result']
            and decision['approved_tasks'] == list(PAIR) and decision['maximum_new_DFT_calls'] == 2
            and decision['allow_retry'] is False, 'Actual central continuation decision differs')
    require(decision['native_receipt'] == central['receipt']
            and isinstance(decision['scheduler_receipts'],list) and decision['scheduler_receipts'],
            'Actual central cost evidence required')
    for record in [decision['native_receipt'],*decision['scheduler_receipts']]:
        verify(record)
    return central


def run(w, out, result, tid):
    require(tid in IDS, 'One named cell required'); pf, full = preflight_gate(w)
    center = continuation_gate(w) if tid in PAIR else None
    task = stage_task(w,tid,out/'native')
    task['prepared_task'] = dict(full_manifest=pf['full_manifest'],task_id=tid)
    validate_task(w,task,tid)
    mp = write(out/'MANIFEST.json',manifest(w,[task]))
    result.update(task_id=tid, selected_tasks=[tid], task=write(out/'TASK.json',task), manifest=mp,
                  preflight_result=w['cfg']['preflight_result'],
                  central_result=w['cfg'].get('central_result') if tid in PAIR else None,
                  continuation_decision=w['cfg'].get('continuation_decision') if tid in PAIR else None)
    runtime = module('consistent_grid_native_allocation',w['qc']['helpers']['reference_runtime'])
    allocation = runtime.allocation(w['native']['render_orca_runtime_input'],out)
    result['allocation'] = allocation
    if center:
        for key in ('allocated_cpus','allocated_memory_MiB','maxcore_MB','workers','ranks_per_worker',
                    'local_mpi_slots','reserve_fraction','memory_source','per_rank_decimal_MB_before_rounding'):
            require(allocation[key] == center['runtime_memory_policy'][key], 'Actual allocation differs from central: '+key)
    write(out/'START.json',result)
    try:
        result['native_runner_invocations'] = 1; result['new_DFT_calls_requested'] = 1
        w['native']['run_orca_task_manifest'].run_manifest(verify(mp),orca_path=verify(w['qc']['orca']),
            workers=1,nprocs=12,selected_tasks={tid},
            expected_runner_sha256=w['qc']['native_implementation']['run_orca_task_manifest.py']['sha256'])
    except Exception as exc:
        result['runner_failure'] = repr(exc)
    finally:
        result['new_DFT_calls_attempted'] = int(Path(task['output_path']).exists())
        for key in ('input','xyz','pointcharges'):
            verify(task[key])
    dft = admit(w,mp,tid); result['DFT'] = dft
    result['rows'][tid] = dict(status=dft['status'], task=result['task'], manifest=mp)
    if dft['status'] == 'available':
        if center:
            require(dft['runtime_memory_policy'] == center['runtime_memory_policy'], 'Executed allocation differs from central')
        if tid == CENTER:
            result['assembled'] = write(out/'ASSEMBLED.json',central_assembly(w,dft))
        result['status'] = 'complete'
    else:
        result['status'] = 'unavailable'


def collect(w, out, result):
    preflight_gate(w); cfg = w['cfg']; records = cfg['endpoint_results']
    require(set(records) == set(IDS), 'All three endpoint slots must be retained')
    result.update(preflight_result=cfg['preflight_result'], endpoint_results=records,
                  rows={tid:dict(status='unavailable',result=records[tid]) for tid in IDS}, checks=None)
    if records[CENTER] is not None:
        require(records[CENTER] == cfg.get('central_result'), 'Collection central identity differs')
    admitted = {}; jobs = set(); costs = []
    for tid in IDS:
        if records[tid] is None:
            continue
        try:
            if tid in PAIR:
                continuation_gate(w)
            parent, row = saved_result(w,records[tid],tid)
            require(parent['slurm_job_id'] not in jobs, 'Duplicate native allocation in three-cell panel')
            jobs.add(parent['slurm_job_id'])
            costs.append(dict(task_id=tid,result=records[tid],slurm_job_id=parent['slurm_job_id'],
                              allocation=parent.get('allocation'),native_receipt=row.get('receipt'),
                              molecular_wall_seconds=row.get('molecular_wall_seconds'),
                              application_wall_seconds=parent.get('application_wall_seconds'),
                              native_runner_invocations=parent['native_runner_invocations'],
                              new_DFT_calls_attempted=parent['new_DFT_calls_attempted']))
            result['rows'][tid].update(status=row['status'], DFT=row)
            if row['status'] == 'available':
                admitted[tid] = row
        except Exception as exc:
            result['rows'][tid].update(error=repr(exc),traceback=traceback.format_exc())
    result['native_allocation_evidence'] = costs
    if set(admitted) != set(IDS):
        result.update(status='unavailable',directional_gate_passed=False)
        return
    central = admitted[CENTER]
    require(all(admitted[t]['runtime_memory_policy'] == central['runtime_memory_policy'] for t in PAIR),
            'Matched executed CPU/memory settings required')
    new_origin = central_assembly(w,central); pair = {}
    for tid in PAIR:
        old = w['failed']['rows'][tid]
        native = dict(old['native_energies'], PBE0_embedded_hartree=admitted[tid]['energy_hartree'])
        pair[tid] = dict(status='complete',DFT=admitted[tid],native_energies=native,
                         MM=old['MM'],ALPB=old['ALPB'])
    checks = w['qual'].compare_pair(w['qc'],{'La':new_origin},w['fixtures'],pair,'MM_O97','La')
    result.update(status='complete' if checks['passed'] else 'failed_derivative_gate',
                  checks=checks,directional_gate_passed=checks['passed'],
                  fresh_central_assembly=write(out/'ASSEMBLED_CENTRAL.json',new_origin),
                  shared_displacement_components={t:dict(MM=pair[t]['MM'],ALPB=pair[t]['ALPB']) for t in PAIR},
                  comparison_rule='Unchanged original compare_pair; new central gradient and new final native +/- energies only.')


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--stage',choices=('preflight','run','collect'),required=True)
    ap.add_argument('--task',choices=IDS)
    ap.add_argument('--config',type=Path,required=True); ap.add_argument('--output',type=Path,required=True)
    args = ap.parse_args()
    require(os.environ.get('SLURM_JOB_ID','').isdigit(), 'All stages require Slurm')
    out = args.output.resolve(); out.mkdir(parents=True,exist_ok=False); started = time.monotonic()
    result = dict(protocol=PROTOCOL,stage=args.stage,status='execution_failure',implementation=pin(__file__),
                  config=pin(args.config),slurm_job_id=os.environ['SLURM_JOB_ID'],settings=SETTINGS,
                  declared_tasks=list(IDS),rows={tid:dict(status='not_requested') for tid in IDS},
                  new_DFT_calls_requested=0,new_DFT_calls_attempted=0,native_runner_invocations=0,
                  new_MM_calls=0,new_ALPB_calls=0,new_MACE_calls=0,new_GFN_calls=0,optimization_steps=0,
                  physical_state_qualified=False,global_force_qualified=False,SCF_error_bound_available=False,
                  original_failed_qualification_preserved=True,original_six_and_response_eight_remain_stopped=True,
                  single_step_truncation_qualified=False)
    try:
        cfg = read(result['config']); cfg['_config_pin'] = result['config']
        result['scientific_binding'] = binding(cfg)
        w = inputs(cfg,out)
        if args.stage == 'preflight':
            require(args.task is None, 'Preflight prepares all three only'); preflight(w,out,result)
        elif args.stage == 'run':
            run(w,out,result,args.task)
        else:
            require(args.task is None, 'Collection retains all three'); collect(w,out,result)
    except Exception as exc:
        result.update(status='execution_failure',error=repr(exc),traceback=traceback.format_exc())
    result['application_wall_seconds'] = time.monotonic()-started
    usage = resource.getrusage(resource.RUSAGE_SELF)
    result['measured_parent_CPU_seconds'] = usage.ru_utime+usage.ru_stime
    result['parent_peak_RSS_KiB'] = usage.ru_maxrss
    result['parent_usage_excludes_MPI_children'] = True
    result['interpretation'] = 'One La external-charge direction under a separate numerical recipe; no physical accuracy, stable-state or broader force qualification. Original failed protocol stays failed.'
    write(out/'RESULT.json',result)
    return 0 if result['status'] in ('prepared','complete') else 1


if __name__ == '__main__':
    raise SystemExit(main())
