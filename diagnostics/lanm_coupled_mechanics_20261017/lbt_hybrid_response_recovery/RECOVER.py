#!/usr/bin/env python3
"""Complete one certified embedded trial, then resume the original finite search."""
import argparse
import copy
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import re
import resource
import time
import traceback

PROTOCOL = 'lbt_response_single_embedded_solver_recovery_v1'
RESPONSE_SHA = '05b1c412feebfc72048405c8723d426168872325ce32d956732f5c898c028981'
SCOUT_SHA = '16b6dd2785f80b5cfdff6edfa2d72c37cb71c1e409312c4e1a582710a3b2e65b'
EVALUATOR_SHA = '686eba84fd3dfcf28724d2d4e7957528c987445d7137b18c88a07dfec055f661'
FAILED_SHA = 'ddcaef5b8396815844164ab4a6a0b0d3cb898ad103efadfacbe9863a217ab6ca'
PRIMARY_SHA = '5110d1300133ed7e071ed39720141f61ce528af3bf18dba9a41d58f71c759bd9'
STATE = ('7CCN', 'La')
REGION = 'baseline'
LIMITS = {'accepted_steps': 80, 'evaluations': 200}
POLICY = {'source_id': '7CCN', 'region': REGION, 'metal': 'La', 'primary_case': 'failed_005',
          'original_accepted_steps': 2, 'original_distinct_evaluations': 3,
          'total_limits': LIMITS, 'maximum_new_composites': 197,
          'missing_trial_ALPB_calls': 1, 'future_requests_use_original_solver': True,
          'stop_at_first_additional_backend_failure': True}


def require(ok, message):
    if not ok:
        raise ValueError(message)


def pin(path):
    p = Path(path).resolve()
    with p.open('rb') as f:
        sha = hashlib.file_digest(f, 'sha256').hexdigest()
    return {'path': str(p), 'sha256': sha}


def verify(p):
    require(pin(p['path'])['sha256'] == p['sha256'], 'Changed pin: ' + p['path'])
    return Path(p['path'])


def read(p):
    return json.loads(verify(p).read_text())


def write(path, value):
    with Path(path).open('x') as f:
        json.dump(value, f, indent=2, allow_nan=False)
        f.write('\n')
    return pin(path)


def module(name, p):
    spec = importlib.util.spec_from_file_location(name, verify(p))
    obj = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(obj)
    return obj


def normalized(value):
    """Verified byte identity permits immutable snapshots of the same inputs."""
    if isinstance(value, dict):
        if set(value) == {'path', 'sha256'}:
            verify(value)
            return {'sha256': value['sha256']}
        return {k: normalized(v) for k, v in value.items()}
    if isinstance(value, list):
        return [normalized(v) for v in value]
    return value


def binding(cfg):
    value = normalized({k: v for k, v in cfg.items() if k != 'preflight_result'})
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False).encode()).hexdigest()


def certify_scc(cfg, out, old):
    """Rerun saved native admission and the three frozen agreement checks; no calls."""
    scout = module('recovery_saved_scc', cfg['scc_implementation'])
    require(cfg['scc_implementation']['sha256'] == SCOUT_SHA, 'SCC diagnostic source differs')
    result = read(cfg['scc_result'])
    scfg = read(result['config'])
    require(result['protocol'] == scout.PROTOCOL and result['stage'] == 'run'
            and result['status'] == 'complete' and result['qualification_passed']
            and result['implementation']['sha256'] == scfg['implementation']['sha256'] == SCOUT_SHA
            and result['slurm_job_id'] == '13036' and result['declared_calls'] == 4
            and result['actual_attempted_GFN_calls'] == result['completed_processes'] == 4
            and result['new_MACE_MM_ALPB_DFT_calls'] == 0
            and normalized(scfg['original_failed_response']) == normalized(cfg['original_failed_result']),
            'Actual four-call SCC qualification unavailable')
    require(scfg['fixtures']['failed']['point'] == old['requests'][-1]
            and scfg['fixtures']['failed']['point']['sha256'] == FAILED_SHA
            and scfg['fixtures']['preceding']['point'] == old['requests'][-2], 'SCC fixtures not linked to exact partial pool')
    d = out / 'scc_saved_replay'
    d.mkdir()
    scratch = {}
    helpers = scout.prepare(scfg, d, scratch)
    backend, gradients, admission, _, np, fixtures = helpers
    require(len(result['cells']) == 4, 'Four SCC diagnostic cells required')
    replayed = {}
    for row, (case_id, fixture_name, damping) in zip(result['cells'], scout.CASES):
        require((row['case_id'], row['fixture'], row['broydamp']) == (case_id, fixture_name, damping)
                and row['status'] == 'accepted' and row['attempted'] is True
                and row['printed_initial_broydamp'] == damping, 'SCC case state differs')
        scout.verify_tree(row)
        require(read(row['receipt']) == {k: v for k, v in row.items() if k != 'receipt'}, 'SCC row receipt differs')
        fixture = fixtures[fixture_name]
        task, execution, native = read(row['task']), read(row['execution']), read(row['native_receipt'])
        staged = task['staged_files']
        directory = verify(row['task']).parent
        require(set(staged) == {'core.xyz', 'parameters', 'external.pc', 'control.inp'}
                and staged['core.xyz']['sha256'] == fixture['task']['xyz']['sha256']
                and staged['parameters']['sha256'] == scfg['runtime']['parameters']['sha256']
                and staged['external.pc']['sha256'] == fixture['task']['staged_files']['external.pc']['sha256']
                and verify(staged['control.inp']).read_text() == scout.controls_text(damping),
                'Changed SCC input bytes/control')
        command = [str(verify(scfg['runtime']['executable']))] + scout.ARGS
        expected = copy.deepcopy(fixture['task'])
        expected.update(xyz=staged['core.xyz'], command=command, staged_files=staged,
                        original_task=fixture['native']['task'], solver_diagnostic_protocol=scout.PROTOCOL,
                        actual_external_pointcharges=staged['external.pc'],
                        solver_control={'broydamp': damping, 'control': staged['control.inp'],
                                        'unchanged_embedding_control': scout.EMBEDDING,
                                        'native_adaptive_behavior_unchanged': True})
        require(task == expected and execution['command'] == command, 'SCC changed settings beyond declared damping')
        values = re.findall(r'Broyden damping\s+([-+\d.Ee]+)', (directory / 'stdout.txt').read_text())
        require(values and float(values[-1]) == damping, 'Actual requested damping not replayed')
        admitted = admission.admit(directory, task, execution, backend=backend, gradients=gradients, np=np)
        require(admitted == row['native'] and all(native[k] == v for k, v in admitted.items())
                and native['status'] == 'accepted' and native['attempted'] is True
                and native['task'] == row['task'] and native['execution'] == row['execution']
                and native['solver_diagnostic_protocol'] == scout.PROTOCOL
                and native['solver_control'] == task['solver_control']
                and native['original_native_result'] == fixture['pins']['native_result']
                and native['original_failure_not_relabelled'] is True, 'Flat native certificate differs')
        scout.verify_tree(native)
        replayed[case_id] = admitted
    checks = []
    for name, left, right, fixture_name in (
            ('target_branch_consistency', 'failed_005', 'failed_010', 'failed'),
            ('preceding_control', 'preceding_005', None, 'preceding'),
            ('origin_control', 'origin_005', None, 'origin')):
        reference = fixtures[fixture_name]['reference'] if right is None else replayed[right]
        check = {'name': name, 'left': left,
                 'right': right or scfg['fixtures'][fixture_name]['native_result'], 'status': 'complete'}
        check.update(scout.compare(replayed[left], reference, fixtures[fixture_name]['mapping'], helpers))
        checks.append(check)
    require(checks == result['checks'] and all(c['pass'] for c in checks), 'Frozen SCC agreement checks do not replay')
    primary = result['cells'][0]['native_receipt']
    require(primary['sha256'] == PRIMARY_SHA and result['cells'][0]['case_id'] == POLICY['primary_case'],
            'Predeclared primary embedded result changed')
    evidence = {'result': cfg['scc_result'], 'implementation': cfg['scc_implementation'],
                'config': result['config'], 'specification': scfg['specification'],
                'primary': primary, 'secondary': result['cells'][1]['native_receipt'],
                'all_four_native_receipts': [r['native_receipt'] for r in result['cells']],
                'checks': checks, 'agreement_gates': scfg['agreement_gates'],
                'new_molecular_calls': 0, 'no_ground_state_or_global_solver_stability_claim': True}
    return evidence, fixtures, helpers


def pending_trial(loaded, checkpoint, previous, failed, old_search):
    """Replay the unchanged direction before any molecular call; exact X required."""
    np, core = loaded['ev'].np, loaded['core']
    x = np.asarray(previous['positions_A'], float)
    g = np.asarray(previous['gradient_kcal_mol_A'], float)
    history = [(np.asarray(h['s']), np.asarray(h['y']), h['rho']) for h in checkpoint['history']]
    pg, _ = core.remove_rigid(x, g)
    d = loaded['geometry'].lbfgs_direction(pg.ravel(), history).reshape(x.shape)
    d, _ = core.remove_rigid(x, d)
    require(float(np.sum(pg*d)) < -1e-14, 'Pending original direction would reset history')
    maximum = float(np.linalg.norm(d, axis=1).max())
    require(maximum >= 1e-14, 'Pending original direction is zero')
    d *= min(1., checkpoint['trust_A']/maximum)
    slope = float(np.sum(g*d))
    require(slope < 0, 'Pending direction is not descending')
    trial = x+d
    positions = read(failed['positions'])
    require(np.array_equal(trial, positions['positions_A']), 'Pending first trial does not exactly replay failed coordinates')
    entry = old_search['trace'][-1]
    require(entry['point'] == checkpoint['point'] and entry['accepted_steps'] == 2
            and entry['distinct_evaluations'] == 2 and entry['trust_A'] == checkpoint['trust_A']
            and len(entry['trials']) == 1 and entry['trials'][0]['multiplier'] == 1.
            and entry['trials'][0]['evaluation_requested'] is True
            and entry['trials'][0]['accepted'] is False and 'point' not in entry['trials'][0],
            'Original interrupted first-trial trace differs')
    geometry = loaded['guards'][STATE].check(trial)
    require(geometry == entry['trials'][0]['geometry'] and geometry['pass'] and geometry['state_compatible'],
            'Pending trial original physical guard does not replay')
    return {'status': 'passed', 'exact_failed_positions': failed['positions'], 'geometry': geometry,
            'multiplier': 1., 'slope_kcal_mol': slope,
            'Armijo_limit_kcal_mol': core.SETTINGS['armijo']*slope,
            'history_entries': len(history), 'history_reset': False, 'new_molecular_calls': 0}


def saved_components(loaded, failed, evidence, helpers):
    """Admit the three successful old components and the certified embedded one."""
    ev, base, np = loaded['ev'], loaded['base'], loaded['ev'].np
    require(failed['status'] == 'unavailable' and failed['identity_sha256'] == ev.identity_hash
            and failed['state'] == ev.states[('7CCN', REGION, 'La')]
            and set(failed['component_receipts']) == {'MACE','MM','MM_execution','GFN_embedded','GFN_isolated'},
            'Failed component inventory/state differs')
    base.verify_tree(failed['component_receipts'])
    mp, task, positions = read(failed['mapping']), read(failed['task']), read(failed['positions'])
    x = np.asarray(positions['positions_A'])
    require(x.shape == (256,3) and np.isfinite(x).all()
            and failed['cache_key'] == ev.key('7CCN', REGION, 'La', x)
            and positions['source_ids'] == [a['id'] for a in mp['physical_atoms']]
            and np.array_equal(x, [a['xyz_A'] for a in mp['physical_atoms']])
            and read(task['physical_source']) == mp['physical_atoms']
            and failed['source_task'] == ev.tasks[('7CCN',REGION,'La')]
            and task['origin_task'] == failed['source_task'], 'Failed physical task differs')
    sy, xyz = ev.admit.xyz(task['xyz'])
    require(sy == ['La']+[a['element'] for a in mp['core'][1:]]
            and np.array_equal(xyz,[a['xyz_A'] for a in mp['core']]), 'Failed core differs')
    receipts = failed['component_receipts']
    ml, mm, isolated = [read(receipts[k]) for k in ('MACE','MM','GFN_isolated')]
    backend, gradients, admission, _, _, _ = helpers
    itask, execution = read(isolated['task']), read(isolated['execution'])
    replay = admission.admit(verify(isolated['task']).parent,itask,execution,backend=backend,gradients=gradients,np=np)
    require(isolated['status'] == 'accepted' and all(isolated[k] == v for k,v in replay.items())
            and itask['kind'] == 'isolated' and itask['origin_task'] == failed['source_task']
            and itask['xyz']['sha256'] == task['xyz']['sha256'], 'Saved isolated component differs')
    require(ml['task'] == failed['task'] and ml['xyz'] == task['xyz'] and ml['physical_source'] == task['physical_source']
            and ml['input_state_check']['charge'] == task['charge'] == -2
            and ml['input_state_check']['spin_multiplicity'] == task['multiplicity'] == 1
            and ml['input_state_check']['atoms'] == task['atom_count'] == 172
            and ml['force_convention'] == 'negative vacuum energy derivative', 'Saved MACE state differs')
    mex = read(receipts['MM_execution'])
    require(mex['returncode'] == 0 and mex['owned_process_group_empty']
            and mm['status'] == 'complete' and mm['positions'] == failed['positions']
            and mm['system'] == ev.partitions[('7CCN',REGION)][0]['systems']['La']['system']
            and mm['implementation']['sha256'] == EVALUATOR_SHA
            and mm['configurations_attempted'] == 1 and mm['group_queries_attempted'] == 4,
            'Saved MM ledger differs')
    embedded = read(evidence['primary'])
    etask = read(embedded['task'])
    require(etask['kind'] == 'embedded' and etask['origin_task'] == failed['source_task']
            and etask['xyz']['sha256'] == task['xyz']['sha256'], 'Certified embedded component differs')
    q = read(task['global_solvent_charges'])['charges_e']
    return {'mapping':mp, 'task':task, 'positions':x, 'charges':q, 'MACE':ml, 'MM':mm,
            'isolated':isolated, 'embedded':embedded, 'solvent_key':ev.solvent_key('7CCN','La',x,q)}


def prepare(cfg, out, result, loaded=None, saved_replay=False):
    require(cfg['protocol'] == PROTOCOL and cfg['implementation']['sha256'] == pin(__file__)['sha256']
            and cfg['recovery_policy'] == POLICY, 'Recovery declaration differs')
    verify(cfg['implementation']); verify(cfg['specification'])
    old = read(cfg['original_failed_result']); rcfg = read(cfg['original_response_config'])
    require(normalized(old['config']) == normalized(cfg['original_response_config'])
            and rcfg['implementation']['sha256'] == old['implementation']['sha256'] == RESPONSE_SHA
            and rcfg['region'] == REGION and old['slurm_job_id'] == '13019'
            and old['status'] == 'unavailable' and old['termination'] == 'execution_failure'
            and old['accepted_steps'] == 2 and old['distinct_evaluations'] == 3
            and old['requests'][-1]['sha256'] == FAILED_SHA, 'Exact failed response lineage differs')
    response = module('embedded_recovery_response',rcfg['implementation'])
    if loaded is None:
        d=out/'original_response_replay';d.mkdir()
        scratch={}
        try: loaded=response.prepare(rcfg,d,scratch)
        finally:
            if '_evaluator' in scratch:result['_evaluator']=scratch['_evaluator']
    ev,core=loaded['ev'],loaded['core'];result['_evaluator']=ev
    require(loaded['base'].pin(loaded['base'].__file__)['sha256'] == EVALUATOR_SHA
            and ev.identity_hash == old['identity_sha256'], 'Original evaluator differs')
    evidence,fixtures,helpers=certify_scc(cfg,out,old)
    failed=read(old['requests'][-1]);search=read(old['search']);checkpoint=read(search['checkpoint'])
    require(search['status']=='unavailable' and search['termination']=='execution_failure'
            and checkpoint['failure']==search['failure'] and checkpoint['failure'] is not None
            and checkpoint['point']==search['terminal_point']==old['terminal_point']
            and checkpoint['initial_point']==search['start_point']==loaded['starts'][STATE]['provenance']
            and checkpoint['accepted_steps']==search['accepted_steps']==2
            and checkpoint['distinct_evaluations']==search['distinct_evaluations']==3
            and checkpoint['trust_A']==.05 and len(checkpoint['history'])==2
            and checkpoint['settings']==search['settings']==core.SETTINGS
            and checkpoint['limits']==search['limits']==LIMITS
            and checkpoint['evaluation_requests']==old['requests']
            and checkpoint['trial_trace']==search['trace'] and len(checkpoint['evaluated_points'])==2,
            'Exact failed checkpoint does not replay')
    require(read(old['binding'])==response.binding(rcfg,loaded,STATE),'Original response binding differs')
    evaluator=response.Evaluator(loaded,STATE,REGION)
    for request,provenance in zip(old['requests'][:-1],checkpoint['evaluated_points']):
        require(provenance['point']==request,'Prior accepted request linkage differs')
        evaluator.read_point(provenance)
    previous=evaluator.read_point(checkpoint['point'])
    require(core.force_metrics(ev.np.asarray(previous['gradient_kcal_mol_A']))==search['final_full_force']
            and evaluator.difference(previous,loaded['starts'][STATE])==search['work_from_start'],
            'Original last accepted full force/work differs')
    pending=pending_trial(loaded,checkpoint,previous,failed,search)
    parts=saved_components(loaded,failed,evidence,helpers)
    if not saved_replay:
        require(parts['solvent_key'] not in ev.solvent_cache,'Pending missing ALPB unexpectedly already cached')
    require(all(ev.cost[k]==0 for k in ('composites_attempted','MACE_calls_attempted','GFN_calls_attempted',
            'MM_configurations_attempted','MM_group_queries_attempted','ALPB_calls_attempted'))
            and ev.native is None and not ev.requests,'Recovery preflight made or retained new requests')
    preflight={'status':'passed','binding_sha256':binding(cfg),'original_failed_result':cfg['original_failed_result'],
        'old_search':old['search'],'old_checkpoint':search['checkpoint'],'old_failed_request':old['requests'][-1],
        'old_last_accepted_point':checkpoint['point'],'saved_scc_qualification':evidence,
        'pending_trial':pending,'limits':LIMITS,'accepted_steps':2,'distinct_evaluations':3,
        'one_missing_ALPB_required':True,'same_Hamiltonian_identity_sha256':ev.identity_hash,
        'new_molecular_calls':0,'new_model_loads':0}
    result['preflight']=write(out/'RECOVERY_PREFLIGHT.json',preflight)
    result.update(original_failed_result=cfg['original_failed_result'],original_cost=old['actual_cost'],
        certification_result=cfg['scc_result'],certification_cost={k:read(cfg['scc_result'])[k] for k in
            ('actual_attempted_GFN_calls','native_child_CPU_seconds','native_wall_seconds')},
        recovery_binding_sha256=binding(cfg),declared_search_states=rcfg['states'],declared_search_denominator=4)
    return {'loaded':loaded,'response':response,'response_config':rcfg,'old':old,'old_search':search,
            'checkpoint':checkpoint,'failed':failed,'evidence':evidence,'parts':parts,'preflight':preflight}


def read_solvent(ev, record, parts):
    """Same saved solvent gates as the qualifier, independent of an assembled POINT."""
    raw=read(record)
    for p in raw['artifacts'].values():verify(p)
    start=read(raw['artifacts']['START.json'])
    require(start['command']==[str(verify(ev.cfg['ALPB_binary'])),str(verify(start['xyz'])),str(verify(start['charges']))],
            'Completion ALPB runtime differs')
    sy,x=ev.admit.xyz(start['xyz'])
    require(sy==[a['element'] for a in ev.partitions[('7CCN',REGION)][1]['physical_atoms'][:-1]]+['La']
            and ev.np.array_equal(x,parts['positions'])
            and [float(v) for v in verify(start['charges']).read_text().split()]==parts['charges']
            and raw['rounding_bounds_hartree']==[0.,0.] and abs(raw['charge_e']+2)<=1e-5,
            'Completion ALPB physical inputs differ')
    g=ev.np.loadtxt(verify(raw['artifacts']['gradient.dat']))
    require(g.shape==(256,3) and ev.np.isfinite(g).all(),'Completion ALPB full gradient absent')
    return {'energy_hartree':raw['total_hartree'],'gradient_hartree_per_bohr':g.tolist(),
            'components_hartree':raw['components_hartree'],'record':record}


def assemble(cfg, bundle, sol):
    """Exact frozen evaluator arithmetic, with the one explicitly certified receipt."""
    loaded,failed,parts=bundle['loaded'],bundle['failed'],bundle['parts']
    ev,base,np=loaded['ev'],loaded['base'],loaded['ev'].np
    mp,ml,mm,embedded,isolated=(parts[k] for k in ('mapping','MACE','MM','embedded','isolated'))
    zero=np.zeros((len(mp['field']),3))
    eg=ev.mapping.scatter(mp,np.asarray(embedded['qm_gradient']['hartree_per_bohr']),
                          np.asarray(embedded['pointcharge_gradient']['hartree_per_bohr']))*base.HA/base.BOHR
    ig=ev.mapping.scatter(mp,np.asarray(isolated['qm_gradient']['hartree_per_bohr']),zero)*base.HA/base.BOHR
    gradients={'MACE_vacuum':ev.mapping.scatter(mp,-np.asarray(ml['forces_eV_A'])*base.EV,zero),
        'GFN2_field_difference':eg-ig,'MM_ledger':np.asarray(mm['gradient_kcal_mol_A']),
        'fixed_FF_charge_ALPB':np.asarray(sol['gradient_hartree_per_bohr'])*base.HA/base.BOHR}
    native={'MACE_eV':ml['energy_eV'],'GFN_embedded_Eh':embedded['energy_hartree'],
        'GFN_isolated_Eh':isolated['energy_hartree'],'MM_groups_kcal_mol':mm['component_energies_kcal_mol'],
        'ALPB_Eh':sol['energy_hartree']}
    components={'MACE_vacuum':native['MACE_eV']*base.EV,
        'GFN2_field_difference':(native['GFN_embedded_Eh']-native['GFN_isolated_Eh'])*base.HA,
        'MM_ledger':mm['energy_kcal_mol'],'fixed_FF_charge_ALPB':native['ALPB_Eh']*base.HA}
    gradient=sum(gradients[k] for k in base.COMPONENTS)
    require(gradient.shape==(256,3) and np.isfinite(gradient).all()
            and all(math.isfinite(v) for v in components.values()),'Nonfinite completed trial')
    point={'positions_A':parts['positions'].tolist(),'gradient_kcal_mol_A':gradient.tolist(),
        'total_energy_kcal_mol':math.fsum(components.values()),'component_energies_kcal_mol':components,
        'component_gradients_kcal_mol_A':{k:gradients[k].tolist() for k in base.COMPONENTS},
        'native_energies':native,'source_id':'7CCN','region':REGION,'metal':'La','state':failed['state']}
    receipts=copy.deepcopy(failed['component_receipts'])
    receipts.update(GFN_embedded=bundle['evidence']['primary'],ALPB=sol['record'])
    certificate={'protocol':PROTOCOL,'binding_sha256':binding(cfg),'policy':POLICY,
        'Hamiltonian_identity_sha256':ev.identity_hash,'Hamiltonian_cache_key':failed['cache_key'],
        'same_Hamiltonian':True,'identical_solver_protocol':False,
        'original_failed_result':cfg['original_failed_result'],'original_failed_POINT':bundle['old']['requests'][-1],
        'original_failed_SEARCH':bundle['old']['search'],'original_failed_checkpoint':bundle['old_search']['checkpoint'],
        'old_failure_remains_unavailable':True,'original_embedded_component':failed['component_receipts']['GFN_embedded'],
        'saved_scc_qualification':bundle['evidence'],'primary_selection':'predeclared failed_005, not energy or speed selection',
        'state':failed['state'],'positions':failed['positions'],'mapping':failed['mapping'],'task':failed['task'],
        'component_receipts':receipts,'one_missing_ALPB_receipt':sol['record'],
        'new_completion_component_calls':{'MACE':0,'GFN2':0,'MM':0,'ALPB':1,'DFT':0},
        'assembly_implementation':cfg['implementation'],'unchanged_Hamiltonian_implementation':ev.cfg['implementation'],
        'solver_exception_scope':'This exact embedded component only; every future request uses the original solver'}
    raw={k:copy.deepcopy(failed[k]) for k in ('cache_key','identity_sha256','state','positions','mapping','source_task','task')}
    raw.update(status='complete',point=point,component_receipts=receipts,assembly_protocol=PROTOCOL,
               original_failed_POINT=bundle['old']['requests'][-1],new_completion_component_calls=certificate['new_completion_component_calls'])
    return certificate,raw


def versioned_checkpoint(bundle, completed):
    cp=copy.deepcopy(bundle['checkpoint'])
    cp['failure']=None
    cp['evaluation_requests'][-1]=completed['provenance']['point']
    # The pending trial was never returned to the optimizer. Its provenance is
    # appended to evaluated_points only by the original loop's first cache hit.
    return cp


def complete_trial(cfg,bundle,out,result):
    previous=read(cfg['preflight_result']);pcfg=read(previous['config'])
    require(previous['protocol']==PROTOCOL and previous['stage']=='preflight' and previous['status']=='complete'
            and previous['implementation']['sha256']==cfg['implementation']['sha256']
            and pcfg.get('preflight_result') is None and binding(pcfg)==binding(cfg)
            and read(previous['preflight'])==bundle['preflight']
            and all(previous['actual_cost'][k]==0 for k in ('composites_attempted','MACE_calls_attempted',
                'GFN_calls_attempted','MM_configurations_attempted','MM_group_queries_attempted','ALPB_calls_attempted')),
            'Actual exact zero-call preflight required')
    result['actual_preflight_result']=cfg['preflight_result']
    loaded,parts=bundle['loaded'],bundle['parts'];ev=loaded['ev']
    d=out/'trial_completion';d.mkdir()
    sol=ev.solvent('7CCN','La',parts['positions'],parts['charges'],parts['task'],d/'ALPB')
    result['completion_ALPB_receipt']=sol['record']
    sol=read_solvent(ev,sol['record'],parts)
    require(ev.cost['ALPB_calls_attempted']==1 and ev.cost['composites_attempted']==0
            and all(ev.cost[k]==0 for k in ('MACE_calls_attempted','GFN_calls_attempted',
                'MM_configurations_attempted','MM_group_queries_attempted')) and ev.native is None,
            'Trial completion repeated a successful component')
    certificate,raw=assemble(cfg,bundle,sol)
    result['solver_exception_certificate']=write(d/'SOLVER_EXCEPTION_CERTIFICATE.json',certificate)
    raw['solver_exception_certificate']=result['solver_exception_certificate']
    raw['point']['solver_execution_provenance']=result['solver_exception_certificate']
    result['assembled_point']=write(d/'ASSEMBLED_POINT.json',raw)
    completed=loaded['qualifier'].restore(ev,result['assembled_point'])
    require(completed['solver_execution_provenance']==result['solver_exception_certificate'], 'Lost exception provenance')
    cp=versioned_checkpoint(bundle,completed)
    result['resume_checkpoint']=write(out/'CERTIFIED_RESUME_CHECKPOINT.json',cp)
    result['original_failed_request']=bundle['old']['requests'][-1]
    result['completion_actual_cost']=copy.deepcopy(ev.cost)
    return completed,cp


def run(cfg,bundle,out,result):
    loaded=bundle['loaded'];ev=loaded['ev'];response=bundle['response']
    completed,cp=complete_trial(cfg,bundle,out,result)
    class PendingEvaluator(response.Evaluator):
        pending_seen=False
        def evaluate(self,x):
            if not self.pending_seen:
                require(ev.np.array_equal(x,completed['positions_A']), 'First resumed trial differs from certified request')
                before=(self.new_evaluations,len(self.evaluation_requests))
                value=super().evaluate(x)
                require(value['provenance']==completed['provenance']
                        and before==(self.new_evaluations,len(self.evaluation_requests)), 'Certified trial was not an exact free cache hit')
                self.pending_seen=True
                return value
            return super().evaluate(x)
    evaluator=PendingEvaluator(loaded,STATE,REGION)
    result['binding']=write(out/'RECOVERY_BINDING.json',{'protocol':PROTOCOL,'binding_sha256':binding(cfg),
        'original_response_binding':read(bundle['old']['binding']),'preflight':result['preflight'],
        'original_failed_request':result['original_failed_request'],'resume_checkpoint':result['resume_checkpoint'],
        'solver_exception_certificate':result['solver_exception_certificate']})
    search=loaded['core'].search(evaluator,loaded['guards'][STATE],loaded['starts'][STATE],out/'search',LIMITS,
                                  loaded['geometry'].lbfgs_direction,resume=cp)
    search.update(source_id=STATE[0],metal=STATE[1],region=REGION,protocol=PROTOCOL,binding=result['binding'],
        original_failed_result=cfg['original_failed_result'],original_failed_request=result['original_failed_request'],
        solver_exception_certificate=result['solver_exception_certificate'],resume_checkpoint=result['resume_checkpoint'],
        start_policy='unchanged original same-region full-parent minimum')
    result['search']=write(out/'SEARCH.json',search)
    result.update(status=search['status'],source_id=STATE[0],metal=STATE[1],region=REGION,
        continuation_execution_and_physical_admission=search['status']=='complete' and search['physical_state_admitted'],
        terminal_point=search['terminal_point'],termination=search['termination'],converged=search['converged'],
        final_full_force=search['final_full_force'],accepted_steps=search['accepted_steps'],
        distinct_evaluations=search['distinct_evaluations'],work_from_start=search['work_from_start'],
        certified_pending_cache_hit=evaluator.pending_seen)
    require(evaluator.pending_seen,'Certified pending trial did not reach original Armijo test')
    require(search['accepted_steps']<=80 and search['distinct_evaluations']==3+ev.cost['composites_attempted']<=200
            and len(ev.requests)==ev.cost['composites_attempted']<=197
            and ev.cost['MACE_calls_attempted']<=197 and ev.cost['GFN_calls_attempted']<=394
            and ev.cost['MM_configurations_attempted']<=197 and ev.cost['MM_group_queries_attempted']<=788
            and ev.cost['ALPB_calls_attempted']<=198,'Unchanged finite component/evaluation budget exceeded')
    checkpoint=read(search['checkpoint'])
    require(checkpoint['evaluation_requests']==cp['evaluation_requests']+ev.requests
            and len(checkpoint['evaluation_requests'])==search['distinct_evaluations'], 'Final cumulative requests differ')
    result['continuation_actual_cost']={k:ev.cost[k]-result['completion_actual_cost'].get(k,0) for k in ev.cost}


def admit_recovered_search(record, response_cfg, loaded, out):
    """Explicit saved-only admission for this new lineage; caller gates all four slots."""
    r=read(record);cfg=read(r['config']);verify(r['implementation'])
    require(r['protocol']==PROTOCOL and r['stage']=='run' and r['status']=='complete'
            and r['continuation_execution_and_physical_admission'] and r['certified_pending_cache_hit']
            and r['implementation']['sha256']==cfg['implementation']['sha256']==pin(__file__)['sha256'],
            'Recovered response execution unavailable')
    rcfg=read(cfg['original_response_config'])
    require(all(normalized(rcfg[k])==normalized(response_cfg[k]) for k in
            ('implementation','specification','plan','optimizer','geometry_helper','parent_pool_result','parent_pool_implementation'))
            and rcfg['region']==response_cfg['region']==REGION and rcfg['states']==response_cfg['states'],
            'Recovered search and paired-pool original model differ')
    out=Path(out);out.mkdir(exist_ok=False);scratch={}
    bundle=prepare(cfg,out,scratch,loaded=loaded,saved_replay=True)
    previous=read(cfg['preflight_result']);pcfg=read(previous['config'])
    require(previous['status']=='complete' and previous['protocol']==PROTOCOL and previous['stage']=='preflight'
            and previous['implementation']['sha256']==cfg['implementation']['sha256']
            and pcfg.get('preflight_result') is None and binding(pcfg)==binding(cfg)
            and read(previous['preflight'])==read(r['preflight'])==bundle['preflight'],
            'Saved zero-call preflight does not replay')
    ev=loaded['ev'];raw=read(r['assembled_point'])
    require(raw['component_receipts']['ALPB']==r['completion_ALPB_receipt'],'Completion solvent receipt differs')
    sol=read_solvent(ev,r['completion_ALPB_receipt'],bundle['parts'])
    certificate,expected=assemble(cfg,bundle,sol)
    require(read(r['solver_exception_certificate'])==certificate,'Explicit solver certificate differs')
    expected['solver_exception_certificate']=r['solver_exception_certificate']
    expected['point']['solver_execution_provenance']=r['solver_exception_certificate']
    require(raw==expected,'Recovered raw component assembly differs')
    completed=loaded['qualifier'].restore(ev,r['assembled_point'])
    resume=versioned_checkpoint(bundle,completed)
    require(read(r['resume_checkpoint'])==resume and r['original_failed_request']==bundle['old']['requests'][-1],
            'Versioned request substitution differs')
    s=read(r['search']);cp=read(s['checkpoint'])
    require(s['status']=='complete' and s['physical_state_admitted'] and s['failure'] is None
            and (s['source_id'],s['metal'],s['region'])==('7CCN','La',REGION)
            and s['protocol']==PROTOCOL and s['termination']!='execution_failure'
            and s['solver_exception_certificate']==r['solver_exception_certificate']
            and cp['failure'] is None and cp['point']==s['terminal_point']==r['terminal_point']
            and cp['initial_point']==s['start_point']==resume['initial_point']
            and cp['settings']==s['settings']==loaded['core'].SETTINGS
            and cp['limits']==s['limits']==LIMITS
            and cp['accepted_steps']==s['accepted_steps']==r['accepted_steps']<=80
            and cp['distinct_evaluations']==s['distinct_evaluations']==r['distinct_evaluations']<=200,
            'Recovered terminal checkpoint differs')
    require(cp['evaluation_requests']==resume['evaluation_requests']+r['requests']
            and len({p['sha256'] for p in cp['evaluation_requests']})==len(cp['evaluation_requests'])==s['distinct_evaluations']
            and len(r['requests'])==r['actual_cost']['composites_attempted']==s['distinct_evaluations']-3,
            'Recovered cumulative work accounting differs')
    evaluator=bundle['response'].Evaluator(loaded,STATE,REGION)
    requested={p['sha256'] for p in cp['evaluation_requests']}
    for request in cp['evaluation_requests']:
        require(read(request)['status']=='complete','Additional failed request prevents pool admission')
        p=loaded['qualifier'].restore(ev,request)
        evaluator.read_point(p['provenance'])
    for provenance in cp['evaluated_points']:
        require(provenance['point']['sha256'] in requested or provenance['cache_key'] in ev.cache,
                'Untracked evaluated point')
        evaluator.read_point(provenance)
    require(cp['evaluated_points'][:3]==resume['evaluated_points']+[completed['provenance']],
            'Pending trial was not first returned as the certified cache hit')
    for request in r['requests']:
        new=read(request)
        require('assembly_protocol' not in new and 'solver_exception_certificate' not in new,
                'Undeclared further solver recovery')
        for kind in ('embedded','isolated'):
            native=read(new['component_receipts']['GFN_'+kind]);task=read(native['task'])
            require('solver_diagnostic_protocol' not in native and 'solver_control' not in task,
                    'Future request changed solver policy')
            if kind=='embedded':
                require(verify(task['staged_files']['control.inp']).read_text()==loaded['base'].CONTROL,
                        'Future embedded control differs from original solver')
    terminal=evaluator.read_point(cp['point']);gate=loaded['guards'][STATE].check(terminal['positions_A'])
    require(gate==s['final_geometry'] and gate['pass'] and gate['state_compatible']
            and loaded['core'].force_metrics(ev.np.asarray(terminal['gradient_kcal_mol_A']))==s['final_full_force']==r['final_full_force']
            and evaluator.difference(terminal,loaded['starts'][STATE])==s['work_from_start']==r['work_from_start'],
            'Recovered terminal geometry/force/work does not replay')
    return r,s,cp,terminal


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--stage',choices=('preflight','run'),required=True)
    parser.add_argument('--config',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    a=parser.parse_args()
    require(os.environ.get('SLURM_JOB_ID','').isdigit(),'Slurm required')
    if a.stage=='run':
        require(os.environ.get('CUDA_VISIBLE_DEVICES') not in (None,'','NoDevFiles'),'Allocated GPU required for unchanged continuation')
    out=a.output.resolve();out.mkdir(parents=True,exist_ok=False)
    tick=time.monotonic();before=resource.getrusage(resource.RUSAGE_SELF)
    result={'protocol':PROTOCOL,'stage':a.stage,'status':'execution_failure','implementation':pin(__file__),
        'config':pin(a.config),'slurm_job_id':os.environ['SLURM_JOB_ID'],'source_id':STATE[0],'metal':STATE[1],
        'region':REGION,'policy':POLICY,'new_DFT_calls':0,'new_derivative_directions':0,
        'old_failed_search_is_not_relabelled':True,'continuation_execution_and_physical_admission':False}
    try:
        cfg=json.loads(a.config.read_text())
        require((cfg.get('preflight_result') is not None)==(a.stage=='run'),'Actual preflight binding required only for run')
        bundle=prepare(cfg,out,result)
        if a.stage=='preflight':result['status']='complete'
        else:run(cfg,bundle,out,result)
        normalized(cfg);verify(result['config'])
    except Exception as exc:
        result.update(status='execution_failure',failure=repr(exc),traceback=traceback.format_exc(),
                      continuation_execution_and_physical_admission=False)
    ev=result.pop('_evaluator',None)
    if ev is not None:
        result.update(actual_cost=ev.cost,requests=ev.requests,points=ev.points)
        if ev.native is not None:
            result.update(native_MACE_calls=ev.native.calls,native_MACE_device=ev.native.device,native_MACE_versions=ev.native.versions)
    after=resource.getrusage(resource.RUSAGE_SELF)
    result.update(elapsed_seconds=time.monotonic()-tick,
        parent_CPU_seconds=(after.ru_utime+after.ru_stime)-(before.ru_utime+before.ru_stime),
        allocation={k:os.environ.get(k) for k in ('SLURM_JOB_ID','SLURM_CPUS_PER_TASK','CUDA_VISIBLE_DEVICES')})
    write(out/'RESULT.json',result)
    print(json.dumps({'status':result['status'],'result':pin(out/'RESULT.json')}))
    return 0 if result['status']=='complete' else 1


if __name__=='__main__':
    raise SystemExit(main())
