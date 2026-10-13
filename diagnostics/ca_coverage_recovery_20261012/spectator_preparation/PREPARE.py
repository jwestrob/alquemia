"""Finite spectator-aware source preparation; scientific work only inside Slurm."""
import argparse
import copy
import hashlib
import importlib.util
import json
import multiprocessing
import os
from pathlib import Path
import resource
import socket
import subprocess
import sys
import time
import traceback
import types
from concurrent.futures import ProcessPoolExecutor, as_completed

PROTOCOL = 'Nikasha_two_Ca_source_explicit_spectator_frozen_envelope4p3_v1'
HERE = Path(__file__).resolve().parent
ROOT = Path('/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs')
THREADS = ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS', 'BLIS_NUM_THREADS', 'OPENMM_CPU_THREADS')
REGRESSION = 'p12293-pqq-la_model__conditioned_Ca__seed-1_sample-0'


def pin(path):
    p = Path(path).resolve()
    return {'path': str(p), 'sha256': hashlib.sha256(p.read_bytes()).hexdigest()}


def verify(p):
    assert pin(p['path']) == p, ('pin changed', p['path'])
    return Path(p['path'])


def read(path):
    return json.loads((verify(path) if isinstance(path, dict) else Path(path)).read_text())


def write(path, value):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('x') as h:
        json.dump(value, h, indent=2, allow_nan=False); h.write('\n')


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(spec); sys.modules[name] = m; spec.loader.exec_module(m)
    return m


def freeze(audits, output, regression_result=None):
    """Path/hash/JSON binding only; does not import a scientific executable."""
    out = Path(output).resolve(); out.mkdir(parents=True, exist_ok=False)
    scope_path = ROOT/'workspaces/orca_parameter_audit_20261012/ca_conditioned_execution_v1/SCIENTIFIC_SCOPE.json'
    scope = read(scope_path); assert len(scope['rows']) == 125
    old_path = ROOT/'workspaces/accommodation_goal_20260920/folds_v1/preparation_reconciled_v1.json'
    old = read(old_path); sources = {r['case_id']: r for r in old['cases']}
    selection_path = ROOT/'workspaces/motion_envelope_20260923/run_v2/SELECTION.json'
    selection = read(selection_path); groups = {g['selection_id']: g for g in selection['groups']}
    anchors_path = ROOT/'workspaces/motion_envelope_20260923/run_v2/PREPARATION.json'
    anchors = {(r['case_id'], r['selection_id']): r for r in read(anchors_path)['cases']}
    primary = {r['source_case_id']: r for r in scope['rows']}
    prior = read(scope['inputs']['preparation']); prior_rows = {r['case_id']: r for r in prior['rows']}

    def request(old_id, source, audit_row=None):
        oldrow = primary[old_id]; group = groups[oldrow['selection_id']]
        anchor = anchors[(group['members'][0], group['selection_id'])]
        assert anchor['status'] == 'prepared' and anchor['mapping_eligible']
        verify(group['union']); verify(anchor['state_signature'])
        return {'original_source_id': old_id, 'source': source, 'group': group,
                'anchor': anchor, 'state_anchor': pin(anchors_path), 'anchor_case_id': anchor['case_id'],
                'frozen_state_key': anchor['state_key'], 'audit_row': audit_row}

    requests = []; audit_pins = []
    for ap in audits:
        audit = read(ap); audit_pins.append(pin(ap))
        assert audit['terminal'] and audit['denominator'] == 5 and audit['unique_PQQ_Ca_sources'] == 5
        assert sorted(r['sample'] for r in audit['rows']) == list(range(5))
        for r in sorted(audit['rows'], key=lambda x: x['sample']):
            assert r['status'] == 'unique_PQQ_Ca_source_geometry' and r['all_source_atoms_preserved_exactly']
            verify(r['full_source']); verify(r['model'])
            s = copy.deepcopy(sources[r['original_source_id']]['source'])
            for key in list(s):
                if key.startswith('archived_') or key in ('source_adapter', 'original_source_structure', 'required_metal_normalization', 'source_readiness'):
                    s.pop(key)
            active = copy.deepcopy(r['active_Ca'])
            s.update(case_id=r['new_source_id'], source_structure=r['full_source'], raw_source_metal=active,
                     metal={**active, 'resname': 'LA'}, original_source_id=r['original_source_id'],
                     source_occupancy_state=r['state_id'], source_conditioning_metal='Ca',
                     canonical_coordinate_match=False, primary_evaluation_pool=False,
                     source_kind='new_two_Ca_occupancy_hypothesis', source_geometry_audit=pin(ap))
            requests.append(request(r['original_source_id'], s, r))
    assert len({r['source']['case_id'] for r in requests}) == len(requests)
    regression_source = copy.deepcopy(sources[REGRESSION]['source'])
    regression_source['source_structure'] = regression_source['original_source_structure']
    for key in list(regression_source):
        if key.startswith('archived_') or key in ('source_adapter', 'original_source_structure'):
            regression_source.pop(key)
    regression = request(REGRESSION, regression_source)
    regression['old_source_preparation'] = sources[REGRESSION]
    regression['old_expanded_preparation'] = prior_rows[REGRESSION]
    regression['old_expanded_preparation_pin'] = scope['inputs']['preparation']
    # Pins are checked for all actually imported repository modules before and after work.
    module_pins = {str(p.resolve()): pin(p) for p in sorted((ROOT/'scripts').glob('*.py'))}
    files = {n: pin(HERE/n) for n in ('PREPARE.py', 'PROTONATE.py', 'EXECUTION_PLAN.md')}
    helpers = {'regression_comparator': pin(ROOT/'diagnostics/orca_parameter_audit_20261012/ca_conditioned_preparation/PREPARE.py')}
    manifest = {'protocol_id': PROTOCOL, 'implementation': files, 'helpers': helpers,
        'root': str(ROOT), 'original_scope': pin(scope_path), 'source_preparation': pin(old_path),
        'original_selection': pin(selection_path), 'original_anchors': pin(anchors_path),
        'canonical_reference': scope['inputs']['reference'], 'config': selection['config'],
        'audits': audit_pins, 'requests': requests, 'denominator': len(requests),
        'regression': regression, 'regression_reuse': pin(regression_result) if regression_result else None,
        'repository_script_pins': module_pins,
        'settings': {'spectator_boundary_A': 4.3, 'PQQ_association_A': 4.0,
          'PDB_max_coordinate_error_A': 0.001, 'regression_tolerance_A': 1e-12,
          'pH': 7.0, 'python_random_seed': 20260914, 'threads_per_preparation': 1,
          'PQQ_state': 'C14H3N2O8_charge_minus3', 'waters': 'dry_exclude_all_source_and_synthetic_waters',
          'context': 'exact_original_canonical_La_union_no_expansion', 'spectator_local_energy': 'not_represented'},
        'energy_calls': 0, 'optimizations': 0, 'optimization_scope': 'no new accommodation optimization; unchanged protonator internally places hydrogens', 'process_start_method': 'spawn', 'production_changed': False,
        'python': pin('/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python')}
    verify(manifest['canonical_reference'])
    write(out/'manifest.json', manifest)
    write(out/'READY.json', {'status': 'metadata_frozen_scientific_preparation_UNRUN', 'manifest': pin(out/'manifest.json'),
        'new_sources': len(requests), 'regression_source': REGRESSION, 'new_energy_calls': 0})
    print(json.dumps({'manifest': pin(out/'manifest.json'), 'new_sources': len(requests), 'status': 'ready_for_Slurm'}))


def allocation(workers):
    job = os.environ.get('SLURM_JOB_ID'); assert job and job.isdigit(), 'actual Slurm allocation required'
    text = subprocess.run(['scontrol', 'show', 'job', '-o', job], capture_output=True, text=True, check=True).stdout
    fields = dict(t.split('=', 1) for t in text.split() if '=' in t)
    count = int(fields['NumCPUs'])
    assert fields['JobState'] == 'RUNNING' and int(fields['NumNodes']) == 1
    assert 1 <= workers <= count <= 12 and len(os.sched_getaffinity(0)) >= workers
    assert int(os.environ['SLURM_CPUS_ON_NODE']) == count
    assert all(os.environ.get(k) == '1' for k in THREADS)
    return {'job_id': job, 'hostname': socket.gethostname(), 'scontrol': text, 'allocated_CPUs': count,
            'workers': workers, 'affinity': sorted(os.sched_getaffinity(0)), 'thread_environment': {k: os.environ[k] for k in THREADS}}


def imported_pins(m):
    used = {}
    for name, module in list(sys.modules.items()):
        path = getattr(module, '__file__', None)
        if path and str(ROOT/'scripts') in str(Path(path).resolve()):
            p = str(Path(path).resolve()); assert p in m['repository_script_pins'], ('unfrozen import', p)
            verify(m['repository_script_pins'][p]); used[name] = m['repository_script_pins'][p]
    return used


def atom_map(path, heavy=False):
    import gemmi
    st = gemmi.read_structure(str(path)); assert len(st) == 1
    records = {}
    for c in st[0]:
        for r in c:
            for a in r:
                if heavy and a.element.name in ('H', 'D'): continue
                key = (c.name, r.seqid.num, r.seqid.icode.strip(), a.name.strip(), a.element.name, a.altloc.strip('\x00 '))
                assert key not in records, ('duplicate source atom', key)
                records[key] = (a.pos.x, a.pos.y, a.pos.z)
    return records


def conserved(before, after, tolerance, heavy=False):
    a, b = atom_map(before, heavy), atom_map(after, heavy)
    assert a.keys() == b.keys(), 'source heavy-atom membership/identity changed'
    maximum = max((abs(x-y) for key in a for x, y in zip(a[key], b[key])), default=0.0)
    assert maximum <= tolerance + 1e-12, ('source coordinates changed', maximum)
    return {'before': pin(before), 'after': pin(after), 'atom_count': len(a), 'maximum_coordinate_error_A': maximum,
            'coordinate_tolerance_A': tolerance, 'identity_element_altloc_preserved': True, 'heavy_only': heavy}


def exact_atom(path, selector):
    import gemmi
    st = gemmi.read_structure(str(path)); assert len(st) == 1
    atoms = [(c, r, a) for c in st[0] for r in c for a in r if
        (c.name, r.seqid.num, r.seqid.icode.strip(), r.name, a.name.strip(), a.element.name) ==
        (selector['chain'], selector['resnum'], selector.get('icode', ''), selector['resname'], selector['atom'], selector['element'])]
    assert len(atoms) == 1, ('selector not unique', selector)
    c, r, atom = atoms[0]; assert len(r) == 1
    return st, r, atom


def alias(source, audit_row, directory):
    import gemmi
    from accommodation_folds import adapt_ca, atoms
    raw = verify(source['source_structure'])
    adapted = adapt_ca(source, directory/'active_alias')
    active_path = verify(adapted['source_structure']); checks = [conserved(raw, active_path, 0)]
    spectator = None
    if audit_row:
        spectator = copy.deepcopy(audit_row['spectator_Ca'])
        st, residue, _ = exact_atom(active_path, spectator)
        before = atoms(st); residue.name = 'CA'
        path = directory/'full_source_two_Ca_aliased.cif'; st.make_mmcif_document().write_file(str(path))
        assert atoms(gemmi.read_structure(str(path))) == before
        checks.append(conserved(active_path, path, 0)); adapted['source_structure'] = pin(path)
        spectator['resname'] = 'CA'
    path = directory/'SOURCE_ALIAS.json'
    write(path, {'protocol_id': PROTOCOL, 'input': source['source_structure'], 'output': adapted['source_structure'],
        'active_before': source['raw_source_metal'], 'active_after': source['metal'],
        'spectator_before': audit_row['spectator_Ca'] if audit_row else None, 'spectator_after': spectator,
        'changes': 'residue aliases only; native elements, atom identity, coordinates and all atoms retained', 'checks': checks})
    adapted['spectator_preparation_alias'] = pin(path)
    return adapted, spectator


def spectator_check(request, prepared, row, spectator, preservation, path):
    import numpy as np
    import gemmi
    from affordable_common import xyz
    core = read(prepared['core']['parent']); full = verify(core['source_structure'])
    st, _, ion = exact_atom(full, spectator); pos = np.asarray([ion.pos.x, ion.pos.y, ion.pos.z])
    ca = [(c, r, a) for c in st[0] for r in c for a in r if a.element.name == 'Ca']
    assert len(ca) == 2
    active = request['source']['metal']; _, _, act = exact_atom(full, active)
    pqq = request['source']['pqq']
    pqq_xyz = [[a.pos.x, a.pos.y, a.pos.z] for c in st[0] for r in c if
        (c.name, r.seqid.num, r.seqid.icode.strip(), r.name) ==
        (pqq['chain'], pqq['resnum'], pqq.get('icode', ''), pqq['resname']) for a in r if a.element.name not in ('H', 'D')]
    assert len(pqq_xyz) == 24
    pqq_array = np.asarray(pqq_xyz)
    distances = [{'chain': c.name, 'resname': r.name, 'atom': a.name,
        'nearest_PQQ_heavy_A': float(np.min(np.linalg.norm(pqq_array-np.asarray([a.pos.x,a.pos.y,a.pos.z]), axis=1)))} for c,r,a in ca]
    assert sum(d['nearest_PQQ_heavy_A'] < 4.0 for d in distances) == 1
    assert float(np.min(np.linalg.norm(pqq_array-np.asarray([act.pos.x,act.pos.y,act.pos.z]),axis=1))) < 4.0
    endpoints = row['representations']['context']['endpoints']
    coordinates = xyz(verify(endpoints['Ca']['xyz']))
    heavy = [(i, x) for i,x in enumerate(coordinates) if x[0] not in ('H','D')]
    all_dist = [(float(np.linalg.norm(pos-np.asarray(x[1:]))), i, x[0]) for i,x in heavy]
    nearest = min(all_dist)
    audit = read(row['representations']['context']['preparation'])
    mapping = audit['mapping']
    def is_spectator(value):
        return (value.get('chain'), value.get('resnum'), value.get('insertion_code',value.get('icode','')),
                value.get('atom'), value.get('element')) == (spectator['chain'],spectator['resnum'],spectator.get('icode',''),spectator['atom'],'Ca')
    def contains(value):
        if isinstance(value,dict): return is_spectator(value) or any(contains(v) for v in value.values())
        if isinstance(value,list): return any(contains(v) for v in value)
        return False
    absent = not contains(mapping)
    anchors = {}
    for name, sel in request['source']['roles'].items():
        points = [[a.pos.x,a.pos.y,a.pos.z] for c in st[0] for r in c if
            (c.name,r.seqid.num,r.seqid.icode.strip(),r.name)==(sel['chain'],sel['resnum'],sel.get('icode',''),sel['resname'])
            for a in r if a.element.name not in ('H','D')]
        assert points
        anchors[name] = float(np.min(np.linalg.norm(np.asarray(points)-pos,axis=1)))
    passed = nearest[0] > 4.3 and absent
    receipt = {'status': 'passed' if passed else 'unsupported', 'passed': passed, 'protocol_id': PROTOCOL,
        'active_selector': active, 'spectator_selector': spectator, 'full_protonated_source': pin(full),
        'Ca_count': len(ca), 'PQQ_association_cutoff_A': 4.0, 'PQQ_association_comparison': 'strictly_less_than', 'metal_PQQ_distances': distances,
        'spectator_to_active_A': float(np.linalg.norm(pos-np.asarray([act.pos.x,act.pos.y,act.pos.z]))),
        'spectator_to_PQQ_A': float(np.min(np.linalg.norm(pqq_array-pos,axis=1))),
        'spectator_to_anchor_heavy_A': anchors,
        'nearest_retained_physical_heavy_A': nearest[0], 'nearest_context_index': nearest[1],
        'nearest_context_element': nearest[2], 'retained_heavy_atom_count': len(heavy),
        'spectator_boundary_cutoff_A': 4.3, 'spectator_absent_from_source_to_QM_caps_bonds': absent,
        'full_source_preservation': preservation, 'spectator_local_energy': 'not_represented',
        'limitation': 'Distant spectator remains in source only; no global electrostatic or mechanical spectator response is calculated.',
        'context_preparation': row['representations']['context']['preparation'], 'endpoint': endpoints['Ca']['xyz']}
    write(path, receipt)
    assert passed, 'spectator at/below 4.3 A from retained physical heavy atoms or represented in local map'
    return pin(path)


def run_request(manifest, request, output, regression=False):
    assert os.environ.get('SLURM_JOB_ID')
    start = time.monotonic(); directory = Path(output); directory.mkdir(parents=True, exist_ok=False)
    m = read(manifest); source = request['source']; cid = source['case_id']
    result = {'case_id': cid, 'original_source_id': request['original_source_id'], 'selection_id': request['group']['selection_id'],
        'source': source, 'status': 'unsupported', 'mapping_eligible': False, 'paired_maps_equal': False,
        'state_anchor': request['state_anchor'], 'anchor_case_id': request['anchor_case_id'],
        'frozen_state_key': request['frozen_state_key'], 'spectator_admission': None}
    try:
        sys.path.insert(0,str(ROOT/'scripts'))
        import pqq_fast_prepare as fast
        import motion_envelope as envelope
        import gemmi
        import numpy as np
        from affordable_common import paired, xyz
        for p in m['implementation'].values(): verify(p)
        fixed, legacy, crystal, normalizer, cp = fast.helpers(m['config'])
        wrapper = load('spectator_protonation', verify(m['implementation']['PROTONATE.py']))
        narrow = types.SimpleNamespace(protonate_standard_only=lambda source, output, *, pins_path, ph:
            wrapper.protonate(source, output, pins_path=pins_path, ph=ph, legacy=legacy, spectator_count=0 if regression else 1))
        scope = dict(fast.prepare_case.__globals__)
        scope['helpers'] = lambda config: (fixed, narrow, crystal, normalizer, cp)
        prepare_case = types.FunctionType(fast.prepare_case.__code__, scope, 'prepare_spectator_source')
        imported_pins(m)
        used, spectator = alias(source, request['audit_row'], directory)
        prepared = prepare_case(used, m['config'], directory/'source_preparation')
        assert prepared['status'] == 'prepared'
        core = read(prepared['core']['parent']); raw = verify(source['source_structure'])
        preservation = [conserved(raw, verify(core['normalized_source_structure']),0.001,True),
                        conserved(raw, verify(core['source_structure']),0.001,True)]
        for endpoints in (prepared['core']['endpoints'],prepared['representations']['context']['endpoints']):
            paired(verify(endpoints['La']['xyz']),verify(endpoints['Ca']['xyz']),endpoints['La']['charge'],endpoints['Ca']['charge'])
        row = envelope.prepare_task(prepared,request['group'],m['config'],directory/'expanded')
        result.update(row, original_source_id=request['original_source_id'],selection_id=request['group']['selection_id'],
            state_anchor=request['state_anchor'],anchor_case_id=request['anchor_case_id'],frozen_state_key=request['frozen_state_key'],
            source_preparation=pin(directory/'source_preparation'/'preparation_result.json'),
            source_preservation=preservation,spectator_admission=None,protocol_id=PROTOCOL)
        if row['status'] != 'prepared_awaiting_group_check': raise ValueError(row.get('reason',row['status']))
        assert row['mapping_eligible'] and row['paired_maps_equal']
        anchor=request['anchor']
        assert row['state_key']==anchor['state_key'] and row['protein_key']==anchor['protein_key'], 'original fixed union state/protein anchor mismatch'
        if not regression:
            result['spectator_admission']=spectator_check(request,prepared,row,spectator,preservation,directory/'SPECTATOR_ADMISSION.json')
        else:
            old_source=request['old_source_preparation']; checks=[]
            for name in ('core','context'):
                neweps=prepared['core']['endpoints'] if name=='core' else prepared['representations']['context']['endpoints']
                oldeps=old_source['core']['endpoints'] if name=='core' else old_source['representations']['context']['endpoints']
                for z in ('Ca','La'):
                    check=fast.coordinate_comparison(verify(neweps[z]['xyz']),verify(oldeps[z]['xyz']))
                    assert check['max_displacement_A'] is not None and check['max_displacement_A']<=1e-12
                    assert (neweps[z]['charge'],neweps[z]['multiplicity'])==(oldeps[z]['charge'],oldeps[z]['multiplicity'])
                    checks.append({'representation':name,'metal':z,**check})
            comparator=load('spectator_original_regression',verify(m['helpers']['regression_comparator']))
            compared=copy.deepcopy(result);compared['status']='prepared'
            checks.append(comparator.compare_probe(compared,request['old_expanded_preparation']))
            result['regression']={'status':'passed','passed':True,'checks':checks,'old_preparation':request['old_expanded_preparation_pin']}
        result['status']='prepared'; result['preparation_profile']=PROTOCOL
        result['imported_sources']=imported_pins(m)
    except Exception as exc:
        result.update(status='unsupported',reason=str(exc),exception_type=type(exc).__name__)
        (directory/'FAILURE.txt').write_text(traceback.format_exc())
        if (directory/'SPECTATOR_ADMISSION.json').exists(): result['spectator_admission']=pin(directory/'SPECTATOR_ADMISSION.json')
    result.update(preparation_wall_seconds=time.monotonic()-start, scientific_energy_calls=0, new_optimizations=0,
                  admission_manifest=pin(manifest))
    # Place the final admission beside the actual maps as required by the existing native executor.
    admission=(Path(result['maps']['Ca']['path']).parent if result.get('maps') else directory)/'ADMISSION.json'
    write(admission,result)
    write(directory/'RESULT.json',{'admission':pin(admission),'status':result['status'],'case_id':cid})
    return result


def execute(path,workers):
    placement=allocation(workers);m=read(path);assert m['protocol_id']==PROTOCOL and m['root']==str(ROOT)
    assert pin(__file__)==m['implementation']['PREPARE.py']
    assert pin(sys.executable)==m['python']
    for p in m['implementation'].values():verify(p)
    for key in ('original_scope','source_preparation','original_selection','original_anchors','canonical_reference'):verify(m[key])
    for p in m['audits']:verify(p)
    output=Path(path).resolve().parent/'execution';output.mkdir(parents=True,exist_ok=False)
    write(output/'PLACEMENT.json',placement);start=time.monotonic();results={};regression=None
    if m['regression_reuse']:
        previous=read(m['regression_reuse']);assert previous['status']=='passed' and previous['regression_source']==REGRESSION
        previous_manifest=read(previous['manifest'])
        assert previous_manifest['implementation']==m['implementation'] and previous_manifest['config']==m['config']
        regression={**previous,'reused_from':m['regression_reuse']}
    else:
        rr=run_request(path,m['regression'],output/'regression',True)
        regression={'status':'passed' if rr['status']=='prepared' and rr.get('regression',{}).get('passed') else 'failed',
            'regression_source':REGRESSION,'manifest':pin(path),'result':rr}
    write(output/'REGRESSION.json',regression)
    if regression['status']=='passed':
        with ProcessPoolExecutor(max_workers=workers, mp_context=multiprocessing.get_context('spawn')) as pool:
            jobs={pool.submit(run_request,path,r,output/'sources'/r['source']['case_id']):r for r in m['requests']}
            for future in as_completed(jobs):
                r=jobs[future];cid=r['source']['case_id']
                try:result=future.result()
                except Exception as exc:
                    result={'case_id':cid,'original_source_id':r['original_source_id'],'selection_id':r['group']['selection_id'],
                        'source':r['source'],'status':'worker_failure','reason':str(exc),'mapping_eligible':False,'paired_maps_equal':False,
                        'frozen_state_key':r['frozen_state_key'],'state_anchor':r['state_anchor'],'spectator_admission':None}
                results[cid]=result;print(json.dumps({'case_id':cid,'status':result['status'],'reason':result.get('reason')}),flush=True)
    rows=[results.get(r['source']['case_id'],{'case_id':r['source']['case_id'],'original_source_id':r['original_source_id'],
        'selection_id':r['group']['selection_id'],'source':r['source'],'status':'not_run_regression_failed',
        'mapping_eligible':False,'paired_maps_equal':False,'spectator_admission':None,
        'frozen_state_key':r['frozen_state_key'],'state_anchor':r['state_anchor']}) for r in m['requests']]
    used=resource.getrusage(resource.RUSAGE_SELF);children=resource.getrusage(resource.RUSAGE_CHILDREN)
    receipt={'placement':placement,'wall_seconds':time.monotonic()-start,
        'CPU_seconds':used.ru_utime+used.ru_stime+children.ru_utime+children.ru_stime,
        'self_peak_RSS_KiB':used.ru_maxrss,'child_peak_RSS_KiB':children.ru_maxrss,'energy_calls':0,'new_optimizations':0}
    final={'protocol_id':PROTOCOL,'manifest':pin(path),'terminal':True,'rows':rows,
        'denominator':len(rows),'supported':sum(r['status']=='prepared' for r in rows),
        'regression':pin(output/'REGRESSION.json'),'canonical_reference':m['canonical_reference'],'config':m['config'],
        'receipt':receipt,'status':'complete' if regression['status']=='passed' else 'regression_failed',
        'primary125_unchanged':True,'new_source_occupancy_hypothesis':'two native Ca; local endpoint substitutes only PQQ-associated Ca',
        'production_changed':False,'classification_results':None}
    write(output/'FINAL_PREPARATION.json',final)
    print(json.dumps({'final_preparation':pin(output/'FINAL_PREPARATION.json'),'supported':final['supported'],'denominator':len(rows)}))
    if regression['status']!='passed':raise SystemExit(2)


def main():
    p=argparse.ArgumentParser();sub=p.add_subparsers(dest='operation',required=True)
    f=sub.add_parser('freeze');f.add_argument('--audit',action='append',required=True);f.add_argument('--output',required=True);f.add_argument('--regression-result')
    e=sub.add_parser('execute');e.add_argument('--manifest',required=True);e.add_argument('--workers',type=int,required=True)
    a=p.parse_args()
    if a.operation=='freeze':freeze(a.audit,a.output,a.regression_result)
    else:execute(a.manifest,a.workers)


if __name__=='__main__':main()
