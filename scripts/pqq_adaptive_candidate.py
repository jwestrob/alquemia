"""Opt-in PQQ source -> minimal adaptive pool; released defaults stay unchanged."""
from __future__ import annotations
import argparse
import copy
import fcntl
import json
import os
from pathlib import Path
import re
import shutil
import time
import numpy as np
import scipy
import scipy.optimize._slsqp_py as scipy_slsqp
import scipy.optimize._slsqplib as scipy_kernel
from affordable_common import InvalidArtifact, cache_key, paired, read_json, record, verify, write_new, xyz
import pqq_standard as standard
import pqq_fast_prepare as source_prepare
import adaptive_origin_recovery as recovery
import adaptive_completion as completion
import adaptive_angular_proposals as angular
from adaptive_force_diagnostic import preview, project
from accommodation_nonlinear import WarmGPU
from accommodation_reference_geometry import mapping_case
from accommodation_torsion_profiles import geometry_check
from mace_site_kinematics import Kinematics
from mace_hybrid import check_atoms, write_xyz
from compact_solvation_scanner import decision, pqq_decision_compatible
from nikasha_pool_compare import call as adaptive_decision
import nikasha_pool as pool

PROTOCOL = 'Nikasha_opt_in_source_minimal_adaptive_PQQ_v1'
REPLAY = 'Nikasha_completed_origin_recovery_replay_v1'
CANDIDATES = list(recovery.CANDIDATES)
ROOT = Path(__file__).resolve().parents[1]
DEFAULT_REFERENCE = ROOT / 'workspaces/adaptive_minimal_pool_20260923/REFERENCE_v1.json'
DEFAULT_AGREEMENT = ROOT / 'diagnostics/pqq_adaptive_candidate_20260923/PLAN.md'


def put(path, data):
    """Immutable restart: only exactly identical records can be reused."""
    p = Path(path)
    if p.exists():
        if read_json(p) != data: raise InvalidArtifact('partial-stage record differs: ' + str(p))
    else: write_new(p, data)
    return record(p)


def identities(req):
    ids = req.get('case_ids') if req.get('schema_version') == REPLAY else [c['case_id'] for c in req['cases']]
    if not ids or len(set(ids)) != len(ids) or any(not re.fullmatch(r'[A-Za-z0-9_.-]+', c) for c in ids):
        raise InvalidArtifact('explicit nonempty safe unique case IDs required')
    return ids


def check_request(req, release, reference):
    ids = identities(req)
    method = {k: reference['signature'][k] for k in ('model', 'software', 'orca')}
    if req.get('schema_version') == REPLAY:
        result = read_json(verify(req['collection'])); mp = verify(result['manifest'])
        recovery.validate_pool(mp); m = read_json(mp)
        if set(ids) - {c['case_id'] for c in result['cases']}: raise InvalidArtifact('unknown archived source')
        if any(m[k] != method[k] for k in method): raise InvalidArtifact('archive method differs')
        if m['reference']['sha256'] != reference['_pin']['sha256']: raise InvalidArtifact('archive reference differs')
        for c in result['cases']:
            if c['case_id'] not in ids: continue
            if pool.choose_rows(c['matrix'], CANDIDATES) != c['pool']: raise InvalidArtifact('archived pool algebra changed')
        return 'exact_archive_replay'
    if req.get('schema_version') != source_prepare.SOURCE_PROTOCOL or req['config'] != release['source_configuration']:
        raise InvalidArtifact('unsupported source protocol/configuration')
    cp = read_json(verify(req['config']['calibration_implementation_pins']))
    fresh_method = {k: req['config'][k] for k in ('model', 'software')} | {'orca': cp['orca_runtime']['executable']}
    if fresh_method != method: raise InvalidArtifact('source electronic method incompatible with reference')
    for c in req['cases']:
        verify(c['source_structure'])
        if c['normalization'] not in ('protenix_generic_PQQ', 'selected_crystal_chain'):
            raise InvalidArtifact('unsupported selector/source normalization policy')
        if c['normalization'] == 'selected_crystal_chain':
            if not {'model','chain','raw_metal','raw_pqq'} <= set(c['assembly']): raise InvalidArtifact('explicit crystal selection required')
        if not {'metal','pqq','roles','assembly'} <= set(c): raise InvalidArtifact('explicit canonical selectors required')
    return 'fresh_source_request'


def reference_data(path):
    ref = read_json(path)
    recovery.check_reference(path, ref['signature'])
    return {**ref, '_pin': record(path)}


def prepare(request, output, reference=DEFAULT_REFERENCE, release=standard.DEFAULT_RELEASE, agreement=DEFAULT_AGREEMENT):
    r = standard.release(release); ref = reference_data(reference); req = read_json(request)
    mode = check_request(req, r, ref); ids = identities(req)
    out = Path(output).resolve()
    if 'workspaces' not in out.parts: raise InvalidArtifact('candidate products must be under workspaces/')
    out.mkdir(parents=True, exist_ok=False)
    impl = recovery.snapshot(Path(__file__).parent, out / 'implementation')
    domains = {c: 'consumed_recovery_geometry' for c in ids} if mode == 'exact_archive_replay' else {
        c['case_id']: 'consumed_reference_geometry' if c['source_structure']['sha256'] in r['reference_source_sha256']
        else 'unvalidated_input_domain' for c in req['cases']}
    p = {'protocol_id': PROTOCOL, 'request': record(request), 'mode': mode, 'case_ids': ids,
         'reference': record(reference), 'release': record(release), 'agreement': record(agreement),
         'implementation': impl, 'settings': completion.SETTINGS, 'candidate_ids': CANDIDATES,
         'source_domain_status': domains, 'production_changed': False, 'automatic_fallback': False,
         'resources': {'CPUs': 32, 'host_memory_MiB': 200000, 'GPUs': 1, 'GFN2_workers': 4, 'GFN2_MPI_ranks': 8},
         'optimizer_software': {'version': scipy.__version__, 'wrapper': record(scipy_slsqp.__file__), 'kernel': record(scipy_kernel.__file__)},
         'declared_maximum_calls': {'origin_MACE': 0 if mode == 'exact_archive_replay' else 2*len(ids),
             'bounded_MACE_searches': 0 if mode == 'exact_archive_replay' else 2*len(ids),
             'cross_MACE': 0 if mode == 'exact_archive_replay' else 2*len(ids),
             'GFN2': 0 if mode == 'exact_archive_replay' else 12*len(ids), 'DFT': 0},
         'fresh_molecular_integration_status_at_prepare': 'unexecuted'}
    put(out/'plan.json', p)
    return dry_run(out/'plan.json')


def checked(plan):
    p = read_json(plan)
    if p['protocol_id'] != PROTOCOL or p['settings'] != completion.SETTINGS or p['candidate_ids'] != CANDIDATES:
        raise InvalidArtifact('candidate method changed')
    if p['production_changed'] or p['automatic_fallback']: raise InvalidArtifact('default/fallback unsupported')
    for pin in p['implementation'].values(): verify(pin)
    for name in ('wrapper', 'kernel'): verify(p['optimizer_software'][name])
    if p['optimizer_software']['version'] != scipy.__version__: raise InvalidArtifact('optimizer version differs')
    verify(p['agreement']); r = standard.release(verify(p['release'])); ref = reference_data(verify(p['reference']))
    req = read_json(verify(p['request']))
    if check_request(req, r, ref) != p['mode'] or identities(req) != p['case_ids']: raise InvalidArtifact('request scope changed')
    if p['resources'] != {'CPUs':32,'host_memory_MiB':200000,'GPUs':1,'GFN2_workers':4,'GFN2_MPI_ranks':8}:
        raise InvalidArtifact('resource layout changed')
    return p, req, r, ref


def dry_run(plan):
    p, _, _, _ = checked(plan)
    return {'status': 'ready', 'plan': record(plan), 'mode': p['mode'], 'case_ids': p['case_ids'],
            'declared_maximum_calls': p['declared_maximum_calls'], 'new_molecular_calls': 0,
            'source_domain_status': p['source_domain_status'],
            'fresh_preparation_and_selector_status': 'not_executed' if p['mode']=='fresh_source_request' else 'archive_verified'}


def common_manifest(plan):
    p, _, r, ref = checked(plan)
    qualification = read_json(verify(ref['numerical_provenance'][0]['qualification']))
    if not qualification['both_controls_pass'] or not qualification['formerly_failed_cell_complete']:
        raise InvalidArtifact('native MaxIter500 qualification missing')
    return {**{k:ref['signature'][k] for k in ('model','software','orca')},
            'protocol_id':PROTOCOL,'settings':{**pool.SETTINGS,'candidate_order':CANDIDATES},
            'proposal_settings':completion.SETTINGS,'reference':p['reference'],'agreement':p['agreement'],
            'plan':record(plan),'cpu_python':r['cpu_python'],'gpu_python':r['gpu_python'],
            'implementation':p['implementation'],'GFN2_maxiter':500,
            'numerical_policy_id':'native_GFN2_MaxIter500_unchanged_convergence_v1',
            'numerical_qualification':ref['numerical_provenance'][0]['qualification'],
            'optimizer_software':p['optimizer_software']}


def low_prepare(manifest):
    """Same recipe/executor as the pool, four 8-rank tasks within 32 CPUs."""
    from affordable_workflow import dry_run as existing_dry_run
    m = read_json(manifest); root = Path(manifest).parent/'solvent'; sd = root/'shard_0'; tasks=[]
    for cell in m['tasks']:
        for medium in ('vacuum','alpb'):
            td = sd/'tasks'/(cell['task_id']+'__'+medium); td.mkdir(parents=True,exist_ok=True)
            xp = td/'core.xyz'; ip = td/'endpoint.inp'; body=pool.input_text(cell['charge'],cell['multiplicity'],medium,'native').replace('%scf\n','%scf\n MaxIter 500\n')
            if xp.exists() and xyz(xp)!=xyz(verify(cell['xyz'])): raise InvalidArtifact('partial coordinate changed')
            if not xp.exists(): shutil.copyfile(verify(cell['xyz']),xp)
            if ip.exists() and ip.read_text()!=body: raise InvalidArtifact('partial input changed')
            if not ip.exists(): ip.write_text(body)
            tasks.append({**cell,'task_id':td.name,'cell_id':cell['task_id'],'case':cell['case_id'],
                          'medium':medium,'xyz':record(xp),'input':record(ip),'output_path':str(td/'endpoint.out')})
    lm={'protocol_id':PROTOCOL,'pool_manifest':record(manifest),'agreement':m['agreement'],'orca':m['orca'],
        'tasks':tasks,'all_tasks':tasks,'execution_resources':{'mpi_ranks':8,'concurrent_tasks':4},
        'execution_policy':{'task_runner':m['implementation']['run_orca_task_manifest.py'],
                            'runtime_renderer':m['implementation']['render_orca_runtime_input.py']}}
    lp=sd/'manifest.json';put(lp,lm);put(root/'INDEX.json',{'pool_manifest':record(manifest),'shards':[record(lp)]})
    if tasks: existing_dry_run(lp)
    return lp


def origin_manifest(plan, preparation, output):
    """Real prepared coordinates/maps only; no force/energy call in this builder."""
    p,_,_,_=checked(plan); prep=read_json(preparation); m=common_manifest(plan)
    if prep['source_manifest'] != p['request']: raise InvalidArtifact('prepared source request differs')
    if [c['case_id'] for c in prep['cases']] != p['case_ids']: raise InvalidArtifact('preparation membership changed')
    out=Path(output);out.mkdir(parents=True,exist_ok=False); cases=[]; tasks=[]
    for src in prep['cases']:
        cid=src['case_id'];c={'case_id':cid,'status':'unavailable','reason':src.get('reason'),
            'candidates':[{'id':'origin'}],'matrix':{'Ca':{},'La':{}},'source':src}
        pending=[]
        try:
            if src['status']!='prepared' or not pqq_decision_compatible(src): raise InvalidArtifact('source chemistry/preparation unsupported')
            rep=src['representations']['context'];ca,la=[rep['endpoints'][z] for z in ('Ca','La')]
            paired(verify(la['xyz']),verify(ca['xyz']),la['charge'],ca['charge'])
            mp=mapping_case(src,prep['config']['topology'],out/'maps')
            if mp['status']!='supported': raise InvalidArtifact(mp['reason'])
            kin=Kinematics(read_json(verify(mp['mapping']))['context']); q=np.zeros(len(kin.modes))
            for z in ('Ca','La'):
                ep=rep['endpoints'][z]; atoms=xyz(verify(ep['xyz'])); check_atoms(atoms,ep['charge'])
                if ep['multiplicity']!=1 or not np.allclose(kin.evaluate(q)[1],[a[1:] for a in atoms],atol=1e-12,rtol=0):
                    raise InvalidArtifact('paired origin state or exact map differs')
                if not geometry_check(kin,q,[a[0] for a in atoms])['pass']: raise InvalidArtifact('origin geometry unsupported')
                tid=cid+'__'+z
                t={'task_id':tid,'case_id':cid,'metal':z,'candidate':'origin','xyz':ep['xyz'],
                   'charge':ep['charge'],'multiplicity':1,'mapping':mp['mapping'],'source_preparation':rep['preparation'],
                   'mode_count':len(kin.modes)}
                pending.append(t);c['matrix'][z]['origin']={'status':'pending','task_id':tid,'xyz':ep['xyz'],'reused':False}
            c.update(status='prepared',reason=None);tasks.extend(pending)
        except (InvalidArtifact,KeyError,ValueError,OSError) as exc: c['reason']=str(exc)
        cases.append(c)
    m.update(cases=cases,tasks=tasks,new_MACE_cells=len(tasks),preparation=record(preparation))
    path=out/'manifest.json';put(path,m);low_prepare(path);return path


def validate_stage(manifest, plan):
    m=read_json(manifest); base=common_manifest(plan)
    for key in ('protocol_id','model','software','orca','reference','agreement','plan','implementation','GFN2_maxiter','numerical_qualification','optimizer_software'):
        if m[key]!=base[key]:raise InvalidArtifact('stage scientific identity changed: '+key)
    p=read_json(plan)
    if [c['case_id'] for c in m['cases']]!=p['case_ids']:raise InvalidArtifact('stage population changed')
    if len({t['task_id'] for t in m['tasks']})!=len(m['tasks']):raise InvalidArtifact('duplicate stage task')
    if 'origin_collection' in m:
        original=read_json(verify(m['origin_collection'])); om=read_json(verify(original['manifest']))
        if om['plan']!=record(plan) or m['settings']!=completion.SETTINGS:raise InvalidArtifact('proposal source/settings changed')
        expected=[]
        for c in m['cases']:
            if c['status']!='prepared':continue
            src=next(x for x in original['cases'] if x['case_id']==c['case_id'])
            pair=[next(t for t in om['tasks'] if t['case_id']==c['case_id'] and t['metal']==z) for z in ('Ca','La')]
            expected.extend(selected_tasks(pair,{z:src['matrix'][z]['origin'] for z in ('Ca','La')},m))
        if expected!=m['tasks']:raise InvalidArtifact('actual-force selector/source changed')
    elif 'preparation' in m:
        prep=read_json(verify(m['preparation']))
        if prep['source_manifest']!=p['request']:raise InvalidArtifact('stage preparation source changed')
        expected=[]
        for c in m['cases']:
            src=next(x for x in prep['cases'] if x['case_id']==c['case_id'])
            if c['source']!=src:raise InvalidArtifact('prepared source row changed')
            if c['status']!='prepared':continue
            pair=[t for t in m['tasks'] if t['case_id']==c['case_id']]
            if [t['metal'] for t in pair]!=['Ca','La']:raise InvalidArtifact('unpaired origins')
            for t in pair:
                ep=src['representations']['context']['endpoints'][t['metal']]
                if (t['xyz'],t['charge'],t['multiplicity'])!=(ep['xyz'],ep['charge'],ep['multiplicity']):raise InvalidArtifact('origin state changed')
                kin=Kinematics(read_json(verify(t['mapping']))['context'])
                if not np.allclose(kin.evaluate(np.zeros(len(kin.modes)))[1],[a[1:] for a in xyz(verify(t['xyz']))],atol=1e-12,rtol=0):raise InvalidArtifact('origin map changed')
                expected.append(t['task_id'])
        if expected!=[t['task_id'] for t in m['tasks']]:raise InvalidArtifact('undeclared origin task')
    elif 'proposal_manifest' in m:
        parent=read_json(verify(m['proposal_manifest']))
        if parent['plan']!=record(plan):raise InvalidArtifact('candidate source differs')
        expected=[]
        for c in m['cases']:
            if c['status']!='prepared':continue
            prior=next(x for x in parent['cases'] if x['case_id']==c['case_id'])
            for z in ('Ca','La'):
                expected_origin={**prior['matrix'][z]['origin'],'reused':True}
                if c['matrix'][z]['origin']!=expected_origin:raise InvalidArtifact('origin component changed')
            if [q['id'] for q in c['candidates']]!=CANDIDATES:raise InvalidArtifact('candidate set changed')
            for candidate in c['candidates'][1:]:
                r=read_json(verify(candidate['proposal_receipt']))
                if r['manifest']!=m['proposal_manifest'] or r['status']!='proposal_available' or r['proposal']['coordinate']!=candidate['xyz']:raise InvalidArtifact('candidate receipt mismatch')
                for z in ('Ca','La'):
                    tid=c['case_id']+'__'+z+'__at_'+candidate['id'];t=next(t for t in m['tasks'] if t['task_id']==tid)
                    src=next(t for t in parent['tasks'] if t['case_id']==c['case_id'] and t['metal']==z)
                    atoms=xyz(verify(candidate['xyz']))
                    if xyz(verify(t['xyz']))!=[(z,*atoms[0][1:]),*atoms[1:]] or (t['charge'],t['multiplicity'])!=(src['charge'],src['multiplicity']):raise InvalidArtifact('cross endpoint state changed')
                    expected.append(tid)
        if expected!=[t['task_id'] for t in m['tasks']]:raise InvalidArtifact('undeclared candidate task')
    else:raise InvalidArtifact('unknown finite stage')
    return m


def native_calls(manifest):
    """One existing warm worker per stage; failed receipts are not silently retried."""
    m=read_json(manifest);root=Path(manifest).parent;gpu=None
    try:
        for t in m['tasks']:
            if t.get('native_reuse'): pool.native_reuse(t,m);continue
            d=root/'cells'/t['task_id'];d.mkdir(parents=True,exist_ok=True)
            request={k:t[k] for k in ('task_id','xyz','charge','multiplicity')}|{'manifest':record(manifest)}
            rp=d/'request.json';put(rp,request);dest=d/'mace_result.json'
            if dest.exists():
                result=read_json(dest)
                if result['request']!=record(rp) or result['model']!=m['model']: raise InvalidArtifact('native cache mismatch')
                continue
            if gpu is None: gpu=WarmGPU(manifest)
            gpu.evaluate(record(rp))
    finally:
        if gpu is not None: gpu.close()


def selected_tasks(pair, cells, method):
    """Project actual native Cartesian forces; select one common physical basis."""
    data=[]
    for t in pair:
        cell=cells[t['metal']]
        if cell['status']!='complete': raise InvalidArtifact('required origin energy/force unavailable')
        native=read_json(verify(cell['MACE']));req=read_json(verify(native['request']))
        if native['status']!='complete' or native['model']!=method['model'] or (req['charge'],req['multiplicity'])!=(t['charge'],t['multiplicity']) or xyz(verify(req['xyz']))!=xyz(verify(t['xyz'])):
            raise InvalidArtifact('origin native model/state/coordinate differs')
        mapping=read_json(verify(t['mapping']))['context'];q=np.zeros(len(mapping['modes']))
        force=np.load(verify(native['forces']),allow_pickle=False);coord,vectors,raw,norm,_=project(mapping,q,force)
        if not np.allclose(coord,[a[1:] for a in xyz(verify(t['xyz']))],atol=1e-12,rtol=0): raise InvalidArtifact('origin force map differs')
        data.append((mapping,vectors,raw,norm,native))
    if data[0][0]!=data[1][0]: raise InvalidArtifact('paired physical maps differ')
    paired(verify(pair[1]['xyz']),verify(pair[0]['xyz']),pair[1]['charge'],pair[0]['charge'])
    choice=preview(data[0][0]['modes'],data[0][1],data[0][3],data[1][3])
    if len(choice['selected'])!=4: raise InvalidArtifact('fewer than four independent supported angular modes')
    ids=[v['id'] for v in choice['selected']]; names=[v['id'] for v in data[0][0]['modes']];indices=[names.index(i) for i in ids]
    result=[]
    for t,(_,_,raw,_,native) in zip(pair,data):
        cell=cells[t['metal']]
        point={'status':'complete','coordinate':t['xyz'],'MACE':cell['MACE'],'MACE_eV':native['energy_eV'],
               'forces':native['forces'],'full_q':[0.]*len(names),'active_q_radian':[0.]*4,
               'gradient_kcal_mol_rad':raw[indices].tolist()}
        task={**t,'active_indices':indices,'active_mode_ids':ids,'active_roles':['adaptive_physical_angular']*4,
              'selector':choice,'q0_status':'available','q0':{'components':cell['components']},
              'origin_reuse':{'point':point,'proposal_receipt':cell['MACE']}}
        angular.final_geometry(Kinematics(data[0][0]),task,np.zeros(4),[a[0] for a in xyz(verify(t['xyz']))])
        result.append(task)
    return result


def proposal_manifest(plan, origins, output):
    data=read_json(origins);om=read_json(verify(data['manifest']));m=common_manifest(plan);out=Path(output)
    out.mkdir(parents=True,exist_ok=False);cases=[];tasks=[]
    for c in data['cases']:
        row={**c,'origin_pool':c['pool']}
        try:
            if c['pool']['status']!='available': raise InvalidArtifact(c.get('reason') or 'required origin cell unavailable')
            pair=[next(t for t in om['tasks'] if t['case_id']==c['case_id'] and t['metal']==z) for z in ('Ca','La')]
            selected=selected_tasks(pair,{z:c['matrix'][z]['origin'] for z in ('Ca','La')},m)
            tasks.extend(selected);row.update(status='prepared',reason=None)
        except (InvalidArtifact,KeyError,ValueError,OSError) as exc:row.update(status='unavailable',reason=str(exc))
        cases.append(row)
    m.update(settings=completion.SETTINGS,cases=cases,tasks=tasks,origin_collection=record(origins))
    path=out/'manifest.json';put(path,m);return path


def candidates_manifest(plan, proposals, output):
    parent=read_json(proposals);m=common_manifest(plan);out=Path(output);out.mkdir(parents=True,exist_ok=False)
    cases=[];tasks=[]
    for src in parent['cases']:
        c=copy.deepcopy(src);c['candidates']=[{'id':'origin'}]
        for z in ('Ca','La'):
            if 'origin' in c['matrix'][z]: c['matrix'][z]['origin']['reused']=True
        pending=[]
        try:
            if src['status']!='prepared': raise InvalidArtifact(src['reason'])
            ts={z:next(t for t in parent['tasks'] if t['case_id']==c['case_id'] and t['metal']==z) for z in ('Ca','La')}
            for maker in ('Ca','La'):
                rp=Path(proposals).parent/'proposals'/ts[maker]['task_id']/'result.json';r=read_json(rp)
                if r['manifest']!=record(proposals) or r['status']!='proposal_available': raise InvalidArtifact('adaptive '+maker+' candidate unavailable: '+str(r.get('reason')))
                point=r['proposal'];q=np.array(point['full_q']);kin=Kinematics(read_json(verify(ts[maker]['mapping']))['context'])
                if any(q[i] for i in set(range(len(q)))-set(ts[maker]['active_indices'])): raise InvalidArtifact('unselected physical mode moved')
                atoms=xyz(verify(point['coordinate']));angular.final_geometry(kin,ts[maker],q[ts[maker]['active_indices']],[a[0] for a in atoms])
                if not np.allclose(kin.evaluate(q)[1],[a[1:] for a in atoms],atol=1e-12,rtol=0): raise InvalidArtifact('proposal coordinate mismatch')
                name='adaptive_'+maker;c['candidates'].append({'id':name,'xyz':point['coordinate'],'proposal_receipt':record(rp),
                    'boundary_flag':r.get('boundary_flag'),'final_geometry':r.get('final_geometry'),'optimizer':r['optimizer']})
                for z in ('Ca','La'):
                    tid=c['case_id']+'__'+z+'__at_'+name;d=out/'cells'/tid;d.mkdir(parents=True);xp=d/'context.xyz'
                    values=[(z,*atoms[0][1:]),*atoms[1:]];write_xyz(xp,values);check_atoms(values,ts[z]['charge'])
                    if xyz(xp)!=values: raise InvalidArtifact('cross coordinate serialization changed')
                    t={'task_id':tid,'case_id':c['case_id'],'metal':z,'candidate':name,'xyz':record(xp),
                       'charge':ts[z]['charge'],'multiplicity':1,'source_mapping':ts[z]['mapping'],'source_preparation':ts[z]['source_preparation']}
                    if z==maker:t['native_reuse']=point['MACE'];pool.native_reuse(t,m)
                    pending.append(t);c['matrix'][z][name]={'status':'pending','task_id':tid,'xyz':record(xp),'reused':False}
            c.update(status='prepared',reason=None);tasks.extend(pending)
        except (InvalidArtifact,KeyError,ValueError,OSError) as exc:c.update(status='unavailable',reason=str(exc))
        cases.append(c)
    m.update(cases=cases,tasks=tasks,proposal_manifest=record(proposals),new_MACE_cells=sum(not t.get('native_reuse') for t in tasks))
    path=out/'manifest.json';put(path,m);low_prepare(path);return path


def execute(plan):
    p,req,_,_=checked(plan);root=Path(plan).parent
    if p['mode']=='exact_archive_replay': return {'status':'archive_replay_requires_no_execution','new_molecular_calls':0}
    recovery.allocation()
    if Path(__file__).resolve()!=verify(p['implementation']['pqq_adaptive_candidate.py']):raise InvalidArtifact('execute the pinned implementation')
    start=time.monotonic();error=None
    with (root/'execute.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        try:
            prep=root/'source_preparation/preparation.json'
            if not prep.exists(): source_prepare.prepare(verify(p['request']),prep.parent)
            op=root/'origins/manifest.json'
            if not op.exists(): origin_manifest(plan,prep,op.parent)
            validate_stage(op,plan);native_calls(op)
            from affordable_workflow import execute as low_execute
            for stage in ('origins',):
                lp=root/stage/'solvent/shard_0/manifest.json'
                if read_json(lp)['tasks']:
                    try: low_execute(lp)
                    except Exception as exc: put(root/stage/'solver_failure.json',{'reason':str(exc)})
            oc=root/'origins/collection.json'
            if not oc.exists(): pool.collect(op,oc)
            pp=root/'proposals/manifest.json'
            if not pp.exists(): proposal_manifest(plan,oc,pp.parent)
            pm=validate_stage(pp,plan);gpu=None
            try:
                if pm['tasks']:gpu=WarmGPU(pp)
                for t in pm['tasks']: completion.optimize(pp,t,gpu)
            finally:
                if gpu is not None:gpu.close()
            cp=root/'pool/manifest.json'
            if not cp.exists(): candidates_manifest(plan,pp,cp.parent)
            validate_stage(cp,plan);native_calls(cp);lp=cp.parent/'solvent/shard_0/manifest.json'
            if read_json(lp)['tasks']:
                try:low_execute(lp)
                except Exception as exc:put(cp.parent/'solver_failure.json',{'reason':str(exc)})
            final=cp.parent/'collection.json'
            if not final.exists():pool.collect(cp,final)
        except Exception as exc:error=str(exc)
        finally:
            elapsed=time.monotonic()-start
            put(root/('execution_'+os.environ['SLURM_JOB_ID']+'.json'),{'plan':record(plan),'error':error,
                'job_id':os.environ['SLURM_JOB_ID'],'wall_seconds':elapsed,'allocated_core_seconds':32*elapsed,
                'allocated_GPU_seconds':elapsed,'source_preparation_included':True})
    if error:raise InvalidArtifact(error)
    return {'status':'execution_finished','collection':record(root/'pool/collection.json')}


def collect(plan, output):
    p,req,r,ref=checked(plan);root=Path(plan).parent
    src=verify(req['collection']) if p['mode']=='exact_archive_replay' else root/'pool/collection.json'
    data=read_json(src) if src.exists() else None
    if data and p['mode']!='exact_archive_replay' and read_json(verify(data['manifest']))['plan']!=record(plan):raise InvalidArtifact('collection belongs to another plan')
    fallback_origin=root/'origins/collection.json'
    origin_data=read_json(fallback_origin) if not data and fallback_origin.exists() else None
    if origin_data and read_json(verify(origin_data['manifest']))['plan']!=record(plan):raise InvalidArtifact('origin collection belongs to another plan')
    static=read_json(verify(r['artifacts']['calibration']))['calibration']['context']['bands'];rows=[]
    for cid in p['case_ids']:
        c=next((c for c in data['cases'] if c['case_id']==cid),None) if data else None
        selected=c['pool'] if c else {'status':'unavailable','reason':'execution_or_collection_incomplete'}
        original_case=c or (next((x for x in origin_data['cases'] if x['case_id']==cid),None) if origin_data else None)
        origin=None
        if original_case and all(original_case['matrix'].get(z,{}).get('origin',{}).get('status')=='complete' for z in ('Ca','La')):
            origin=pool.score(original_case['matrix']['Ca']['origin']['components'],original_case['matrix']['La']['origin']['components'])['composite_R_model_kcal_mol']
        if c and selected['status']=='available':
            if [v['id'] for v in c['candidates']]!=CANDIDATES or pool.choose_rows(c['matrix'],CANDIDATES)!=selected:raise InvalidArtifact('incomplete/changed three-candidate matrix')
        variants={}
        for mode in ('mathematical','operational'):
            value=selected[mode]['composite_R_model_kcal_mol'] if selected['status']=='available' else None
            variants[mode]={'R_model_kcal_mol':value,'decision':adaptive_decision(value,ref['variants'][mode]['bands']),
                            'delta_R_from_origin':None if value is None or origin is None else value-origin}
        rows.append({'case_id':cid,'status':selected['status'],'reason':c.get('reason') or selected.get('reason') if c else selected['reason'],
            'source_domain_status':p['source_domain_status'][cid],'original_R_model_kcal_mol':origin,
            'original_released_decision':decision(origin,static),'candidate':variants,'pool':selected,
            'candidate_geometry':[q for q in c['candidates'] if q['id']!='origin'] if c else [],
            'reference_compatibility':'compatible_protocol_only_not_independent_biological_validation'})
    result={'protocol_id':PROTOCOL,'plan':record(plan),'mode':p['mode'],'source_collection':record(src) if data else None,
        'reference':p['reference'],'original_reference':r['artifacts']['calibration'],'rows':rows,
        'denominator':len(rows),'available':sum(c['status']=='available' for c in rows),'production_changed':False,
        'aquo_referenced_score':None,'affinity_probability':None,'new_reference_fitted':False,
        'execution_receipts':[record(q) for q in root.glob('execution_*.json')],
        'interpretation':'finite electronic geometry candidate; no equilibrium population, affinity or minimum claim'}
    put(output,result);return {k:v for k,v in result.items() if k!='rows'}


def report(result, output):
    d=read_json(result);checked(verify(d['plan']))
    lines=['# Opt-in minimal adaptive PQQ candidate','',f"Mode: {d['mode']}. Available {d['available']}/{d['denominator']}; defaults unchanged.",'',
           '| Source | Original R / call | Candidate R / call | Delta R | Status |','|---|---|---|---|---|']
    for r in d['rows']:
        a=r['candidate']['operational'];lines.append(f"| {r['case_id']} | {r['original_R_model_kcal_mol']} / {r['original_released_decision']} | {a['R_model_kcal_mol']} / {a['decision']} | {a['delta_R_from_origin']} | {r['status']}: {r['reason'] or ''} |")
    lines+=['','R is model kcal/mol, Ca minus La. Each method uses its own frozen bands.',
        'This is an opt-in development candidate. Protocol compatibility does not establish',
        'new biological validation. Missing candidate scores never inherit an original score.',
        'Complete matrix, works, selections, source-domain labels and geometry flags are in '+str(Path(result).resolve())+'.']
    with Path(output).open('x') as f:f.write('\n'.join(lines)+'\n')
    return {'report':record(output),'result':record(result)}


def main():
    parser=argparse.ArgumentParser(description=__doc__);sub=parser.add_subparsers(dest='op',required=True)
    q=sub.add_parser('prepare');q.add_argument('--request',required=True,type=Path);q.add_argument('--output',required=True,type=Path)
    q.add_argument('--reference',type=Path,default=DEFAULT_REFERENCE);q.add_argument('--release',type=Path,default=standard.DEFAULT_RELEASE);q.add_argument('--agreement',type=Path,default=DEFAULT_AGREEMENT)
    for name in ('dry-run','execute','collect'):
        q=sub.add_parser(name);q.add_argument('--plan',required=True,type=Path)
        if name=='collect':q.add_argument('--output',required=True,type=Path)
    q=sub.add_parser('report');q.add_argument('--result',required=True,type=Path);q.add_argument('--output',required=True,type=Path)
    args=vars(parser.parse_args());op=args.pop('op').replace('-','_');print(json.dumps(globals()[op](**args),indent=2))

if __name__=='__main__':main()
