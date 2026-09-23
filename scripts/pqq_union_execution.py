"""Fresh matched execution of released static and the opt-in union precision profile.

Source eligibility stays in pqq_union_candidate. Existing molecular kernels and
executors remain unchanged; this adapter names their approved composition.
"""
from __future__ import annotations
import argparse
import copy
import fcntl
import importlib.util
import itertools
import json
import os
from pathlib import Path
import re
import shutil
import sys
import time
from types import SimpleNamespace
import numpy as np
import scipy
import scipy.optimize._slsqp_py as slsqp_wrapper
import scipy.optimize._slsqplib as slsqp_kernel

from affordable_common import InvalidArtifact, cache_key, read_json, record, verify, write_new, xyz
from adaptive_origin_recovery import snapshot
import pqq_union_candidate as source_interface
import pqq_adaptive_candidate as existing_candidate
import pqq_standard as standard
import pqq_fast_release
import compact_solvation_scanner as static_scanner
import slsqp_precision as precision
import nikasha_pool as shared_pool
import union_adaptive
from accommodation_folds_compare import decision as union_decision

PROTOCOL='Nikasha_fresh_union_ftol1e8_rank1_execution_v1'
PROFILE='native_OMOL_union_SLSQP200_ftol1e8_GFN2_MaxIter500_rank1_v1'
STATIC_PROFILE='released_local_context_native_OMOL_GFN2_rank8_v1'
ARMS=('released-static','union-candidate')
GROUP='a0a3f2yly8-pqq-la_model'
CANDIDATES=['origin','adaptive_Ca','adaptive_La']
RESOURCES={'CPUs':32,'GPUs':1,'host_memory_MiB':200000}
REFERENCE_SHA='1f8470bbdf7a056098ca261aa10b72bdc940cf71ab860fbbc54090bda1e73250'
ROOT=Path(__file__).resolve().parents[1]


def declared_calls(arm):
    return {'native_MACE_scalar':20 if arm=='released-static' else 0,
            'native_MACE_q0_gradient':20 if arm=='union-candidate' else 0,
            'bounded_native_MACE_searches':20 if arm=='union-candidate' else 0,
            'maximum_cross_MACE':20 if arm=='union-candidate' else 0,
            'maximum_GFN2':40 if arm=='released-static' else 120,'DFT':0}


def private(name):
    path=Path(__file__).with_name(name+'.py');key='_fresh_union_execution_'+name
    spec=importlib.util.spec_from_file_location(key,path);module=importlib.util.module_from_spec(spec)
    sys.modules[key]=module;spec.loader.exec_module(module);return module


def check_request(req):
    identity=source_interface.check_request(req)
    if len(req['groups'])!=1 or req['groups'][0]['protein_id']!=GROUP or len(req['cases'])!=10:
        raise InvalidArtifact('this finite integration requires all ten A0A3 reference sources')
    if req['groups'][0]['kind']!='ten_folds':raise InvalidArtifact('ten-fold union required')
    return identity


def prepare(request,arm,reference,rank_qualification,maxiter_qualification,agreement,output):
    if arm not in ARMS:raise InvalidArtifact('unknown matched arm')
    req=read_json(request);identities=check_request(req)
    out=Path(output).resolve()
    if 'workspaces' not in out.parts:raise InvalidArtifact('execution products belong in workspaces')
    out.mkdir(parents=True,exist_ok=False)
    impl=snapshot(Path(__file__).parent,out/'implementation')
    p={'protocol_id':PROTOCOL,'profile':PROFILE if arm=='union-candidate' else STATIC_PROFILE,
       'arm':arm,'request':record(request),'release':req['release'],'reference':record(reference),
       'rank_qualification':record(rank_qualification),'maxiter_qualification':record(maxiter_qualification),
       'agreement':record(agreement),'implementation':impl,'case_ids':[c['case_id'] for c in req['cases']],
       'source_identities':identities,'settings':precision.SETTINGS,'resources':RESOURCES,
       'optimizer_software':{'version':scipy.__version__,'wrapper':record(slsqp_wrapper.__file__),'kernel':record(slsqp_kernel.__file__)},
       'declared_calls':declared_calls(arm),
       'mode':'fresh_source_request','source_domain_status':{c['case_id']:'consumed_reference_structural_repeat' for c in req['cases']},
       'cache_policy':'no_archived_energy_or_force_reuse; exact within-run artifacts only',
       'production_changed':False,'automatic_fallback':False,'molecular_submission_status':'not_submitted'}
    existing_candidate.put(out/'plan.json',p);return dry_run(out/'plan.json')


def checked(plan):
    p=read_json(plan);req=read_json(verify(p['request']))
    if (p['protocol_id']!=PROTOCOL or p['arm'] not in ARMS or
        p['profile']!=(PROFILE if p['arm']=='union-candidate' else STATIC_PROFILE) or
        p['settings']!=precision.SETTINGS or p['resources']!=RESOURCES or
        p['production_changed'] or p['automatic_fallback']):raise InvalidArtifact('fresh execution profile differs')
    if check_request(req)!=p['source_identities'] or p['case_ids']!=[c['case_id'] for c in req['cases']]:
        raise InvalidArtifact('declared source identity differs')
    if (p['declared_calls']!=declared_calls(p['arm']) or p['mode']!='fresh_source_request' or
        p['cache_policy']!='no_archived_energy_or_force_reuse; exact within-run artifacts only'):
        raise InvalidArtifact('fresh execution/call/reuse policy changed')
    for pin in p['implementation'].values():verify(pin)
    for k in ('wrapper','kernel'):verify(p['optimizer_software'][k])
    if p['optimizer_software']['version']!=scipy.__version__:raise InvalidArtifact('optimizer software changed')
    verify(p['agreement']);release=standard.release(verify(p['release']))
    rp=verify(p['reference']);ref=read_json(rp)
    if (record(rp)['sha256']!=REFERENCE_SHA or ref['protocol_id']!=precision.POOL_PROTOCOL or
        ref['optimizer_settings']!=precision.SETTINGS or ref['settings']!=union_adaptive.POOL_SETTINGS):
        raise InvalidArtifact('frozen precision reference differs')
    if ref['model']!=req['config']['model'] or ref['crystals_or_noncanonical_used_for_fit']:
        raise InvalidArtifact('reference method or calibration membership differs')
    if any(v['status']!='available' or v['available_calibration']!=25 for v in ref['variants'].values()):
        raise InvalidArtifact('complete fixed canonical reference required')
    q=read_json(verify(p['rank_qualification']));n=read_json(verify(p['maxiter_qualification']))
    if not(q['protocol_id']=='Nikasha_union_reference28_native_GFN2_rank1_scalar_qualification_v1' and
           q['all_numerical_gates_pass'] and q['complete_cells']==336 and q['complete']==28):
        raise InvalidArtifact('rank-one scalar qualification absent')
    if not(n['both_controls_pass'] and n['formerly_failed_cell_complete']):
        raise InvalidArtifact('native MaxIter500 qualification absent')
    return p,req,release,ref


def dry_run(plan):
    p,req,_,_=checked(plan)
    root=Path(plan).parent;prep=root/'source_preparation/preparation.json'
    preparation_status=('prepared' if read_json(prep)['supported']==len(req['cases']) else 'partial') if prep.exists() else 'not_executed'
    return {'plan':record(plan),'profile':p['profile'],'arm':p['arm'],'sources':len(req['cases']),
            'declared_calls':p['declared_calls'],'resources':p['resources'],'new_molecular_calls':0,
            'fresh_source_preparation_status':preparation_status,
            'submission_status':'submitted' if (root/'SUBMISSION.json').exists() else 'not_submitted'}


def fresh_preparation(plan,output):
    """Independent fresh source preparation, using the existing Ca adapter."""
    p,req,_,_=checked(plan);out=Path(output);start=time.monotonic()
    if p['arm']=='union-candidate':
        target=out.parent/'union_preparation'
        source_interface.prepare(verify(p['request']),target,source_mode='fresh',agreement=verify(p['agreement']))
        prepared=read_json(target/'PREPARATION.json')
        data={'source_manifest':p['request'],'config':req['config'],'cases':prepared['cases'],
              'denominator':10,'supported':prepared['supported'],'union_preparation':record(target/'PREPARATION.json'),
              'source_preparation':record(target/'SOURCE_PREPARATION.json')}
    else:
        import accommodation_folds as folds
        rows=[folds.prepare_one(c,req['config'],out) for c in req['cases']]
        data={'source_manifest':p['request'],'config':req['config'],'cases':rows,'denominator':10,
              'supported':sum(c['status']=='prepared' for c in rows)}
    data.update(execution_plan=record(plan),preparation_reused=False,
                fresh_preparation_seconds=time.monotonic()-start,new_molecular_energy_calls=0)
    existing_candidate.put(out/'preparation.json',data);return record(out/'preparation.json')


def same_run_preparation(plan, preparation):
    p,req,_,_=checked(plan);old=read_json(preparation)
    if (old.get('execution_plan')!=record(plan) or old.get('preparation_reused') is not False or
            old['source_manifest']!=p['request'] or old['config']!=req['config'] or
            old['denominator']!=10 or [c['case_id'] for c in old['cases']]!=p['case_ids']):
        raise InvalidArtifact('fresh run requires its own exact completed preparation')
    for src,row in zip(req['cases'],old['cases']):
        if not source_interface.matching_source(src,row.get('source',{})) or row.get('source_preparation_reuse'):
            raise InvalidArtifact('same-run preparation source/reuse differs')
    if p['arm']=='union-candidate':
        actual=read_json(verify(old['source_preparation']))
        if actual['source_mode']!='fresh':raise InvalidArtifact('fresh union source preparation required')
        source_interface.checked(Path(verify(old['union_preparation'])).parent/'plan.json')
    return old


def static_manifest_check(plan, preparation, manifest):
    p,req,release,_=checked(plan);same_run_preparation(plan,preparation)
    m=read_json(manifest);inv=read_json(verify(m['inventory']));prep=read_json(preparation)
    if (m['cases']!=p['case_ids'] or inv['source_preparation']!=record(preparation) or
            inv['cases']!=prep['cases'] or inv['model']!=req['config']['model'] or
            inv['software']!=req['config']['software'] or m['agreement']!=p['agreement'] or
            m['archived_comparison']!=release['artifacts']['calibration'] or
            m['cpu_python_invocation']!=release['cpu_python'] or
            m['gpu_python_invocation']!=release['gpu_python'] or
            {k:v['sha256'] for k,v in m['implementation'].items()}!=
            {k:v['sha256'] for k,v in p['implementation'].items()}):
        raise InvalidArtifact('same-run static scoring manifest differs')
    return static_scanner.validate(manifest)


def static_execution_state(manifest):
    """The released scanner has exclusive-create stage paths, not a resume API."""
    root=Path(manifest).parent;m=read_json(manifest)
    markers=[p for p in (root/'execution_started.json',root/'execution_complete.json',
                        root/'timing',root/'mace') if p.exists()]
    for pin in m['solvent_manifests'].values():
        for task in read_json(verify(pin))['tasks']:
            path=Path(task['output_path']);inp=Path(task['input']['path'])
            markers.extend(p for p in (path,inp.with_name(inp.stem+'.runtime.inp'),
                                       path.with_name(path.name+'.execution.json')) if p.exists())
    return {'status':'already_started_collect_only' if markers else 'not_started',
            'evidence_paths':sorted({str(p) for p in markers}),
            'automatic_molecular_retry':False}


def union_mapping(source,topology,destination):
    from second_shell_context import parent_state
    from coordination_preparation_context import geometry
    state=parent_state(source['original_core'],topology,require_endpoint_receipts=False)
    rep=source['representations']['context'];context=read_json(verify(rep['preparation']));maps=[]
    for z in ('Ca','La'):
        maps.append(geometry(state,context,xyz(verify(rep['endpoints'][z]['xyz'])),
                             xyz(verify(source['original_core']['endpoints'][z]['xyz']))))
    if maps[0]!=maps[1]:raise InvalidArtifact('paired physical union maps differ')
    path=Path(destination)/(source['case_id']+'.json');existing_candidate.put(path,maps[0])
    return {'status':'supported','mapping':record(path)}


def common_manifest(plan):
    p,req,r,_=checked(plan);cp=read_json(verify(req['config']['calibration_implementation_pins']))
    return {'protocol_id':PROTOCOL,'settings':{**shared_pool.SETTINGS,'candidate_order':CANDIDATES},
        'proposal_settings':precision.SETTINGS,'reference':p['reference'],'agreement':p['agreement'],
        'plan':record(plan),'model':req['config']['model'],'software':req['config']['software'],
        'orca':cp['orca_runtime']['executable'],'cpu_python':r['cpu_python'],'gpu_python':r['gpu_python'],
        'implementation':p['implementation'],'GFN2_maxiter':500,'numerical_policy_id':PROFILE,
        'numerical_qualification':p['maxiter_qualification'],'rank_qualification':p['rank_qualification'],
        'optimizer_software':p['optimizer_software']}


def low_prepare(manifest):
    from affordable_workflow import dry_run as old_dry_run
    m=read_json(manifest);root=Path(manifest).parent/'solvent';d=root/'shard_0';tasks=[]
    for cell in m['tasks']:
        for medium in ('vacuum','alpb'):
            td=d/'tasks'/(cell['task_id']+'__'+medium);td.mkdir(parents=True,exist_ok=False)
            xp=td/'core.xyz';shutil.copyfile(verify(cell['xyz']),xp)
            ip=td/'endpoint.inp';ip.write_text(shared_pool.input_text(cell['charge'],cell['multiplicity'],medium,'native').replace('%scf\n','%scf\n MaxIter 500\n'))
            tasks.append({**cell,'task_id':td.name,'cell_id':cell['task_id'],'case':cell['case_id'],
                          'medium':medium,'xyz':record(xp),'input':record(ip),'output_path':str(td/'endpoint.out')})
    lm={'protocol_id':PROTOCOL,'pool_manifest':record(manifest),'agreement':m['agreement'],'orca':m['orca'],
        'tasks':tasks,'all_tasks':tasks,'execution_resources':{'mpi_ranks':1,'concurrent_tasks':32},
        'execution_policy':{'task_runner':m['implementation']['run_orca_task_manifest.py'],
                            'runtime_renderer':m['implementation']['render_orca_runtime_input.py']}}
    path=d/'manifest.json';existing_candidate.put(path,lm)
    existing_candidate.put(root/'INDEX.json',{'pool_manifest':record(manifest),'shards':[record(path)]})
    if tasks:old_dry_run(path)
    return path


def candidate_engine(plan):
    """Private configured instances: no mutation of imported/default executors."""
    engine=private('pqq_adaptive_candidate');pool=private('nikasha_pool');completion=private('adaptive_completion')
    completion.SETTINGS=copy.deepcopy(precision.SETTINGS)
    engine.PROTOCOL=PROTOCOL;engine.completion=completion;engine.pool=pool
    engine.checked=checked;engine.common_manifest=common_manifest;engine.low_prepare=low_prepare
    engine.mapping_case=union_mapping
    engine.pqq_decision_compatible=lambda src: src['status']=='prepared' and src['representations']['context']['protocol_id']=='native_OMOL_GFN2_fixed_group_source_fragment_union_v1'
    engine.source_prepare=SimpleNamespace(prepare=lambda request,output:fresh_preparation(plan,output))
    base_validate=engine.validate_stage
    def validate_stage(manifest,expected_plan):
        m=base_validate(manifest,expected_plan)
        if 'origin_collection' not in m:
            low=read_json(Path(manifest).parent/'solvent/shard_0/manifest.json')
            if (low['pool_manifest']!=record(manifest) or
                low['execution_resources']!={'mpi_ranks':1,'concurrent_tasks':32} or
                len(low['tasks'])!=2*len(m['tasks'])):
                raise InvalidArtifact('qualified rank-one stage policy differs')
            source={t['task_id']:t for t in m['tasks']}
            expected={(t['task_id'],s) for t in m['tasks'] for s in ('vacuum','alpb')}
            if {(t['cell_id'],t['medium']) for t in low['tasks']}!=expected:
                raise InvalidArtifact('scalar stage membership differs')
            for t in low['tasks']:
                old=source[t['cell_id']]
                body=shared_pool.input_text(old['charge'],old['multiplicity'],t['medium'],'native').replace('%scf\n','%scf\n MaxIter 500\n')
                if (verify(t['input']).read_text()!=body or verify(t['xyz']).read_bytes()!=verify(old['xyz']).read_bytes() or
                    (t['charge'],t['multiplicity'])!=(old['charge'],old['multiplicity'])):
                    raise InvalidArtifact('scalar stage geometry/state/recipe changed')
        return m
    engine.validate_stage=validate_stage
    pool.validate=lambda manifest:engine.validate_stage(manifest,plan)
    return engine


def execute(plan):
    p,_,release,_=checked(plan);out=Path(plan).parent
    if Path(__file__).resolve()!=verify(p['implementation']['pqq_union_execution.py']):
        raise InvalidArtifact('execute the immutable implementation snapshot')
    if not os.environ.get('SLURM_JOB_ID') or int(os.environ.get('SLURM_CPUS_ON_NODE','0'))!=32 or int(os.environ.get('SLURM_MEM_PER_NODE','0'))!=200000:
        raise InvalidArtifact('declared32CPU/200000MiB/oneH200 allocation required')
    existing_prep=out/'source_preparation/preparation.json'
    if existing_prep.exists():
        same_run_preparation(plan,existing_prep)
    if p['arm']=='union-candidate':return candidate_engine(plan).execute(plan)
    start=time.monotonic();error=None
    with (out/'execute.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        try:
            if not existing_prep.exists():
                if existing_prep.parent.exists():raise InvalidArtifact('partial source preparation retained; automatic restaging unsupported')
                fresh_preparation(plan,existing_prep.parent)
            same_run_preparation(plan,existing_prep)
            manifest=out/'prepared_score/scoring/manifest.json'
            if not manifest.exists():
                if (out/'prepared_score').exists():raise InvalidArtifact('partial scoring preparation retained; automatic restaging unsupported')
                pqq_fast_release.prepare_score(existing_prep,verify(release['artifacts']['calibration']),verify(p['agreement']),
                                              out/'prepared_score',release['cpu_python'],release['gpu_python'])
            static_manifest_check(plan,existing_prep,manifest)
            state=static_execution_state(manifest)
            result=(static_scanner.execute(manifest) if state['status']=='not_started' else
                    {**state,'manifest':record(manifest),'new_molecular_calls':0})
        except Exception as exc:error=repr(exc)
        finally:
            elapsed=time.monotonic()-start
            existing_candidate.put(out/('execution_'+os.environ['SLURM_JOB_ID']+'.json'),{'plan':record(plan),'error':error,
                'job_id':os.environ['SLURM_JOB_ID'],'wall_seconds':elapsed,'allocated_core_seconds':32*elapsed,
                'allocated_GPU_seconds':elapsed,'source_preparation_included':True})
    if error:raise InvalidArtifact(error)
    return result


def collect(plan,output):
    p,req,release,ref=checked(plan);root=Path(plan).parent;rows=[]
    pp=root/'source_preparation/preparation.json'
    prepared={x['case_id']:x for x in read_json(pp)['cases']} if pp.exists() else {}
    if p['arm']=='released-static':
        manifest=root/'prepared_score/scoring/manifest.json';result=None
        if manifest.exists():
            generated=[path for path in manifest.parent.glob('result_*.json')
                       if re.fullmatch(r'result_[0-9]+\.json',path.name)]
            if len(generated)==1:
                result=read_json(generated[0])
                if result['manifest']!=record(manifest):raise InvalidArtifact('static collection manifest differs')
            elif not generated:
                artifact=Path(output).with_name(Path(output).stem+'__components.json')
                result=static_scanner.collect(manifest,artifact)
            else:raise InvalidArtifact('ambiguous existing static collections')
        by={x['case_id']:x for x in result['rows']} if result else {}
        for source in req['cases']:
            row=by.get(source['case_id']);rows.append({'case_id':source['case_id'],
                'status':row['status'] if row else 'unavailable','reason':row.get('reason') if row else 'fresh_preparation_or_scoring_unavailable',
                'R':row['composite_R_model_kcal_mol'] if row else None,'decision':row['published_PQQ_decision'] if row else 'unavailable',
                'components':row,'source_conditioning_metal':source['source_conditioning_metal']})
    else:
        path=root/'pool/collection.json';result=read_json(path) if path.exists() else None
        by={x['case_id']:x for x in result['cases']} if result else {}
        origin_path=root/'origins/collection.json'
        origins={x['case_id']:x for x in read_json(origin_path)['cases']} if origin_path.exists() else {}
        for source in req['cases']:
            case=by.get(source['case_id']);pool=case['pool'] if case else {'status':'unavailable','reason':'fresh_execution_incomplete'}
            variants={}
            for mode in ('mathematical','operational'):
                value=pool[mode]['composite_R_model_kcal_mol'] if pool['status']=='available' else None
                variants[mode]={'R':value,'decision':union_decision(value,ref['variants'][mode]['bands'])}
            origin=None
            origin_case=case or origins.get(source['case_id'])
            if origin_case and all(origin_case['matrix'][z].get('origin',{}).get('status')=='complete' for z in ('Ca','La')):
                origin=shared_pool.score(origin_case['matrix']['Ca']['origin']['components'],origin_case['matrix']['La']['origin']['components'])
            rows.append({'case_id':source['case_id'],'status':pool['status'],'reason':pool.get('reason') or (case.get('reason') if case else None),
                'R':variants['operational']['R'],'decision':variants['operational']['decision'],'variants':variants,
                'union_origin':origin,'pool':pool,'matrix':case['matrix'] if case else None,
                'candidates':case['candidates'] if case else [],'source_conditioning_metal':source['source_conditioning_metal']})
    for row in rows:
        src=prepared.get(row['case_id'])
        row['preparation_status']=src['status'] if src else 'not_executed'
        row['preparation_reason']=src.get('reason') if src else 'fresh_source_preparation_not_available'
    value={'protocol_id':PROTOCOL,'profile':p['profile'],'arm':p['arm'],'plan':record(plan),
        'reference':p['reference'] if p['arm']=='union-candidate' else release['artifacts']['calibration'],
        'rows':rows,'denominator':10,'available':sum(x['status']=='available' for x in rows),
        'new_calls_in_collection':0,'production_changed':False,'archive_energy_reuse':False,
        'execution_receipts':[record(x) for x in root.glob('execution_*.json')]}
    existing_candidate.put(output,value);return {'collection':record(output),'available':value['available'],'denominator':10}


def compare(static,candidate,output):
    from accommodation_fold_proposals import strict_summary
    a,b=read_json(static),read_json(candidate)
    pa,req,release,_=checked(verify(a['plan']));pb,_,_,ref=checked(verify(b['plan']))
    if (a['arm']!='released-static' or b['arm']!='union-candidate' or pa['request']!=pb['request'] or
        a['archive_energy_reuse'] or b['archive_energy_reuse']):
        raise InvalidArtifact('fresh matched arms/source identity required')
    ids=[c['case_id'] for c in req['cases']]
    if any([x['case_id'] for x in r['rows']]!=ids for r in (a,b)):
        raise InvalidArtifact('all ten sources must remain in declared order')
    sb=read_json(verify(release['artifacts']['calibration']))['calibration']['context']['bands']
    bands={'released_static':{'Ca_max':sb['Ca_supported_max_R_model_kcal_mol'],
                             'La_min':sb['La_supported_min_R_model_kcal_mol']},
           'union_candidate':ref['variants']['operational']['bands']}
    rows=[]
    for src,x,y in zip(req['cases'],a['rows'],b['rows']):
        expected=src['expected_class'];methods={}
        for name,row in (('released_static',x),('union_candidate',y)):
            call=row['decision'];status=('unavailable' if row['R'] is None else 'inconclusive' if call=='inconclusive' else
                                        'correct' if call==expected+'-supported' else 'wrong')
            methods[name]={'R':row['R'],'decision':call,'outcome':status,'status':row['status']}
        rows.append({'case_id':src['case_id'],'source_conditioning_metal':src['source_conditioning_metal'],
                     'canonical_coordinate_match':src['canonical_coordinate_match'],'expected_class':expected,
                     'methods':methods,'delta_R':y['R']-x['R'] if x['R'] is not None and y['R'] is not None else None,
                     'interpretation':'method contrast includes local-to-union context and finite accommodation'})
    normalized={r['case_id']:r for r in rows};g=req['groups'][0]
    la=[c['case_id'] for c in req['cases'] if c['source_conditioning_metal']=='La' and c['case_id']!=g['canonical_case_id']]
    ca=[c['case_id'] for c in req['cases'] if c['source_conditioning_metal']=='Ca']
    sets={'La4':la,'Ca5':ca,**{'La3_'+str(i+1):list(v) for i,v in enumerate(itertools.combinations(la,3))}}
    aggregates=[{'descriptor':name,'members':members,'methods':{method:strict_summary(members,normalized,method,band,'La')
                for method,band in bands.items()}} for name,members in sets.items()]
    balanced={}
    for method,band in bands.items():
        arms=[next(x for x in aggregates if x['descriptor']==name)['methods'][method] for name in ('La4','Ca5')]
        value=None if any(x['R'] is None for x in arms) else sum(x['R'] for x in arms)/2
        balanced[method]={'R':value,'decision':union_decision(value,band),
                          'missing_members':sorted({v for a in arms for v in a['missing_members']})}
    aggregates.append({'descriptor':'balanced','members':la+ca,'methods':balanced})
    from collections import Counter
    common=[r for r in rows if all(v['R'] is not None for v in r['methods'].values())]
    result={'protocol_id':PROTOCOL,'static':record(static),'candidate':record(candidate),'request':pa['request'],
        'rows':rows,'aggregates':aggregates,'source_denominator':10,'independent_protein_groups':1,
        'common_coverage':len(common),'all_counts':{m:dict(Counter(r['methods'][m]['outcome'] for r in rows)) for m in bands},
        'common_counts':{m:dict(Counter(r['methods'][m]['outcome'] for r in common)) for m in bands},
        'calibration_refitted':False,'new_molecular_calls':0,'production_changed':False,
        'timing_caution':'Different released/candidate batching and model loads; not an isolated speedup of accommodation or MPI ranks.'}
    existing_candidate.put(output,result);return {'comparison':record(output),'common_coverage':len(common),'all_counts':result['all_counts']}


def report(result,output):
    r=read_json(result);checked(verify(r['plan']))
    lines=['# Fresh PQQ execution: '+r['arm'],'',f"Available {r['available']}/{r['denominator']}; reference and defaults unchanged.",'',
           '| Source | R (model kcal/mol) | Decision | Status |','|---|---:|---|---|']
    for row in r['rows']:lines.append(f"| {row['case_id']} | {row['R']} | {row['decision']} | {row['status']}: {row.get('reason') or ''} |")
    lines+=['','These are consumed structural repeats of one reference protein. Candidate geometry selection is a finite energy comparison, not an affinity probability or stationary-minimum claim.']
    with Path(output).open('x') as f:f.write('\n'.join(lines)+'\n')
    return {'report':record(output)}


def main():
    p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest='op',required=True)
    for op,fields in {'prepare':('request','arm','reference','rank_qualification','maxiter_qualification','agreement','output'),
                     'dry_run':('plan',),'execute':('plan',),'collect':('plan','output'),'report':('result','output'),
                     'compare':('static','candidate','output')}.items():
        q=sub.add_parser(op.replace('_','-'))
        for f in fields:q.add_argument('--'+f.replace('_','-'),required=True)
    args=vars(p.parse_args());op=args.pop('op').replace('-','_');print(json.dumps(globals()[op](**args),indent=2))

if __name__=='__main__':main()
