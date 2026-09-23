"""Opt-in execution adapter for explicit prepared three-La-source PQQ groups."""
from __future__ import annotations
import argparse
import copy
import json
import os
from pathlib import Path
import statistics
import time
import scipy
import scipy.optimize._slsqp_py as slsqp_wrapper
import scipy.optimize._slsqplib as slsqp_kernel
from affordable_common import InvalidArtifact,read_json,record,verify
from adaptive_origin_recovery import snapshot
import pqq_three_source as preparation
import pqq_union_execution as union_execution
import pqq_adaptive_candidate as candidate
import pqq_standard as standard
import nikasha_pool as pool
from accommodation_folds_compare import decision

PROTOCOL='Nikasha_explicit_three_La_source_union_execution_v1'
PROFILE=union_execution.PROFILE
RESOURCES=union_execution.RESOURCES
CANDIDATES=union_execution.CANDIDATES
CACHE='reuse_exact_source_preparation_only; no_archive_energy_or_force; exact_within_run_only'


def declared_calls(n):
    return {'native_MACE_q0_gradient':2*n,'bounded_native_MACE_searches':2*n,
            'maximum_cross_MACE':2*n,'maximum_GFN2':12*n,'DFT':0,'new_protonation':0}


def request_from_preparations(paths):
    groups=[];cases=[];config=None;release=None;reference=None
    for path in paths:
        preparation.dry_run(path);p=read_json(path);r=preparation.check_request(read_json(verify(p['request'])))
        if config is None:config=r['config'];release=r['release'];reference=r['reference']
        if (r['config'],r['release'],r['reference'])!=(config,release,reference):
            raise InvalidArtifact('groups require the same qualified source configuration and reference')
        groups.append({'protein_id':r['protein_id'],'members':[c['case_id'] for c in r['cases']],
            'preparation':record(path),'request':p['request'],'physical_state_anchor':p['physical_state_anchor'],
            'preparation_status':p['status'],'source_evidence':r['source_evidence']})
        cases.extend(r['cases'])
    if not groups or len({g['protein_id'] for g in groups})!=len(groups) or len({c['case_id'] for c in cases})!=len(cases):
        raise InvalidArtifact('unique explicit groups and source identifiers required')
    return {'protocol_id':PROTOCOL,'groups':groups,'cases':cases,'config':config,'release':release,'reference':reference}


def prepare(preparations,rank_qualification,maxiter_qualification,agreement,output):
    req=request_from_preparations(preparations);out=Path(output).resolve()
    if 'workspaces' not in out.parts:raise InvalidArtifact('candidate products belong under workspaces')
    out.mkdir(parents=True,exist_ok=False);rp=out/'REQUEST.json';candidate.put(rp,req)
    p={'protocol_id':PROTOCOL,'profile':PROFILE,'request':record(rp),'release':req['release'],'reference':req['reference'],
       'rank_qualification':record(rank_qualification),'maxiter_qualification':record(maxiter_qualification),
       'agreement':record(agreement),'implementation':snapshot(Path(__file__).parent,out/'implementation'),
       'settings':preparation.precision.SETTINGS,'resources':RESOURCES,'case_ids':[c['case_id'] for c in req['cases']],
       'mode':'explicit_prepared_three_source_request','cache_policy':CACHE,'automatic_fallback':False,
       'production_changed':False,'declared_calls':declared_calls(len(req['cases'])),
       'optimizer_software':{'version':scipy.__version__,'wrapper':record(slsqp_wrapper.__file__),'kernel':record(slsqp_kernel.__file__)},
       'source_domain_status':{c['case_id']:'explicit_source_protocol_compatible_biological_domain_unvalidated' for c in req['cases']},
       'molecular_submission_status':'not_submitted','new_molecular_calls_in_preparation':0}
    path=out/'plan.json';candidate.put(path,p);adopt_preparation(path,out/'source_preparation')
    engine(path).origin_manifest(path,out/'source_preparation/preparation.json',out/'origins')
    return dry_run(path)


def checked(plan):
    p=read_json(plan);req=read_json(verify(p['request']))
    if (p['protocol_id']!=PROTOCOL or p['profile']!=PROFILE or p['resources']!=RESOURCES or
        p['settings']!=preparation.precision.SETTINGS or p['production_changed'] or p['automatic_fallback'] or
        p['cache_policy']!=CACHE or p['mode']!='explicit_prepared_three_source_request' or p['new_molecular_calls_in_preparation']):
        raise InvalidArtifact('explicit three-source execution policy differs')
    expected=request_from_preparations([verify(g['preparation']) for g in req['groups']])
    if expected!=req or p['case_ids']!=[c['case_id'] for c in req['cases']] or p['declared_calls']!=declared_calls(len(req['cases'])):
        raise InvalidArtifact('explicit preparation/source membership or call count differs')
    if (p['release'],p['reference'])!=(req['release'],req['reference']):raise InvalidArtifact('execution source/reference differs')
    ref=preparation.check_reference(p['reference'],req['config']);release=standard.release(verify(p['release']))
    for pin in p['implementation'].values():verify(pin)
    for k in ('wrapper','kernel'):verify(p['optimizer_software'][k])
    if p['optimizer_software']['version']!=scipy.__version__:raise InvalidArtifact('optimizer software differs')
    verify(p['agreement'])
    rank=read_json(verify(p['rank_qualification']));numerical=read_json(verify(p['maxiter_qualification']))
    if not(rank['protocol_id']=='Nikasha_union_reference28_native_GFN2_rank1_scalar_qualification_v1' and
           rank['all_numerical_gates_pass'] and rank['complete_cells']==336 and rank['complete']==28):
        raise InvalidArtifact('rank-one scalar qualification absent')
    if not(numerical['both_controls_pass'] and numerical['formerly_failed_cell_complete']):
        raise InvalidArtifact('native MaxIter500 qualification absent')
    return p,req,release,ref


def source_envelope(plan):
    p,req,_,_=checked(plan);rows=[]
    for group in req['groups']:
        original=read_json(verify(group['preparation']))
        for source in original['cases']:
            row=copy.deepcopy(source)
            if original['status']!='prepared':
                row.update(source_preparation_status=source['status'],status='group_preparation_unavailable',
                           reason='required three-source preparation unavailable: '+', '.join(original['missing_members']))
            rows.append(row)
    return {'source_manifest':p['request'],'config':req['config'],'cases':rows,'denominator':len(rows),
            'supported':sum(c['status']=='prepared' for c in rows),'execution_plan':record(plan),
            'source_preparations':[g['preparation'] for g in req['groups']],
            'preparation_reused':True,'new_protonation':False,'new_molecular_energy_calls':0}


def adopt_preparation(plan,output):
    path=Path(output)/'preparation.json';candidate.put(path,source_envelope(plan));return record(path)


def engine(plan):
    # Configure private copies only; the tested orchestration and physical search are unchanged.
    own=union_execution.private('pqq_union_execution');own.PROTOCOL=PROTOCOL;own.PROFILE=PROFILE
    own.checked=checked;own.fresh_preparation=adopt_preparation
    own.union_mapping=lambda source,topology,destination:{'status':'supported','mapping':source['physical_mapping']}
    return own.candidate_engine(plan)


def dry_run(plan):
    p,req,_,_=checked(plan);root=Path(plan).parent
    if read_json(root/'source_preparation/preparation.json')!=source_envelope(plan):raise InvalidArtifact('adopted source state differs')
    e=engine(plan);m=e.validate_stage(root/'origins/manifest.json',plan)
    if len(m['tasks'])>p['declared_calls']['native_MACE_q0_gradient']:raise InvalidArtifact('origin count exceeds declaration')
    if any(t.get('native_reuse') for t in m['tasks']):raise InvalidArtifact('archived origin energy/force reuse prohibited')
    low=read_json(root/'origins/solvent/shard_0/manifest.json')
    from affordable_workflow import dry_run as solver_dry_run
    if low['tasks']:solver_dry_run(root/'origins/solvent/shard_0/manifest.json')
    return {'plan':record(plan),'status':'preflight_passed','groups':len(req['groups']),'sources':len(req['cases']),
            'prepared_origin_sources':len(m['tasks'])//2,'origin_MACE_tasks':len(m['tasks']),
            'origin_GFN2_tasks':len(low['tasks']),'declared_calls':p['declared_calls'],'resources':p['resources'],
            'source_preparation_reused':True,'new_molecular_calls':0,
            'submission_status':'submitted' if (root/'SUBMISSION.json').exists() else 'not_submitted'}


def execute(plan):
    p,_,_,_=checked(plan);root=Path(plan).parent
    if Path(__file__).resolve()!=verify(p['implementation']['pqq_three_source_execution.py']):
        raise InvalidArtifact('execute the immutable implementation snapshot')
    if not os.environ.get('SLURM_JOB_ID') or int(os.environ.get('SLURM_CPUS_ON_NODE','0'))!=32 or int(os.environ.get('SLURM_MEM_PER_NODE','0'))!=200000:
        raise InvalidArtifact('32CPU/200000MiB/oneH200 allocation required')
    dry_run(plan)
    if list(root.glob('execution_*.json')) or (root/'EXECUTION_STARTED.json').exists():
        raise InvalidArtifact('execution already started; collect existing outputs, no automatic scientific retry')
    candidate.put(root/'EXECUTION_STARTED.json',{'plan':record(plan),'job_id':os.environ['SLURM_JOB_ID'],
        'time_unix':time.time(),'source_preparation_reused':True,'archive_energy_or_force_reused':False})
    return engine(plan).execute(plan)


def source_result(source,case,origin_case,ref):
    selected=case['pool'] if case else {'status':'unavailable','reason':'execution_or_collection_incomplete'}
    if case and selected['status']=='available':
        if [q['id'] for q in case['candidates']]!=CANDIDATES or pool.choose_rows(case['matrix'],CANDIDATES)!=selected:
            raise InvalidArtifact('required common pool changed or incomplete')
    original=case or origin_case;origin=None
    if original and all(original['matrix'].get(z,{}).get('origin',{}).get('status')=='complete' for z in ('Ca','La')):
        origin=pool.score(original['matrix']['Ca']['origin']['components'],original['matrix']['La']['origin']['components'])
    origin_r=origin['composite_R_model_kcal_mol'] if origin else None
    variants={}
    for mode in ('mathematical','operational'):
        value=selected[mode]['composite_R_model_kcal_mol'] if selected['status']=='available' else None
        variants[mode]={'R_model_kcal_mol':value,'decision':decision(value,ref['variants'][mode]['bands']),
                        'delta_R_from_union_origin':None if value is None or origin_r is None else value-origin_r}
    return {'case_id':source['case_id'],'status':selected['status'],
            'reason':selected.get('reason') or (case.get('reason') if case else 'execution_or_collection_incomplete'),
            'union_origin_R_model_kcal_mol':origin_r,'union_origin_components':origin,
            'variants':variants,'pool':selected,'matrix':case['matrix'] if case else None,
            'candidates':case['candidates'] if case else [],'affinity_probability':None}


def group_results(groups,rows,ref):
    by={r['case_id']:r for r in rows};result=[]
    for group in groups:
        members=[by[c] for c in group['members']];variants={}
        for mode in ('mathematical','operational'):
            values=[r['variants'][mode]['R_model_kcal_mol'] for r in members]
            value=statistics.median(values) if len(values)==3 and all(v is not None for v in values) else None
            variants[mode]={'R_model_kcal_mol':value,'decision':decision(value,ref['variants'][mode]['bands']),
                'member_R_model_kcal_mol':values,'range_kcal_mol':max(values)-min(values) if value is not None else None,
                'rule':'strict_median_all_three_declared_sources'}
        result.append({'protein_id':group['protein_id'],'members':group['members'],'denominator':3,
            'available':sum(r['status']=='available' for r in members),'variants':variants,
            'source_evidence':group['source_evidence'],'biological_accuracy_evaluated':False,
            'physical_state_anchor':group['physical_state_anchor'],'classification_is_prediction':True})
    return result


def collect(plan,output):
    p,req,_,ref=checked(plan);root=Path(plan).parent;e=engine(plan);stages={};pins={}
    for name in ('origins','pool'):
        mp=root/name/'manifest.json'
        if not mp.exists():continue
        e.validate_stage(mp,plan)
        # Report-only snapshots never install a partial stage's authoritative collection.json.
        dest=Path(output).with_name(Path(output).stem+'__'+name+'.json')
        e.pool.collect(mp,dest);stages[name]=read_json(dest);pins[name]=record(dest)
    origins={x['case_id']:x for x in stages.get('origins',{}).get('cases',[])}
    cases={x['case_id']:x for x in stages.get('pool',{}).get('cases',[])}
    rows=[source_result(c,cases.get(c['case_id']),origins.get(c['case_id']),ref) for c in req['cases']]
    for group in req['groups']:
        for row in rows:
            if row['case_id'] in group['members']:
                row['protein_id']=group['protein_id'];row['source_evidence']=group['source_evidence'][row['case_id']]
    value={'protocol_id':PROTOCOL,'profile':PROFILE,'plan':record(plan),'reference':p['reference'],
        'rows':rows,'groups':group_results(req['groups'],rows,ref),'source_denominator':len(rows),
        'source_available':sum(r['status']=='available' for r in rows),'group_denominator':len(req['groups']),
        'source_collections':pins,'production_changed':False,'new_calls_in_collection':0,'new_reference_fitted':False,
        'archive_energy_or_force_reused':False,'source_preparation_reused':True,'biological_accuracy_evaluated':False,
        'execution_receipts':[record(q) for q in sorted(root.glob('execution_*.json'))],
        'receipt_timing_qualification':'source_preparation_included in reused engine receipt refers only to exact prepared-envelope handoff, not fresh protonation'}
    candidate.put(output,value);return {'result':record(output),'sources':len(rows),'available':value['source_available']}


def report(result,output):
    d=read_json(result);checked(verify(d['plan']))
    lines=['# Explicit three-source opt-in PQQ candidate','',f"Available sources: {d['source_available']}/{d['source_denominator']}. Biological labels were not evaluated.",'',
        '| Source | Union origin R | Accommodated R | Prediction | Status |','|---|---:|---:|---|---|']
    for r in d['rows']:
        v=r['variants']['operational'];lines.append(f"| {r['case_id']} | {r['union_origin_R_model_kcal_mol']} | {v['R_model_kcal_mol']} | {v['decision']} | {r['status']}: {r['reason'] or ''} |")
    lines+=['','## Strict three-source summaries','','| Protein | Available | Median R | Prediction |','|---|---:|---:|---|']
    for g in d['groups']:
        v=g['variants']['operational'];lines.append(f"| {g['protein_id']} | {g['available']}/3 | {v['R_model_kcal_mol']} | {v['decision']} |")
    lines+=['','R is the raw Ca-minus-La model contrast, in kcal/mol, not an aquo-referenced binding free energy.',
        'The original value is for the same fixed union; it is not the released local-context score.',
        'The candidate uses its pinned threefold reference. Missing members remain unavailable.',
        'Source labels, both numerical variants, full matrices, selections, work and geometry flags remain in the result JSON.',
        'Preparation reused exact archived physical states; endpoint energies and forces are fresh for this run. Defaults are unchanged.']
    with Path(output).open('x') as f:f.write('\n'.join(lines)+'\n')
    return {'report':record(output)}


def main():
    ap=argparse.ArgumentParser(description=__doc__);sub=ap.add_subparsers(dest='op',required=True)
    q=sub.add_parser('prepare');q.add_argument('--preparations',nargs='+',required=True)
    for f in ('rank-qualification','maxiter-qualification','agreement','output'):q.add_argument('--'+f,required=True)
    for op in ('dry-run','execute'):
        q=sub.add_parser(op);q.add_argument('--plan',required=True)
    q=sub.add_parser('collect');q.add_argument('--plan',required=True);q.add_argument('--output',required=True)
    q=sub.add_parser('report');q.add_argument('--result',required=True);q.add_argument('--output',required=True)
    args=vars(ap.parse_args());op=args.pop('op').replace('-','_');print(json.dumps(globals()[op](**args),indent=2))

if __name__=='__main__':main()
