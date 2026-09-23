"""Thin strict-native execution profile for explicit prepared three-source envelopes."""
from __future__ import annotations
import argparse
import copy
import json
import os
from pathlib import Path
import statistics
import scipy
import scipy.optimize._slsqp_py as slsqp_wrapper
import scipy.optimize._slsqplib as slsqp_kernel
from affordable_common import InvalidArtifact,read_json,record,verify
from adaptive_origin_recovery import snapshot
import pqq_three_source_envelope as preparation
import pqq_union_execution as union
import pqq_adaptive_candidate as candidate
import motion_envelope_scalar as scalar

PROTOCOL='Nikasha_explicit_three_La_envelope4p3_strict_native_execution_v1'
PROFILE=scalar.PROFILE
RESOURCES=union.RESOURCES
CACHE='exact_archived_protonation_only; fresh_energy_force; same_run_reuse_only'
PROTEINS=('PQQSEQ_07ab500e3df76b30d71c','PQQSEQ_83440678cbbd658047c9')


def interface():
    ex=union.private('pqq_three_source_execution')
    ex.PROTOCOL=PROTOCOL;ex.PROFILE=PROFILE;ex.CACHE=CACHE
    ex.preparation=preparation;ex.checked=checked;ex.engine=engine
    return ex


def checked(plan):
    p=read_json(plan);req=read_json(verify(p['request']))
    if (p['protocol_id']!=PROTOCOL or p['profile']!=PROFILE or p['resources']!=RESOURCES or
        p['settings']!=preparation.precision.SETTINGS or p['production_changed'] or p['automatic_fallback'] or
        p['cache_policy']!=CACHE or p['mode']!='explicit_prepared_three_source_request' or p['new_molecular_calls_in_preparation']):
        raise InvalidArtifact('explicit envelope scientific execution profile differs')
    ex=interface();actual=ex.request_from_preparations([verify(g['preparation']) for g in req['groups']])
    if actual!=req or [g['protein_id'] for g in req['groups']]!=list(PROTEINS):
        raise InvalidArtifact('fixed two previously declared PLM triples required')
    if p['case_ids']!=[c['case_id'] for c in req['cases']] or p['declared_calls']!=ex.declared_calls(6):
        raise InvalidArtifact('exact finite six-source scope differs')
    if (p['reference'],p['release'])!=(req['reference'],req['release']):raise InvalidArtifact('source/reference differs')
    if preparation.reference_status(p['reference'],req['config'])['status']!='available':
        raise InvalidArtifact('no separating actual envelope reference')
    scalar.qualified(verify(p['strict_qualification']))
    for pin in p['implementation'].values():verify(pin)
    for key in ('wrapper','kernel'):verify(p['optimizer_software'][key])
    if p['optimizer_software']['version']!=scipy.__version__:raise InvalidArtifact('optimizer runtime changed')
    for key in ('agreement','rank_qualification','maxiter_qualification'):verify(p[key])
    return p,req,preparation.original.standard.release(verify(p['release'])),read_json(verify(p['reference']))


def strict_base_input(charge,multiplicity,medium,solver):
    if multiplicity!=1 or solver!='native':raise InvalidArtifact('native closed-shell scalar required')
    # The existing union runner inserts MaxIter once; expose the exact remaining
    # strict recipe so its builder and validator reproduce the same final text.
    return scalar.recipe(charge,medium,'fresh').replace(' MaxIter 500\n','')


def engine(plan):
    ex=interface();own=union.private('pqq_union_execution')
    own.PROTOCOL=PROTOCOL;own.PROFILE=PROFILE;own.checked=checked
    own.fresh_preparation=ex.adopt_preparation
    own.union_mapping=lambda row,topology,out:{'status':'supported','mapping':row['physical_mapping']}
    own.shared_pool=union.private('nikasha_pool');own.shared_pool.input_text=strict_base_input
    e=own.candidate_engine(plan)
    def completed(manifest,task_id):
        task=next(t for t in read_json(manifest)['tasks'] if t['task_id']==task_id)
        result=scalar.endpoint(manifest,task)
        if result['status']!='complete':
            raise InvalidArtifact(result.get('reason','strict native scalar unavailable'))
        return result
    e.pool.completed=completed
    return e


def prepare(preparations,strict_qualification,rank_qualification,maxiter_qualification,agreement,output):
    ex=interface();req=ex.request_from_preparations(preparations);out=Path(output).resolve()
    if 'workspaces' not in out.parts:raise InvalidArtifact('candidate products belong under workspaces')
    out.mkdir(parents=True,exist_ok=False);candidate.put(out/'REQUEST.json',req)
    p={'protocol_id':PROTOCOL,'profile':PROFILE,'request':record(out/'REQUEST.json'),'release':req['release'],'reference':req['reference'],
       'strict_qualification':scalar.qualified(strict_qualification),'rank_qualification':record(rank_qualification),
       'maxiter_qualification':record(maxiter_qualification),'agreement':record(agreement),
       'implementation':snapshot(Path(__file__).parent,out/'implementation'),'settings':preparation.precision.SETTINGS,
       'resources':RESOURCES,'case_ids':[c['case_id'] for c in req['cases']],
       'mode':'explicit_prepared_three_source_request','cache_policy':CACHE,'automatic_fallback':False,'production_changed':False,
       'declared_calls':ex.declared_calls(6),'optimizer_software':{'version':scipy.__version__,
            'wrapper':record(slsqp_wrapper.__file__),'kernel':record(slsqp_kernel.__file__)},
       'source_domain_status':{c['case_id']:'experimental_PLM_prediction_unknown_label' for c in req['cases']},
       'molecular_submission_status':'not_submitted','new_molecular_calls_in_preparation':0}
    plan=out/'plan.json';candidate.put(plan,p);ex.adopt_preparation(plan,out/'source_preparation')
    engine(plan).origin_manifest(plan,out/'source_preparation/preparation.json',out/'origins')
    result=dry_run(plan);candidate.put(out/'PREFLIGHT.json',result);return result


def dry_run(plan):
    result=interface().dry_run(plan);p,_,_,_=checked(plan)
    if result['origin_MACE_tasks']!=12 or result['origin_GFN2_tasks']!=24:
        raise InvalidArtifact('six complete paired origins required for this integration')
    root=Path(plan).parent;low=read_json(root/'origins/solvent/shard_0/manifest.json')
    for task in low['tasks']:
        if verify(task['input']).read_text()!=scalar.recipe(task['charge'],task['medium'],'fresh'):
            raise InvalidArtifact('staged recipe differs from qualified strict native')
    result.update(strict_qualification=p['strict_qualification'],biological_accuracy_evaluated=False,
                  scalar_profile=PROFILE,actual_new_molecular_calls=0)
    return result


def execute(plan):
    p,_,_,_=checked(plan)
    if Path(__file__).resolve()!=verify(p['implementation']['pqq_three_source_envelope_execution.py']):
        raise InvalidArtifact('execute the immutable prepared implementation')
    return interface().execute(plan)


def collect(plan,output):
    ex=interface();raw=Path(output).with_name(Path(output).stem+'__engine.json')
    if not raw.exists():ex.collect(plan,raw)
    d=read_json(raw)
    if d['plan']!=record(plan):raise InvalidArtifact('collection belongs to another plan')
    rows={r['case_id']:r for r in d['rows']}
    for group in d['groups']:
        values=[rows[c]['union_origin_R_model_kcal_mol'] for c in group['members']]
        complete=all(v is not None for v in values)
        group['union_origin']={'member_R_model_kcal_mol':values,
            'median_R_model_kcal_mol':statistics.median(values) if complete else None,
            'range_kcal_mol':max(values)-min(values) if complete else None,'decision':None,
            'reason':'origin-only contrast shown without substituting adaptive reference'}
    d.update(engine_collection=record(raw),candidate_status='experimental_not_promoted',
             calibration_limitation='canonical25 and3 crystals retained; A0A3Ca3 probe inconclusive; full transfer separate',
             scope='exact_two_consumed_unlabeled_PLM_triples_usability_not_accuracy',
             scalar_profile=PROFILE,original_DFT_and_production_results_changed=False)
    candidate.put(output,d)
    return {'result':record(output),'sources':d['source_denominator'],'available':d['source_available'],
            'groups':d['group_denominator']}


def report(result,output):
    d=read_json(result);checked(verify(d['plan']))
    ex=interface();r=ex.report(result,output)
    # A separate footer keeps the existing report renderer unchanged.
    with Path(output).open('a') as f:
        f.write('\n\nStrict native scalar profile: '+PROFILE+'.\n'
                'This is a usability experiment on two unlabeled, consumed PLM triples, not accuracy validation.\n'
                'The actual envelope calibration preserves25canonical+3crystal calls, but A0A3Ca3 is inconclusive.\n'
                'Source protonation is archived reuse; all molecular work is fresh within this run.\n')
    return {'report':record(output)}


def main():
    ap=argparse.ArgumentParser(description=__doc__);sub=ap.add_subparsers(dest='op',required=True)
    q=sub.add_parser('prepare');q.add_argument('--preparations',nargs='+',required=True)
    for k in ('strict-qualification','rank-qualification','maxiter-qualification','agreement','output'):q.add_argument('--'+k,required=True)
    for op,keys in {'dry-run':('plan',),'execute':('plan',),'collect':('plan','output'),'report':('result','output')}.items():
        q=sub.add_parser(op)
        for k in keys:q.add_argument('--'+k,required=True)
    args=vars(ap.parse_args());op=args.pop('op').replace('-','_');print(json.dumps(globals()[op](**args),indent=2))

if __name__=='__main__':main()
