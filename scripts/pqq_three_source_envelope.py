"""Prepare explicit three-source PQQ motion-envelope inputs; no molecular execution."""
from __future__ import annotations
import argparse
import copy
import json
from pathlib import Path
import re
import time
import types
import numpy as np

from affordable_common import InvalidArtifact,cache_key,paired,read_json,record,verify,write_new,xyz
import pqq_three_source as original
import motion_envelope as envelope
import motion_envelope_pool as pool
import motion_envelope_scalar as scalar
import slsqp_precision as precision
from mace_site_kinematics import Kinematics
from nikasha_pool_compare import extrema_reference

ROOT=Path(__file__).resolve().parents[1]
PROTOCOL='Nikasha_explicit_three_La_motion_envelope4p3_preparation_v1'
REFERENCE_ID='Nikasha_motion_envelope4p3_adaptive_strict_native_canonical25_v1'
REFERENCE_SHA='4a1d3ae6e83b5b6ace408f729c39390d44f673554e338df932fd38fc971f67d9'
POLICY={**original.POLICY,'union':'all_three_declared_fixed_4p3_A_source_fragment_lists',
        'polar_neighbor_cutoff_A':4.3,'original_anchors':'exact_saved_anchor_replay',
        'source_preparation':'reuse_exact_archived_protonation_only'}
IMPLEMENTATIONS=('pqq_three_source_envelope.py','pqq_three_source.py','pqq_union_candidate.py',
 'motion_envelope.py','consistent_context.py','second_shell_context.py','environment_context_chemistry.py',
 'union_triple_transfer.py','coordination_preparation_context.py','mace_site_kinematics.py',
 'motion_envelope_run.py','motion_envelope_scalar.py','motion_envelope_pool.py','strict_native_pool.py')


def reference_status(pin,config):
    if pin is None:return {'status':'unsupported_reference','reason':'actual envelope canonical reference is absent','bands':None}
    r=read_json(verify(pin))
    if (pin['sha256']!=REFERENCE_SHA or r.get('protocol_id')!=pool.POOL_PROTOCOL or r.get('reference_id')!=REFERENCE_ID or
        r.get('canonical_denominator')!=25 or r.get('crystals_or_probes_used_for_fit') is not False or
        r.get('model')!=config['model'] or r.get('pool_settings')!=pool.POOL_SETTINGS or
        r.get('optimizer_settings')!=precision.SETTINGS or r.get('scalar_profile')!=scalar.PROFILE):
        raise InvalidArtifact('incompatible actual envelope reference/profile; historical bands are not substitutes')
    verify(r['collection']);verify(r['inventory']);verify(r['implementation'])
    members=[v for v in r['rows'] if v['role']=='canonical_calibration']
    if len(members)!=25 or len({v['root_case_id'] for v in members})!=25:
        raise InvalidArtifact('exact canonical25 reference membership required')
    for variant in ('mathematical','operational'):
        expected=extrema_reference([{**v,'R_model_kcal_mol':v['scores'][variant]} for v in members],variant,REFERENCE_ID)
        fields=('status','bands','gap_model_kcal_mol','available_calibration','calibration_denominator','class_extrema')
        if any(r['variants'][variant].get(k)!=expected.get(k) for k in fields):
            raise InvalidArtifact('actual canonical reference algebra differs')
    if any(v['status']!='available' or v['available_calibration']!=25 or v['bands'] is None for v in r['variants'].values()):
        return {'status':'unsupported_reference','reason':'envelope canonical25 does not supply a complete separating reference','bands':None}
    return {'status':'available','reference_id':REFERENCE_ID,
            'bands':{v:r['variants'][v]['bands'] for v in ('mathematical','operational')}}


def check_request(req):
    # Reuse the exact physical source validator with private globals. The older
    # public API and its reference/policy constants are never changed.
    gl=dict(original.check_request.__globals__)
    gl.update(PROTOCOL=PROTOCOL,POLICY=POLICY,check_reference=reference_status)
    checked=types.FunctionType(original.check_request.__code__,gl,'check_envelope_sources')(req)
    actual=reference_status(req['reference'],req['config'])
    if req['reference_status']!=actual:raise InvalidArtifact('frozen reference readiness differs')
    if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*',req['protein_id']):raise InvalidArtifact('safe protein ID required')
    if req['anchor_is_biological_label_or_calibration_member'] or req['new_molecular_calls']:
        raise InvalidArtifact('physical source anchor is not a label or an executed calculation')
    return checked


def request(sources,config,reference,protein_id,agreement,output,release=original.standard.DEFAULT_RELEASE):
    cfg=original.configuration(config);rows=original.descriptors(read_json(sources));cases=[]
    for row in rows:
        c={k:copy.deepcopy(row[k]) for k in original.SOURCE_FIELDS}
        if isinstance(c['source_structure'],str):c['source_structure']=record(c['source_structure'])
        cases.append(c)
    cases.sort(key=lambda c:c['case_id'])
    if len(cases)!=3:raise InvalidArtifact('exactly three declared source structures required')
    ref=record(reference) if Path(reference).is_file() else None
    req={'protocol_id':PROTOCOL,'policy':POLICY,'protein_id':protein_id,'cases':cases,'config':cfg,
         'source_identity':cache_key(original.source_identity(cases[0],cfg)),
         'physical_state_anchor':cases[0]['case_id'],'reference':ref,'reference_requested_path':str(Path(reference).resolve()),
         'reference_status':reference_status(ref,cfg),'release':record(release),'agreement':record(agreement),
         'source_descriptors':record(sources),'configuration_source':record(config),
         'source_evidence':{r['case_id']:{k:copy.deepcopy(r.get(k)) for k in ('expected_class','evidence_stratum','label_scope','source_provenance')} for r in rows},
         'anchor_is_biological_label_or_calibration_member':False,'production_changed':False,'new_molecular_calls':0}
    check_request(req);write_new(output,req)
    return {'request':record(output),'status':req['reference_status']['status'],'sources':3,'new_molecular_calls':0}


def finite_calls(available):
    n=3 if available else 0
    return {'origin_MACE_force_calls':2*n,'maximum_bounded_MACE_searches':2*n,
            'maximum_cross_MACE_calls':2*n,'maximum_strict_GFN2_scalar_calls':12*n,
            'DFT_calls':0,'protonation_calls':0,'new_folds':0}


def prepare(request,source_preparation,output):
    start=time.monotonic();req=check_request(read_json(request));out=Path(output).resolve()
    if 'workspaces' not in out.parts:raise InvalidArtifact('candidate products belong under workspaces')
    out.mkdir(parents=True,exist_ok=False)
    old=read_json(source_preparation)
    if old['config']!=req['config']:raise InvalidArtifact('archived source preparation configuration differs')
    sources=[]
    for src in req['cases']:
        matching=[r for r in old['cases'] if original.matching_source(src,r.get('source',{}))]
        if len(matching)!=1:raise InvalidArtifact('one exact archived source preparation required for '+src['case_id'])
        sources.append(copy.deepcopy(matching[0]))
    write_new(out/'SOURCE_PREPARATION.json',{'request':record(request),'config':req['config'],'cases':sources,
        'source_archive':record(source_preparation),'source_mode':'reuse-exact','new_protonation':False})
    missing=[r['case_id'] for r in sources if r['status']!='prepared'];discovered=[];rows=[];tasks=[];scalar_tasks=[];group=None
    ready_ref=req['reference_status']['status']=='available'
    if ready_ref and not missing:
        discovered=[envelope.discover_source(r,req['config']) for r in sources]
        missing=[r['case_id'] for r in discovered if r['status']!='discovered']
        write_new(out/'DISCOVERY.json',{'policy':envelope.POLICY,'cases':discovered,'new_molecular_calls':0})
        if not missing:
            keys=set().union(*({original.context.fragment_id(f) for f in r['fragments']} for r in discovered))
            group={'root_case_id':req['protein_id'],'members':[r['case_id'] for r in sources],
                   'physical_state_anchor':req['physical_state_anchor'],'fragments':[original.context.fragment_record(k) for k in sorted(keys)],
                   'source_request':record(request),'discovery':record(out/'DISCOVERY.json'),'selection_rule':POLICY['union'],
                   'anchor_has_no_label_or_calibration_role':True}
            write_new(out/'UNION.json',group);group['union']=record(out/'UNION.json')
            rows=[envelope.prepare_task(r,group,req['config'],out/'union') for r in sources]
            anchor=next(r for r in rows if r['case_id']==req['physical_state_anchor'])
            for r in rows:
                if r['status']!='prepared_awaiting_group_check':continue
                if (not anchor.get('mapping_eligible') or r['state_key']!=anchor.get('state_key') or r['protein_key']!=anchor.get('protein_key')):
                    r.update(status='source_state_mismatch',reason='exact three-source atom/proton/cap/bond/charge/water/protein state differs')
                else:
                    r.update(status='prepared',state_anchor_signature=anchor['state_signature'],
                        physical_state_anchor=req['physical_state_anchor'],physical_mapping=r['maps']['Ca'])
            missing=[r['case_id'] for r in rows if r['status']!='prepared']
            if not missing:
                for r in rows:
                    for z in ('Ca','La'):
                        ep=r['representations']['context']['endpoints'][z]
                        task={'task_id':r['case_id']+'__origin__'+z,'case_id':r['case_id'],'protein_id':req['protein_id'],
                              'metal':z,'candidate':'origin',**ep,'mapping':r['maps'][z],
                              'source_preparation':r['representations']['context']['preparation'],'model':req['config']['model'],
                              'forces_requested':True}
                        tasks.append(task)
                        for medium in ('vacuum','alpb'):
                            ip=out/'scalar_inputs'/(task['task_id']+'__'+medium+'.inp')
                            ip.parent.mkdir(parents=True,exist_ok=True);ip.write_text(scalar.recipe(ep['charge'],medium,'fresh'))
                            scalar_tasks.append({**{k:task[k] for k in ('case_id','protein_id','metal','candidate','xyz','charge','multiplicity')},
                                'task_id':task['task_id']+'__'+medium,'medium':medium,'input':record(ip),'profile':scalar.PROFILE})
    if not rows:
        rows=[{'case_id':s['case_id'],'source':s['source'],'status':'unavailable','original_status':s['status'],
               'reason':req['reference_status'].get('reason') or 'one or more declared sources/discoveries unavailable'} for s in sources]
    status='prepared' if ready_ref and not missing else 'unsupported_reference' if not ready_ref else 'unavailable'
    result={'protocol_id':PROTOCOL,'request':record(request),'source_preparation':record(out/'SOURCE_PREPARATION.json'),
            'source_archive':record(source_preparation),'config':req['config'],'reference':req['reference'],
            'reference_status':req['reference_status'],'envelope_policy':envelope.POLICY,'union':group['union'] if group else None,
            'discovery':record(out/'DISCOVERY.json') if discovered else None,'cases':rows,'denominator':3,
            'supported_sources':sum(r['status']=='prepared' for r in rows),'missing_members':missing,'status':status,
            'physical_state_anchor':req['physical_state_anchor'],'origin_MACE_tasks':tasks,'origin_scalar_tasks':scalar_tasks,
            'prospective_calls':finite_calls(status=='prepared'),'scalar_profile':scalar.PROFILE,'optimizer_settings':precision.SETTINGS,
            'future_execution':{'status':'not_launched_requires_separate_integration','group_summary':'strict_median_of_all_three_sources',
                 'pool':'origin_adaptive_Ca_adaptive_La_shared_between_both_metals','all_required_cells_must_complete':True,
                 'selector':'motion_envelope_run.build_tasks_after_fresh_paired_origin_forces',
                 'molecular_cache_reuse':False,'archived_protonation_reused':True},
            'new_molecular_calls':0,'new_protonation':False,'production_changed':False,
            'preparation_wall_seconds':time.monotonic()-start,'implementation':{n:record(ROOT/'scripts'/n) for n in IMPLEMENTATIONS}}
    write_new(out/'PREPARATION.json',result);result=dry_run(out/'PREPARATION.json');write_new(out/'DRY_RUN.json',result);return result


def dry_run(preparation):
    p=read_json(preparation);req=check_request(read_json(verify(p['request'])))
    if (p['protocol_id']!=PROTOCOL or p['envelope_policy']!=envelope.POLICY or p['scalar_profile']!=scalar.PROFILE or
        p['optimizer_settings']!=precision.SETTINGS or p['new_molecular_calls'] or p['new_protonation'] or p['production_changed']):
        raise InvalidArtifact('frozen preparation-only scientific profile differs')
    for pin in p['implementation'].values():verify(pin)
    src=read_json(verify(p['source_preparation']));arch=read_json(verify(p['source_archive']))
    ids=[c['case_id'] for c in req['cases']]
    if (p['denominator']!=3 or p['reference']!=req['reference'] or p['reference_status']!=req['reference_status'] or
        p['physical_state_anchor']!=req['physical_state_anchor'] or [r['case_id'] for r in p['cases']]!=ids or
        [r['case_id'] for r in src['cases']]!=ids or p['config']!=req['config'] or src['config']!=req['config'] or arch['config']!=req['config']):
        raise InvalidArtifact('source membership/configuration/reference/anchor differs')
    for s,c in zip(src['cases'],req['cases']):
        exact=[r for r in arch['cases'] if original.matching_source(c,r.get('source',{}))]
        if len(exact)!=1 or s!=exact[0]:raise InvalidArtifact('archived source preparation changed')
    if p['supported_sources']!=sum(r['status']=='prepared' for r in p['cases']):raise InvalidArtifact('supported source count differs')
    if p['prospective_calls']!=finite_calls(p['status']=='prepared'):raise InvalidArtifact('finite prospective call count differs')
    if p['status']=='prepared':
        if p['reference_status']['status']!='available' or p['missing_members'] or p['supported_sources']!=3:
            raise InvalidArtifact('complete prepared triple and actual reference required')
        discovered=read_json(verify(p['discovery']));group=read_json(verify(p['union']))
        if discovered['policy']!=envelope.POLICY or [r['case_id'] for r in discovered['cases']]!=ids or any(r['status']!='discovered' for r in discovered['cases']):
            raise InvalidArtifact('exact complete source discovery required')
        keys=set().union(*({original.context.fragment_id(f) for f in r['fragments']} for r in discovered['cases']))
        if group['members']!=ids or group['fragments']!=[original.context.fragment_record(k) for k in sorted(keys)]:
            raise InvalidArtifact('envelope union differs from all three declared fragment lists')
        if len(p['origin_MACE_tasks'])!=6 or len(p['origin_scalar_tasks'])!=12:raise InvalidArtifact('origin task denominator differs')
        if len({r['state_key'] for r in p['cases']})!=1 or len({r['protein_key'] for r in p['cases']})!=1:
            raise InvalidArtifact('whole-three physical state/protein mismatch')
        for r in p['cases']:
            signature=read_json(verify(r['state_signature']));identity=read_json(verify(r['protein_identity']))
            if signature!=read_json(verify(r['state_anchor_signature'])) or cache_key(signature)!=r['state_key'] or cache_key(identity)!=r['protein_key']:
                raise InvalidArtifact('source atom/proton/cap/bond/water/protein state changed')
            maps={z:read_json(verify(r['maps'][z])) for z in ('Ca','La')}
            if maps['Ca']!=maps['La']:raise InvalidArtifact('paired physical map differs')
            kin=Kinematics(maps['Ca']['context']);q0=kin.evaluate(np.zeros(len(kin.modes)))[1]
            ep=r['representations']['context']['endpoints'];paired(verify(ep['La']['xyz']),verify(ep['Ca']['xyz']),ep['La']['charge'],ep['Ca']['charge'])
            for z in ('Ca','La'):
                t=next(x for x in p['origin_MACE_tasks'] if (x['case_id'],x['metal'])==(r['case_id'],z))
                if any(t[k]!=ep[z][k] for k in ('xyz','charge','multiplicity')) or t['mapping']!=r['maps'][z] or t['model']!=req['config']['model'] or t['forces_requested'] is not True:
                    raise InvalidArtifact('native origin task scientific identity differs')
                if not np.allclose(q0,[a[1:] for a in xyz(verify(ep[z]['xyz']))],rtol=0,atol=1e-12):raise InvalidArtifact('paired q0 mapping differs')
                for medium in ('vacuum','alpb'):
                    low=next(x for x in p['origin_scalar_tasks'] if (x['case_id'],x['metal'],x['medium'])==(r['case_id'],z,medium))
                    if any(low[k]!=ep[z][k] for k in ('xyz','charge','multiplicity')) or verify(low['input']).read_text()!=scalar.recipe(ep[z]['charge'],medium,'fresh'):
                        raise InvalidArtifact('strict origin scalar scientific identity differs')
    elif p['origin_MACE_tasks'] or p['origin_scalar_tasks']:
        raise InvalidArtifact('incomplete reference or triple cannot emit molecular tasks')
    return {'preparation':record(preparation),'protein_id':req['protein_id'],'status':p['status'],
            'source_ids':ids,'supported_sources':p['supported_sources'],'missing_members':p['missing_members'],
            'reference_status':p['reference_status'],'reference':p['reference'],'prospective_calls':p['prospective_calls'],
            'execution_implemented':False,'new_molecular_calls':0,'new_protonation':False}


def report(preparation,output):
    result=dry_run(preparation);p=read_json(preparation)
    result.update(preparation_wall_seconds=p['preparation_wall_seconds'],
                  atoms_by_source={r['case_id']:r.get('atom_count') for r in p['cases']},
                  cost_status='molecular_cost_unmeasured_no_execution; archived protonation cost excluded',
                  requested_execution={'warm_MACE':'one H200 /32 CPU /200000 MiB, existing warm worker',
                      'strict_scalar':'rank1, up to32 workers; CPU-only32 CPU/65536 MiB existing runner'},
                  eligibility='preparation-only; no new source score or affinity classification')
    write_new(output,result);return result


def main():
    ap=argparse.ArgumentParser(description=__doc__);sub=ap.add_subparsers(dest='op',required=True)
    for op,fields in {'request':('sources','config','reference','protein-id','agreement','output'),
        'prepare':('request','source-preparation','output'),'dry-run':('preparation',),'report':('preparation','output')}.items():
        q=sub.add_parser(op)
        for f in fields:q.add_argument('--'+f,required=True)
    args=vars(ap.parse_args());op=args.pop('op').replace('-','_');print(json.dumps(globals()[op](**args),indent=2))

if __name__=='__main__':main()
