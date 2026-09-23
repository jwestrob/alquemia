"""Prepare an explicit three-La-source PQQ union; no historical-label requirement."""
from __future__ import annotations
import argparse
import copy
import json
from pathlib import Path
import re
import time
import numpy as np

from affordable_common import InvalidArtifact,cache_key,paired,read_json,record,verify,write_new,xyz
import consistent_context as context
import pqq_fast_prepare as source_prepare
import pqq_standard as standard
from pqq_union_candidate import source_identity,matching_source
from pqq_union_execution import union_mapping
from adaptive_origin_recovery import snapshot
from mace_site_kinematics import Kinematics
import slsqp_precision as precision
import union_adaptive

PROTOCOL='Nikasha_explicit_three_La_source_union_preparation_v1'
REFERENCE_PROTOCOL='Nikasha_three_La_union_minimal_native_OMOL_GFN2_rank1_v1'
REFERENCE_ID='Nikasha_three_La_membership_ftol1e8_canonical25_v1'
REFERENCE_SHA='e725ea836d4a76bff8476cdedbc12be5d6a8bc5595da092d312f5dcec5c4ac9c'
SOURCE_FIELDS=('case_id','source_structure','metal','pqq','roles','assembly','normalization','source_conditioning_metal','raw_source_metal')
POLICY={'source_count':3,'conditioning':'La','state_anchor':'lexicographically_first_declared_source_id_before_energies',
        'union':'all_three_declared_local_fragment_lists','missing_member':'whole_group_unavailable_no_reduced_union',
        'paired_coordinate_tolerance_A':1e-12,'historical_calibration_member_required':False}


def configuration(path):
    d=read_json(path);return d.get('config',d)


def descriptors(data):
    keys=[k for k in ('sources','cases') if k in data]
    if len(keys)!=1:raise InvalidArtifact('one explicit sources or existing cases list required')
    return data[keys[0]]


def check_reference(reference,config):
    r=read_json(verify(reference))
    if (reference['sha256']!=REFERENCE_SHA or r['protocol_id']!=REFERENCE_PROTOCOL or r['reference_id']!=REFERENCE_ID or
        r['canonical_denominator']!=25 or r['crystals_or_stress_used_for_fit'] or
        r['optimizer_settings']!=precision.SETTINGS or r['settings']!=union_adaptive.POOL_SETTINGS or r['model']!=config['model']):
        raise InvalidArtifact('incompatible frozen three-source reference')
    if any(v['status']!='available' or v['available']!=25 for v in r['variants'].values()):
        raise InvalidArtifact('complete frozen canonical reference required')
    return r


def check_request(req):
    if req.get('protocol_id')!=PROTOCOL or req.get('policy')!=POLICY or req['production_changed']:
        raise InvalidArtifact('unsupported explicit-three-source policy')
    release=standard.release(verify(req['release']))
    if req['config']!=release['source_configuration']:raise InvalidArtifact('unqualified source preparation configuration')
    if configuration(verify(req['configuration_source']))!=req['config']:
        raise InvalidArtifact('explicit source configuration pin differs')
    rows=descriptors(read_json(verify(req['source_descriptors'])));declared=[]
    for row in rows:
        clean={k:copy.deepcopy(row[k]) for k in SOURCE_FIELDS}
        if isinstance(clean['source_structure'],str):clean['source_structure']=record(clean['source_structure'])
        declared.append(clean)
    if sorted(declared,key=lambda c:c['case_id'])!=req['cases']:
        raise InvalidArtifact('explicit source descriptor pin differs')
    check_reference(req['reference'],req['config']);verify(req['agreement'])
    cases=req['cases'];ids=[c['case_id'] for c in cases]
    if len(cases)!=3 or len(set(ids))!=3 or ids!=sorted(ids) or not all(re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*',c) for c in ids):
        raise InvalidArtifact('exactly three unique safe source IDs in deterministic order required')
    if req['physical_state_anchor']!=ids[0]:raise InvalidArtifact('predeclared physical-state anchor differs')
    if len({verify(c['source_structure']).read_bytes() for c in cases})!=3:
        raise InvalidArtifact('three distinct real source structures required')
    identities=[]
    for c in cases:
        if (set(c)!=set(SOURCE_FIELDS) or c['source_conditioning_metal']!='La' or
                c['normalization']!='protenix_generic_PQQ' or c['raw_source_metal']['element']!='La'):
            raise InvalidArtifact('explicit La-conditioned Protenix PQQ source/selector required')
        identities.append(source_identity(c,req['config']))
    if any(i!=identities[0] for i in identities[1:]):raise InvalidArtifact('source protein sequence/numbering/site roles/assembly differ')
    if req['source_identity']!=cache_key(identities[0]):raise InvalidArtifact('source identity digest differs')
    return req


def request(sources,config,reference,protein_id,agreement,output,release=standard.DEFAULT_RELEASE):
    data=read_json(sources);rows=descriptors(data);cfg=configuration(config)
    if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*',protein_id):raise InvalidArtifact('safe explicit protein ID required')
    cases=[]
    for row in rows:
        c={k:copy.deepcopy(row[k]) for k in SOURCE_FIELDS}
        if isinstance(c['source_structure'],str):c['source_structure']=record(c['source_structure'])
        cases.append(c)
    cases.sort(key=lambda c:c['case_id'])
    if len(cases)!=3:raise InvalidArtifact('exactly three sources required')
    req={'protocol_id':PROTOCOL,'policy':POLICY,'protein_id':protein_id,'cases':cases,'config':cfg,
         'source_identity':cache_key(source_identity(cases[0],cfg)),'physical_state_anchor':cases[0]['case_id'],
         'reference':record(reference),'release':record(release),'agreement':record(agreement),
         'source_descriptors':record(sources),'configuration_source':record(config),
         'source_evidence':{r['case_id']:{k:copy.deepcopy(r.get(k)) for k in ('expected_class','evidence_stratum','label_scope','source_provenance')} for r in rows},
         'anchor_is_biological_label_or_calibration_member':False,'production_changed':False,'new_molecular_calls':0}
    check_request(req);write_new(output,req)
    return {'request':record(output),'sources':3,'physical_state_anchor':req['physical_state_anchor'],'new_molecular_calls':0}


def prepare(request,output,source_mode='reuse-exact',source_preparation=None):
    start=time.monotonic();req=check_request(read_json(request));out=Path(output).resolve()
    if 'workspaces' not in out.parts:raise InvalidArtifact('candidate products belong under workspaces')
    if source_mode not in ('reuse-exact','fresh'):raise InvalidArtifact('explicit source preparation mode required')
    if (source_mode=='reuse-exact') != bool(source_preparation):raise InvalidArtifact('reuse requires an explicit source-preparation pin; fresh prohibits it')
    out.mkdir(parents=True,exist_ok=False);sources=[]
    if source_mode=='reuse-exact':
        old=read_json(source_preparation)
        if old['config']!=req['config']:raise InvalidArtifact('source preparation configuration differs')
        catalog=old['cases']
    for src in req['cases']:
        try:
            if source_mode=='reuse-exact':
                candidates=[c for c in catalog if matching_source(src,c.get('source',{}))]
                if len(candidates)!=1:raise InvalidArtifact('one exact source preparation required')
                row=copy.deepcopy(candidates[0]);row['source_preparation_reuse']=record(source_preparation)
            else:row=source_prepare.prepare_case(src,req['config'],out/'source'/src['case_id'])
        except Exception as exc:row={'case_id':src['case_id'],'source':src,'status':'unsupported','reason':str(exc)}
        sources.append(row)
    source_path=out/'SOURCE_PREPARATION.json'
    write_new(source_path,{'request':record(request),'config':req['config'],'cases':sources,'source_mode':source_mode,
                           'source_archive':record(source_preparation) if source_preparation else None,'new_molecular_calls':0})
    missing=[r['case_id'] for r in sources if r['status']!='prepared'];selection=None;rows=[];tasks=[]
    if missing:
        for row in sources:rows.append({'case_id':row['case_id'],'source':row['source'],'status':'group_source_unavailable',
            'source_status':row['status'],'reason':row.get('reason') or 'another declared source is unavailable','missing_members':missing})
    else:
        fragments=set().union(*({context.fragment_id(f) for f in read_json(verify(r['representations']['context']['preparation']))['added_fragments']} for r in sources))
        group={'root_case_id':req['protein_id'],'members':[c['case_id'] for c in req['cases']],
            'physical_state_anchor':req['physical_state_anchor'],'fragments':[context.fragment_record(f) for f in sorted(fragments)],
            'source_request':record(request),'selection_rule':POLICY['union'],'anchor_has_no_label_or_calibration_role':True}
        selection=out/'UNION.json';write_new(selection,group);group['union']=record(selection)
        rows=[context.prepare_one(row,group,req['config'],out/'union') for row in sources]
        anchor=next(r for r in rows if r['case_id']==req['physical_state_anchor'])
        for row in rows:
            if row['status']!='prepared_awaiting_group_check':continue
            row['physical_state_anchor']=req['physical_state_anchor']
            if ('state_key' not in anchor or row['state_key']!=anchor['state_key'] or row['protein_key']!=anchor['protein_key']):
                row.update(status='source_state_mismatch',reason='source atom/cofactor/H-parent/cap/bond/charge/water/protein state differs from declared anchor')
                continue
            row['state_anchor_signature']=anchor['state_signature']
            try:
                mapping=union_mapping(row,req['config']['topology'],out/'maps');row['physical_mapping']=mapping['mapping']
                kin=Kinematics(read_json(verify(mapping['mapping']))['context']);coord=kin.evaluate(np.zeros(len(kin.modes)))[1]
                pair=row['representations']['context']['endpoints']
                paired(verify(pair['La']['xyz']),verify(pair['Ca']['xyz']),pair['La']['charge'],pair['Ca']['charge'])
                if not all(np.allclose(coord,[a[1:] for a in xyz(verify(pair[z]['xyz']))],atol=1e-12,rtol=0) for z in ('Ca','La')):
                    raise InvalidArtifact('physical paired q0 mapping differs')
                row.update(status='prepared',physical_mode_ids=[m['id'] for m in kin.modes],
                           selected_four_modes_status='requires_fresh_paired_native_origin_forces')
            except Exception as exc:row.update(status='physical_mapping_unsupported',reason=str(exc))
        missing=[r['case_id'] for r in rows if r['status']!='prepared']
        if not missing:
            for row in rows:
                for z in ('Ca','La'):
                    ep=row['representations']['context']['endpoints'][z]
                    tasks.append({'task_id':row['case_id']+'__'+z,'case_id':row['case_id'],'metal':z,**ep,
                        'mapping':row['physical_mapping'],'source_preparation':row['representations']['context']['preparation']})
    status='prepared' if not missing else 'unavailable'
    data={'protocol_id':PROTOCOL,'request':record(request),'config':req['config'],'source_manifest':record(request),
          'source_preparation':record(source_path),'union':record(selection) if selection else None,'cases':rows,
          'status':status,'denominator':3,'supported_sources':sum(r['status']=='prepared' for r in rows),
          'missing_members':missing,'physical_state_anchor':req['physical_state_anchor'],'reference':req['reference'],
          'tasks':tasks,'future_execution':{'status':'deferred_pending_threefold_transfer_qualification','launched':False,
              'maximum_origin_MACE':6 if not missing else 0,'maximum_searches':6 if not missing else 0,
              'maximum_cross_MACE':6 if not missing else 0,'maximum_GFN2':36 if not missing else 0},
          'new_molecular_calls':0,'new_protonation':source_mode=='fresh','wall_seconds':time.monotonic()-start,
          'production_changed':False,'implementation':snapshot(Path(__file__).parent,out/'implementation')}
    path=out/'PREPARATION.json';write_new(path,data);return dry_run(path)


def dry_run(preparation):
    p=read_json(preparation);req=check_request(read_json(verify(p['request'])))
    if p['protocol_id']!=PROTOCOL or p['new_molecular_calls'] or p['production_changed'] or p['future_execution']['launched']:
        raise InvalidArtifact('preparation-only scope differs')
    for v in p['implementation'].values():verify(v)
    if p['reference']!=req['reference'] or p['physical_state_anchor']!=req['physical_state_anchor'] or [r['case_id'] for r in p['cases']]!=[c['case_id'] for c in req['cases']]:
        raise InvalidArtifact('prepared request membership/reference/anchor differs')
    if p['denominator']!=3 or p['supported_sources']!=sum(r['status']=='prepared' for r in p['cases']):
        raise InvalidArtifact('prepared source denominator differs')
    if p['status']=='prepared':
        src=read_json(verify(p['source_preparation']))
        if src['request']!=p['request'] or len(src['cases'])!=3 or any(r['status']!='prepared' for r in src['cases']):
            raise InvalidArtifact('complete exact local source preparation required')
        selection=read_json(verify(p['union']))
        if selection['members']!=[c['case_id'] for c in req['cases']]:raise InvalidArtifact('fixed three-member selection differs')
        expected=set().union(*({context.fragment_id(f) for f in read_json(verify(r['representations']['context']['preparation']))['added_fragments']} for r in src['cases']))
        if selection['fragments']!=[context.fragment_record(f) for f in sorted(expected)]:
            raise InvalidArtifact('union is not exactly the three source fragment lists')
        if p['missing_members'] or len(p['tasks'])!=6:raise InvalidArtifact('complete paired triple required')
        if len({r['protein_key'] for r in p['cases']})!=1:
            raise InvalidArtifact('common prepared protein identity differs')
        for row in p['cases']:
            own=read_json(verify(row['state_signature']));anchor=read_json(verify(row['state_anchor_signature']))
            if own!=anchor or cache_key(own)!=row['state_key']:raise InvalidArtifact('prepared common state differs')
            actual_protein=read_json(verify(row['protein_identity']))
            if cache_key(actual_protein)!=row['protein_key']:raise InvalidArtifact('prepared protein identity differs')
            pair=row['representations']['context']['endpoints']
            paired(verify(pair['La']['xyz']),verify(pair['Ca']['xyz']),pair['La']['charge'],pair['Ca']['charge'])
            kin=Kinematics(read_json(verify(row['physical_mapping']))['context']);coord=kin.evaluate(np.zeros(len(kin.modes)))[1]
            for z in ('Ca','La'):
                ep=row['representations']['context']['endpoints'][z];t=next(t for t in p['tasks'] if t['case_id']==row['case_id'] and t['metal']==z)
                if (t['xyz'],t['charge'],t['multiplicity'],t['mapping'])!=(ep['xyz'],ep['charge'],ep['multiplicity'],row['physical_mapping']):
                    raise InvalidArtifact('prepared endpoint/task identity differs')
                if not np.allclose(coord,[a[1:] for a in xyz(verify(ep['xyz']))],atol=1e-12,rtol=0):raise InvalidArtifact('q0 physical map differs')
    elif p['tasks']:raise InvalidArtifact('incomplete triple cannot emit executable tasks')
    return {'preparation':record(preparation),'status':p['status'],'sources':3,'supported_sources':p['supported_sources'],
            'missing_members':p['missing_members'],'physical_state_anchor':p['physical_state_anchor'],
            'historical_canonical_member_required':False,'reference':p['reference'],'future_execution':p['future_execution'],'new_molecular_calls':0}


def main():
    ap=argparse.ArgumentParser(description=__doc__);sub=ap.add_subparsers(dest='op',required=True)
    q=sub.add_parser('request')
    for f in ('sources','config','reference','protein-id','agreement','output'):q.add_argument('--'+f,required=True)
    q=sub.add_parser('prepare');q.add_argument('--request',required=True);q.add_argument('--output',required=True)
    q.add_argument('--source-mode',choices=('reuse-exact','fresh'),required=True);q.add_argument('--source-preparation')
    q=sub.add_parser('dry-run');q.add_argument('--preparation',required=True)
    args=vars(ap.parse_args());op=args.pop('op').replace('-','_');print(json.dumps(globals()[op](**args),indent=2))

if __name__=='__main__':main()
