"""Frozen union/adaptive225 source partition and comparison; no new scorer."""
from __future__ import annotations
import argparse
from collections import Counter
from datetime import datetime,timezone
import json
from pathlib import Path
import statistics
from affordable_common import InvalidArtifact,read_json,record,verify,write_new
import union_adaptive as union
from accommodation_folds_compare import decision
from accommodation_fold_proposals import outcome,summarize_rows


def declared_rows(transfer,sources):
    result=read_json(transfer);src=read_json(sources)
    primary=[s for s in src['cases'] if s['primary_evaluation_pool']]
    idx={r['case_id']:r for r in result['rows']}
    if len(primary)!=225 or set(idx)!={s['case_id'] for s in primary}:raise InvalidArtifact('original225 population differs')
    rows=[]
    for s in primary:
        old=idx[s['case_id']]
        if old['expected_class']!=s['expected_class'] or old['canonical_coordinate_match']:raise InvalidArtifact('transfer label/canonical exclusion differs')
        rows.append({'case_id':s['case_id'],'actual_union_case_id':s['case_id'],'known_class':s['expected_class'],
            'biological_group':s['biological_group'],'root_case_id':s['root_case_id'],
            'source_conditioning_metal':s['source_conditioning_metal'],'role':'consumed_noncanonical_transfer',
            'union_origin_status':old['status'],'reason':old.get('reason')})
    return rows


def selection(transfer,sources,calibration,reference,pilot,agreement,output):
    rows=declared_rows(transfer,sources);p=read_json(pilot)
    reuse=[r['case_id'] for r in p['cases'] if r['case_id'] not in union.CANONICAL_REUSE]
    ready=[r for r in rows if r['union_origin_status']=='complete'];new=[r for r in ready if r['case_id'] not in reuse]
    if len(ready)!=208 or len(reuse)!=4 or len(new)!=204:raise InvalidArtifact('declared208/4/204 source count differs')
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False)
    top={'population':'union_adaptive_primary225','transfer':record(transfer),'sources':record(sources),
        'calibration':record(calibration),'reference':record(reference),'pilot_collection':record(pilot),
        'agreement':record(agreement),'all_cases':rows,'reuse_case_ids':reuse,'shard_count':4,'new_sources':204,
        'total_denominator':225,'origin_available':208,'preparation_exclusions':17,'new_searches':408,
        'maximum_cross_MACE_calls':408,'maximum_GFN2_calls':1632,'new_q0_calls':0,'new_DFT_calls':0,
        'selection_UTC':datetime.now(timezone.utc).isoformat()}
    path=out/'SELECTION.json';write_new(path,top);shards=[]
    for i in range(4):
        m={'population':'transfer225_shard','selection':record(path),'shard':i,'shard_count':4,
            'reference':record(reference),'cases':new[i::4]}
        q=out/('INPUTS_shard_'+str(i)+'.json');write_new(q,m);validate_selection(q,calibration);shards.append(record(q))
    ready={'selection':record(path),'shards':shards,'new_sources':204,'sources_per_shard':51,'new_searches':408,'molecular_calls':0}
    write_new(out/'READY.json',ready);return ready


def validate_selection(inputs,calibration):
    data=read_json(inputs);top=read_json(verify(data['selection']))
    if data['population']!='transfer225_shard' or data['shard_count']!=4 or data['shard'] not in range(4):raise InvalidArtifact('undeclared shard')
    if top['calibration']!=record(calibration) or data['reference']!=top['reference']:raise InvalidArtifact('frozen origin/reference changed')
    verify(top['agreement']);ref=read_json(verify(top['reference']))
    if (ref['protocol_id']!=union.POOL_PROTOCOL or ref['calibration_denominator']!=25 or ref['crystals_or_noncanonical_used_for_fit'] or
        any(r['status']!='available' for r in ref['variants'].values())):raise InvalidArtifact('complete frozen canonical reference required')
    rows=declared_rows(verify(top['transfer']),verify(top['sources']))
    if rows!=top['all_cases']:raise InvalidArtifact('original source roles/membership changed')
    pilot=read_json(verify(top['pilot_collection']));pm=read_json(verify(pilot['manifest']))
    if pm['protocol_id']!=union.POOL_PROTOCOL or pm['settings']!=union.POOL_SETTINGS or pm['model']!=ref['model']:raise InvalidArtifact('reuse method differs')
    expected_reuse=[c['case_id'] for c in pilot['cases'] if c['case_id'] not in union.CANONICAL_REUSE]
    if top['reuse_case_ids']!=expected_reuse or len(expected_reuse)!=4:raise InvalidArtifact('four prior fold reuse differs')
    prior={r['case_id']:r for r in read_json(verify(top['transfer']))['rows']}
    for cid in expected_reuse:
        c=next(c for c in pilot['cases'] if c['case_id']==cid)
        if c['source']['origin_row']!=prior[cid] or c['pool']['status']!='available':raise InvalidArtifact('exact prior source pool unavailable')
    new=[r for r in rows if r['union_origin_status']=='complete' and r['case_id'] not in expected_reuse]
    if len(new)!=204 or data['cases']!=new[data['shard']::4] or len(data['cases'])!=51:raise InvalidArtifact('finite disjoint shard differs')
    if datetime.fromisoformat(ref['frozen_UTC'])>=datetime.fromisoformat(top['selection_UTC']):raise InvalidArtifact('transfer selection preceded reference freeze')
    return {'status':'validated','new_sources':len(data['cases']),'total_denominator':225,'preparation_exclusions':17,'molecular_calls':0}


def main():
    p=argparse.ArgumentParser(description=__doc__);s=p.add_subparsers(dest='operation',required=True)
    q=s.add_parser('selection')
    for k in ('transfer','sources','calibration','reference','pilot','agreement','output'):q.add_argument('--'+k,required=True)
    q=s.add_parser('validate_selection')
    for k in ('inputs','calibration'):q.add_argument('--'+k,required=True)
    args=vars(p.parse_args());op=args.pop('operation');print(json.dumps(globals()[op](**args),indent=2))
if __name__=='__main__':main()
