"""Deterministic canonical geometry preparation under triple-specific membership."""
from __future__ import annotations
import argparse
import copy
from concurrent.futures import ProcessPoolExecutor,as_completed
import json
from pathlib import Path
import time
from affordable_common import InvalidArtifact,read_json,record,verify,write_new
import consistent_context as context
import union_triple_preparation as triples

ROOT=Path(__file__).resolve().parents[1]
PLAN=ROOT/'diagnostics/union_triple_preparation_20260923/CANONICAL_PLAN.md'
PROTOCOL='Nikasha_canonical25_three_La_source_fragment_union_preparation_v1'


def prepare(primary,crystals,output,agreement=PLAN,workers=4):
    start=time.monotonic();p=read_json(primary);m=triples.validate_selection(verify(p['selection']))
    source=read_json(verify(m['source_preparation']));si={r['case_id']:r for r in source['cases']}
    full=read_json(verify(m['tenfold_preparation']));old={r['case_id']:r for r in full['cases']}
    groups={g['selection_id']:g for g in m['groups']};prepared={(r['selection_id'],r['case_id']):r for r in p['cases']}
    out=Path(output).resolve()
    if 'workspaces' not in out.parts:raise InvalidArtifact('candidate products belong under workspaces/')
    out.mkdir(parents=True,exist_ok=False);chosen=[]
    for gid in sorted({t['protein_id'] for t in m['triples']}):
        options=sorted((t for t in m['triples'] if t['protein_id']==gid),key=lambda t:tuple(sorted(t['members'])))
        complete=[t for t in options if not t['missing_members']]
        if not complete:raise InvalidArtifact('no complete triple for designated canonical source: '+gid)
        selected=complete[0];canon=[r for r in source['cases'] if r['source']['root_case_id']==gid and r['source']['canonical_coordinate_match']]
        if len(canon)!=1:raise InvalidArtifact('actual canonical membership unavailable')
        skipped=[{'triple_id':t['triple_id'],'members':t['members'],'missing_members':t['missing_members']}
            for t in options[:options.index(selected)] if t['missing_members']]
        chosen.append({'protein_id':gid,'case_id':canon[0]['case_id'],'selected_triple':selected,
            'skipped_incomplete_before_selected':skipped,'selection_id':selected['selection_id'],
            'source':canon[0]['source'],'expected_class':canon[0]['expected_class']})
    if len(chosen)!=25:raise InvalidArtifact('all25 canonical sources required')
    crystal=read_json(crystals)
    if [c['case_id'] for c in crystal['cases']]!=['1H4I','4MAE','1KB0'] or crystal['config']!=m['config']:
        raise InvalidArtifact('existing three crystal selectors required')
    selection={'protocol_id':PROTOCOL,'primary_preparation':record(primary),'primary_selection':p['selection'],
        'agreement':record(agreement),'canonical_sources':chosen,'crystals':record(crystals),'denominator':28,
        'calibration_denominator':25,'crystal_transfer_denominator':3,'new_molecular_calls':0,
        'prospective_numerical_profile':'native_OMOL_four_mode_SLSQP200_ftol1e-8_Hartree_native_GFN2_1rank_MaxIter500',
        'implementation':record(__file__),'rule':'lexicographically_first_originally_complete_noncanonical_La_triple'}
    write_new(out/'SELECTION.json',selection)
    results={}
    with ProcessPoolExecutor(max_workers=workers) as pool:
        pending={pool.submit(context.prepare_one,si[x['case_id']],groups[x['selection_id']],m['config'],out/'prepared'):x for x in chosen}
        for f in as_completed(pending):
            x=pending[f];r=f.result();g=groups[x['selection_id']];anchor=prepared[x['selection_id'],g['members'][0]]
            r.update(selection_id=x['selection_id'],protein_id=x['protein_id'],is_separate_stress_probe=False,
                role='designated_canonical_calibration_preparation',expected_class=x['expected_class'])
            if r['status']=='prepared_awaiting_group_check':
                if anchor['status']!='prepared' or (r['state_key'],r['protein_key'])!=(anchor['state_key'],anchor['protein_key']):
                    r.update(status='selected_triple_state_mismatch',reason='canonical graph/state differs from selected triple')
                else:
                    r['status']='prepared';prev=old[r['case_id']]
                    same=r['state_key']==prev['state_key'] and all(context.reusable_state(r['representations']['context']['endpoints'][z],prev['representations']['context']['endpoints'][z]) for z in ('Ca','La'))
                    r['tenfold_comparison']={'exact_reuse_eligible':same,'tenfold_state_signature':prev['state_signature'],
                        'atom_count':prev['atom_count'],'atom_count_difference':r['atom_count']-prev['atom_count'],
                        'La_charge':prev['new_La_charge'],'cap_count':prev['cap_count']}
            results[r['case_id']]=r;print(json.dumps({'case_id':r['case_id'],'status':r['status']}),flush=True)
    rows=[results[x['case_id']] for x in chosen]
    for c in crystal['cases']:
        r=copy.deepcopy(c);r.update(selection_id='singleton_'+r['case_id'],protein_id=r['case_id'],
            role='consumed_crystal_transfer_preparation',is_separate_stress_probe=False,
            tenfold_comparison={'exact_reuse_eligible':True,'tenfold_state_signature':c['state_signature'],
                'atom_count':c['atom_count'],'atom_count_difference':0,'La_charge':c['new_La_charge'],'cap_count':c['cap_count']})
        for z in ('Ca','La'):verify(r['representations']['context']['endpoints'][z]['xyz'])
        rows.append(r)
    result={'protocol_id':PROTOCOL,'selection':p['selection'],'canonical_selection':record(out/'SELECTION.json'),
        'cases':rows,'denominator':28,'calibration_denominator':25,'crystal_transfer_denominator':3,
        'supported':sum(r['status']=='prepared' for r in rows),'new_molecular_calls':0,'new_reference':None,
        'wall_seconds':time.monotonic()-start,'workers':workers,'production_changed':False}
    write_new(out/'PREPARATION.json',result);return {k:v for k,v in result.items() if k!='cases'}


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--primary',required=True);p.add_argument('--crystals',required=True)
    p.add_argument('--output',required=True);p.add_argument('--agreement',default=str(PLAN));p.add_argument('--workers',type=int,default=4)
    print(json.dumps(prepare(**vars(p.parse_args())),indent=2))
