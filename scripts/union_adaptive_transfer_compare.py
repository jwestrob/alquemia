"""Join fixed union/adaptive225 results to existing frozen-reference ledgers."""
from __future__ import annotations
import argparse
from collections import Counter
from pathlib import Path
import json
from affordable_common import InvalidArtifact,read_json,record,verify,write_new
from accommodation_folds_compare import decision
from accommodation_fold_proposals import outcome,strict_summary,summarize_rows
from accommodation_nonlinear import relative_components
import nikasha_pool as pool
import union_adaptive as union
from union_adaptive_transfer import validate_selection


def key(row,kind):
    return row['case_id'] if kind=='rows' else (row['root_case_id'],row.get('pool'),tuple(row['members']))


def add_ledger(base,other,kind,names):
    idx={key(r,kind):r for r in other}
    if len(idx)!=len(other) or set(idx)!={key(r,kind) for r in base}:raise InvalidArtifact('ledger source/group membership differs')
    result=[]
    for r in base:
        q=idx[key(r,kind)]
        if r['root_case_id']!=q['root_case_id'] or r['expected_class']!=q['expected_class']:raise InvalidArtifact('ledger biological label differs')
        if kind=='rows' and r['source_conditioning_metal']!=q['source_conditioning_metal']:raise InvalidArtifact('conditioning differs')
        result.append({**r,'methods':{**r['methods'],**{new:q['methods'][old] for old,new in names.items()}}})
    return result


def compare(selection,collections,baseline,union_static,native,standalone,output):
    top=read_json(selection);ref=read_json(verify(top['reference']));bands_new=ref['variants']['operational']['bands']
    if not bands_new:raise InvalidArtifact('frozen canonical reference unavailable')
    old=read_json(baseline);us=read_json(union_static);nr=read_json(native);sr=read_json(standalone)
    if len(collections)!=4:raise InvalidArtifact('all four terminal shard collections required')
    available={};shards=set();pins=[];optimizer={}
    for path in collections:
        c=read_json(path);m=read_json(verify(c['manifest']));sm=read_json(verify(m['source_manifest']));inputs=read_json(verify(sm['inputs']))
        validate_selection(verify(sm['inputs']),verify(sm['calibration']))
        if inputs['selection']!=record(selection) or inputs['shard'] in shards:raise InvalidArtifact('different or duplicate frozen shard')
        if m['reference']!=top['reference'] or m['model']!=ref['model'] or m['settings']!=ref['settings'] or sm['settings']!=ref['optimizer_settings']:raise InvalidArtifact('frozen method/reference changed')
        if c['protocol_id']!=union.POOL_PROTOCOL or c['denominator']!=51 or [q['case_id'] for q in c['cases']]!=m['declared_case_ids']:raise InvalidArtifact('terminal shard denominator differs')
        shards.add(inputs['shard']);pins.append(record(path))
        proposals=read_json(verify(m['source']))
        for e in proposals['endpoints']:
            r=read_json(verify(e['proposal_receipt'])) if e['proposal_receipt'] else None
            optimizer[(e['case_id'],e['metal'])]={'status':e['status'],'reason':e.get('reason'),'boundary_flag':e.get('boundary_flag'),
                'receipt':e['proposal_receipt'],'optimizer':r.get('optimizer') if r else None,'wall_seconds':r.get('wall_seconds') if r else None}
        for q in c['cases']:
            if q['case_id'] in available:raise InvalidArtifact('duplicate physical source across shards')
            available[q['case_id']]={'case':q,'collection':record(path)}
    if len(available)!=204:raise InvalidArtifact('missing requested fresh source')
    pilot=read_json(verify(top['pilot_collection']))
    pm=read_json(verify(pilot['manifest']));pr=read_json(verify(pm['source']))
    for cid in top['reuse_case_ids']:
        c=next(c for c in pilot['cases'] if c['case_id']==cid);available[cid]={'case':c,'collection':top['pilot_collection']}
        for e in pr['endpoints']:
            if e['case_id']!=cid:continue
            r=read_json(verify(e['proposal_receipt']))
            optimizer[(cid,e['metal'])]={'status':e['status'],'reason':e.get('reason'),'boundary_flag':e.get('boundary_flag'),
                'receipt':e['proposal_receipt'],'optimizer':r['optimizer'],'wall_seconds':r['wall_seconds']}
    if len(available)!=208:raise InvalidArtifact('prepared source union changed')
    methods_new=('union_adaptive','union_adaptive_mathematical')
    bands={**old['bands'],'context_union':read_json(verify(us['reference']))['bands'],
        'native_minimal_recovered':nr['bands']['minimal_recovered'],'full_adaptive':nr['bands']['full_adaptive'],
        'standalone_static':sr['bands']['standalone_static'],'union_adaptive':bands_new,
        'union_adaptive_mathematical':ref['variants']['mathematical']['bands']}
    joined={}
    for kind,count in (('rows',225),('pools',75),('triples',100)):
        if len(old[kind])!=count:raise InvalidArtifact('original strict denominator differs')
        rows=add_ledger(old[kind],us[kind],kind,{'context_union':'context_union'})
        rows=add_ledger(rows,nr[kind],kind,{'minimal_recovered':'native_minimal_recovered','full_adaptive':'full_adaptive'})
        joined[kind]=add_ledger(rows,sr[kind],kind,{'standalone_static':'standalone_static'})
    expected={x['case_id']:x for x in top['all_cases']};rows=[]
    origins={x['case_id']:x for x in read_json(verify(top['transfer']))['rows']}
    for r in joined['rows']:
        src=expected[r['case_id']];entry=available.get(r['case_id']);c=entry['case'] if entry else None
        if r['expected_class']!=src['known_class'] or r['root_case_id']!=src['root_case_id']:raise InvalidArtifact('frozen biological label differs')
        chosen=c['pool'] if c else None;status=chosen['status'] if chosen else src['union_origin_status'];work={};new={}
        if c and c['status']=='prepared' and pool.choose_rows(c['matrix'],[q['id'] for q in c['candidates']])!=chosen:raise InvalidArtifact('parsed pool/component minima differ')
        for variant,name in (('operational','union_adaptive'),('mathematical','union_adaptive_mathematical')):
            value=chosen[variant]['composite_R_model_kcal_mol'] if chosen and chosen['status']=='available' else None
            call=decision(value,bands[name]);new[name]={'R':value,'decision':call,'outcome':outcome(call,r['expected_class'])}
            if value is not None:
                work[variant]={z:relative_components(c['matrix'][z][chosen['rows'][z][variant+'_candidate']]['components'],c['matrix'][z]['origin']['components']) for z in ('Ca','La')}
                raw0=origins[r['case_id']]['composite_R_model_kcal_mol']
                if abs((value-raw0)-(work[variant]['Ca']['composite_kcal_mol']-work[variant]['La']['composite_kcal_mol']))>1e-6:raise InvalidArtifact('work/sign/unit accounting differs')
        rows.append({**r,'methods':{**r['methods'],**new},'union_adaptive_status':status,
            'failure_reason':chosen.get('reason') if chosen else src['reason'],'union_selected_work':work or None,
            'union_pool_components':{v:chosen[v] for v in ('mathematical','operational')} if chosen else None,
            'selected_candidates':{v:{z:chosen['rows'][z][v+'_candidate'] for z in ('Ca','La')} for v in ('mathematical','operational')} if chosen and chosen['status']=='available' else None,
            'optimizers':{z:optimizer.get((r['case_id'],z)) for z in ('Ca','La')},'source_collection':entry['collection'] if entry else top['transfer'],
            'reused_common8':r['case_id'] in top['reuse_case_ids']})
    index={r['case_id']:r for r in rows};pools=[]
    for r in joined['pools']:
        new={}
        for name in methods_new:
            if r['pool'] in ('La4','Ca5'):value=strict_summary(r['members'],index,name,bands[name],r['expected_class'])
            elif r['pool']=='balanced':
                arms=[next(p for p in pools if p['root_case_id']==r['root_case_id'] and p['pool']==arm)['methods'][name] for arm in ('La4','Ca5')]
                number=None if any(a['R'] is None for a in arms) else sum(a['R'] for a in arms)/2
                call=decision(number,bands[name]);value={'R':number,'decision':call,'outcome':outcome(call,r['expected_class']),
                    'required_members':len(r['members']),'missing_members':sorted({x for a in arms for x in a['missing_members']}),
                    'within_protein_range':None,'equal_mean_of_complete_arm_medians':True}
            else:raise InvalidArtifact('undeclared fold aggregate')
            new[name]=value
        pools.append({**r,'methods':{**r['methods'],**new}})
    triples=[{**r,'methods':{**r['methods'],**{name:strict_summary(r['members'],index,name,bands[name],r['expected_class']) for name in methods_new}}} for r in joined['triples']]
    methods=list(rows[0]['methods']);subsets={'all225':rows,'La100':[r for r in rows if r['source_conditioning_metal']=='La'],
        'Ca125':[r for r in rows if r['source_conditioning_metal']=='Ca'],**{name:[r for r in pools if r['pool']==name] for name in ('La4','Ca5','balanced')},'La100_triples':triples}
    counts={name:summarize_rows(rr,methods) for name,rr in subsets.items()};matched={}
    for name,rr in subsets.items():
        matched[name]={}
        for other in methods:
            if other in methods_new:continue
            common=[r for r in rr if all(r['methods'][x]['R'] is not None for x in ('union_adaptive',other))]
            matched[name][other]={'declared':len(rr),'common':len(common),'counts':summarize_rows(common,('union_adaptive',other)),
                'transitions':dict(Counter(r['methods'][other]['outcome']+'->'+r['methods']['union_adaptive']['outcome'] for r in common))}
    result={'protocol_id':union.POOL_PROTOCOL,'selection':record(selection),'reference':top['reference'],'collections':pins,
        'comparators':{k:record(v) for k,v in dict(baseline=baseline,union_static=union_static,native=native,standalone=standalone).items()},
        'rows':rows,'pools':pools,'triples':triples,'bands':bands,'counts':counts,'matched':matched,
        'new_thresholds_fitted':False,'new_molecular_calls':0,'production_changed':False,'implementation':record(__file__),
        'interpretation':'Consumed structures; all225 plus strict25arm/balanced groups and100correlated triples; no new biological validation.'}
    write_new(output,result);return counts['all225']

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('selection','baseline','union_static','native','standalone','output'):p.add_argument('--'+name,required=True)
    p.add_argument('--collections',nargs='+',required=True)
    print(json.dumps(compare(**vars(p.parse_args())),indent=2))
