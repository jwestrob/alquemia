"""Frozen threefold-union transfer comparison; no fitting or molecular calls."""
import argparse
from collections import Counter
import json
from pathlib import Path
import statistics
from affordable_common import InvalidArtifact,read_json,record,verify,write_new
from accommodation_folds_compare import decision
from accommodation_fold_proposals import outcome
from accommodation_nonlinear import contrast_components,relative_components
import nikasha_pool as pool
import union_triple_transfer as transfer
import union_triple_transfer_run as runner

METHODS=('context_composite','union_adaptive','union_precision','triple_union','triple_union_mathematical')


def tagged(value,bands,label):
    call=decision(value,bands)
    return {'R':value,'decision':call,'outcome':outcome(call,label)}


def counts(rows):return {m:dict(Counter(r['methods'][m]['outcome'] for r in rows)) for m in METHODS}


def compare(audit,collections,output):
    a=read_json(audit);ref=read_json(verify(a['reference']));prior=read_json(verify(a['precision_comparison']))
    if a['protocol_id']!=transfer.PROTOCOL or len(collections)!=2:raise InvalidArtifact('fixed all100 scope/two collections required')
    bands={m:prior['bands'][m] for m in METHODS[:3]}
    bands.update(triple_union=ref['variants']['operational']['bands'],triple_union_mathematical=ref['variants']['mathematical']['bands'])
    old={r['case_id']:r for r in prior['rows']};fresh={};sources={};shards=set();unavailable={}
    for path in collections:
        d=read_json(path);pm=read_json(verify(d['manifest']));sm=read_json(verify(pm['source_manifest']));runner.validate_pool(verify(d['manifest']))
        if sm['audit']!=record(audit) or sm['reference']!=a['reference'] or sm['shard'] in shards:raise InvalidArtifact('collection scope/reference/shard differs')
        shards.add(sm['shard'])
        for c in d['cases']:
            if c['case_id'] in fresh:raise InvalidArtifact('duplicate fresh source-membership pair')
            fresh[c['case_id']]=c;sources[c['case_id']]=record(path)
        for c in sm['unavailable']:unavailable[c['case_id']]=c
    if shards!={0,1} or set(fresh)|set(unavailable)!={r['pair_id'] for r in a['rows'] if not r['pool_reuse']}:
        raise InvalidArtifact('terminal55 new-pair coverage differs')
    rows=[];by_pair={}
    for r in a['rows']:
        pid=r['pair_id'];label=r['source']['source']['expected_class'];saved=r['pool_reuse']
        c=saved['pool'] if saved else fresh.get(pid);methods={m:old[r['case_id']]['methods'][m] for m in METHODS[:3]}
        if c and c['status']=='prepared' and pool.choose_rows(c['matrix'],[x['id'] for x in c['candidates']])!=c['pool']:
            raise InvalidArtifact('actual pool recomputation differs')
        if c and c['status']!='prepared' and c['pool']['status']=='available':raise InvalidArtifact('failed required proposal gained an origin-only result')
        available=bool(c and c['pool']['status']=='available')
        for name,variant in [('triple_union','operational'),('triple_union_mathematical','mathematical')]:
            value=c['pool'][variant]['composite_R_model_kcal_mol'] if available else None
            methods[name]=tagged(value,bands[name],label)
        works=None;static=None
        if available:
            static=contrast_components(c['matrix']['Ca']['origin']['components'],c['matrix']['La']['origin']['components'])
            works={z:relative_components(c['matrix'][z][c['pool']['rows'][z]['operational_candidate']]['components'],c['matrix'][z]['origin']['components']) for z in ('Ca','La')}
            delta=works['Ca']['composite_kcal_mol']-works['La']['composite_kcal_mol']
            if abs(delta-(methods['triple_union']['R']-static['composite_R_model_kcal_mol']))>1e-6:raise InvalidArtifact('selected work sign/unit differs')
        row={'pair_id':pid,'case_id':r['case_id'],'selection_id':r['selection_id'],'root_case_id':r['source']['protein_id'],
            'expected_class':label,'source_conditioning_metal':'La','status':'available' if available else 'unavailable',
            'methods':methods,'pool_reused':bool(saved),'source_collection':saved['collection'] if saved else sources.get(pid),
            'pool':c['pool'] if c else None,'matrix':c['matrix'] if c else None,'reason':c.get('reason') if c else unavailable.get(pid),
            'static_components':static,'selected_work':works,'atom_count':r['source']['atom_count'],
            'tenfold_atom_count':r['source']['tenfold_comparison']['atom_count'],'La_charge':r['source']['new_La_charge']}
        rows.append(row);by_pair[r['selection_id'],r['case_id']]=row
    triples=[]
    for t in a['triples']:
        label=old[t['members'][0]]['expected_class'];methods={}
        if any(old[cid]['expected_class']!=label for cid in t['members']):raise InvalidArtifact('source group label differs')
        for m in METHODS:
            cells=[old[cid]['methods'][m] if m in METHODS[:3] else by_pair.get((t['selection_id'],cid),{}).get('methods',{}).get(m,{'R':None}) for cid in t['members']]
            missing=[cid for cid,x in zip(t['members'],cells) if x['R'] is None]
            if m.startswith('triple_union') and t['status']!='prepared':missing=list(t['members'])
            values=[x['R'] for x in cells if x['R'] is not None];v=statistics.median(values) if not missing else None
            methods[m]={**tagged(v,bands[m],label),'missing_members':missing,'required_members':3,
                'range':max(values)-min(values) if not missing else None,'member_R':[x['R'] for x in cells]}
        triples.append({**t,'expected_class':label,'methods':methods})
    matched={}
    for m in METHODS[:3]:
        rr=[r for r in rows if r['methods'][m]['R'] is not None and r['methods']['triple_union']['R'] is not None]
        tt=[r for r in triples if r['methods'][m]['R'] is not None and r['methods']['triple_union']['R'] is not None]
        matched[m]={'pairs':len(rr),'pair_counts':counts(rr),'triples':len(tt),'triple_counts':counts(tt)}
    result={'protocol_id':transfer.PROTOCOL,'audit':record(audit),'reference':a['reference'],'prior':a['precision_comparison'],
        'collections':[record(p) for p in collections],'bands':bands,'rows':rows,'triples':triples,
        'pair_denominator':125,'triple_denominator':100,'counts':{'pairs':counts(rows),'triples':counts(triples)},'matched':matched,
        'new_thresholds_fitted':False,'new_molecular_calls':0,'production_changed':False,
        'interpretation':'100 overlapping developmental triples from25 consumed groups, not100 independent proteins; strict three-member median',
        'implementation':record(__file__)}
    write_new(output,result);return {k:v for k,v in result.items() if k not in ('rows','triples')}


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('audit','collections','output'):p.add_argument('--'+k,required=True,**({'nargs':'+'} if k=='collections' else {}))
    print(json.dumps(compare(**vars(p.parse_args())),indent=2))
