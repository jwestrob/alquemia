"""Exact full100 envelope comparison; preserve source/context occurrences and missingness."""
import argparse
from collections import Counter
import json
import statistics
from affordable_common import InvalidArtifact,read_json,record,verify,write_new
from accommodation_folds_compare import decision
from accommodation_fold_proposals import outcome
from accommodation_nonlinear import contrast_components,relative_components
from nikasha_pool import choose_rows
from motion_envelope_transfer_inventory import PROTOCOL,key
from motion_envelope_scalar import PROFILE

BASELINES={'released_static':'context_composite','tenfold_precision':'union_precision','tenfold_strict':'union_precision_strict'}
NEW=('envelope_origin','envelope_operational','envelope_mathematical')
METHODS=tuple(BASELINES)+('threefold3p5',)+NEW


def call(value,bands,label):
    d=decision(value,bands);return {'R':value,'decision':d,'outcome':outcome(d,label)}


def summarize(rows,methods=METHODS):
    return {m:dict(Counter(r['methods'][m]['outcome']for r in rows))for m in methods}


def aggregate(cells,bands,label,force_missing=False):
    values=[r['R']for r in cells];missing=force_missing or any(v is None for v in values)
    v=None if missing else statistics.median(values)
    return {**call(v,bands,label),'member_R':values,'range':None if missing else max(values)-min(values),'complete':not missing,'required_members':3}


def compare(inventory,collections,strict,threefold,output):
    a=read_json(inventory);ref=read_json(verify(a['reference']));old=read_json(strict);prior=read_json(threefold)
    if a['protocol_id']!=PROTOCOL or len(collections)!=2:raise InvalidArtifact('finite100/two shard comparison required')
    oldrows={r['case_id']:r for r in old['rows']};priorpairs={(r['selection_id'],r['case_id']):r for r in prior['rows']}
    priort={t['triple_id']:t for t in prior['triples']}
    bands={m:old['bands'][v]for m,v in BASELINES.items()};bands['threefold3p5']=prior['bands']['triple_union']
    bands.update(envelope_origin=ref['origin_only_reference']['bands'],envelope_operational=ref['variants']['operational']['bands'],envelope_mathematical=ref['variants']['mathematical']['bands'])
    fresh={};where={};unavailable={};shards=set()
    for path in collections:
        d=read_json(path);pm=read_json(verify(d['manifest']));sm=read_json(verify(pm['source_manifest']))
        if sm['inventory']!=record(inventory)or sm['reference']!=a['reference']or sm['shard']in shards or pm['numerical_policy_id']!=PROFILE:
            raise InvalidArtifact('actual shard/reference/profile differs')
        shards.add(sm['shard'])
        for c in d['cases']:
            k=key(c['source']['union'])
            if k in fresh:raise InvalidArtifact('duplicate source/context')
            fresh[k]=c;where[k]=record(path)
        by={r['case_id']:r for r in read_json(verify(sm['origin_collection']))['rows']}
        for r in sm['unavailable']:unavailable[key(by[r['case_id']]['source'])]=r
    expected={key(r)for r in a['rows']if not r['reuse']}
    if shards!={0,1}or set(fresh)|set(unavailable)!=expected:raise InvalidArtifact('terminal101 source coverage differs')
    pilot=read_json(verify(a['pilot_collection']));pilotrows={c['case_id']:c for c in pilot['cases']};rows=[];by_pair={}
    for r in a['rows']:
        k=key(r);cid=r['case_id'];label=oldrows[cid]['expected_class'];saved=r['reuse']
        c=pilotrows[saved['case_id']]if saved else fresh.get(k)
        if c and c['source']['union']!=r['prepared']:raise InvalidArtifact('actual prepared source/context changed')
        if c and c['status']=='prepared' and choose_rows(c['matrix'],[x['id']for x in c['candidates']])!=c['pool']:
            raise InvalidArtifact('actual matrix selection differs')
        available=bool(c and c['pool']['status']=='available');methods={m:oldrows[cid]['methods'][v]for m,v in BASELINES.items()}
        origin=None;work=None
        if c and all(c['matrix'][z]['origin']['status']=='complete'for z in ('Ca','La')):
            origin=contrast_components(c['matrix']['Ca']['origin']['components'],c['matrix']['La']['origin']['components'])
        methods['envelope_origin']=call(origin['composite_R_model_kcal_mol']if origin else None,bands['envelope_origin'],label)
        for variant in ('operational','mathematical'):
            methods['envelope_'+variant]=call(c['pool'][variant]['composite_R_model_kcal_mol']if available else None,bands['envelope_'+variant],label)
        if available:
            work={z:relative_components(c['matrix'][z][c['pool']['rows'][z]['operational_candidate']]['components'],c['matrix'][z]['origin']['components'])for z in ('Ca','La')}
            if abs(work['Ca']['composite_kcal_mol']-work['La']['composite_kcal_mol']-(methods['envelope_operational']['R']-methods['envelope_origin']['R']))>1e-6:
                raise InvalidArtifact('origin/accommodation unit or sign differs')
            for cells in c['matrix'].values():
                for cell in cells.values():
                    if cell['status']!='complete'or any(low.get('observed_TolE_hartree')!=1e-10 for low in cell['low'].values()):
                        raise InvalidArtifact('available pool lacks complete strict scalar cells')
        row={'case_id':cid,'selection_id':r['selection_id'],'root_case_id':r['prepared']['source']['root_case_id'],
            'expected_class':label,'status':'available'if available else'unavailable','methods':methods,'pool':c['pool']if c else None,
            'matrix':c['matrix']if c else None,'origin_components':origin,'selected_work':work,'atom_count':r['atom_count'],'La_charge':r['La_charge'],
            'pool_reused':bool(saved),'source_collection':saved['collection']if saved else where.get(k),'reason':c.get('reason')if c else unavailable.get(k),
            'old_threefold_mapping':'per-triple occurrence; merged new contexts can map to different old contexts'}
        rows.append(row);by_pair[k]=row
    triples=[];occurrences=[]
    for t in a['triples']:
        ot=priort[t['triple_id']]
        if t['members']!=ot['members']or t['protein_id']!=ot['protein_id']:raise InvalidArtifact('original triple membership changed')
        label=oldrows[t['members'][0]]['expected_class'];members=[]
        for cid in t['members']:
            if oldrows[cid]['expected_class']!=label:raise InvalidArtifact('same-group label differs')
            new=by_pair.get((t['selection_id'],cid));previous=priorpairs.get((ot['selection_id'],cid))
            methods={m:oldrows[cid]['methods'][v]for m,v in BASELINES.items()}
            methods['threefold3p5']=previous['methods']['triple_union']if previous else call(None,bands['threefold3p5'],label)
            for m in NEW:methods[m]=new['methods'][m]if new and t['status']=='prepared'else call(None,bands[m],label)
            r={'triple_id':t['triple_id'],'case_id':cid,'protein_id':t['protein_id'],'expected_class':label,
                'selection_id':t['selection_id'],'old_threefold_selection_id':ot['selection_id'],'methods':methods}
            occurrences.append(r);members.append(r)
        result={m:aggregate([r['methods'][m]for r in members],bands[m],label,t['status']!='prepared'if m in NEW else ot['status']!='prepared'if m=='threefold3p5'else False)for m in METHODS}
        triples.append({**t,'expected_class':label,'methods':result})
    if(len(rows),len(triples),len(occurrences))!=(104,100,300):raise InvalidArtifact('104/100/300 denominator differs')
    matched={}
    for m in METHODS:
        chosen=[r for r in triples if r['methods'][m]['R']is not None and r['methods']['envelope_operational']['R']is not None]
        matched[m]={'triple_count':len(chosen),'counts':summarize(chosen)}
    groups=[]
    for gid in dict.fromkeys(t['protein_id']for t in triples):
        group=[t for t in triples if t['protein_id']==gid]
        if len(group)!=4:raise InvalidArtifact('four declared triples per protein required')
        groups.append({'protein_id':gid,'triple_ids':[t['triple_id']for t in group],
            'counts':summarize(group),'all_four_complete':{m:all(t['methods'][m]['R']is not None for t in group)for m in METHODS}})
    result={'protocol_id':PROTOCOL,'inventory':record(inventory),'reference':a['reference'],'collections':[record(p)for p in collections],
        'strict_baseline':record(strict),'old_threefold_baseline':record(threefold),'bands':bands,'rows':rows,'triples':triples,
        'source_context_occurrences':occurrences,'groups':groups,'denominators':{'prepared_unique_pairs':104,'triples':100,'occurrences':300,'proteins':25},
        'counts':{'pairs':summarize(rows,tuple(BASELINES)+NEW),'triples':summarize(triples),'occurrences':summarize(occurrences)},
        'matched_triples':matched,'new_thresholds_fitted':False,'new_molecular_calls':0,'production_changed':False,
        'interpretation':'Overlapping developmental subsets of25 consumed proteins; old3.5 uses its actual old scalar profile. Ca-conditioned pilot regression remains separate.',
        'implementation':record(__file__)}
    write_new(output,result);return {k:v for k,v in result.items()if k not in ('rows','triples','source_context_occurrences','groups')}

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('inventory','collections','strict','threefold','output'):p.add_argument('--'+k,required=True,**({'nargs':'+'}if k=='collections'else{}))
    print(json.dumps(compare(**vars(p.parse_args())),indent=2))
