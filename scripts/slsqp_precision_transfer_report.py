"""Join uniform-precision transfer with immutable 225-source prior ledgers."""
from __future__ import annotations
import argparse
from collections import Counter
import json
from pathlib import Path
from affordable_common import InvalidArtifact,read_json,record,verify,write_new,HA_TO_KCAL
from accommodation_folds_compare import decision
from accommodation_fold_proposals import outcome,strict_summary,summarize_rows
from accommodation_nonlinear import relative_components
from mace_hybrid import EV_TO_KCAL
from slsqp_precision_expansion_report import endpoint
import nikasha_pool as pool
import slsqp_precision as pilot
import slsqp_precision_transfer as transfer

METHODS=('union_precision','union_precision_mathematical')
VARIANTS=('operational','mathematical')


def numerical(old,new,old_manifest,new_manifest):
    """Compare actual successful cells only; failed-old endpoints stay distinct."""
    op=verify(old_manifest).parent;np=verify(new_manifest).parent
    cid=new['case_id'];ends={}
    for z in ('Ca','La'):
        a=op/'proposals'/(cid+'__'+z)/'result.json';b=np/'proposals'/(cid+'__'+z)/'result.json'
        ends[z]=endpoint(record(a),record(b) if b.exists() else None)
    both=old['pool']['status']=='available' and new['pool']['status']=='available';cells=[]
    if both:
        for z in ('Ca','La'):
            for name in ('adaptive_Ca','adaptive_La'):
                a=old['matrix'][z][old['aliases'][name]['representative']]['components']
                b=new['matrix'][z][new['aliases'][name]['representative']]['components']
                dv=(b['GFN2_vacuum_hartree']-a['GFN2_vacuum_hartree'])*HA_TO_KCAL
                da=(b['GFN2_ALPB_hartree']-a['GFN2_ALPB_hartree'])*HA_TO_KCAL
                cells.append({'metal':z,'candidate':name,'MACE_delta_kcal_mol':(b['MACE_eV']-a['MACE_eV'])*EV_TO_KCAL,
                    'vacuum_delta_kcal_mol':dv,'ALPB_delta_kcal_mol':da,'transfer_delta_kcal_mol':da-dv,
                    'vacuum_pass':abs(dv)<=pilot.GATES['GFN2_cell_work_kcal_mol'],
                    'ALPB_pass':abs(da)<=pilot.GATES['GFN2_cell_work_kcal_mol']})
    delta={v:new['pool'][v]['composite_R_model_kcal_mol']-old['pool'][v]['composite_R_model_kcal_mol'] if both else None for v in VARIANTS}
    eligible=both and all(e['successful_old_comparison'] for e in ends.values())
    passed=all(e['native_energy_pass'] for e in ends.values()) and len(cells)==4 and all(c['vacuum_pass'] and c['ALPB_pass'] for c in cells) and abs(delta['operational'])<=pilot.GATES['pooled_R_kcal_mol'] if eligible else None
    return {'endpoints':ends,'components':cells,'delta_R':delta,'numerical_pass':passed,
            'successful_pool_comparison':both,'old_status':old['pool']['status'],'new_status':new['pool']['status']}


def compare(selection,collections,baseline,output):
    top=read_json(selection);inv,declared,fresh=transfer.population(verify(top['inventory']),verify(top['rank_qualification']))
    ref=read_json(verify(top['reference']));oldref=read_json(verify(ref['old_reference']));old=read_json(baseline)
    if len(collections)!=4 or top['all_cases']!=declared or top['reference']!=inv['new_reference']:
        raise InvalidArtifact('full frozen-reference four-shard denominator required')
    available={};shards=set();pins=[];cache={}
    def load(pin):
        path=verify(pin)
        if str(path) not in cache:cache[str(path)]=read_json(path)
        return cache[str(path)]
    for path in collections:
        coll=read_json(path);m=load(coll['manifest']);sm=load(m['source_manifest']);i=sm['shard']
        expected=[r['case_id'] for r in fresh[i::4]]
        if i in shards or sm['selection']!=record(selection) or [c['case_id'] for c in coll['cases']]!=expected or coll['denominator']!=len(expected):
            raise InvalidArtifact('duplicate/changed transfer shard')
        if (coll['protocol_id']!=pilot.POOL_PROTOCOL or m['reference']!=top['reference'] or m['settings']!=ref['settings'] or
            m['model']!=ref['model'] or sm['settings']!=ref['optimizer_settings'] or sm['rank_qualification']!=top['rank_qualification']):
            raise InvalidArtifact('frozen candidate method differs')
        low=read_json(Path(path).parent/'solvent/shard_0/manifest.json')
        if low['pool_manifest']!=coll['manifest'] or low['execution_resources']!=transfer.LOW_RESOURCES or low['rank_qualification']!=top['rank_qualification']:
            raise InvalidArtifact('scalar profile changed')
        shards.add(i);pins.append(record(path))
        sources={s['case_id']:s for s in sm['sources']}
        for c in coll['cases']:
            cid=c['case_id']
            if cid in available:raise InvalidArtifact('duplicate source')
            src=sources[cid];oc=load(src['old_pool']);prior=next(x for x in oc['cases'] if x['case_id']==cid)
            available[cid]={'case':c,'collection':record(path),'manifest':m,'source_manifest':m['source_manifest'],
                'old_case':prior,'old_manifest':src['parent'],'reuse':False}
    if len(available)!=202 or shards!={0,1,2,3}:raise InvalidArtifact('202 fresh sources required')
    prior34=load(inv['precision34_comparison']);priorby={r['case_id']:r for r in prior34['rows']}
    for source in declared:
        if not source['precision34_reuse']:continue
        cid=source['case_id'];coll=load(source['precision34_reuse']);m=load(coll['manifest']);sm=load(m['source_manifest'])
        c=next(c for c in coll['cases'] if c['case_id']==cid);r=priorby[cid]
        if m['settings']!=ref['settings'] or sm['settings']!=ref['optimizer_settings'] or m['model']!=ref['model']:
            raise InvalidArtifact('reused precision model differs')
        available[cid]={'case':c,'collection':source['precision34_reuse'],'manifest':m,'source_manifest':m['source_manifest'],
            'reuse':True,'prior_numeric':{k:r[k] for k in ('endpoints','components','delta_R','numerical_pass')},'old_manifest':r['parent']}
    if len(available)!=208:raise InvalidArtifact('208 physical sources required')
    by={r['case_id']:r for r in declared};bands={**old['bands'],**{name:ref['variants'][v]['bands'] for v,name in zip(VARIANTS,METHODS)}}
    if len(old['rows'])!=225 or set(by)!={r['case_id'] for r in old['rows']}:raise InvalidArtifact('original225 ledger differs')
    rows=[]
    for r in old['rows']:
        cid=r['case_id'];source=by[cid];entry=available.get(cid);c=entry['case'] if entry else None
        if (r['expected_class'],r['root_case_id'],r['source_conditioning_metal'])!=(source['known_class'],source['root_case_id'],source['source_conditioning_metal']):
            raise InvalidArtifact('biological/source grouping changed')
        chosen=c['pool'] if c else None;ok=bool(chosen and chosen['status']=='available');new={};work={}
        if c and c['status']=='prepared' and pool.choose_rows(c['matrix'],[x['id'] for x in c['candidates']])!=chosen:
            raise InvalidArtifact('pool selection/component replay changed')
        if c and c['source']!=next(x for x in load(entry['old_manifest'])['cases'] if x['case_id']==cid):
            raise InvalidArtifact('original source context changed')
        for v,name in zip(VARIANTS,METHODS):
            value=chosen[v]['composite_R_model_kcal_mol'] if ok else None;call=decision(value,bands[name])
            new[name]={'R':value,'decision':call,'outcome':outcome(call,r['expected_class'])}
            if ok:
                work[v]={z:relative_components(c['matrix'][z][chosen['rows'][z][v+'_candidate']]['components'],c['matrix'][z]['origin']['components']) for z in ('Ca','La')}
                origin=pool.score(c['matrix']['Ca']['origin']['components'],c['matrix']['La']['origin']['components'])['composite_R_model_kcal_mol']
                if abs(value-origin-(work[v]['Ca']['composite_kcal_mol']-work[v]['La']['composite_kcal_mol']))>1e-6:
                    raise InvalidArtifact('energy/sign/unit accounting differs')
        numbers=(entry['prior_numeric'] if entry['reuse'] else numerical(entry['old_case'],c,entry['old_manifest'],entry['source_manifest'])) if entry else None
        rows.append({**r,'methods':{**r['methods'],**new},'precision_status':chosen['status'] if chosen else source['union_origin_status'],
            'precision_reason':chosen.get('reason') if chosen else source.get('reason'),'precision_work':work or None,
            'precision_candidates':{v:{z:chosen['rows'][z][v+'_candidate'] for z in ('Ca','La')} for v in VARIANTS} if ok else None,
            'precision_numerical_comparison':numbers,'precision_collection':entry['collection'] if entry else None,
            'precision_reused_from34':bool(entry and entry['reuse']),
            'precision_old_reference_transfer':{v:decision(new[name]['R'],oldref['variants'][v]['bands']) for v,name in zip(VARIANTS,METHODS)}})
    idx={r['case_id']:r for r in rows};groups=[]
    if (len(old['pools']),len(old['triples']))!=(75,100):raise InvalidArtifact('original strict aggregation denominator differs')
    for r in old['pools']:
        new={}
        for name in METHODS:
            if r['pool'] in ('La4','Ca5'):value=strict_summary(r['members'],idx,name,bands[name],r['expected_class'])
            elif r['pool']=='balanced':
                arms=[next(p for p in groups if p['root_case_id']==r['root_case_id'] and p['pool']==arm)['methods'][name] for arm in ('La4','Ca5')]
                number=None if any(a['R'] is None for a in arms) else sum(a['R'] for a in arms)/2
                call=decision(number,bands[name]);value={'R':number,'decision':call,'outcome':outcome(call,r['expected_class']),
                    'required_members':len(r['members']),'missing_members':sorted({x for a in arms for x in a['missing_members']}),
                    'within_protein_range':None,'equal_mean_of_complete_arm_medians':True}
            else:raise InvalidArtifact('unsupported aggregate')
            new[name]=value
        groups.append({**r,'methods':{**r['methods'],**new}})
    triples=[{**r,'methods':{**r['methods'],**{name:strict_summary(r['members'],idx,name,bands[name],r['expected_class']) for name in METHODS}}} for r in old['triples']]
    subsets={'all225':rows,'La100':[r for r in rows if r['source_conditioning_metal']=='La'],
        'Ca125':[r for r in rows if r['source_conditioning_metal']=='Ca'],**{p:[r for r in groups if r['pool']==p] for p in ('La4','Ca5','balanced')},'La100_triples':triples}
    names=list(rows[0]['methods']);counts={key:summarize_rows(rs,names) for key,rs in subsets.items()};matched={}
    for key,rs in subsets.items():
        matched[key]={}
        for other in names:
            if other in METHODS:continue
            common=[r for r in rs if all(r['methods'][x]['R'] is not None for x in (METHODS[0],other))]
            matched[key][other]={'declared':len(rs),'common':len(common),'counts':summarize_rows(common,(METHODS[0],other)),
                'transitions':dict(Counter(r['methods'][other]['outcome']+'->'+r['methods'][METHODS[0]]['outcome'] for r in common))}
    result={'protocol_id':pilot.POOL_PROTOCOL,'selection':record(selection),'reference':top['reference'],'rank_qualification':top['rank_qualification'],
        'baseline':record(baseline),'collections':pins,'six_reuses_from_precision34':inv['precision34_comparison'],
        'rows':rows,'pools':groups,'triples':triples,'bands':bands,'counts':counts,'matched':matched,
        'numerical_gates':pilot.GATES,'numerical_comparison_denominator':sum(bool(r['precision_numerical_comparison'] is not None and r['precision_numerical_comparison']['numerical_pass'] is not None) for r in rows),
        'numerical_passes':sum(bool(r['precision_numerical_comparison'] is not None and r['precision_numerical_comparison']['numerical_pass'] is True) for r in rows),
        'new_thresholds_fitted':False,'new_molecular_calls':0,'production_changed':False,'implementation':record(__file__),
        'interpretation':'Separately calibrated candidate; consumed source folds and correlated aggregates, not new independent biological validation. Failed numerical equivalence gates retained.'}
    write_new(output,result);return counts['all225']


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for arg in ('selection','baseline','output'):p.add_argument('--'+arg,required=True,type=Path)
    p.add_argument('--collections',required=True,type=Path,nargs='+')
    print(json.dumps(compare(**vars(p.parse_args())),indent=2))
