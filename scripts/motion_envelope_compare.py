"""Canonical-only envelope reference and same-source strict-tenfold comparison."""
import argparse
from collections import Counter
from datetime import datetime,timezone
import json
from pathlib import Path
import numpy as np
from affordable_common import InvalidArtifact,read_json,record,verify,write_new,xyz
from accommodation_folds_compare import decision
from accommodation_fold_proposals import outcome
from accommodation_nonlinear import contrast_components,relative_components
from nikasha_pool_compare import extrema_reference
from nikasha_pool import choose_rows
import motion_envelope_pool as envelope

VARIANTS=('mathematical','operational')
REFERENCE_ID='Nikasha_motion_envelope4p3_adaptive_strict_native_canonical25_v1'


def checked_pool(c):
    if c['status']=='prepared':
        actual=choose_rows(c['matrix'],[r['id']for r in c['candidates']])
        if actual!=c['pool']:raise InvalidArtifact('actual pool algebra differs')
    else:actual=c['pool']
    return actual


def reference(collection,output):
    d=read_json(collection);m=read_json(verify(d['manifest']));envelope.validate_pool(verify(d['manifest']))
    source=read_json(verify(m['source_manifest']));native=read_json(verify(source['origin_collection']));inventory=read_json(verify(native['inventory']))
    declared=inventory['pilot'];rows=[]
    if len(declared)!=34 or len(d['cases'])!=34:raise InvalidArtifact('fixed34 source denominator differs')
    for c,s in zip(d['cases'],declared):
        origin=c['source']['origin_row']
        if(origin['source_case_id'],origin['selection_id'])!=(s['case_id'],s['selection_id']):raise InvalidArtifact('declared reference/probe identity differs')
        selected=checked_pool(c);actual_source=origin['source']['source']
        canonical=s['role']=='canonical_calibration'
        if canonical and not actual_source.get('canonical_coordinate_match'):raise InvalidArtifact('noncanonical geometry entered calibration')
        root=actual_source.get('root_case_id',s['case_id'])
        rows.append({'case_id':c['case_id'],'source_case_id':s['case_id'],'root_case_id':root,'role':s['role'],
            'expected_class':s['expected_class'],'declared_triple_id':s.get('declared_triple_id'),'selection_id':s['selection_id'],
            'origin_R_model_kcal_mol':contrast_components(c['matrix']['Ca']['origin']['components'],c['matrix']['La']['origin']['components'])['composite_R_model_kcal_mol']if all(c['matrix'][z]['origin']['status']=='complete'for z in ('Ca','La'))else None,
            'status':selected['status'],'scores':{v:selected[v]['composite_R_model_kcal_mol']if selected['status']=='available'else None for v in VARIANTS},'pool':selected})
    canonical=[r for r in rows if r['role']=='canonical_calibration']
    if len(canonical)!=25 or len({r['root_case_id']for r in canonical})!=25:raise InvalidArtifact('25 original canonical proteins required')
    variants={v:extrema_reference([{**r,'R_model_kcal_mol':r['scores'][v]}for r in canonical],v,REFERENCE_ID)for v in VARIANTS}
    origin_reference=extrema_reference([{**r,'R_model_kcal_mol':r['origin_R_model_kcal_mol']}for r in canonical],'origin',REFERENCE_ID+'_origin_only')
    for r in rows:
        r['origin_own_reference']=call(r['origin_R_model_kcal_mol'],origin_reference['bands'],r['expected_class'])
        r['own_reference']={v:call(r['scores'][v],variants[v]['bands'],r['expected_class'])for v in VARIANTS}
    result={'protocol_id':envelope.POOL_PROTOCOL,'reference_id':REFERENCE_ID,'collection':record(collection),'inventory':native['inventory'],
        'variants':variants,'origin_only_reference':origin_reference,'origin_ablation_agreement':record(Path(__file__).resolve().parents[1]/'diagnostics/motion_envelope_20260923/ORIGIN_ABLATION_PLAN.md'),'rows':rows,'denominator':34,'canonical_denominator':25,'crystal_denominator':3,'probe_denominator':6,
        'crystals_or_probes_used_for_fit':False,'frozen_UTC':datetime.now(timezone.utc).isoformat(),'implementation':record(__file__),
        'model':m['model'],'pool_settings':m['settings'],'optimizer_settings':source['settings'],'scalar_profile':m['numerical_policy_id'],
        'production_changed':False,'new_molecular_calls':0}
    write_new(output,result);return {'reference':record(output),'variants':{v:{k:r[k]for k in ('status','bands','gap_model_kcal_mol')}for v,r in variants.items()}}


def call(value,bands,expected):
    result=decision(value,bands)if bands is not None else'unavailable'
    return {'decision':result,'outcome':outcome(result,expected)}


def strict_baselines(primary,supplement):
    p=read_json(primary);s=read_json(supplement);rows={}
    for r in p['rows']:
        fresh=r['branches']['fresh'];rows[r['case_id']]={'matrix':fresh['matrix'],'pool':fresh['pool'],'expected_class':r['expected_class'],
            'source':record(primary),'old_collection':r['collection'],'old_pool':r['cold_pool']}
    for r in s['cases']:
        if r['case_id']in rows:raise InvalidArtifact('duplicate strict comparison source')
        rows[r['case_id']]={'matrix':r['matrix'],'pool':r['pool'],'expected_class':r['expected_class'],'source':record(supplement),'old_collection':r['collection'],'old_pool':r['old_pool']}
    for cid,row in rows.items():
        archived=next(c for c in read_json(verify(row['old_collection']))['cases']if c['case_id']==cid)
        row['original_core']=archived['source']['union']['original_core']
        row['source_preparation']=archived['source']['origin_row']['source_preparation']
        row['selected_modes']=[v['id']for v in archived['source']['selection']['selected']]
        row['La_charge']=archived['source']['union']['new_La_charge']
    return rows


def core_match(source,baseline):
    matrix=baseline['matrix'];mapping=read_json(verify(baseline['source_preparation']))['core_to_context']
    differences={}
    for z in ('Ca','La'):
        core=xyz(verify(source['original_core']['endpoints'][z]['xyz']));old_core=xyz(verify(baseline['original_core']['endpoints'][z]['xyz']))
        old=xyz(verify(matrix[z]['origin']['xyz']))
        if [a[0]for a in old_core]!=[a[0]for a in core]:raise InvalidArtifact('matched baseline original core membership differs')
        delta=float(np.max(np.abs(np.array([a[1:]for a in old_core])-np.array([a[1:]for a in core]))))
        for ci,xi in mapping.items():
            if xi is None:continue
            if old[xi][0]!=core[int(ci)][0]:raise InvalidArtifact('baseline source-to-context atom differs')
            delta=max(delta,float(np.max(np.abs(np.array(old[xi][1:])-np.array(core[int(ci)][1:])))))
        if delta>1e-12:raise InvalidArtifact('strict comparator is not the same source core')
        differences[z]=delta
    return differences


def shifts(new,old):
    aliases={'composite_R_model_kcal_mol':'composite_R_model_kcal_mol','native_R_model_kcal_mol':'native_MACE_R_model_kcal_mol','solvation_delta_R_kcal_mol':'solvent_R_kcal_mol'}
    return {k:new[k if k in new else alias]-old[k if k in old else alias]for k,alias in aliases.items()}


def history(case_id,root_case_id,canonical,triple_reference,triple_transfer):
    rows=[]
    for pin,data in ((record(triple_reference),read_json(triple_reference)),(record(triple_transfer),read_json(triple_transfer))):
        for r in data['rows']:
            if r['case_id']==case_id or(canonical and r.get('root_case_id')==root_case_id and r.get('role')=='canonical_calibration'):
                rows.append({'source':pin,'case_id':r['case_id'],'selection_id':r.get('selection_id'),'role':r.get('role'),
                    'pool':r['pool'],'historical_loose_scalar_only':True})
    return rows


def compare(reference,strict_primary,strict_supplement,strict_reference,strict_origin_reference,triple_reference,triple_transfer,output):
    ref=read_json(reference);data=read_json(verify(ref['collection']));baselines=strict_baselines(strict_primary,strict_supplement)
    sr=read_json(strict_reference);strict_bands={v:sr['branches']['fresh']['variants'][v]['bands']for v in VARIANTS};rows=[]
    static=read_json(strict_origin_reference);static_bands=static['static_reference']['bands']
    expected=extrema_reference(static['static_reference']['rows'],'origin','checked_strict_origin')
    if expected['bands']!=static_bands or {r['case_id']for r in static['static_reference']['rows']}!={r['root_case_id']for r in ref['rows']if r['role']=='canonical_calibration'}:raise InvalidArtifact('strict origin reference source or algebra differs')
    for c,r in zip(data['cases'],ref['rows']):
        canonical=r['role']=='canonical_calibration';bid=r['root_case_id']if canonical else r['source_case_id']
        if bid not in baselines:raise InvalidArtifact('missing strict same-source comparator '+bid)
        old=baselines[bid];source=c['source']['origin_row']['source'];identity=core_match(source,old)
        if old['expected_class']!=r['expected_class']:raise InvalidArtifact('matched label differs')
        before=old['pool'];new=c['pool'];origin={label:contrast_components(mx['Ca']['origin']['components'],mx['La']['origin']['components'])
            if all(mx[z]['origin']['status']=='complete'for z in ('Ca','La'))else None for label,mx in [('envelope',c['matrix']),('strict_tenfold',old['matrix'])]}
        origin_shift=shifts(origin['envelope'],origin['strict_tenfold'])if all(origin.values())else None
        variants={}
        for v in VARIANTS:
            available=new['status']=='available'and before['status']=='available';work={}
            if available:
                for label,pool,mx in [('envelope',new,c['matrix']),('strict_tenfold',before,old['matrix'])]:
                    work[label]={z:relative_components(mx[z][pool['rows'][z][v+'_candidate']]['components'],mx[z]['origin']['components'])for z in ('Ca','La')}
            pool_shift=shifts(new[v],before[v])if available else None
            variants[v]={'pool_shift':pool_shift,'origin_shift':origin_shift,'differential_accommodation_shift':{k:pool_shift[k]-origin_shift[k]for k in pool_shift}if pool_shift and origin_shift else None,
                'selected_endpoint_work':work,'envelope_own_reference':r['own_reference'][v],
                'strict_tenfold_own_reference':call(before[v]['composite_R_model_kcal_mol']if before['status']=='available'else None,strict_bands[v],r['expected_class']),
                'envelope_through_strict_tenfold_bands':call(r['scores'][v],strict_bands[v],r['expected_class'])}
        rows.append({**r,'strict_tenfold_case_id':bid,'strict_baseline_source':old['source'],'same_original_core_max_A':identity,
            'preparation_changes':{'new_atom_count':source['atom_count'],'strict_tenfold_atom_count':len(xyz(verify(old['matrix']['La']['origin']['xyz']))),
                'new_La_charge':source['new_La_charge'],'strict_tenfold_La_charge':old['La_charge'],
                'envelope_selected_modes':[v['id']for v in c['source']['selection']['selected']],'tenfold_selected_modes':old['selected_modes']},
            'origin_only_calls':{'envelope':r['origin_own_reference'],'strict_tenfold':call(origin['strict_tenfold']['composite_R_model_kcal_mol']if origin['strict_tenfold']else None,static_bands,r['expected_class'])},
            'origin_components':origin,'matrix':c['matrix'],'strict_tenfold_matrix':old['matrix'],'strict_tenfold_pool':before,'variants':variants,
            'historical_tenfold_loose_pool':old['old_pool'],'historical_threefold_loose':history(r['source_case_id'],r['root_case_id'],canonical,triple_reference,triple_transfer)})
    counts={}
    for role in ['all34',*dict.fromkeys(r['role']for r in rows)]:
        subset=rows if role=='all34'else[r for r in rows if r['role']==role]
        counts[role]={'denominator':len(subset),'origin_only':{name:dict(Counter(r['origin_only_calls'][name]['outcome']for r in subset))for name in ('envelope','strict_tenfold')},**{v:{name:dict(Counter(r['variants'][v][name]['outcome']for r in subset))for name in ('envelope_own_reference','strict_tenfold_own_reference','envelope_through_strict_tenfold_bands')}for v in VARIANTS}}
    spreads=[]
    for root in dict.fromkeys(r['root_case_id']for r in rows):
        group=[r for r in rows if r['root_case_id']==root and r['role']!='canonical_calibration']
        if len(group)<2:continue
        vals={v:{'envelope':[r['scores'][v]for r in group],'strict_tenfold':[r['strict_tenfold_pool'][v]['composite_R_model_kcal_mol']if r['strict_tenfold_pool']['status']=='available'else None for r in group]}for v in VARIANTS}
        spreads.append({'root_case_id':root,'source_ids':[r['source_case_id']for r in group],
            'range_model_kcal_mol':{v:{k:max(x)-min(x)if all(y is not None for y in x)else None for k,x in a.items()}for v,a in vals.items()}})
    result={'protocol_id':envelope.POOL_PROTOCOL,'reference':record(reference),'strict_primary':record(strict_primary),'strict_supplement':record(strict_supplement),
        'strict_reference':record(strict_reference),'strict_origin_reference':record(strict_origin_reference),'historical_threefold_inputs':[record(triple_reference),record(triple_transfer)],'rows':rows,'counts':counts,'prescribed_probe_pair_spreads':spreads,
        'denominator':34,'new_fit_in_comparison':False,'new_molecular_calls':0,'production_changed':False,'implementation':record(__file__),
        'interpretation':'Consumed development calibration/probes; context composition, geometry search and cavity effects change together. No independent biological accuracy claim.'}
    write_new(output,result);return {'comparison':record(output),'counts':counts,'spreads':spreads}


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);s=p.add_subparsers(dest='command',required=True)
    for op,keys in {'reference':('collection','output'),'compare':('reference','strict_primary','strict_supplement','strict_reference','strict_origin_reference','triple_reference','triple_transfer','output')}.items():
        q=s.add_parser(op)
        for k in keys:q.add_argument('--'+k.replace('_','-'),required=True)
    a=vars(p.parse_args());cmd=a.pop('command');print(json.dumps(globals()[cmd](**a),indent=2))
