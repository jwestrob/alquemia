"""Report the fixed eight-source finite-probe outcome without recalibration."""
import argparse
from collections import Counter, defaultdict
import json
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
from affordable_common import InvalidArtifact, read_json, record, verify, write_new
from accommodation_nonlinear import relative_components
from nikasha_pool import choose_rows


def classify(value,bands):
    if value is None:return 'unavailable'
    if value<=bands['Ca_max']:return 'Ca'
    if value>=bands['La_min']:return 'La'
    return 'inconclusive'


def outcome(call,known):
    return call if call in ('unavailable','inconclusive') else 'correct' if call==known else 'wrong'


def report(collection,specification,output):
    col=read_json(collection);spec=read_json(specification);m=read_json(verify(col['manifest']))
    if m['specification']!=record(specification) or m['branch']!='solvent_guided':raise InvalidArtifact('wrong finite branch/specification')
    inputs=read_json(verify(spec['inputs']));ref=read_json(verify(inputs['adaptive_reference']))
    if [c['case_id'] for c in col['cases']]!=[c['case_id'] for c in inputs['cases']]:raise InvalidArtifact('common8 source order differs')
    bands={'released':ref['old_frozen_bands'],'adaptive':ref['variants']['operational']['bands']}
    rows=[];parents={}
    for c,src,proposal in zip(col['cases'],inputs['cases'],spec['cases']):
        pin=src['pool_collection'];path=pin['path']
        if path not in parents:parents[path]=read_json(verify(pin))
        base=next(b for b in parents[path]['cases'] if b['case_id']==c['case_id'])
        if c['prior_pool']!=base['pool']:raise InvalidArtifact('archived adaptive result changed')
        if c['status']=='prepared' and choose_rows(c['matrix'],[p['id'] for p in c['candidates']])!=c['pool']:
            raise InvalidArtifact('actual finite selection algebra differs')
        original=base['old_result']['R0']['composite_R_model_kcal_mol']
        prior=base['pool']['operational']['composite_R_model_kcal_mol']
        good=c['pool']['status']=='available'
        value=c['pool']['operational']['composite_R_model_kcal_mol'] if good else None
        selected={};works={}
        if good:
            for z in ('Ca','La'):
                a=c['pool']['rows'][z]['operational_candidate'];b=base['pool']['rows'][z]['operational_candidate']
                selected[z]={'prior':b,'new':a}
                works[z]=relative_components(c['matrix'][z][a]['components'],base['matrix'][z][b]['components'])
            difference=works['Ca']['composite_kcal_mol']-works['La']['composite_kcal_mol']
            if abs((value-prior)-difference)>1e-7:raise InvalidArtifact('incremental contrast sign differs')
        transfer={name:{model:classify(R,b) for model,R in [('released_static',original),('adaptive_prior',prior),('new_probes',value)]} for name,b in bands.items()}
        rows.append({'case_id':c['case_id'],'known_class_report_only':src['known_class'],
            'biological_group':base['old_result'].get('biological_group',c['case_id']),
            'conditioning':base['old_result'].get('source_conditioning_metal','reference'),
            'status':'available' if good else 'unavailable','reason':c['pool'].get('reason'),
            'search_probe_coverage':proposal['probe_coverage'],'source_probe_failures':[p for p in proposal['probes'] if p['status']!='available'],
            'released_R':original,'prior_adaptive_R':prior,'probe_R':value,
            'delta_R_vs_adaptive':value-prior if good else None,
            'endpoint_selection':selected,'endpoint_incremental_work_kcal_mol':works,
            'reference_transfer_only':transfer,'new_protocol_reference':None,'validated_new_protocol_decision':None})
    summaries={}
    for band in bands:
        summaries[band]={model:dict(Counter(outcome(r['reference_transfer_only'][band][model],r['known_class_report_only']) for r in rows))
                         for model in ('released_static','adaptive_prior','new_probes')}
    groups=defaultdict(list)
    for r in rows:groups[r['biological_group']].append(r)
    spreads=[]
    for group,rs in groups.items():
        if len(rs)<2:continue
        available=all(r['probe_R'] is not None for r in rs)
        spread=lambda key:max(r[key] for r in rs)-min(r[key] for r in rs)
        spreads.append({'biological_group':group,'source_ids':[r['case_id'] for r in rs],
            'selected_source_count':len(rs),'complete_new_coverage':available,
            'prior_adaptive_range_kcal_mol':spread('prior_adaptive_R'),
            'probe_range_kcal_mol':spread('probe_R') if available else None,
            'scope':'only the predeclared selected source pair, not a full protein conformer distribution'})
    result={'collection':record(collection),'specification':record(specification),'adaptive_reference':inputs['adaptive_reference'],
        'rows':rows,'denominator':8,'available':sum(r['status']=='available' for r in rows),
        'reference_transfer_counts':summaries,'selected_source_group_spreads':spreads,
        'new_MACE_complete':col['MACE_calls_complete'],'new_MACE_required':m['new_MACE_cells'],
        'new_GFN2_complete':col['GFN2_complete'],'new_GFN2_required':m['maximum_new_GFN2_calls'],
        'reference_status':'new protocol uncalibrated; old bands are development transfer only',
        'numerical_stability_status':'predictive interpretation held pending parallel near-identical Q88JH5 native GFN2 reproducibility audit',
        'new_fits_or_molecular_calls_by_report':0,'production_changed':False}
    write_new(output,result)
    lines=['# Finite solvent-guided common8 result','',
        '**Predictive interpretation is on hold:** a parallel near-identical Q88JH5 geometry comparison exposed a large native GFN2 solvent-contrast discrepancy. Until that numerical issue is resolved, the following energies, band transitions and spread changes are descriptive outcomes only. The approved finite pilot was not changed in response.','',
        f'Complete new composite scores: {result["available"]}/8. The existing adaptive pool remains the comparator; all32declared probe slots, including5geometrically unsupported directions, remain in the preparation record.','',
        '| Source | Known class | Previous adaptive R | Probe R | ΔR | Adaptive-band transfer |',
        '|---|---|---:|---:|---:|---|']
    fmt=lambda x:'unavailable' if x is None else f'{x:.6f}'
    for r in rows:
        t=r['reference_transfer_only']['adaptive']
        lines.append('| '+r['case_id']+' | '+r['known_class_report_only']+' | '+fmt(r['prior_adaptive_R'])+' | '+fmt(r['probe_R'])+' | '+fmt(r['delta_R_vs_adaptive'])+' | '+t['adaptive_prior']+' → '+t['new_probes']+' |')
    lines+=['','These are consumed, deliberately selected development sources. Frozen-band transfer is not an independently calibrated decision for the expanded search. No threshold was fitted to these outcomes. R=Ca−La, in model kcal/mol.','',
            '## Actual execution and incremental work','',
            f'Native MACE: {result["new_MACE_complete"]}/{result["new_MACE_required"]}; native GFN2: {result["new_GFN2_complete"]}/{result["new_GFN2_required"]}. No newDFT, optimization, new fold or CPCM call.','',
            '| Source | Ca work | La work | Native contribution to ΔR | Solvent contribution to ΔR |',
            '|---|---:|---:|---:|---:|']
    for r in rows:
        w=r['endpoint_incremental_work_kcal_mol']
        vals=[w[z]['composite_kcal_mol'] for z in ('Ca','La')]+[w['Ca'][k]-w['La'][k] for k in ('native_MACE_kcal_mol','solvent_transfer_kcal_mol')] if w else [None]*4
        lines.append('| '+r['case_id']+' | '+' | '.join(map(fmt,vals))+' |')
    lines+=['','Both metals compete on the same admitted candidate union. Mathematical minima, every component and selected candidate are retained in the raw collection. A required energy failure makes the source unavailable rather than silently shrinking the pool.','',
            '## Frozen-reference transfer counts','',json.dumps(summaries,indent=2),'',
            '## Predeclared selected source pairs','',json.dumps(spreads,indent=2),'',
            'No full225 expansion or default promotion is authorized by this result.','']
    Path(output).with_suffix('.md').write_text('\n'.join(lines))
    return {'output':record(output),'available':result['available'],'reference_transfer_counts':summaries,'spreads':spreads}


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--collection',required=True);p.add_argument('--specification',required=True);p.add_argument('--output',required=True)
    print(json.dumps(report(**vars(p.parse_args()))))
