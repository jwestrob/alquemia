"""Conditional canonical extension of the declared charge-feature ablation."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
from affordable_common import InvalidArtifact,cache_key,read_json,record,verify,write_new,xyz
from mace_file_checks import cached_file_checks
from mace_hybrid import check_atoms,write_xyz,EV_TO_KCAL
from mace_omol_ablation import ADAPTER,COMPONENT
from mace_omol_ablation_run import PROTOCOL,SEMANTICS,descriptor_model,validate as development_validate,collect as development_collect
from mace_omol_intact import EVALUATION,geometry,collect as collect_endpoints
from mace_omol_panel import DECISION,expected_tasks as native_tasks
from mace_omol_panel_prepare import validate as prepared,source_rows
from mace_omol_backbone_audit import inspect as inspect_backbone

PREPARATION_SHA='fe91ea6fa3090f01897e828c434ac24076ffad54dcb5d0269659ddf5fa10fee8'


def inputs(preparation,development_collection):
    if record(preparation)['sha256']!=PREPARATION_SHA:
        raise InvalidArtifact('canonical ablation requires the declared strict v2 preparation')
    data=prepared(preparation)
    saved=read_json(development_collection);mp=verify(saved['manifest']);parent=read_json(mp)
    if parent.get('stage')!='ablation_development' or parent.get('protocol_id')!=PROTOCOL:
        raise InvalidArtifact('canonical ablation requires its own completed development protocol')
    development_validate(mp);actual=development_collect(mp)
    if saved!=actual or not actual['canonical_extension_permitted']:
        raise InvalidArtifact('all actual ablation development gates must pass before canonical extension')
    if parent['inventory']!=data['source_inventory']:
        raise InvalidArtifact('canonical source inventory differs')
    tasks=[];reuse={};values={}
    for t in native_tasks(data):
        t.update(energy_component=COMPONENT,charge_feature_adapter=ADAPTER,output_semantics=SEMANTICS)
        if t['case_id'] in ('1H4I','4MAE'):
            name='PQQ_'+t['case_id'];r=actual['scores'][name]['endpoints'][t['metal']][t['position']]
            source=next(s for s in parent['tasks'] if s['task_id']==r['task_id'])
            for field in ('source_xyz','preparation','charge','spin_multiplicity','state','metal_index','position',
                          'variant','assembly','microstate','explicit_waters','energy_component','charge_feature_adapter',
                          'edge_adapter','require_isolated_metal','energy_only','capture_native_readout'):
                if t[field]!=source[field]:raise InvalidArtifact('canonical descriptor reuse differs: '+field)
            if xyz(verify(source['xyz']))!=geometry(t):raise InvalidArtifact('canonical reuse geometry differs')
            reuse[t['task_id']]={'source_collection':record(development_collection),'source_task_id':source['task_id']}
            values[t['task_id']]=r
        else:tasks.append(t)
    if len(tasks)!=100 or len(reuse)!=8 or len(data['rows'])!=28:
        raise InvalidArtifact('declared 100-new/8-reuse/28-case extension changed')
    if [(r['case_id'],r['status']) for r in data['rows'] if r['status']!='prepared']!=[('1KB0','unsupported')]:
        raise InvalidArtifact('unsupported case inventory changed')
    return data,parent,tasks,reuse,values


@cached_file_checks
def prepare(preparation,development_collection,agreement,output):
    from mace_omol import common,seal
    data,parent,tasks,reused,_=inputs(preparation,development_collection)
    _,out,m=common(verify(data['source_inventory']),verify(parent['software']),agreement,output,'ablation_canonical')
    if m['agreement']!=parent['agreement']:raise InvalidArtifact('canonical extension belongs to another plan')
    for name in ('mace_omol_ablation.py','mace_omol_edges.py','mace_omol_products.py','mace_omol_readout.py'):
        if m['implementation'][name]['sha256']!=parent['implementation'][name]['sha256']:
            raise InvalidArtifact('qualified descriptor adapter changed')
    m.update(protocol_id=PROTOCOL,model=descriptor_model(verify(parent['software'])),energy_evaluation=EVALUATION,
             preparation=record(preparation),development_collection=record(development_collection),
             decision_policy=DECISION,output_semantics=SEMANTICS,tasks=tasks,reused=reused,
             unprepared_cases=[r for r in data['rows'] if r['status']!='prepared'])
    for t in tasks:
        if t['position']=='bound':t['xyz']=t['source_xyz']
        else:
            p=out/(t['task_id']+'.xyz');write_xyz(p,geometry(t));t['xyz']=record(p)
    return seal(out,m)


@cached_file_checks
def validate(manifest):
    from mace_omol import SCHEMA,TOL
    m=read_json(manifest)
    data,parent,tasks,reused,_=inputs(verify(m['preparation']),verify(m['development_collection']))
    if (m['schema_version']!=SCHEMA or m['protocol_id']!=PROTOCOL or m['stage']!='ablation_canonical'
            or m['model']!=descriptor_model(verify(m['software'])) or m['energy_evaluation']!=EVALUATION
            or m['decision_policy']!=DECISION or m['output_semantics']!=SEMANTICS or m['tolerances']!=TOL
            or m['software']!=parent['software'] or m['inventory']!=parent['inventory']
            or m['agreement']!=parent['agreement'] or m['reused']!=reused
            or m['unprepared_cases']!=[r for r in data['rows'] if r['status']!='prepared']):
        raise InvalidArtifact('canonical descriptor model/preparation/decision policy changed')
    for pin in [m['agreement'],*m['implementation'].values()]:verify(pin)
    for name in ('mace_omol_ablation.py','mace_omol_edges.py','mace_omol_products.py','mace_omol_readout.py'):
        if m['implementation'][name]['sha256']!=parent['implementation'][name]['sha256']:
            raise InvalidArtifact('descriptor adapter differs from qualification')
    if len(m['tasks'])!=100:raise InvalidArtifact('canonical descriptor task count changed')
    for t,e in zip(m['tasks'],tasks):
        if {k:v for k,v in t.items() if k not in ('xyz','cache_key')}!=e:
            raise InvalidArtifact('canonical descriptor state differs from frozen preparation')
        if xyz(verify(t['xyz']))!=geometry(t) or check_atoms(xyz(verify(t['xyz'])),t['charge'])!=t['state']:
            raise InvalidArtifact('canonical descriptor coordinates/electrons changed')
        payload={k:v for k,v in t.items() if k!='cache_key'}
        if t['cache_key']!=cache_key({'task':payload,'model':m['model'],'software':m['software'],'implementation':m['implementation']}):
            raise InvalidArtifact('canonical descriptor cache differs')
    return {'status':'pass','tasks':100,'reused_endpoints':8,'unsupported_cases':1,'manifest':record(manifest)}


@cached_file_checks
def collect(manifest):
    m=read_json(manifest)
    _,_,_,refs,values=inputs(verify(m['preparation']),verify(m['development_collection']))
    result=collect_endpoints(manifest)
    result.update(reused_rows={k:{**ref,'result':values[k]} for k,ref in refs.items()},
                  unprepared_cases=m['unprepared_cases'],output_semantics=SEMANTICS,
                  preparation_complete=False)
    # Execution/numerical completeness refer to the finite supported inventory.
    # The report retains the unsupported case in the full denominator.
    return result


@cached_file_checks
def report(collection,output):
    from mace_omol_panel_report import summarize
    saved=read_json(collection);mp=verify(saved['manifest']);m=read_json(mp)
    validate(mp);actual=collect(mp)
    if saved!=actual:raise InvalidArtifact('saved canonical descriptor collection differs from receipts')
    data=prepared(verify(m['preparation']));integrity=inspect_backbone(data)
    rows={**actual['rows'],**{k:r['result'] for k,r in actual['reused_rows'].items()}}
    summary=summarize(data,rows,integrity)
    summary['calibration_gap_model_kcal']=summary.pop('calibration_gap_kcal_mol')
    if summary['bands'] is not None:
        b=summary['bands'];summary['bands']={
            'Ca_max_inclusive_model_kcal':b['Ca_max_inclusive_kcal_mol'],
            'La_min_inclusive_model_kcal':b['La_min_inclusive_kcal_mol'],
            'scale':'R_mask','scope':PROTOCOL+'; canonical PQQ only','evidence_use':b['evidence_use']}
    originals={r['case_id']:r for r in source_rows(verify(data['source_inventory']))}
    for score in summary['scores']:
        score['R_mask_model_kcal']=score.pop('R_coord_kcal_mol')
        score['bound_minus_detached_model_eV']=score.pop('bound_minus_detached_eV')
        baseline=originals[score['case_id']]['baseline']
        score['baseline']={k:baseline[k] for k in ('class','published_R_kcal_mol','published_S_kcal_mol','release')}
    result={'status':actual['status'],'preparation_complete':False,'preparation_integrity':integrity,
            'collection':record(collection),'protocol_id':PROTOCOL,'model':m['model'],'agreement':m['agreement'],
            'preparation':m['preparation'],'development_collection':m['development_collection'],
            'output_semantics':SEMANTICS,'score_unit':'kcal_equivalent_model_units',
            'formula':'(T_bound_Ca-T_detached_Ca)-(T_bound_La-T_detached_La)',
            'model_eV_to_model_kcal':EV_TO_KCAL,'decision_policy':DECISION,
            'numerical_gate_pass':actual['numerical_gate_pass'],**summary,
            'baseline_changed':False,'production_promotion':False,'broad_affinity_validated':False,
            'reference':None,'binding_free_energy_kcal_mol':None,'report_implementation':record(__file__)}
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False);write_new(out/'result.json',result)
    lines=['# Canonical charge-feature ablation descriptor','',
           'Modified model units; not quantum energies or binding free energies. All cases consumed.',
           f"Calibration {result['calibration_valid_count']}/25; gap {result['calibration_gap_model_kcal']} model kcal.",
           f"Transfer scores {result['transfer_valid_count']}/3; decisions in expected region {result['transfer_correct_count']}/3.",
           '1KB0 remains unsupported. An unavailable decision is not a measured wrong prediction.','',
           '|Case|Expected|R_mask, model kcal|Decision|','|---|---|---:|---|']
    for s in result['scores']:lines.append(f"|{s['case_id']}|{s['expected_class']}|{s['R_mask_model_kcal']}|{s['calibrated_class']}|")
    lines+=['','Canonical composition already separates these classes. This is not independent proof of broad affinity discrimination.',
            'Baseline and old references unchanged. Missing corrections remain unavailable; no automatic promotion.']
    (out/'REPORT.md').write_text('\n'.join(lines)+'\n')
    return {'status':result['status'],'result':record(out/'result.json')}


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest='command',required=True)
    a=sub.add_parser('prepare')
    for k in ('preparation','development-collection','agreement','output'):a.add_argument('--'+k,required=True)
    a=sub.add_parser('report')
    for k in ('collection','output'):a.add_argument('--'+k,required=True)
    args=vars(p.parse_args());command=args.pop('command')
    print(json.dumps({'prepare':prepare,'report':report}[command](**args),indent=2))
