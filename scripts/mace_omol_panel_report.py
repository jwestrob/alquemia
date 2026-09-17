"""Report frozen whole-chain canonical bands and all three consumed transfers."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import shutil
from affordable_common import InvalidArtifact,read_json,record,verify,write_new
from mace_hybrid import EV_TO_KCAL
from mace_omol_panel import PROTOCOL,DECISION,validate,collect
from mace_omol_panel_prepare import source_rows,validate as prepared
from mace_file_checks import cached_file_checks
from mace_omol_backbone_audit import inspect as inspect_backbone


def summarize(data,rows,integrity):
    """Pure score algebra; missing observations never become numerical zero."""
    scores=[]
    audited={r['case_id']:r for r in integrity['rows']}
    if set(audited)!={r['case_id'] for r in data['rows']}:
        raise InvalidArtifact('preparation audit case inventory differs')
    for case in data['rows']:
        terms={}
        endpoints={}
        for metal in ('La','Ca'):
            endpoints[metal]={p:rows.get(f'{case["case_id"]}_{metal}_{p}_primary',
                                         {'status':'unavailable','energy_eV':None}) for p in ('bound','detached')}
            pair=endpoints[metal]
            terms[metal]=(pair['bound']['energy_eV']-pair['detached']['energy_eV']) if all(
                r['status']=='computed' for r in pair.values()) else None
        supported=audited[case['case_id']]['status']=='pass'
        available=supported and all(v is not None for v in terms.values())
        scores.append({**case,'status':'computed' if available else 'unavailable',
                       'preparation_status':case['status'] if supported else 'unsupported_geometry',
                       'preparation_integrity':audited[case['case_id']],
                       'R_coord_kcal_mol':None if not available else (terms['Ca']-terms['La'])*EV_TO_KCAL,
                       'bound_minus_detached_eV':terms if supported else {'Ca':None,'La':None},
                       'endpoints':endpoints,'endpoint_use':'scoring' if supported else 'invalid_preparation_diagnostic_only',
                       'calibrated_class':None,'decision_status':'calibration_unavailable'})
    calibration=[s for s in scores if s['evaluation_role']=='calibration']
    transfers=[s for s in scores if s['evaluation_role']!='calibration']
    if len(calibration)!=25 or len(transfers)!=3:raise InvalidArtifact('declared canonical roles changed')
    gap=None;bands=None
    if all(s['status']=='computed' for s in calibration):
        upper=max(s['R_coord_kcal_mol'] for s in calibration if s['expected_class']=='Ca')
        lower=min(s['R_coord_kcal_mol'] for s in calibration if s['expected_class']=='La')
        gap=lower-upper
        if gap>DECISION['minimum_gap_kcal_mol']:
            bands={'Ca_max_inclusive_kcal_mol':upper,'La_min_inclusive_kcal_mol':lower,
                   'scale':'R_coord','scope':'this_whole_chain_canonical_protocol_only',
                   'evidence_use':'consumed_calibration; no independent accuracy claim'}
    if bands is not None:
        for score in scores:
            value=score['R_coord_kcal_mol']
            if value is None:continue
            score['calibrated_class']='Ca' if value<=bands['Ca_max_inclusive_kcal_mol'] else (
                'La' if value>=bands['La_min_inclusive_kcal_mol'] else 'inconclusive')
            score['decision_status']='research_calibrated'
    for score in scores:
        score['expected_direction_reached']=score['calibrated_class']==score['expected_class']
    return {'scores':scores,'calibration_gap_kcal_mol':gap,'bands':bands,
            'calibration_valid_count':sum(s['status']=='computed' for s in calibration),
            'calibration_total_count':25,'transfer_valid_count':sum(s['status']=='computed' for s in transfers),
            'transfer_total_count':3,'transfer_correct_count':sum(s['expected_direction_reached'] for s in transfers),
            'canonical_decision_gate_pass':bands is not None and all(s['expected_direction_reached'] for s in transfers),
            'unavailable_score_count':sum(s['status']!='computed' for s in scores),
            'outside_training_charge_cases':[s['case_id'] for s in scores if s.get('outside_reported_training_charge_range')],
            'outside_training_size_cases':[s['case_id'] for s in scores if s.get('outside_reported_training_size_range')]}


@cached_file_checks
def report(collection,preparation_audit,output):
    saved=read_json(collection);mp=verify(saved['manifest']);m=read_json(mp)
    validate(mp);actual=collect(mp)
    if saved!=actual:raise InvalidArtifact('saved panel collection differs from actual receipts')
    data=prepared(verify(m['preparation']))
    audit=read_json(preparation_audit)
    if audit['preparation_manifest']!=m['preparation']:
        raise InvalidArtifact('preparation audit belongs to another input manifest')
    verify(audit['implementation']);verify(audit['agreement'])
    current=inspect_backbone(data)
    if {k:audit[k] for k in current}!=current:
        raise InvalidArtifact('preparation audit differs from actual geometry')
    rows={**actual['rows'],**{k:r['result'] for k,r in actual['reused_rows'].items()}}
    summary=summarize(data,rows,audit)
    originals={r['case_id']:r for r in source_rows(verify(data['source_inventory']))}
    for score in summary['scores']:
        baseline=originals[score['case_id']]['baseline']
        score['baseline']={k:baseline[k] for k in ('class','published_R_kcal_mol','published_S_kcal_mol','release')}
    result={'status':actual['status'] if audit['pass'] else 'incomplete_preparation',
            'execution_status':actual['status'],'preparation_integrity_pass':audit['pass'],
            'preparation_audit':record(preparation_audit),'protocol_id':PROTOCOL,'collection':record(collection),
            'preparation':m['preparation'],'agreement':m['agreement'],'model':m['model'],
            'decision_policy':DECISION,'product_equivalence':m['product_equivalence'],
            'formula':'(E_bound_Ca-E_detached_Ca)-(E_bound_La-E_detached_La)',
            'energy_units':'eV','score_units':'kcal/mol','eV_to_kcal_mol':EV_TO_KCAL,
            'numerical_gate_pass':actual['numerical_gate_pass'],**summary,
            'canonical_operational_gate_pass':actual['numerical_gate_pass'] and summary['canonical_decision_gate_pass'],
            'baseline_changed':False,'production_promotion':False,'broad_affinity_validated':False,
            'evidence_use':'consumed_retrospective_development','binding_free_energy_kcal_mol':None,
            'solvent_correction_kcal_mol':None,'reference':None,'relaxation_correction_kcal_mol':None,
            'force_validation_status':'not_requested_energy_only'}
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False)
    shutil.copyfile(__file__,out/Path(__file__).name);result['report_implementation']=record(out/Path(__file__).name)
    write_new(out/'result.json',result)
    lines=['# Intact OMOL canonical extension','',
           f'Numerical gate: {result["numerical_gate_pass"]}. Canonical operational gate: {result["canonical_operational_gate_pass"]}.',
           f'Calibration: {result["calibration_valid_count"]}/25 computed; gap {result["calibration_gap_kcal_mol"]} kcal/mol.',
           f'Transfers: {result["transfer_valid_count"]}/3 computed, {result["transfer_correct_count"]}/3 reach their expected region.',
           '', '| Case | Role | Expected | R_coord, kcal/mol | Decision | Charge-range extrapolation |',
           '|---|---|---|---:|---|---|']
    for s in summary['scores']:
        lines.append(f'| {s["case_id"]} | {s["evaluation_role"]} | {s["expected_class"]} | {s["R_coord_kcal_mol"]} | {s["calibrated_class"]} | {s.get("outside_reported_training_charge_range")} |')
    lines+=['','All sources and charges were fixed before scoring. No outcome-based exclusions or separate charge-specific bands.',
            'All whole proteins exceed reported training sizes; four also exceed the training charge range. Numerical success does not remove these extrapolations.',
            'These are consumed PQQ functional-class references. Calibration separation is not an independent accuracy estimate, and structural replicates/homologues are not independent biological observations.',
            'Do not apply these bands to non-PQQ affinity labels. The prior alpha/GGR comparison remains a separate, qualified development result.',
            'No absolute binding-free-energy interpretation, aquo reference, solvent, relaxation, entropy or production promotion. Baseline and historical references remain unchanged.']
    (out/'REPORT.md').write_text('\n'.join(lines)+'\n')
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    for key in ('collection','preparation-audit','output'):parser.add_argument('--'+key,required=True)
    result=report(**vars(parser.parse_args()))
    print(json.dumps({k:result[k] for k in ('status','numerical_gate_pass','canonical_operational_gate_pass',
                     'calibration_gap_kcal_mol','transfer_correct_count','unavailable_score_count')},indent=2))
