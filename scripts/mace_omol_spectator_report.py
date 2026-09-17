"""Report the actual frozen spectator test with explicit baseline/new fields."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import shutil
from affordable_common import InvalidArtifact,read_json,record,verify,write_new
from mace_file_checks import cached_file_checks
from mace_omol_spectator import PROTOCOL,SETTINGS,collect


@cached_file_checks
def report(collection,output):
    saved=read_json(collection);mp=verify(saved['manifest']);manifest=read_json(mp)
    if manifest['protocol_id']!=PROTOCOL or manifest['settings']!=SETTINGS:
        raise InvalidArtifact('unexpected spectator method')
    for pin in [manifest['agreement'],*manifest['implementation'].values()]:verify(pin)
    actual=collect(mp)
    if saved!=actual:raise InvalidArtifact('saved spectator collection differs from actual receipts')
    contrasts=[]
    for old in actual['contrasts']:
        # The immutable collector retained the original contrast's `pass` field.
        # Name that field explicitly here so it cannot be mistaken for the new
        # direction or the separately tested spectator-consistency criterion.
        row={k:v for k,v in old.items() if k!='pass'}
        row['original_direction_pass']=old['pass']
        value=row['modified_delta_R_kcal_mol']
        row['modified_direction_pass']=value is not None and value>row['required_min_kcal_mol']
        contrasts.append(row)
    result={'status':actual['status'],'protocol_id':PROTOCOL,'collection':record(collection),
            'manifest':saved['manifest'],'source_report':manifest['source_report'],
            'agreement':manifest['agreement'],'settings':SETTINGS,'scores':actual['scores'],
            'contrasts':contrasts,'numerical_checks':actual['checks'],
            'numerical_gate_pass':actual['numerical_gate_pass'],
            'consistency_checks':actual['consistency_checks'],
            'spectator_consistency_pass':actual['spectator_consistency_pass'],
            'modified_development_directions_pass':all(c['modified_direction_pass'] for c in contrasts),
            'computed_endpoints':sum(r['status']=='computed' for r in actual['rows'].values()),
            'planned_endpoints':20,'successful_attempts':sum(a['accepted'] for a in actual['attempts']),
            'failed_attempts':sum(not a['accepted'] for a in actual['attempts']),
            'evaluation_seconds':sum(r['evaluation_seconds'] for r in actual['rows'].values() if r['status']=='computed'),
            'maximum_direct_Coulomb_bound_kcal_mol':max(t['spectator']['direct_Coulomb_score_bound_kcal_mol'] for t in manifest['tasks']),
            'fragment_charge_status':'intended_formal_fragment_state_not_constrained_or_predicted_by_checkpoint',
            'mechanistic_interpretation':'disconnected spectator changes global charge conditioning; not measured sodium affinity',
            'calibrated_class':None,'baseline_changed':False,'broad_affinity_validated':False,
            'production_promotion':False}
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False)
    shutil.copyfile(__file__,out/Path(__file__).name)
    result['report_implementation']=record(out/Path(__file__).name)
    write_new(out/'result.json',result)
    lines=['# Disconnected sodium consistency check','',
           f'Execution: {result["status"]}; {result["computed_endpoints"]}/20 endpoints.',
           f'Numerical gate: {result["numerical_gate_pass"]}. Spectator consistency: {result["spectator_consistency_pass"]}.',
           f'Modified development directions: {result["modified_development_directions_pass"]}.',
           '', '| Case | Original R_coord | With spectator | Change, kcal/mol | Consistent within0.1 |',
           '|---|---:|---:|---:|---|']
    for name,s in result['scores'].items():
        lines.append(f'| {name} | {s["original_R_coord_kcal_mol"]} | {s["R_coord_kcal_mol"]} | {s["change_kcal_mol"]} | {s["pass"]} |')
    lines+=['','| Contrast | Original | With spectator | Change, kcal/mol | Direction retained |',
            '|---|---:|---:|---:|---|']
    for c in contrasts:
        lines.append(f'| {c["positive_case"]} minus {c["negative_case"]} | {c["delta_R_kcal_mol"]} | {c["modified_delta_R_kcal_mol"]} | {c["change_kcal_mol"]} | {c["modified_direction_pass"]} |')
    lines+=['','The constructed Na+ is at least10000A from every original atom and has zero model graph edges. Original coordinates/protonation/waters remain identical. Total charge increases by one; electron parity remains valid.',
            'The checkpoint receives total charge and cannot constrain independent fragment charges. Interpret sensitivity as a representation issue for the intended separated ionic state, not as measured sodium effects or an exact ground-state error.',
            'All cases are consumed development. No threshold adjustment, new biological label, absolute affinity claim or production promotion.']
    (out/'REPORT.md').write_text('\n'.join(lines)+'\n')
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    for key in ('collection','output'):parser.add_argument('--'+key,required=True)
    result=report(**vars(parser.parse_args()))
    print(json.dumps({k:result[k] for k in ('status','numerical_gate_pass','spectator_consistency_pass',
                     'modified_development_directions_pass','scores')},indent=2))
