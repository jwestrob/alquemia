"""Fixed canonical calibration and retrospective transfer; no score refitting."""
from __future__ import annotations
import argparse
import json
import math
from pathlib import Path
import shutil
from affordable_common import InvalidArtifact,read_json,record,verify,write_new
from mace_hybrid import EV_TO_KCAL
from mace_canonical_run import inventory,validate,collect,MACE,GB


def fit_bands(rows):
    calibration=[r for r in rows if r['evaluation_role']=='calibration']
    if len(calibration)!=25 or sum(r['expected_class']=='La' for r in calibration)!=11:
        raise InvalidArtifact('fixed25 calibration roles/classes required')
    valid=[r for r in calibration if r['status']=='computed' and r['R_kcal_mol'] is not None]
    result={'status':'unavailable','calibration_count':25,'valid_calibration_count':len(valid),
            'U_max_Ca_kcal_mol':None,'L_min_La_kcal_mol':None,'gap_kcal_mol':None,
            'minimum_gap_kcal_mol':.02,'decision_bands':None,'transfer_rows_used':False}
    if len(valid)!=25:return result
    if not all(math.isfinite(r['R_kcal_mol']) for r in valid):raise InvalidArtifact('nonfinite canonical contrast')
    upper=max(r['R_kcal_mol'] for r in valid if r['expected_class']=='Ca')
    lower=min(r['R_kcal_mol'] for r in valid if r['expected_class']=='La')
    result.update(U_max_Ca_kcal_mol=upper,L_min_La_kcal_mol=lower,gap_kcal_mol=lower-upper,
                  status='calibratable' if lower-upper>.02 else 'calibration_separation_failed')
    if result['status']=='calibratable':result['decision_bands']={'Ca_supported_max_R_kcal_mol':upper,'La_supported_min_R_kcal_mol':lower}
    return result


def decision(value,bands):
    if value is None:return 'unavailable_endpoint'
    if not math.isfinite(value):raise InvalidArtifact('nonfinite decision input')
    if bands is None:return 'unavailable_calibration'
    if value<=bands['Ca_supported_max_R_kcal_mol']:return 'Ca-supported'
    if value>=bands['La_supported_min_R_kcal_mol']:return 'La-supported'
    return 'inconclusive'


def contrast_components(endpoints):
    ca,la=endpoints['Ca'],endpoints['La']
    vacuum=(ca['MACE_energy_eV']-la['MACE_energy_eV'])*EV_TO_KCAL
    solvent=ca['GB_reaction_kcal_mol']-la['GB_reaction_kcal_mol']
    short=(ca['MACE_components_eV']['interaction_energy']-la['MACE_components_eV']['interaction_energy'])*EV_TO_KCAL
    return {'R_vacuum_kcal_mol':vacuum,'R_GB_kcal_mol':solvent,
            'R_short_kcal_mol':short,'R_kcal_mol':vacuum+solvent}


def checked(path):
    saved=read_json(path);mp=verify(saved['manifest']);validate(mp);actual=collect(mp)
    if saved!=actual:raise InvalidArtifact('collection differs from current successful/failed execution receipts; recollect explicitly')
    rows=dict(saved['rows'])
    for key,reused in saved['reused_rows'].items():rows[key]=reused['result']
    if len(rows)!=56:raise InvalidArtifact('all56 logical endpoints, including unavailable ones, must remain visible')
    return saved,read_json(mp),rows


def checkpoint_report(mace,gb):
    mc,mm,mr=checked(mace);gc,gm,gr=checked(gb)
    if (mm['schema_version']!=MACE or gm['schema_version']!=GB or gm['source_mace_collection']!=record(mace) or
            mm['inventory']!=gm['inventory'] or mm['checkpoint_label']!=gm['checkpoint_label']):
        raise InvalidArtifact('MACE/GB source pair differs')
    data=inventory(verify(mm['inventory']));rows=[]
    for source in data['rows']:
        name=source['case_id'];endpoint={};failure=[]
        for metal in ('La','Ca'):
            key=name+'_'+metal;m=mr[key];g=gr[key]
            if m['status']!='computed' or g['status']!='computed':
                failure.append({'metal':metal,'MACE_status':m['status'],'GB_status':g['status']})
                endpoint[metal]={'status':'unavailable','MACE_energy_eV':None,'GB_reaction_kcal_mol':None}
            else:
                endpoint[metal]={'status':'computed','MACE_energy_eV':m['energy_eV'],
                                 'GB_reaction_kcal_mol':g['GB_reaction_kcal_mol'],
                                 'MACE_components_eV':m['energy_components_eV'],
                                 'density_coefficients':m['density_coefficients'],
                                 'MACE_manifest':m['manifest'],'GB_manifest':g['manifest'],
                                 'total_charge_error_e':m['total_charge_error_e']}
        row={'case_id':name,'expected_class':source['expected_class'],'evaluation_role':source['evaluation_role'],
             'evidence_stratum':source['evidence_stratum'],'sequence_accession_group':source['sequence_accession_group'],
             'prospectively_blind':False,'baseline':source['baseline'],'endpoints':endpoint,
             'status':'unavailable' if failure else 'computed','failures':failure,
             'R_kcal_mol':None,'R_vacuum_kcal_mol':None,'R_GB_kcal_mol':None,'R_short_kcal_mol':None,
             'S_kcal_mol':None,'aquo_reference_status':'not_used_own_raw_R_calibration',
             'entropy_correction_kcal_mol':None,'relaxation_correction_kcal_mol':None}
        if not failure:
            row.update(contrast_components(endpoint))
        rows.append(row)
    calibration=fit_bands(rows)
    for row in rows:
        row['decision']=decision(row['R_kcal_mol'],calibration['decision_bands'])
        row['expected_supported_region']=row['expected_class']+'-supported'
        row['expected_region_pass']=row['decision']==row['expected_supported_region']
    transfer=[r for r in rows if r['evaluation_role']!='calibration']
    if {r['case_id'] for r in transfer}!={'1H4I','4MAE','1KB0'}:raise InvalidArtifact('transfer denominator changed')
    return {'checkpoint':mm['checkpoint_label'],'sources':{'MACE':record(mace),'GB':record(gb),'inventory':mm['inventory']},
            'model':mm['model'],'software':mm['software'],'calibration':calibration,'rows':rows,
            'transfer_count':3,'valid_transfer_count':sum(r['status']=='computed' for r in transfer),
            'supported_transfer_count':sum(r['expected_region_pass'] for r in transfer),
            'inconclusive_transfer_count':sum(r['decision']=='inconclusive' for r in transfer),
            'unavailable_calibration_transfer_count':sum(r['decision']=='unavailable_calibration' for r in transfer),
            'incorrect_supported_transfer_count':sum(r['decision'].endswith('-supported') and not r['expected_region_pass'] for r in transfer),
            'operational_gate_pass':calibration['status']=='calibratable' and all(r['expected_region_pass'] for r in transfer),
            'broad_affinity_validated':False,'incremental_information_beyond_composition_established':False}


def report(medium_mace,medium_gb,large_mace,large_gb,output):
    medium=checkpoint_report(medium_mace,medium_gb);large=checkpoint_report(large_mace,large_gb)
    if medium['checkpoint']!='medium' or large['checkpoint']!='large' or medium['sources']['inventory']!=large['sources']['inventory']:
        raise InvalidArtifact('primary/sensitivity checkpoint roles or inventory changed')
    pairs=[]
    for m,l in zip(medium['rows'],large['rows']):
        if m['case_id']!=l['case_id']:raise InvalidArtifact('comparison row alignment differs')
        pairs.append({'case_id':m['case_id'],'large_minus_medium_R_kcal_mol':None if m['R_kcal_mol'] is None or l['R_kcal_mol'] is None else l['R_kcal_mol']-m['R_kcal_mol'],
                      'medium_decision':m['decision'],'large_decision':l['decision']})
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False)
    implementation=out/'mace_canonical_report.py';shutil.copyfile(__file__,implementation)
    result={'status':'complete','primary_checkpoint':'medium','primary_operational_gate_pass':medium['operational_gate_pass'],
            'implementation':record(implementation),'models':{'medium':medium,'large':large},'checkpoint_sensitivity':pairs,
            'baseline_changed':False,'prospective_validation_claimed':False,'broad_affinity_validated':False,
            'incremental_information_beyond_composition_established':False,
            'interpretation':'functional-class calibration plus consumed structural transfer; both1H4I/P16027 and4MAE/I0JWN7 share calibration accessions'}
    write_new(out/'result.json',result)
    lines=['# Direct canonical PQQ MACE scorer','',f'Primary medium operational gate: **{medium["operational_gate_pass"]}**. Baseline unchanged.',
           '', '| Model | Calibration gap, kcal/mol | Valid calibration rows | Supported transfer cases | Transfer cases lacking calibration | Operational gate |',
           '|---|---:|---:|---:|---:|---|']
    for name,r in result['models'].items():
        c=r['calibration'];lines.append(f'| {name} | {c["gap_kcal_mol"]} | {c["valid_calibration_count"]}/25 | {r["supported_transfer_count"]}/3 | {r["unavailable_calibration_transfer_count"]}/3 | {r["operational_gate_pass"]} |')
    lines+=['','Calibration rows determine only their own checkpoint bands; no transfer row was used to fit or shift a threshold. Larger raw R is La-like. No aquo reference, old threshold, mechanical correction or sign reversal is used. If calibration fails, transfer decisions remain unavailable; zero supported transfer cases then does not mean three wrong classifications.','',
            '| Model | Transfer structure | Expected | Raw R, kcal/mol | Decision |','|---|---|---|---:|---|']
    for name,r in result['models'].items():
        for row in r['rows']:
            if row['evaluation_role']!='calibration':lines.append(f'| {name} | {row["case_id"]} | {row["expected_class"]} | {row["R_kcal_mol"]} | {row["decision"]} |')
    lines+=['','## Scope','',
            'All cases are consumed.1H4I shares accessionP16027 and4MAE sharesI0JWN7 with calibration; these are structural-transfer observations, not new independent proteins. Homologous calibration members are not assumed independent biological replicates.1KB0 is a separate retrospective PQQ class control. The25-member panel is perfectly confounded with motif/core charge; separation does not establish added information beyond composition.','',
            'This is PQQ functional-class evidence. Known non-PQQ affinity failures are not repaired by a successful result here. An incomplete or failed endpoint is retained as unavailable; no baseline or other-checkpoint substitution occurs. Original frozen source H, water exclusions, caps and PQQ state remain unchanged;4MAE retains its declared dry noncore15P vacancy.','',
            'Unrounded endpoint energies, component contrasts, classifications, baseline records, implementation/source pins and every denominator are retained in result.json. Execution receipts and measured costs remain in the source collections and their Slurm accounting records. No production promotion or prospective validation is claimed.']
    (out/'REPORT.md').write_text('\n'.join(lines)+'\n')
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for key in ('medium-mace','medium-gb','large-mace','large-gb','output'):p.add_argument('--'+key,required=True)
    r=report(**vars(p.parse_args()));print(json.dumps({'primary_operational_gate_pass':r['primary_operational_gate_pass'],
         'models':{k:{'calibration':v['calibration'],'supported_transfer_count':v['supported_transfer_count']} for k,v in r['models'].items()}},indent=2))
