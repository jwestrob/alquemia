"""Frozen OMOL calibration and distinct retrospective affinity-direction tests."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import shutil
from affordable_common import HA_TO_KCAL, InvalidArtifact, energy, read_json, record, verify, write_new
from mace_hybrid import EV_TO_KCAL
from mace_omol import validate, collect, benchmark_tasks
from mace_canonical_run import inventory
from mace_canonical_report import fit_bands, decision


def paired_energy(ca, la):
    if ca['status']!='computed' or la['status']!='computed':
        return None
    return (ca['energy_eV']-la['energy_eV'])*EV_TO_KCAL


def baseline_nonpqq(assessment, tasks):
    from mace_curvature import dft_source
    from mace_mechanics_assess import verified_new_dft
    a=read_json(assessment); prepared=read_json(verify(a['sources']['prepared']))
    actual,_=verified_new_dft(verify(a['sources']['DFT']))
    old,_=dft_source(verify(prepared['DFT_reference']))
    source_tasks={t['task_id']:t for t in tasks}; out={}
    for name in ('GGR_extended','GGR_connected','ALPHA_1F6S','ALPHA_6IP9'):
        endpoints={}
        for metal in ('La','Ca'):
            key=f'{name}_center_{metal}'; ref=a['gradient_sources'][key]
            for pin in ref.values():verify(pin)
            if ref['xyz']['sha256']!=source_tasks[name+'_'+metal]['xyz']['sha256']:
                raise InvalidArtifact('baseline center coordinates differ from OMOL core')
            value=energy(verify(ref['output']))
            archived=(old['energies'][prepared['reused_GGR'][key]['source_task_id']]
                      if name.startswith('GGR') else actual['rows'][key])
            if archived['energy_hartree']!=value:
                raise InvalidArtifact('baseline energy differs from actual executed center')
            endpoints[metal]={'energy_hartree':value,'source':ref}
        r=(endpoints['Ca']['energy_hartree']-endpoints['La']['energy_hartree'])*HA_TO_KCAL
        if abs(r-a['rows'][name]['DFT_R_kcal_mol'])>1e-8:
            raise InvalidArtifact('saved baseline paired contrast changed')
        out[name]={'R_kcal_mol':r,'endpoints':endpoints,'protocol':'native_r2scan3c_CPCM_TightSCF_center',
                   'calibrated_class':None,'reference_status':'raw_contrast_only'}
    return out


def report(collection, mechanics_assessment, output):
    saved=read_json(collection);mp=verify(saved['manifest']);validate(mp);current=collect(mp);m=read_json(mp)
    if saved!=current or m['stage']!='benchmark':
        raise InvalidArtifact('actual current benchmark collection required')
    data=inventory(verify(m['inventory']));tasks=benchmark_tasks(data,verify(m['mechanics']))
    rows=dict(current['rows']);rows.update({k:v['result'] for k,v in current['reused_rows'].items()})
    if len(rows)!=64:raise InvalidArtifact('all64 logical endpoints must remain in the report')
    task_lookup={t['task_id']:t for t in tasks};baseline=baseline_nonpqq(mechanics_assessment,tasks)
    scores=[]
    for source in data['rows']:
        name=source['case_id'];ca,la=(rows[name+'_'+metal] for metal in ('Ca','La'))
        value=paired_energy(ca,la)
        scores.append({'case_id':name,'expected_class':source['expected_class'],
                       'evaluation_role':source['evaluation_role'],'evidence_stratum':source['evidence_stratum'],
                       'sequence_accession_group':source['sequence_accession_group'],'prospectively_blind':False,
                       'status':'computed' if value is not None else 'unavailable','R_kcal_mol':value,
                       'endpoints':{'Ca':ca,'La':la},'baseline':source['baseline']})
    bands=fit_bands(scores)
    for row in scores:
        row['decision']=decision(row['R_kcal_mol'],bands['decision_bands'])
        row['expected_region_pass']=row['decision']==row['expected_class']+'-supported'
    transfer=[r for r in scores if r['evaluation_role']!='calibration']
    if {r['case_id'] for r in transfer}!={'1H4I','4MAE','1KB0'}:
        raise InvalidArtifact('transfer case selection changed')
    affinity={}
    for name in ('GGR_extended','GGR_connected','ALPHA_1F6S','ALPHA_6IP9'):
        ca,la=(rows[name+'_'+metal] for metal in ('Ca','La'));t=task_lookup[name+'_La']
        value=paired_energy(ca,la)
        affinity[name]={'R_kcal_mol':value,'status':'computed' if value is not None else 'unavailable',
                        'endpoints':{'Ca':ca,'La':la},'evidence':t['evidence'],'preparation':t['preparation'],
                        'baseline':baseline[name],'calibrated_class':None,'prospectively_blind':False}
    contrasts=[]
    for ggr in ('GGR_extended','GGR_connected'):
        for alpha in ('ALPHA_1F6S','ALPHA_6IP9'):
            av,gv=(affinity[x]['R_kcal_mol'] for x in (alpha,ggr))
            value=av-gv if av is not None and gv is not None else None
            contrasts.append({'positive_case':alpha,'negative_case':ggr,
                              'role':'primary' if ggr=='GGR_extended' else 'representation_robustness',
                              'delta_R_kcal_mol':value,'required_min_kcal_mol':.02,
                              'pass':value is not None and value>.02,
                              'baseline_delta_R_kcal_mol':baseline[alpha]['R_kcal_mol']-baseline[ggr]['R_kcal_mol']})
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False)
    implementation=out/Path(__file__).name;shutil.copyfile(__file__,implementation)
    result={'status':'complete','scoring_status':current['status'],'protocol_id':m['protocol_id'],'implementation':record(implementation),
            'sources':{'OMOL':record(collection),'mechanics_baseline':record(mechanics_assessment),
                       'qualification':m['qualification'],'inventory':m['inventory'],'agreement':m['agreement']},
            'model':m['model'],'software':m['software'],'calibration':bands,'canonical_rows':scores,
            'canonical_operational_gate_pass':bands['status']=='calibratable' and all(r['expected_region_pass'] for r in transfer),
            'transfer_count':3,'supported_transfer_count':sum(r['expected_region_pass'] for r in transfer),
            'unavailable_calibration_transfer_count':sum(r['decision']=='unavailable_calibration' for r in transfer),
            'inconclusive_transfer_count':sum(r['decision']=='inconclusive' for r in transfer),
            'nonPQQ_rows':affinity,'nonPQQ_contrasts':contrasts,
            'nonPQQ_primary_gate_pass':all(r['pass'] for r in contrasts if r['role']=='primary'),
            'nonPQQ_robustness_gate_pass':all(r['pass'] for r in contrasts),
            'broad_affinity_validated':False,'incremental_information_beyond_composition_established':False,
            'training_overlap_absence_claimed':False,'baseline_changed':False,'solution_score':None,
            'solvent_status':'not_part_of_vacuum_descriptor','relaxation_correction_kcal_mol':None,
            'entropy_correction_kcal_mol':None}
    write_new(out/'result.json',result)
    text=['# MACE-OMOL descriptor benchmark','',
          f'Canonical gate: **{result["canonical_operational_gate_pass"]}**. Separate non-PQQ primary/robustness gates: **{result["nonPQQ_primary_gate_pass"]}/{result["nonPQQ_robustness_gate_pass"]}**.',
          '',f'Calibration gap: {bands["gap_kcal_mol"]} kcal/mol; {bands["valid_calibration_count"]}/25 valid. Supported transfer: {result["supported_transfer_count"]}/3; calibration unavailable for {result["unavailable_calibration_transfer_count"]}/3.',
          '', '| Transfer | Expected | Raw R, kcal/mol | Decision |','|---|---|---:|---|']
    for r in transfer:text.append(f'| {r["case_id"]} | {r["expected_class"]} | {r["R_kcal_mol"]} | {r["decision"]} |')
    text+=['','| Non-PQQ relative comparison | OMOL delta R, kcal/mol | Baseline delta R | Direction passes |','|---|---:|---:|---|']
    for r in contrasts:text.append(f'| {r["positive_case"]} minus {r["negative_case"]} | {r["delta_R_kcal_mol"]} | {r["baseline_delta_R_kcal_mol"]} | {r["pass"]} |')
    text+=['','All cases are retrospective. Both crystals share calibration accessions; both alpha structures form one biological group and both GGR representations form another. Canonical functional class and qualified non-PQQ affinity evidence remain separate. The canonical panel is confounded with motif/composition. No broad affinity, absent training overlap, or incremental physics claim follows.',
           '', 'These are vacuum electronic descriptors, not solvent-corrected binding energies. OMOL predicts no atomic density here; no POLAR charge distribution, CPCM energy, baseline band, aquo offset, relaxation or entropy has been added. Default baseline is unchanged. Source H, caps, waters and atom mappings are preserved exactly.',
           '', 'Unrounded endpoints, failures, original baseline comparators and execution/source pins remain in result.json. Numerical qualification and measured allocation costs are separate from predictive evidence.']
    (out/'REPORT.md').write_text('\n'.join(text)+'\n')
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('collection','mechanics-assessment','output'):p.add_argument('--'+k,required=True)
    r=report(**vars(p.parse_args()))
    print(json.dumps({k:r[k] for k in ('calibration','canonical_operational_gate_pass','nonPQQ_primary_gate_pass','nonPQQ_robustness_gate_pass','nonPQQ_contrasts')},indent=2))
