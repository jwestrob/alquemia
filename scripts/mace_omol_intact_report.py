"""Report the declared intact-chain descriptor without fitting a decision band."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import shutil
from affordable_common import InvalidArtifact,read_json,record,verify,write_new
from mace_hybrid import EV_TO_KCAL
from mace_global_benchmark import CONTRASTS
from mace_global_prepare import CASES
from mace_omol_intact import PROTOCOL,validate,collect,preparations
from mace_omol_readout import checked_arrays


def report(collection,output):
    saved=read_json(collection);mp=verify(saved['manifest']);m=read_json(mp)
    if m['stage']!='intact_benchmark':raise InvalidArtifact('intact benchmark required')
    validate(mp);actual=collect(mp)
    if saved!=actual:raise InvalidArtifact('collection no longer matches actual receipts')
    physical=preparations(verify(m['physical_preparation']))
    rows=dict(actual['rows']);rows.update({key:r['result'] for key,r in actual['reused_rows'].items()})
    full=read_json(verify(read_json(verify(m['full_qualification']))['manifest']))
    tasks={t['task_id']:t for t in m['tasks']+full['tasks']}
    scores={}
    for case in CASES:
        endpoints={};p=physical[case];i=p['metal_index']
        for metal in ('La','Ca'):
            pair={position:rows[f'{case}_{metal}_{position}_primary'] for position in ('bound','detached')}
            complete=all(r['status']=='computed' for r in pair.values());components=None
            difference=None
            if complete:
                difference=pair['bound']['energy_eV']-pair['detached']['energy_eV']
                arrays={pos:checked_arrays(r,tasks[f'{case}_{metal}_{pos}_primary']) for pos,r in pair.items()}
                # These are native model bookkeeping terms, not unique physical atomic energies.
                dn=arrays['bound']['node_energy_eV']-arrays['detached']['node_energy_eV']
                de=arrays['bound']['embedding_energy_eV']-arrays['detached']['embedding_energy_eV']
                dr=arrays['bound']['atomic_reference_eV']-arrays['detached']['atomic_reference_eV']
                components={'native_node_sum_eV':float(dn.sum()),'native_embedding_sum_eV':float(de.sum()),
                            'atomic_reference_sum_eV':float(dr.sum()),'metal_native_node_eV':float(dn[i]),
                            'other_native_node_eV':float(dn.sum()-dn[i]),
                            'accounting_error_kcal_mol':float((dn.sum()+de.sum()-difference)*EV_TO_KCAL)}
            endpoints[metal]={'status':'computed' if complete else 'unavailable','bound_minus_detached_eV':difference,
                              'bound':pair['bound'],'detached':pair['detached'],'components':components}
        available=all(e['status']=='computed' for e in endpoints.values())
        scores[case]={'case_id':case,'status':'computed' if available else 'unavailable','endpoints':endpoints,
                      'R_coord_kcal_mol':None if not available else
                          (endpoints['Ca']['bound_minus_detached_eV']-endpoints['La']['bound_minus_detached_eV'])*EV_TO_KCAL,
                      'atoms':len(p['data']['physical_atoms']),'metal_index':i,'evidence':p['data']['evidence'],
                      'assembly':p['data']['assembly'],'preparation':p['source'],
                      'biological_group':'ALPHA' if case.startswith('ALPHA') else case,
                      'calibrated_class':None,'decision_status':'no_absolute_calibration'}
    contrasts=[]
    for positive,negative in CONTRASTS:
        a,b=(scores[key]['R_coord_kcal_mol'] for key in (positive,negative));v=None if a is None or b is None else a-b
        contrasts.append({'positive_case':positive,'negative_case':negative,'delta_R_kcal_mol':v,
                          'required_min_kcal_mol':.02,'pass':v is not None and v>.02,
                          'evidence_stratum':'PQQ_functional_class' if positive.startswith('PQQ') else 'qualified_affinity_order'})
    result={'status':actual['status'],'protocol_id':PROTOCOL,'collection':record(collection),
            'agreement':m['agreement'],'physical_preparation':m['physical_preparation'],'model':m['model'],
            'energy_evaluation':m['energy_evaluation'],'formula':'(E_bound_Ca-E_detached_Ca)-(E_bound_La-E_detached_La)',
            'numerical_gate_pass':actual['numerical_gate_pass'],'scores':scores,'contrasts':contrasts,
            'predictive_development_gate_pass':actual['numerical_gate_pass'] and all(c['pass'] for c in contrasts),
            'unavailable_score_count':sum(r['status']!='computed' for r in scores.values()),
            'evidence_use':'consumed_retrospective_development','independent_affinity_pairs':1,
            'baseline_changed':False,'broad_affinity_validated':False,'binding_free_energy_kcal_mol':None,
            'solvent_correction_kcal_mol':None,'reference':None,'calibration':None,
            'force_validation_status':'not_requested_energy_only'}
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False)
    shutil.copyfile(__file__,out/Path(__file__).name);result['report_implementation']=record(out/Path(__file__).name)
    write_new(out/'result.json',result)
    lines=['# Intact-chain OMOL coordination descriptor','',
           f'Numerical gate: {result["numerical_gate_pass"]}. Three-comparison development gate: {result["predictive_development_gate_pass"]}.',
           '', '| Case | Atoms | Coordination contrast, kcal/mol |','|---|---:|---:|']
    for r in scores.values():lines.append(f'| {r["case_id"]} | {r["atoms"]} | {r["R_coord_kcal_mol"]} |')
    lines+=['','| Frozen comparison | Difference, kcal/mol | Pass (>0.02) |','|---|---:|---|']
    for c in contrasts:lines.append(f'| {c["positive_case"]} minus {c["negative_case"]} | {c["delta_R_kcal_mol"]} | {c["pass"]} |')
    lines+=['','Whole prepared chain A, fixed physical preparation, no synthetic carve caps. Native float64 energy-only inference; no forces or relaxation are claimed.',
            'The finite-range model does not provide full long-range electrostatics, solvent or certified separated ionic states. System sizes extrapolate beyond reported training sizes.',
            'All cases were already consumed. Two alpha structures represent one biological group and one qualified alpha/GGR affinity comparison. PQQ functional class is a separate evidence stratum.',
            'No absolute decision bands, aquo reference, physical binding-free-energy claim or production change. Native component values are model bookkeeping, not unique physical atomic energies.']
    (out/'REPORT.md').write_text('\n'.join(lines)+'\n')
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for key in ('collection','output'):p.add_argument('--'+key,required=True)
    r=report(**vars(p.parse_args()))
    print(json.dumps({k:r[k] for k in ('status','numerical_gate_pass','predictive_development_gate_pass','contrasts')},indent=2))
