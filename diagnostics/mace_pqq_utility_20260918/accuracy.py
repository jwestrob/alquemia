#!/usr/bin/env python3
"""Replay the fixed, consumed 25+3 PQQ references without new calculations."""
import argparse
import csv
import json
from pathlib import Path
import re
from benchmark import read, record, verify, write, dft_score, classify


def output_energy(output, receipt):
    path=verify(output);r=read(verify(receipt))
    if not (r['returncode']==0 and r['normal_termination'] and r['scf_converged']
            and r['artifacts']['output']==output):
        raise ValueError('incomplete/mismatched archived DFT endpoint')
    values=re.findall(r'FINAL SINGLE POINT ENERGY\s+([-+0-9.]+)',path.read_text())
    if len(values)!=1:raise ValueError('ambiguous final energy')
    return float(values[0])


def inventory(root, output):
    root=Path(root).resolve();out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False)
    base=root/'diagnostics/pqq_pmdh_fixed_core_calibration_20260914/result.json'
    crystals=root/'diagnostics/pqq_pmdh_fixed_core_calibration_20260914/reserved_crystal_holdout/result/holdout_result.json'
    external=root/'diagnostics/baseline_benchmark_20260915/RESULT.json'
    mace=root/'workspaces/mace_omol_20260917/masked_calibration_v1/reference.json'
    b=read(base);m=read(mace);h=read(crystals)
    dft={}
    for r in b['scores']+h['scores']:
        energies={metal:output_energy(a['output'],a['execution']) for metal,a in r['artifacts'].items()}
        if energies!=r['energies_hartree']:raise ValueError('stored DFT energy differs from output')
        s=dft_score(energies,b)
        if abs(s['score']-r['S_aquo_gauge_kcal_mol'])>1e-9:raise ValueError('DFT reporting gauge differs')
        dft[r.get('panel_id',r.get('pdb_id'))]={**s,'status':'computed'}
    e=next(r for r in read(external)['scores'] if r['case']=='pqq_1kb0' and r['lane']=='fixed_core')
    energies={metal:output_energy(a['output'],a['receipt']) for metal,a in e['endpoints'].items()}
    s=dft_score(energies,b)
    if abs(s['score']-e['score']['S_kcal_mol'])>1e-9:raise ValueError('external reporting gauge differs')
    dft['1KB0']={**s,'status':'computed'}
    rows=[]
    for r in m['scores']:
        case=r['case_id'];d=dft.pop(case);score=r['two_call_model_kcal']
        if score is not None:
            bound,free=r['bound_model_eV'],r['disconnected_atom_model_eV']
            exact=(bound['Ca']-bound['La']-(free['Ca']-free['La']))*23.06054783061903
            if score!=exact:raise ValueError('two-call score differs from archived algebra')
            label=classify(score,m['bands']['Ca_max_inclusive_model_kcal'],m['bands']['La_min_inclusive_model_kcal'])
            if label!=r['calibrated_class']:raise ValueError('frozen MACE class differs')
        else:label=None
        rows.append({'case_id':case,'expected_class':r['expected_class'],
                     'sequence_group':r['sequence_accession_group'],'evaluation_role':r['evaluation_role'],
                     'evidence_stratum':r['evidence_stratum'],'prospectively_blind':False,
                     'DFT_S_kcal_mol':d['score'],'DFT_class':d['class'],
                     'MACE_R_model_kcal':score,'MACE_class':label,'MACE_preparation':r['preparation_status'],
                     'DFT_correct':d['class']==r['expected_class'],
                     'MACE_correct':label==r['expected_class']})
    if dft or len(rows)!=28:raise ValueError('28-reference inventory differs')
    counts={}
    for name,chosen in [('calibration',[r for r in rows if r['evaluation_role']=='calibration']),
                        ('transfer',[r for r in rows if r['evaluation_role']!='calibration'])]:
        counts[name]={'total':len(chosen),'DFT_correct':sum(r['DFT_correct'] for r in chosen),
                      'MACE_correct':sum(r['MACE_correct'] for r in chosen),
                      'MACE_unavailable':sum(r['MACE_class'] is None for r in chosen)}
    result={'status':'complete','source_records':list(map(record,(base,crystals,external,mace))),
            'counts':counts,'rows':rows,'independent_sequence_groups':len({r['sequence_group'] for r in rows}),
            'DFT_protocol':b['protocol_id'],'MACE_protocol':m['protocol_id'],
            'DFT_bands':b['calibration']['released_supported_bands'],'MACE_bands':m['bands'],
            'new_energy_evaluations':0,'baseline_changed':False,
            'limits':['Consumed calibration/transfer replay, not prospective accuracy.',
                      '1H4I and 4MAE repeat sequences in the canonical panel.',
                      'Sequence groups are not independent protein-family observations.',
                      'Motif/core composition already separates the canonical panel.',
                      'Functional PQQ class is not quantitative relative affinity.',
                      'DFT kcal/mol and MACE model kcal are distinct scales.']}
    write(out/'result.json',result)
    with (out/'scores.tsv').open('x') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
    lines=['# Frozen PQQ reference fidelity','',
           '| Set | DFT correct / total | Masked MACE correct / total | MACE unavailable |',
           '|---|---:|---:|---:|']
    for name,c in counts.items():
        lines.append(f"| {name} | {c['DFT_correct']}/{c['total']} | {c['MACE_correct']}/{c['total']} | {c['MACE_unavailable']} |")
    lines+=['','1KB0 remains unsupported by whole-chain MACE. It is valid in the fixed-core DFT protocol.',
            '',*['- '+x for x in result['limits']],'','Measured timing comparison is separate and pending.']
    (out/'REPORT.md').write_text('\n'.join(lines)+'\n')
    return {'counts':counts,'report':record(out/'REPORT.md'),'result':record(out/'result.json')}


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',required=True);p.add_argument('--output',required=True)
    print(json.dumps(inventory(**vars(p.parse_args())),indent=2))
