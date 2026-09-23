"""Inspect actual Q4W6G0 scalar outputs; no molecular calculations."""
from pathlib import Path
import argparse,json,re,sys
import numpy as np
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'scripts'))
from affordable_common import read_json,record,verify,write_new,xyz
from accommodation_folds_compare import decision


def audit(comparison,output):
    result=read_json(comparison);cid='q4w6g0-pqq-la_model__conditioned_Ca__seed-1_sample-2'
    row=next(r for r in result['rows'] if r['case_id']==cid);records={}
    for version,pin in [('old',row['source_collection']),('new',row['precision_collection'])]:
        c=next(c for c in read_json(verify(pin))['cases'] if c['case_id']==cid);cells=[]
        for z in ('Ca','La'):
            for q in ('adaptive_Ca','adaptive_La'):
                cell=c['matrix'][z][c['aliases'][q]['representative']]
                for medium in ('vacuum','alpb'):
                    low=cell['low'][medium];text=verify(low['output']).read_text();receipt=read_json(verify(low['receipt']))
                    lm=read_json(verify(low['manifest']));task=next(t for t in lm['tasks'] if t['task_id']==low['task_id'])
                    cycles=[int(x) for x in re.findall(r'SCF CONVERGED AFTER\s+(\d+)\s+CYCLES',text)]
                    final=[float(x) for x in re.findall(r'FINAL SINGLE POINT ENERGY\s+([-\d.]+)',text)]
                    if len(cycles)!=1 or len(final)!=1 or final[0]!=low['energy_hartree']:raise ValueError('actual convergence/energy marker differs')
                    a=low['audit']
                    cells.append({'metal':z,'candidate':q,'medium':medium,'energy_hartree':low['energy_hartree'],
                        'SCF_cycles':cycles[0],'normal_termination':receipt['normal_termination'],'scf_converged':receipt['scf_converged'],
                        'returncode':receipt['returncode'],'parallelism':receipt['parallelism'],'xyz':cell['xyz'],'input':task['input'],
                        'output':low['output'],'receipt':low['receipt'],'parameter_export':a['parameter_export'],
                        'native_valence_electrons':a['native_valence_electrons'],'charge_sum_e':a['charge_sum_e'],
                        'metal_mulliken_charge_e':a['mulliken_charges_e'][0],'PMIX_diagnostic_count':a['PMIX_diagnostic_count'],
                        'printed_SCF_components_hartree':a['printed_SCF_components_hartree'],
                        'runtime_line':next((line for line in text.splitlines() if line.startswith('TOTAL RUN TIME:')),None)})
        records[version]={'collection':pin,'cells':cells,'selection':c['pool']['rows'],
            'R':c['pool']['operational']['composite_R_model_kcal_mol'],
            'reference_calls':{ref:decision(c['pool']['operational']['composite_R_model_kcal_mol'],result['bands'][ref]) for ref in ('union_adaptive','union_precision')}}
    old={(r['metal'],r['candidate'],r['medium']):r for r in records['old']['cells']};checks=[]
    for r in records['new']['cells']:
        key=(r['metal'],r['candidate'],r['medium']);a=old[key]
        checks.append({'metal':key[0],'candidate':key[1],'medium':key[2],
            'same_input_bytes':verify(a['input']).read_bytes()==verify(r['input']).read_bytes(),
            'same_parameter_export_bytes':verify(a['parameter_export']).read_bytes()==verify(r['parameter_export']).read_bytes(),
            'maximum_coordinate_delta_A':float(np.linalg.norm(np.array([x[1:] for x in xyz(verify(a['xyz']))])-np.array([x[1:] for x in xyz(verify(r['xyz']))]),axis=1).max())})
    out={'comparison':record(comparison),'case_id':cid,'versions':records,'pair_checks':checks,
        'numerical_comparison':row['precision_numerical_comparison'],'production_changed':False,'new_molecular_calls':0,
        'interpretation':'Change localized to La-at-adaptiveCa vacuum; both SCFs converged with identical recipe/parameter exports. Geometry and parallelism both changed, so their separate causal roles remain unresolved. No retry, seed selection or threshold change.',
        'implementation':record(__file__)}
    write_new(output,out);print(json.dumps({v:{'R':r['R'],'calls':r['reference_calls']} for v,r in records.items()},indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('comparison','output'):p.add_argument('--'+k,required=True,type=Path)
    audit(**vars(p.parse_args()))
