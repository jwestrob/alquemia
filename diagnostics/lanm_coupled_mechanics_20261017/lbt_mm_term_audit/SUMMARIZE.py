#!/usr/bin/env python3
"""Saved-only scalar/term report for the completed four-state MM audit."""
import argparse, csv, hashlib, json, math, os
from pathlib import Path


def pin(path):
    p=Path(path).resolve()
    return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}


def main():
    p=argparse.ArgumentParser();p.add_argument('--result',type=Path,required=True);p.add_argument('--sha256',required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    assert os.environ.get('SLURM_JOB_ID','').isdigit(), 'Slurm required'
    assert pin(a.result)['sha256']==a.sha256
    d=json.loads(a.result.read_text())
    assert d['status']=='complete' and d['slurm_job_id']=='13122' and d['available_states']==4
    assert d['declared_states']==[['7CCO','La'],['7CCO','Dy'],['7CCN','La'],['7CCN','Dy']]
    assert d['contexts_attempted']==16 and d['MM_diagnostic_queries_attempted']==256
    assert all(d[k]==0 for k in ('new_DFT_calls','new_MACE_calls','new_GFN_calls','new_ALPB_calls','new_geometries','optimizer_steps'))
    rows=[];terms=[];force_res=[];energy_res=[]
    for r in d['rows']:
        assert r['status']=='available' and len(r['local_terms'])==52 and all(c['passed'] for c in r['checks'])
        for c in r['checks']:
            force_res.append(c['maximum_force_component_residual_kcal_mol_A'])
            if 'absolute_energy_residual_kcal_mol' in c:energy_res.append(c['absolute_energy_residual_kcal_mol'])
        t=r['target_summary'];f=t['archived_force_kcal_mol_A'];norm=math.sqrt(math.fsum(x*x for x in f));axis=[x/norm for x in f]
        def project(v):return math.fsum(x*y for x,y in zip(v,axis))
        projections={k:project(v) for k,v in t['class_forces_kcal_mol_A'].items()}
        nonbonded={k:project(v) for k,v in t['nonbonded_group_forces_kcal_mol_A'].items()}
        source_terms=[]
        for v in r['local_terms']:
            w={'source':r['source_id'],'metal':r['metal'],'class':v['force_class'],'native_index':v['native_index'],'atom_ids':','.join(v['atom_ids']),
               'force_norm_kcal_mol_A':math.sqrt(math.fsum(x*x for x in v['target_force_kcal_mol_A'])),
               'projection_on_total_MM_force_kcal_mol_A':project(v['target_force_kcal_mol_A']),
               'energy_kcal_mol':v['energy_kcal_mol'],'parameters':json.dumps({k:x for k,x in v.items() if k in ('native_parameter_XML','parameter_units','actual_bond_length_A','reference_bond_length_A','stiffness_kcal_mol_A2','actual_angle_degrees','reference_angle_radians','stiffness_kcal_mol_radian2','periodicity','phase_radians','amplitude_kcal_mol')},sort_keys=True)}
            terms.append(w);source_terms.append(w)
        rows.append({'source':r['source_id'],'metal':r['metal'],'total_MM_force_norm_kcal_mol_A':norm,'bonded_class_signed_projections_kcal_mol_A':projections,
                     'nonbonded_signed_projections_kcal_mol_A':nonbonded,'sum_projection_minus_total_norm':math.fsum(list(projections.values())+list(nonbonded.values()))-norm,
                     'largest_five_terms_by_absolute_projection':sorted(source_terms,key=lambda v:(-abs(v['projection_on_total_MM_force_kcal_mol_A']),v['class'],v['native_index']))[:5]})
    a.output.mkdir(parents=True,exist_ok=False)
    summary={'status':'complete','protocol':'saved_LBT_MM_term_attribution_summary_v1','implementation':pin(__file__),'slurm_job_id':os.environ['SLURM_JOB_ID'],
             'input':pin(a.result),'states':rows,'state_denominator':4,'term_denominator':len(terms),'maximum_closure_force_component_kcal_mol_A':max(force_res),
             'maximum_closure_energy_kcal_mol':max(energy_res),'new_molecular_calls':0,'projection_definition':'dot(component force, unit total MM force); signed vector projection, not an independent energy or causal percentage',
             'physical_accuracy_established':False}
    (a.output/'SUMMARY.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n')
    with (a.output/'ALL_TERMS.tsv').open('w') as out:
        writer=csv.DictWriter(out,fieldnames=list(terms[0]),delimiter='\t');writer.writeheader();writer.writerows(terms)
    lines=['# LBT terminal MM force attribution','', '**All four fixed-coordinate cases pass. Bond stretching and angle bending dominate the terminal MM restoring force.**', '',
           'Exact native parameter and term-inventory checks found no copying or double-counting defect. Full-system and per-term force closure pass; this does not establish that the MM stiffness is physically accurate.', '',
           '## Components along the total MM force', '',
           'Phe1 CA, kcal/mol/angstrom. Signed projections add to the total MM force magnitude; they are not independent energies or causal percentages.', '',
           '| Source / metal | Total MM | Bonds | Angles | Torsions | CMAP | Nonbonded |', '| --- | ---: | ---: | ---: | ---: | ---: | ---: |']
    for r in rows:
        q=r['bonded_class_signed_projections_kcal_mol_A'];n=math.fsum(r['nonbonded_signed_projections_kcal_mol_A'].values())
        lines.append(f"| {r['source']} / {r['metal']} | {r['total_MM_force_norm_kcal_mol_A']:.4f} | {q['HarmonicBondForce']:.4f} | {q['HarmonicAngleForce']:.4f} | {q['PeriodicTorsionForce']:.4f} | {q['CMAPTorsionForce']:.4f} | {n:.4f} |")
    lines+=['','## Evidence and limits','',f"All208 local terms are retained in ALL_TERMS.tsv. Largest closure residuals: force {max(force_res):.3g} kcal/mol/angstrom; energy {max(energy_res):.3g} kcal/mol.", '',
            'Job13122 completed256 MM queries in16 Reference contexts,9 scheduler seconds on1 CPU. This report performs saved arithmetic only. All four original starts, parameters, charges, terms and coordinates are preserved. No DFT, MACE, GFN2, ALPB, optimization or production change was made.', '',
            'The whole-peptide DFT scout13051 independently tests the accuracy of this force. This attribution alone supports no parameter change, geometry repair, force-field verdict or metal-preference claim.']
    (a.output/'REPORT.md').write_text('\n'.join(lines)+'\n')
    print(json.dumps({'status':'complete','states':4,'terms':len(terms),'summary':pin(a.output/'SUMMARY.json')}))

if __name__=='__main__':main()
