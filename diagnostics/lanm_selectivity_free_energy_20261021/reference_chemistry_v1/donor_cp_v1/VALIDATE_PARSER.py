#!/usr/bin/env python3
"""Bounded actual-output replays and larger-center failure controls; zero calls."""
import argparse
import copy
import importlib.util
import json
from pathlib import Path
import sys
import time

HERE=Path(__file__).resolve().parent
sp=importlib.util.spec_from_file_location('donor_cp',HERE/'REFERENCE.py');m=importlib.util.module_from_spec(sp);sp.loader.exec_module(m);b=m.b

def fixture(c):
    lines=['SYNTHETIC PARSER CONTROL ONLY; NO MOLECULAR RESULT','****END OF INPUT****','CARTESIAN COORDINATES (A.U.)']
    for i,(s,z,x) in enumerate(zip(c['symbols'],c['expected_effective_nuclear_charges'],c['positions_A'])):lines.append(f'{i} {s} {z:.4f} 0 1.0 '+' '.join(f'{v/b.BOHR:.6f}' for v in x))
    lines+=['INTERNAL COORDINATES','BASIS SET INFORMATION'];groups={}
    for s in c['symbols']:
        if s not in groups:
            groups[s]=len(groups)+1;contraction='6s4p3d2f' if s in m.METALS else '3s1p' if s=='H' else '5s3p2d1f'
            lines.append(f'Group {groups[s]} Type {s} : primitive set contracted to {contraction}')
    for i,s in enumerate(c['symbols']):lines.append(f'Atom {i}{s} basis set group => {groups[s]}')
    lines+=['ECP PARAMETER INFORMATION']
    if c['expected_ecp_core']:lines += [f'Group 1, Type {c["metal"]} ECP (replacing {c["expected_ecp_core"]} core electrons, lmax=4)',f'Atom 0{c["metal"]} ECP group => 1']
    lines+=['SHARK',f'Number of basis functions ... {c["expected_basis_functions"]}',f'Total Charge Charge .... {c["charge"]}','Multiplicity Mult .... 1',f'Number of Electrons NEL .... {c["expected_electrons"]}',
        'Convergence Check Mode ConvCheckMode .... Total+1el-Energy','SCF CONVERGED AFTER 10 CYCLES','**** Energy Check signals convergence ****','Last Energy change ... 1.0e-10 Tolerance : 1.0e-9']
    if c['method']=='cc':
        lines+=['Frozen core treatment ... NO frozen core','Reference Wavefunction ... RHF',f'Number of correlated electrons ... {c["expected_correlated_electrons"]}',
            'Convergence tolerance (max. residuum) ... 1.000e-08','Iter E(tot) E(Corr) Delta-E Residual Time',
            '-100.000000000 -0.200000000 0.000000001 0.000000006 0.12','--- The Coupled-Cluster iterations have converged ---',
            'COUPLED CLUSTER ENERGY','T1 diagnostic ... 0.007','LARGEST AMPLITUDES',' 1-> 6 -1-> -1 0.030000',' !!!',
            'Triples Correction (T) ... -0.005','E(CCSD(T)) ... -100.000000000']
    else:
        lines+=['Hartree-Fock type HFTyp .... RHF','Exchange Functional Exchange .... PBE','Correlation Functional Correlation .... PBE','Fraction HF Exchange ScalHFX .... 0.250000','RIJ-COSX (HFX calculated with COS-X)).... on','Total Energy : -100.000000000 Eh']
    if c['d4']:lines+=['DFT DISPERSION CORRECTION','DFTD4 V3.4.0','Dispersion correction -0.000100000']
    lines+=['FINAL SINGLE POINT ENERGY '+('-100.000100000' if c['d4'] else '-100.000000000'),'ORCA TERMINATED NORMALLY'];return '\n'.join(lines)+'\n'

def main():
    tick=time.monotonic();ap=argparse.ArgumentParser();ap.add_argument('--panel',required=True);ap.add_argument('--out',type=Path,required=True);a=ap.parse_args();p=m.load(a.panel);a.out.mkdir(parents=True,exist_ok=False);records=[]
    for donor in m.DONORS:
        for metal in m.METALS:
            c=copy.deepcopy(next(x for x in p['cells'] if x['cell_id']==f'cc__{donor}__{metal}__complex'));c.update(method='pbe0',d4=True)
            output=b.verify(p['saved_complexes'][donor][metal]['output']);r=m.parse(output.read_text(),c,0);assert r['status']=='complete',(c['cell_id'],r['checks'])
            records.append({'case':'actual_larger_DFT_'+donor+'_'+metal,'source':b.pin(output),'checks':r['checks'],'basis':r['basis_functions_by_center']})
    water=json.loads(b.verify(p['water_baseline']).read_text());coll=json.loads(b.verify(water['canonical_collection']).read_text())
    for rp in coll['result_pins']:
        r0=json.loads(b.verify(rp).read_text());c=copy.deepcopy(r0['cell']);c.update(method='cc',d4=False);out=b.verify(r0['output']);r=m.parse(out.read_text(),c,r0['returncode']);assert r['status']=='complete',(r0['cell_id'],r['checks'])
        records.append({'case':'actual_canonical_water_'+r0['cell_id'],'source':b.pin(out),'checks':r['checks'],'CC_final_residual':r['CC_residual_trace'][-1][3]})
    for c in p['cells']:
        text=fixture(c);r=m.parse(text,c,0);assert r['status']=='complete',(c['cell_id'],r['checks']);(a.out/(c['cell_id']+'.synthetic.txt')).write_text(text);records.append({'case':'synthetic_positive_'+c['cell_id'],'checks':r['checks'],'synthetic':True})
    c=next(x for x in p['cells'] if x['cell_id']=='cc__acetamide__La__ligand_ghost_ion');text=fixture(c)
    cases=[('last_center_missing','\n'.join(line for line in text.splitlines() if not line.startswith('9 C ')),'all_centers_and_coordinates'),
        ('last_basis_missing','\n'.join(line for line in text.splitlines() if not line.startswith('Atom 9C basis')),'per_center_basis'),
        ('ghost_charge',text.replace('0 La 0.0000','0 La 11.0000'),'ghost_nuclear_charges'),
        ('ghost_ECP',text.replace('ECP PARAMETER INFORMATION','ECP PARAMETER INFORMATION\nAtom 0La ECP group => 1\nGroup 1 Type La ECP replacing 46 core electrons'),'ecp_assignment'),
        ('wrong_total_AO',text.replace('Number of basis functions ... 201','Number of basis functions ... 189'),'basis_function_count'),
        ('frozen_correlated_electrons',text.replace('Number of correlated electrons ... 32','Number of correlated electrons ... 24'),'correlated_electrons'),
        ('implicit_frozen_core',text.replace('NO frozen core','chemical core'),'no_frozen_core'),
        ('looser_CC_target',text.replace('1.000e-08','1.000e-06'),'declared_CC_residual'),
        ('CC_residual_exceeds_target',text.replace('0.000000006','0.000000060'),'declared_CC_residual'),
        ('CC_failed',text.replace('iterations have converged','iterations have NOT converged'),'cc_converged')]
    for name,t,key in cases:
        r=m.parse(t,c,0);assert r['status']=='failed_admission' and not r['checks'][key],name;records.append({'case':name,'expected_rejection':key,'synthetic':True})
    ac=next(x for x in p['cells'] if x['cell_id']=='cc__acetate_bi__Dy__ligand_ghost_ion');t=fixture(ac).replace('Total Charge Charge .... -1','Total Charge Charge .... 0');r=m.parse(t,ac,0);assert not r['checks']['charge'];records.append({'case':'acetate_anion_wrong_charge','expected_rejection':'charge','synthetic':True})
    df=next(x for x in p['cells'] if x['cell_id']=='pbe0__acetamide__La__ligand_ghost_ion');t=fixture(df).replace('FINAL SINGLE POINT ENERGY','DFTD4 V3.4.0\nDispersion correction -0.001\nFINAL SINGLE POINT ENERGY');r=m.parse(t,df,0);assert not r['checks']['dispersion_scope'];records.append({'case':'ghost_D4_rejected','expected_rejection':'dispersion_scope','synthetic':True})
    free=next(x for x in p['cells'] if x['cell_id']=='pbe0__acetate_bi__real_ligand_D4');r=m.parse(fixture(free).replace('Dispersion correction -0.000100000',''),free,0);assert not r['checks']['dispersion_scope'];records.append({'case':'real_anion_D4_required','expected_rejection':'dispersion_scope','synthetic':True})
    receipt={'status':'passed','panel':b.pin(a.panel),'source':b.pin(__file__),'records':records,'molecular_calls':0,'local_threads':1,'elapsed_seconds':time.monotonic()-tick,'command':sys.argv,'utc':b.now(),'donor_canonical_native_output_available':False}
    pin=b.write(a.out/'VALIDATION.json',receipt);print(json.dumps({'status':'passed','records':len(records),'validation':pin,'molecular_calls':0}))

if __name__=='__main__':main()
