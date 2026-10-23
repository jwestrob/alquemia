#!/usr/bin/env python3
"""Real-output replays and explicit synthetic failure controls, zero native calls."""
import argparse
import copy
import importlib.util
import json
from pathlib import Path
import sys
import time

HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('dft_cp',HERE/'REFERENCE.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
b=m.b

def fixture(c):
    lines=['SYNTHETIC PARSER CONTROL, NOT A MOLECULAR RESULT','****END OF INPUT****','CARTESIAN COORDINATES (A.U.)']
    for i,(s,z,x) in enumerate(zip(c['symbols'],c['expected_effective_nuclear_charges'],c['positions_A'])):lines.append(f'{i} {s} {z:.4f} 0 1.0 '+' '.join(f'{v/b.BOHR:.6f}' for v in x))
    lines+=['INTERNAL COORDINATES','BASIS SET INFORMATION'];groups={}
    for s in c['symbols']:
        if s not in groups:
            groups[s]=len(groups)+1; contraction={'La':'6s4p3d2f','Dy':'6s4p3d2f','O':'5s3p2d1f','H':'3s1p'}[s]
            lines.append(f'Group {groups[s]} Type {s} : primitive set contracted to {contraction}')
    for i,s in enumerate(c['symbols']):lines.append(f'Atom {i}{s} basis set group => {groups[s]}')
    lines+=['ECP PARAMETER INFORMATION']
    if c['expected_ecp_core']:lines += [f'Group 1, Type {c["metal"]} ECP (replacing {c["expected_ecp_core"]} core electrons, lmax=4)',f'Atom 0{c["metal"]} ECP group => 1']
    lines+=['SHARK',f'Number of basis functions ... {c["expected_basis_functions"]}',f'Total Charge Charge .... {c["charge"]}','Multiplicity Mult .... 1',f'Number of Electrons NEL .... {c["expected_electrons"]}','Hartree-Fock type HFTyp .... RHF','Exchange Functional Exchange .... PBE','Correlation Functional Correlation .... PBE','Fraction HF Exchange ScalHFX .... 0.250000','RIJ-COSX (HFX calculated with COS-X)).... on','SCF CONVERGED AFTER 10 CYCLES','Total Energy : -100.000000000 Eh']
    if c['d4']:lines+=['DFT DISPERSION CORRECTION','DFTD4 V3.4.0','Dispersion correction -0.000100000']
    lines+=['FINAL SINGLE POINT ENERGY '+('-100.000100000' if c['d4'] else '-100.000000000'),'ORCA TERMINATED NORMALLY']
    return '\n'.join(lines)+'\n'

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--panel',required=True);ap.add_argument('--out',type=Path,required=True);a=ap.parse_args();tick=time.monotonic();p=m.load(a.panel);a.out.mkdir(parents=True,exist_ok=False);records=[]
    original=m.cp.load(m.PARENT)
    for metal in ('La','Dy'):
        c=copy.deepcopy(next(x for x in original['cells'] if x['cell_id']==metal+'__complex'));c['d4']=True
        path=b.verify(p['saved_complexes'][metal]['output']);r=m.parse(path.read_text(),c,0);assert r['status']=='complete',r['checks'];records.append({'case':'actual_DFT_complex_'+metal,'source':b.pin(path),'checks':r['checks']})
    for c in p['cells']:
        text=fixture(c);r=m.parse(text,c,0);assert r['status']=='complete',(c['cell_id'],r['checks']);(a.out/(c['cell_id']+'.synthetic.txt')).write_text(text);records.append({'case':'synthetic_positive_'+c['cell_id'],'checks':r['checks'],'synthetic':True})
    c=next(x for x in p['cells'] if x['cell_id']=='La__water_ghost_ion');t=fixture(c)
    mutations=[('ghost_charge',t.replace('0 La 0.0000','0 La 11.0000'),'ghost_nuclear_charges',0),
        ('ghost_ECP',t.replace('ECP PARAMETER INFORMATION','ECP PARAMETER INFORMATION\nAtom 0La ECP group => 1\nGroup 1 Type La ECP replacing 46 core electrons'),'ecp_assignment',0),
        ('ghost_D4',t.replace('FINAL SINGLE POINT ENERGY','DFTD4 V3.4.0\nDispersion correction -0.001\nFINAL SINGLE POINT ENERGY'),'dispersion_scope',0),
        ('basis_deleted',t.replace('Number of basis functions ... 90','Number of basis functions ... 43'),'basis_function_count',0),
        ('wrong_spin',t.replace('Multiplicity Mult .... 1','Multiplicity Mult .... 6'),'multiplicity',0),
        ('wrong_electrons',t.replace('Number of Electrons NEL .... 10','Number of Electrons NEL .... 8'),'total_electrons',0),
        ('wrong_exchange',t.replace('0.250000','0.200000'),'exact_exchange_25_percent',0),
        ('wrong_correlation',t.replace('Correlation .... PBE','Correlation .... LYP'),'PBE_correlation',0),
        ('missing_COSX',t.replace('RIJ-COSX (HFX calculated with COS-X)).... on',''),'RIJCOSX',0),
        ('SCF_failed',t.replace('SCF CONVERGED AFTER','SCF NOT CONVERGED AFTER'),'scf_converged',0),
        ('energy_closure',t.replace('FINAL SINGLE POINT ENERGY -100.000000000','FINAL SINGLE POINT ENERGY -99.000000000'),'electronic_final_closure',0),
        ('abnormal_exit',t.replace('ORCA TERMINATED NORMALLY','ERROR'),'normal_termination',0),('process_error',t,'process_success',1)]
    water=next(x for x in p['cells'] if x['role']=='real_water');wt=fixture(water)
    for name,text,key,rc in mutations:
        r=m.parse(text,c,rc);assert r['status']=='failed_admission' and not r['checks'][key],name;records.append({'case':name,'expected_rejected':key,'synthetic':True})
    for name,text,key in [('real_water_D4_missing',wt.replace('Dispersion correction -0.000100000',''),'dispersion_scope'),('real_water_inventory',wt.replace('Atom 0O basis set group => 1',''),'per_center_basis')]:
        r=m.parse(text,water,0);assert r['status']=='failed_admission' and not r['checks'][key],name;records.append({'case':name,'expected_rejected':key,'synthetic':True})
    receipt={'status':'passed','panel':b.pin(a.panel),'source':b.pin(__file__),'records':records,'elapsed_seconds':time.monotonic()-tick,'local_threads':1,'molecular_calls':0,'command':sys.argv,'utc':b.now(),'real_ghost_DFT_output_available':False}
    pin=b.write(a.out/'VALIDATION.json',receipt);print(json.dumps({'status':'passed','records':len(records),'validation':pin}))

if __name__=='__main__':main()
