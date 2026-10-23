#!/usr/bin/env python3
"""Native-format replay plus deliberately synthetic failure controls; zero calls."""
import argparse
import importlib.util
import json
from pathlib import Path
import re
import sys
import time

HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('cp_reference',HERE/'REFERENCE.py')
cp=importlib.util.module_from_spec(spec);spec.loader.exec_module(cp)
base=cp.base
ROOT=cp.ROOT

def require(ok,message):
    if not ok:raise AssertionError(message)

def control_fixture(c):
    # Every numeric value below is test data, not a molecular result.
    lines=['SYNTHETIC PARSER CONTROL ONLY','****END OF INPUT****','CARTESIAN COORDINATES (A.U.)',
        'NO LB ZA FRAG MASS X Y Z']
    for i,(symbol,z,x) in enumerate(zip(c['symbols'],c['expected_effective_nuclear_charges'],c['positions_A'])):
        lines.append(f'{i} {symbol} {z:.4f} 0 1.000 '+' '.join(f'{v/base.BOHR:.6f}' for v in x))
    lines+=['INTERNAL COORDINATES','BASIS SET INFORMATION',
        f'Group 1 Type {c["metal"]} : 9s8p5d2f contracted to 6s4p3d2f pattern {{411111/4211/311/11}}',
        'Group 2 Type O : 11s6p2d1f contracted to 5s3p2d1f pattern {62111/411/11/1}',
        'Group 3 Type H : 5s1p contracted to 3s1p pattern {311/1}',
        f'Atom 0{c["metal"]} basis set group => 1','Atom 1O basis set group => 2',
        'Atom 2H basis set group => 3','Atom 3H basis set group => 3','ECP PARAMETER INFORMATION']
    if c['expected_ecp_core']:
        lines += [f'Group 1, Type {c["metal"]} ECP (replacing {c["expected_ecp_core"]} core electrons, lmax=4)',f'Atom 0{c["metal"]} ECP group => 1']
    lines+=['SHARK INTEGRAL PACKAGE','Number of basis functions ... 90',
        f'Total Charge Charge .... {c["charge"]}','Multiplicity Mult .... 1',
        f'Number of Electrons NEL .... {c["expected_electrons"]}','SCF CONVERGED AFTER 10 CYCLES',
        'Last Energy change ... 1.0000e-11 Tolerance : 1.0000e-09',
        'Last Orbital Gradient ... 1.0000e-08 Tolerance : 2.0000e-06',
        'Frozen core treatment ... NO frozen core','Reference Wavefunction ... RHF',
        f'Number of correlated electrons ... {c["expected_correlated_electrons"]}',
        '--- The Coupled-Cluster iterations have converged ---',
        'T1 diagnostic ... 0.004924535','E(CCSD(T)) ... -100.000000000',
        'FINAL SINGLE POINT ENERGY -100.000000000','****ORCA TERMINATED NORMALLY****']
    return '\n'.join(lines)+'\n'

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--panel',required=True);ap.add_argument('--out',type=Path,required=True)
    a=ap.parse_args();tick=time.monotonic();p=cp.load(a.panel);a.out.mkdir(parents=True,exist_ok=False)
    records=[]
    # Replay actual same-version native atomic/basis/electron/ECP formats.
    for metal in ('La','Dy'):
        f=ROOT/f'workspaces/lanm_selectivity_free_energy_20261021/reference_chemistry_v1/pilot_v1/tasks/ecp_pbe0__water_central_{metal}/execution/endpoint.out'
        c=next(c for c in p['cells'] if c['cell_id']==metal+'__complex');r=cp.parse(f.read_text(),c,0)
        subset=['normal_termination','scf_converged','charge','multiplicity','total_electrons','basis_function_count',
            'per_center_basis','four_centers_and_coordinates','ghost_nuclear_charges','ecp_assignment']
        require(all(r['checks'][k] for k in subset),'Actual ORCA6 format replay failed '+metal+': '+str(r['checks']))
        require(not r['checks']['correlated_electrons'] and not r['checks']['cc_converged'],'DFT mistaken for CC')
        records.append({'name':'actual_ORCA6_'+metal,'source':base.pin(f),'checks':{k:r['checks'][k] for k in subset},
            'correctly_rejected_as_CC':r['status']=='failed_admission','coordinate_max_error_A':r['coordinate_max_error_A']})
    # Real canonical CC output contributed to the cclib developer test repository.
    ccfile=HERE.with_name('canonical_water_cp_v1')/'validation_assets/orca42_water_ccsd_t.out'
    ccsource=json.loads((ccfile.parent/'SOURCE.json').read_text())
    require(base.pin(ccfile)['sha256']==ccsource['sha256'],'Downloaded real-output fixture changed')
    c=dict(p['cells'][0]);c.update(charge=0,expected_electrons=10,expected_correlated_electrons=8)
    r=cp.parse(ccfile.read_text(),c,0)
    subset=['normal_termination','scf_converged','cc_converged','RHF_reference','charge','multiplicity','total_electrons','correlated_electrons','final_CC_energy']
    require(all(r['checks'][k] for k in subset),'Real CC format replay failed: '+str(r['checks']))
    require(not r['checks']['no_frozen_core'],'Actual chemical-core output mistaken for NoFrozenCore')
    require(abs(r['energy_hartree']+75.013556310)<1e-12,'Wrong CC total extracted')
    records.append({'name':'real_ORCA42_canonical_CC_format','source':base.pin(ccfile),'provenance':ccsource,
        'checks':{k:r['checks'][k] for k in subset},'correctly_rejected_frozen_core':not r['checks']['no_frozen_core'],
        'not_our_molecular_result':True})
    # Match current executable's literal formats without executing the molecule.
    binary=Path(p['orca']['path']).parent/'orca_mdci';data=binary.read_bytes()
    patterns=[b'--- The Coupled-Cluster iterations have converged ---',b'--- The Coupled-Cluster iterations have NOT converged ---',
        b'--- The CCSD iterations have converged ---',b'--- The CCSD iterations have NOT converged ---',
        b'Frozen core treatment                      ... NO frozen core',b'Number of correlated electrons             ... %6d',
        b'E(CCSD(T))                                     ... %16.9lf']
    require(all(x in data for x in patterns),'Current native diagnostic format differs')
    records.append({'name':'ORCA611_mdci_literal_formats','binary':base.pin(binary),'formats':[x.decode() for x in patterns],
        'evidence':'Read-only executable string inspection; no execution'})
    # Positive controls cover all three distinct ghost/electron roles and metals.
    for c in p['cells']:
        text=control_fixture(c);r=cp.parse(text,c,0);require(r['status']=='complete','Synthetic positive rejected '+c['cell_id']+str(r['checks']))
        (a.out/(c['cell_id']+'.synthetic.txt')).write_text(text)
        records.append({'name':'synthetic_positive_'+c['cell_id'],'checks':r['checks'],'synthetic':True})
    c=next(c for c in p['cells'] if c['cell_id']=='La__water_ghost_ion');text=control_fixture(c)
    cases=[
        ('ghost_charge',text.replace('0 La 0.0000','0 La 11.0000'),'ghost_nuclear_charges',0),
        ('ghost_ECP',text.replace('ECP PARAMETER INFORMATION','ECP PARAMETER INFORMATION\nAtom 0La ECP group => 1\nGroup 1, Type La ECP (replacing 46 core electrons, lmax=4)'),'ecp_assignment',0),
        ('basis_deleted',text.replace('Number of basis functions ... 90','Number of basis functions ... 43'),'basis_function_count',0),
        ('center_basis_missing',text.replace('Atom 0La basis set group => 1',''),'per_center_basis',0),
        ('wrong_metal_contraction',text.replace('contracted to 6s4p3d2f','contracted to 5s4p3d2f'),'per_center_basis',0),
        ('wrong_total_electrons',text.replace('Number of Electrons NEL .... 10','Number of Electrons NEL .... 8'),'total_electrons',0),
        ('frozen_correlated_electrons',text.replace('Number of correlated electrons ... 10','Number of correlated electrons ... 8'),'correlated_electrons',0),
        ('frozen_treatment',text.replace('NO frozen core','chemical core (2 el)'),'no_frozen_core',0),
        ('wrong_charge',text.replace('Total Charge Charge .... 0','Total Charge Charge .... 3'),'charge',0),
        ('wrong_spin',text.replace('Multiplicity Mult .... 1','Multiplicity Mult .... 6'),'multiplicity',0),
        ('wrong_reference',text.replace('Reference Wavefunction ... RHF','Reference Wavefunction ... UHF'),'RHF_reference',0),
        ('failed_CC',text.replace('iterations have converged','iterations have NOT converged'),'cc_converged',0),
        ('failed_CCSD',text.replace('Coupled-Cluster iterations have converged','CCSD iterations have NOT converged'),'cc_converged',0),
        ('Lambda_is_not_CC_convergence',text.replace('The Coupled-Cluster iterations have converged','The Lambda iterations for RHF Coupled Cluster have converged'),'cc_converged',0),
        ('SCF_not_converged',text.replace('SCF CONVERGED AFTER 10 CYCLES','SCF NOT CONVERGED'),'scf_converged',0),
        ('energy_mismatch',text.replace('FINAL SINGLE POINT ENERGY -100.000000000','FINAL SINGLE POINT ENERGY -99.000000000'),'final_CC_energy',0),
        ('missing_CC_energy',text.replace('E(CCSD(T)) ... -100.000000000',''),'final_CC_energy',0),
        ('abnormal_terminal',text.replace('ORCA TERMINATED NORMALLY','ORCA FINISHED BY ERROR TERMINATION'),'normal_termination',0),
        ('process_failure',text,'process_success',1)]
    for name,value,expected,rc in cases:
        r=cp.parse(value,c,rc);require(r['status']=='failed_admission' and not r['checks'][expected],'Failure not rejected '+name)
        records.append({'name':'synthetic_negative_'+name,'expected_rejected_check':expected,'checks':r['checks'],'synthetic':True})
    # Native may print an unrelated auxiliary basis after the orbital summary.
    aux=text.replace('ECP PARAMETER INFORMATION','AUXILIARY/J BASIS SET INFORMATION\nGroup 1 Type La : 1s contracted to 100s\nECP PARAMETER INFORMATION')
    require(cp.parse(aux,c,0)['checks']['per_center_basis'],'Auxiliary summary overwrote primary basis')
    records.append({'name':'auxiliary_summary_isolation','passed':True,'synthetic':True})
    # Independently inspect deck bytes and physical bookkeeping, not expected fields alone.
    olddir=ROOT/'workspaces/lanm_selectivity_free_energy_20261021/reference_chemistry_v1/canonical_water_cp_v1/prepared_v1/inputs'
    static=[]
    for c in p['cells']:
        raw=Path(c['template_input']['path']).read_text();old=olddir/(c['cell_id']+'.inp')
        require(raw==old.read_text(),'Scientific input changed from version1')
        header,*coords=raw.split('* xyz ',1)[1].splitlines();coords=coords[:-1]
        q,m=map(int,header.split());labels=[line.split()[0] for line in coords]
        realZ=sum({'La':11,'Dy':11,'O':8,'H':1}[label] for label in labels if not label.endswith(':'))
        expected={'complex':(3,18,0),'ion_ghost_water':(3,8,3),'water_ghost_ion':(0,10,1)}[c['role']]
        require((q,realZ-q,sum(x.endswith(':') for x in labels))==expected and m==1,'Independent input bookkeeping failed')
        require(('NewECP' in raw)==(c['role']!='water_ghost_ion'),'Ghost ECP declaration failed')
        static.append({'cell_id':c['cell_id'],'charge':q,'effective_electrons':realZ-q,'ghost_centers':sum(x.endswith(':') for x in labels),
            'unchanged_scientific_input':True,'input':base.pin(c['template_input']['path'])})
    receipt={'status':'passed_parser_and_static_controls','panel':base.pin(a.panel),'source':base.pin(__file__),
        'adapter':base.pin(HERE/'REFERENCE.py'),'records':records,'independent_deck_checks':static,
        'real_ghost_output_available':False,'real_current_canonical_output_available':False,
        'first_native_admission_still_required':True,'molecular_calls':0,'local_threads':1,'elapsed_seconds':time.monotonic()-tick,
        'command':sys.argv,'utc':base.now()}
    base.write(a.out/'VALIDATION.json',receipt)
    print(json.dumps({'status':receipt['status'],'records':len(records),'decks':len(static),'molecular_calls':0,'seconds':receipt['elapsed_seconds']}))

if __name__=='__main__':main()
