#!/usr/bin/env python3
"""Saved native replays, geometry checks and synthetic failure controls; zero calls."""
import argparse
import copy
import importlib.util
import json
from pathlib import Path
import sys
import time

HERE=Path(__file__).resolve().parent
sp=importlib.util.spec_from_file_location('water_force',HERE/'REFERENCE.py');m=importlib.util.module_from_spec(sp);sp.loader.exec_module(m);b=m.b
sp=importlib.util.spec_from_file_location('donor_controls',HERE.parent/'donor_cp_v1/VALIDATE_PARSER.py');controls=importlib.util.module_from_spec(sp);sp.loader.exec_module(controls)

def main():
    tick=time.monotonic();ap=argparse.ArgumentParser();ap.add_argument('--panel',required=True);ap.add_argument('--out',type=Path,required=True);a=ap.parse_args();p=m.load(a.panel)
    a.out.mkdir(parents=True,exist_ok=False);records=[];geometry=m.validate(a.panel)
    collection=json.loads(b.verify(p['canonical_collection']).read_text())
    for pin in collection['result_pins']:
        old=json.loads(b.verify(pin).read_text());c=copy.deepcopy(old['cell']);c.update(method='cc',d4=False);output=b.verify(old['output']);r=m.parse(output.read_text(),c,old['returncode'])
        assert r['status']=='complete',(old['cell_id'],r['checks']);records.append({'case':'actual_CC_'+old['cell_id'],'result':pin,'output':old['output'],'checks':r['checks'],'final_CC_residual':r['CC_residual_trace'][-1][3]})
    baseline=json.loads(b.verify(p['water_baseline']).read_text())
    for pin in baseline['new_result_pins']:
        old=json.loads(b.verify(pin).read_text());c=copy.deepcopy(old['cell'])
        if c['metal'] is None:continue
        c.update(method='pbe0');output=b.verify(old['output']);r=m.parse(output.read_text(),c,old['returncode']);assert r['status']=='complete',(old['cell_id'],r['checks'])
        records.append({'case':'actual_PBE0_'+old['cell_id'],'result':pin,'output':old['output'],'checks':r['checks'],'energy_stages':r['native_energy_stages']})
    for metal in m.METALS:
        center=p['central'][metal];c=copy.deepcopy(center['origin_cell']);c.update(method='pbe0',d4=True);output=b.verify(center['PBE0_saved']['output']);r=m.parse(output.read_text(),c,0)
        assert r['status']=='complete',(metal,r['checks']);assert r['native_energy_stages']['has_analytic_gradient_stage']
        records.append({'case':'actual_PBE0_central_'+metal,'output':center['PBE0_saved']['output'],'checks':r['checks'],'energy_stages':r['native_energy_stages']})
    for c in p['cells']:
        text=controls.fixture(c);r=m.parse(text,c,0);assert r['status']=='complete',(c['cell_id'],r['checks'])
        records.append({'case':'synthetic_positive_'+c['cell_id'],'synthetic':True,'checks':r['checks']})
    c=next(x for x in p['cells'] if x['cell_id']=='cc__La__plus_010__water_ghost_ion');text=controls.fixture(c)
    cases=[('missing_center','\n'.join(line for line in text.splitlines() if not line.startswith('3 H ')),'all_centers_and_coordinates'),
        ('missing_basis','\n'.join(line for line in text.splitlines() if not line.startswith('Atom 3H basis')),'per_center_basis'),
        ('ghost_charge',text.replace('0 La 0.0000','0 La 11.0000'),'ghost_nuclear_charges'),
        ('ghost_ECP',text.replace('ECP PARAMETER INFORMATION','ECP PARAMETER INFORMATION\nAtom 0La ECP group => 1\nGroup 1 Type La ECP replacing 46 core electrons'),'ecp_assignment'),
        ('wrong_AO',text.replace('Number of basis functions ... 90','Number of basis functions ... 43'),'basis_function_count'),
        ('wrong_charge',text.replace('Total Charge Charge .... 0','Total Charge Charge .... 3'),'charge'),
        ('wrong_electrons',text.replace('Number of Electrons NEL .... 10','Number of Electrons NEL .... 8'),'total_electrons'),
        ('wrong_correlated',text.replace('Number of correlated electrons ... 10','Number of correlated electrons ... 8'),'correlated_electrons'),
        ('frozen_core',text.replace('NO frozen core','chemical core'),'no_frozen_core'),
        ('CC_residual',text.replace('0.000000006','0.000000060'),'declared_CC_residual'),
        ('CC_not_converged',text.replace('iterations have converged','iterations have NOT converged'),'cc_converged'),
        ('no_termination',text.replace('ORCA TERMINATED NORMALLY','ORCA TERMINATED ABNORMALLY'),'normal_termination')]
    for name,raw,key in cases:
        r=m.parse(raw,c,0);assert r['status']=='failed_admission' and not r['checks'][key],name;records.append({'case':name,'expected_rejection':key,'synthetic':True})
    r=m.parse(text,c,1);assert not r['checks']['process_success'];records.append({'case':'nonzero_returncode','expected_rejection':'process_success','synthetic':True})
    ghost=next(x for x in p['cells'] if x['cell_id']=='pbe0__Dy__minus_005__water_ghost_ion');raw=controls.fixture(ghost).replace('FINAL SINGLE POINT ENERGY','DFTD4\nDispersion correction -0.001\nFINAL SINGLE POINT ENERGY')
    r=m.parse(raw,ghost,0);assert not r['checks']['dispersion_scope'];records.append({'case':'D4_on_ghost','expected_rejection':'dispersion_scope','synthetic':True})
    real=next(x for x in p['cells'] if x['cell_id']=='pbe0__La__plus_010__complex');raw=controls.fixture(real).replace('Dispersion correction -0.000100000','');r=m.parse(raw,real,0)
    assert not r['checks']['dispersion_scope'];records.append({'case':'real_D4_missing','expected_rejection':'dispersion_scope','synthetic':True})
    # Known analytic polynomials test units, sign, CP subtraction and constant cancellation.
    rows={}
    for method in m.METHODS:
        for metal in m.METALS:
            for label,q in m.SHIFTS.items():
                for role in m.ROLES:
                    slope={'complex':2. if method=='cc' else 3.,'ion_ghost_water':.25,'water_ghost_ion':-.5}[role]
                    energy=-100.+(slope*q+7.*q*q)/b.HA;rows[f'{method}__{metal}__{label}__{role}']={'energy_hartree':energy,'electronic_energy_hartree':energy,'D4_energy_hartree':(-.007+.5*q/b.HA) if role=='complex' else 0.}
    for metal in m.METALS:
        for h in p['steps_A']:
            for method,expected in (('cc',{'complex':2.,'CP':2.25}),('pbe0',{'complex':3.,'CP':3.25,'D4':.5,'complex_D4':3.5,'CP_D4':3.75})):
                got=m.derivatives(rows,method,metal,h);assert all(abs(got[key]-value)<1e-6 for key,value in expected.items()),(method,got)
    projection=m.analytic_record([[2.,0.,0.],[-1.,0.,0.],[-.5,0.,0.],[-.5,0.,0.]],[1.,0.,0.]);assert projection['coordinate_derivative_kcal_mol_A']==2. and projection['net_force_norm_kcal_mol_A']==0.
    records.append({'case':'known_polynomial_CP_D4_derivative_and_force_sign','synthetic':True,'status':'passed'})
    central=p['central']['La']['PBE0_native_energy_stages'];same=copy.deepcopy(central);same['has_analytic_gradient_stage']=False
    review=m.compare_stages(central,same);assert not review['changed_energy_stage_signature'] and review['expected_EnGrad_to_SP_difference']
    changed=copy.deepcopy(same);changed['COSX_grid_settings_signature']=['different grid'];assert m.compare_stages(central,changed)['changed_energy_stage_signature']
    changed=copy.deepcopy(same);changed['recomputed_final_exchange']=not changed['recomputed_final_exchange'];assert m.compare_stages(central,changed)['changed_energy_stage_signature']
    records.append({'case':'final_grid_signature_change_and_expected_SP_stage','synthetic':True,'status':'passed'})
    receipt={'status':'passed','source':b.pin(__file__),'fixture_helper':b.pin(HERE.parent/'donor_cp_v1/VALIDATE_PARSER.py'),'panel':b.pin(a.panel),'geometry_check':geometry,
        'records':records,'actual_native_output_replays':12,'synthetic_controls':len(records)-12,'molecular_calls':0,'local_threads':1,'command':sys.argv,'elapsed_seconds':time.monotonic()-tick,'utc':b.now()}
    pin=b.write(a.out/'VALIDATION.json',receipt);print(json.dumps({'status':'passed','records':len(records),'validation':pin,'molecular_calls':0}))

if __name__=='__main__':main()
