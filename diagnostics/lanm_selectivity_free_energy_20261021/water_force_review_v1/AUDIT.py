#!/usr/bin/env python3
"""Independent frozen-input review. Reads saved files; invokes no scientific code."""
import hashlib
import itertools
import json
import math
from pathlib import Path
import re
import sys
import time


def pin(p):
    p=Path(p).resolve()
    return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}


def main():
    start=time.monotonic(); file=Path(sys.argv[1]).resolve(); out=Path(sys.argv[2]).resolve()
    p=json.loads(file.read_text()); checks=0; pins={}
    def check(v,why):
        nonlocal checks
        checks+=1
        if not v:raise AssertionError(why)
    def walk(x):
        if isinstance(x,dict):
            if set(x)=={'path','sha256'}:
                check(pin(x['path'])==x,'Changed '+x['path']);pins[x['path']]=x
            else:
                for v in x.values():walk(v)
        elif isinstance(x,list):
            for v in x:walk(v)
    walk(p)
    check(p['adapter']['sha256']=='2871d929c576281fba152de0a7f58069e8c8e88068cd173f970c11fe0952db7c','Reviewed adapter')
    check(p['orca']['sha256']=='38b5f057452fef275c0a1b98d270ad03d0411dab70453820d1c678bd732f6c83','Proven ORCA')
    expected={f'{method}__{metal}__{q}__{role}' for method,metal,q,role in itertools.product(
        ['cc','pbe0'],['La','Dy'],['minus_010','minus_005','plus_005','plus_010'],['complex','ion_ghost_water','water_ghost_ion'])}
    check(len(p['cells'])==48 and {c['cell_id']for c in p['cells']}==expected,'Exact48 inventory, no central')
    check(p['steps_A']==[.01,.005] and p['derivative_resolution_kcal_mol_A']==.1,'Frozen steps/resolution')
    check(p['scout_cell_id']=='cc__La__plus_010__complex','Declared single scout')
    max_translation_error=0.;max_internal_error=0.;derivatives={}
    for g in p['geometries']:
        central=p['central'][g['metal']];x=central['origin_cell']['positions_A'];actual=g['positions_A']
        ray=[x[1][k]-x[0][k]for k in range(3)];r=math.sqrt(sum(v*v for v in ray));u=[v/r for v in ray]
        check(actual[0]==x[0],'Fixed metal');check(g['q_A']in(-.01,-.005,.005,.01),'No extra displacement')
        for i in [1,2,3]:
            err=max(abs(actual[i][k]-x[i][k]-g['q_A']*u[k])for k in range(3));max_translation_error=max(max_translation_error,err)
            check(err<1e-12,'Rigid translation')
        for i,j in [(1,2),(1,3),(2,3)]:
            err=abs(math.dist(actual[i],actual[j])-math.dist(x[i],x[j]));max_internal_error=max(max_internal_error,err);check(err<2e-12,'Preserved internal geometry')
        check(abs(math.dist(actual[0],actual[1])-r-g['q_A'])<1e-12,'Positive coordinate increases metal-O distance')
        saved=json.loads(Path(central['PBE0_saved']['result']['path']).read_text())
        for model,forces,target in [('PBE0',saved['forces_kcal_mol_A'],central['PBE0_analytic']),('AMOEBA',central['AMOEBA_native']['force_kcal_mol_A'],central['AMOEBA_analytic'])]:
            derivative=-sum(forces[i][k]*u[k]for i in [1,2,3]for k in range(3))
            check(abs(derivative-target['coordinate_derivative_kcal_mol_A'])<1e-12,'Derivative equals negative water force sum')
            derivatives[g['metal']+'_'+model]=derivative
    old=json.loads((file.parent.parent/'prepared_v1/PANEL.json').read_text());oldcells={c['cell_id']:c for c in old['cells']}
    roles={'complex':([False]*4,3,18),'ion_ghost_water':([False,True,True,True],3,8),'water_ghost_ion':([True,False,False,False],0,10)}
    for c in p['cells']:
        raw=Path(c['template_input']['path']).read_text();oldraw=Path(oldcells[c['cell_id']]['template_input']['path']).read_text()
        check(raw==oldraw,'Unchanged scientific input after reporting revision')
        ghost,q,ne=roles[c['role']];check(c['ghost_mask']==ghost and c['charge']==q and c['expected_electrons']==ne,'Fragment state')
        check(c['multiplicity']==1 and c['expected_correlated_electrons']==ne,'Effective singlet and correlated inventory')
        check(c['physical_multiplicity_of_parent']==(1 if c['metal']=='La'else 6),'Parent physical state retained')
        check(c['expected_basis_functions']==90 and c['expected_basis_functions_by_center']==[47,31,6,6],'Same90 AO inventory')
        header,*coordinates=raw.split('* xyz ',1)[1].splitlines();coordinates=coordinates[:-1]
        check(list(map(int,header.split()))==[q,1] and len(coordinates)==4,'Written coordinate state')
        for i,line in enumerate(coordinates):
            t=line.split();check(t[0]==c['symbols'][i]+(':'if ghost[i]else ''),'Written real/ghost labels')
            check(list(map(float,t[1:]))==c['positions_A'][i],'Written exact coordinates')
        basis=Path(p['assets'][c['metal']]['basis']['path']).read_text();gto='NewGTO '+basis.split('NewGTO ',1)[1]
        check(gto in raw,'Exact proven metal orbital coefficients')
        if ghost[0]:check('DelECP '+c['metal']in raw and 'NewECP'not in raw and c['expected_ecp_core']==0,'No ghost-metal ECP')
        else:check(basis in raw and c['expected_ecp_core']==(46 if c['metal']=='La'else 55),'Exact real-metal ECP')
        check('%maxcore 3000' in raw and '%pal nprocs 4 end' in raw,'Four-rank template matches proven resource baseline')
        check('VeryTightSCF' in raw and 'Guess PModel' in raw and 'MaxIter 300' in raw,'Unchanged SCF controls')
        if c['method']=='cc':
            check('RHF CCSD(T)'in raw and 'NoFrozenCore'in raw and 'STol 1e-8'in raw and 'MaxIter 100'in raw,'CC recipe')
            check('D4'not in raw and 'RIJCOSX'not in raw and 'EnGrad'not in raw,'Pure canonical energy')
        else:
            check(('D4'in raw)==(c['role']=='complex'),'D4 on real complexes only')
            check('PBE0'in raw and 'RIJCOSX'in raw and 'DefGrid3'in raw and 'EnGrad'not in raw,'PBE0 displacedSP recipe')
            check(Path(p['assets'][c['metal']]['aux']['path']).read_text()in raw,'Exact auxiliary basis retained')
    validation=file.parent.parent/'validation_v2/VALIDATION.json';v=json.loads(validation.read_text())
    check(v['status']=='passed' and v['panel']==pin(file),'Actual saved validation')
    check(v['actual_native_output_replays']==12 and v['synthetic_controls']==65 and len(v['records'])==77,'Actual replay/failure control denominator')
    source=Path(p['adapter']['path']).read_text()
    check('paired_two_step_and_AMOEBA_comparison' in source and 'paired_PBE0_minus_CC_derivatives' in source,'Paired-specific summaries present')
    check('changed_energy_stage_signature' in source and 'prespecified_smaller_step_summary' in source,'Surface/step flags explicit')
    report={'status':'passed_independent_saved_input_review','implementation':pin(__file__),'panel':pin(file),'validation':pin(validation),
       'checks':checks,'unique_verified_pins':len(pins),'source_pins':list(pins.values()),'cells':48,'new_central_calls':0,
       'max_translation_roundoff_A':max_translation_error,'max_internal_distance_roundoff_A':max_internal_error,
       'independently_recomputed_saved_coordinate_derivatives_kcal_mol_A':derivatives,
       'native_calls':0,'local_threads':1,'physical_accuracy_qualified':False,'submission_released_by_this_review':False,
       'resource_note':'Prior canonical water actual grant16000 MiB on4 ranks; MaxCore3000. New runner derives MaxCore=floor(.75*actualMiB/4).',
       'scope':'Read-only existing source/input/native-replay evidence; no molecular computation or reserved outcome access',
       'elapsed_seconds':time.monotonic()-start}
    out.parent.mkdir(parents=True,exist_ok=True)
    with out.open('x')as f:json.dump(report,f,indent=2,allow_nan=False);f.write('\n')
    print(json.dumps({k:report[k]for k in ('status','checks','unique_verified_pins','cells','max_translation_roundoff_A','max_internal_distance_roundoff_A','elapsed_seconds')}))


if __name__=='__main__':main()
