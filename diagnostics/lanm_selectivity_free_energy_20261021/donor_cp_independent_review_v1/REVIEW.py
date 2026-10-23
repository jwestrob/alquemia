#!/usr/bin/env python3
"""Strict saved-only donor CP provenance, raw-energy and balanced-cycle review.

Full mode requires all22 actual results. Explicit dft-subset mode reviews the
already completed10 and never supplies missing CC energies. No native calls.
"""
import argparse
from decimal import Decimal, getcontext
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import time

PANEL_SHA='1a4bca36f2e5f0d0e223b235b0915708a455627444bcb355fc7cbcba2930a701'
ADAPTER_SHA='f4baaeddce09cb26551e845f05e7f3c1bfff604a84fc02c666908ada2420dea2'
HA=Decimal('627.509474')
getcontext().prec=40


def pin(path):
    path=Path(path).resolve()
    return {'path':str(path),'sha256':hashlib.sha256(path.read_bytes()).hexdigest()}


def dec(value):return Decimal(str(value))


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--panel',required=True);ap.add_argument('--out',required=True)
    ap.add_argument('--mode',choices=('full','dft-subset'),default='full');ap.add_argument('results',nargs='+');a=ap.parse_args()
    start=time.monotonic();out=Path(a.out).resolve();out.mkdir(parents=True,exist_ok=False)
    result={'status':'review_failed','implementation':pin(__file__),'mode':a.mode,'panel':pin(a.panel),
       'new_molecular_calls':0,'local_threads':1,'full22_review_passed':False,
       'physical_accuracy_qualified':False,'fit_or_sampling_released':False}
    checks=0;pins={};raw_rows={};reviews=[]
    def check(ok,label):
        nonlocal checks
        checks+=1
        if not ok:raise AssertionError(label)
    def verify(x):
        check(pin(x['path'])==x,'Changed pin '+x['path']);pins[x['path']]=x;return Path(x['path'])
    try:
        for name in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):
            check(os.environ.get(name)=='1','One local numerical thread required')
        check(result['panel']['sha256']==PANEL_SHA,'Exact frozen22 panel')
        panel=json.loads(Path(a.panel).read_text());check(panel['adapter']['sha256']==ADAPTER_SHA,'Frozen runner')
        spec=importlib.util.spec_from_file_location('donor_cp_frozen_review',verify(panel['adapter']));ref=importlib.util.module_from_spec(spec);spec.loader.exec_module(ref)
        panel=ref.load(a.panel);check(ref.b.HA==float(HA),'Unit conversion')
        for key in ('adapter','specification','parent_panel','base_helper','canonical_helper','dft_helper','water_baseline','orca'):verify(panel[key])
        cells={c['cell_id']:c for c in panel['cells']};check(len(cells)==22,'Full22 inventory')
        expected=set(cells)if a.mode=='full'else {k for k,c in cells.items()if c['method']=='pbe0'}
        check(len(expected)==(22 if a.mode=='full'else 10),'Declared review denominator')
        for filename in a.results:
            rp=pin(filename);pins[rp['path']]=rp;r=json.loads(Path(filename).read_text());cid=r['cell_id']
            check(cid in expected and cid not in raw_rows,'Unknown or duplicate result cell')
            c=cells[cid];check(r['status']=='complete'and r['returncode']==0 and all(r['checks'].values()),'Stored native admission')
            check(r['panel']==result['panel'] and r['adapter']==panel['adapter'] and r['cell']==c,'Exact result/cell/source binding')
            stdout=verify(r['output']);verify(r['input']);verify(r['stderr'])
            startup=Path(filename).parent/'START.json';sp=pin(startup);pins[sp['path']]=sp;s=json.loads(startup.read_text())
            check(s['cell']==c and s['panel']==result['panel'] and s['input']==r['input'],'Actual START/input binding')
            check(s['job_id']==r['job_id'] and s['ranks']==4,'Actual job/rank identity')
            check(s['memory_MiB']>0 and s['maxcore_MiB']==math.floor(.75*s['memory_MiB']/4),'Actual memory rule')
            check(Path(r['input']['path']).read_text()==ref.deck(panel,c,4,s['maxcore_MiB']),'Executed full deck matches expected cell')
            replay=ref.parse(stdout.read_text(),c,r['returncode'])
            check(replay['status']=='complete'and all(replay['checks'].values())and replay['checks']==r['checks'],'Raw native replay admission')
            fields=('energy_hartree',)if c['method']=='cc'else ('energy_hartree','electronic_energy_hartree','D4_energy_hartree')
            for k in fields:check(math.isfinite(replay[k])and replay[k]==r[k],'Raw/stored energy '+cid+' '+k)
            raw_rows[cid]=replay;reviews.append({'cell_id':cid,'result':rp,'start':sp,'input':r['input'],'output':r['output'],
                'native_checks':replay['checks'],'raw_energy_fields_hartree':{k:replay[k]for k in fields},
                'printed_scf_criteria':replay['printed_scf_criteria'],'native_scf_mode':replay['native_scf_mode'],
                'native_scf_stop_signals':replay['native_scf_stop_signals']})
        check(set(raw_rows)==expected,'All declared rows required, missing values never zero')
        # Saved real-complex energies are independently re-extracted from raw
        # output, preserving the electronic/D4 split and exact source geometry.
        saved={}
        parent=ref.b.load_panel(verify(panel['parent_panel']))
        for donor in ref.DONORS:
            saved[donor]={}
            for metal in ref.METALS:
                old=panel['saved_complexes'][donor][metal]
                for key in ('result','output','source_geometry'):verify(old[key])
                r=json.loads(Path(old['result']['path']).read_text());check(r['status']=='complete'and r['returncode']==0 and all(r['checks'].values()),'Saved complex admitted')
                pc=next(c for c in parent['cells']if c['geometry_id']==f'matched_dry_{donor}_central_{metal}')
                check(r['panel']==panel['parent_panel']and r['geometry']==pc['geometry']['xyz']==old['source_geometry'],'Saved matched-dry complex lineage')
                es,ed,ef=ref.dft.electronic_parts(Path(old['output']['path']).read_text())
                check(len(es)==len(ed)==len(ef)==1 and abs(es[0]+ed[0]-ef[0])<1e-8,'Raw saved electronic/D4/final closure')
                for k,v in [('electronic_energy_hartree',es[0]),('D4_energy_hartree',ed[0]),('final_energy_hartree',ef[0])]:
                    check(v==old[k],'Raw saved complex component')
                saved[donor][metal]={'electronic':dec(es[0]),'D4':dec(ed[0])}
        water=json.loads(verify(panel['water_baseline']).read_text());check(water['status']=='complete','Frozen complete water comparison')
        check(panel['isolated_ion_D4_hartree']==0 and panel['primary_exchange_resolution_kcal_mol']==.25,'Frozen dispersion and interpretation')
        methods=['PBE0_CP_kcal_mol','PBE0_D4_CP_kcal_mol']
        if a.mode=='full':methods.insert(0,'CCSD_T_CP_kcal_mol')
        donors={};max_arithmetic_error=Decimal(0)
        for donor in ref.DONORS:
            pair={};ligD4=dec(raw_rows[f'pbe0__{donor}__real_ligand_D4']['D4_energy_hartree'])
            for metal in ref.METALS:
                def E(method,role,field):return dec(raw_rows[f'{method}__{donor}__{metal}__{role}'][field])
                electro=(saved[donor][metal]['electronic']-E('pbe0','ion_ghost_ligand','electronic_energy_hartree')-E('pbe0','ligand_ghost_ion','electronic_energy_hartree'))*HA
                disp=(saved[donor][metal]['D4']-ligD4)*HA
                interaction={'PBE0_CP_kcal_mol':electro,'PBE0_D4_CP_kcal_mol':electro+disp}
                if a.mode=='full':interaction['CCSD_T_CP_kcal_mol']=(E('cc','complex','energy_hartree')-E('cc','ion_ghost_ligand','energy_hartree')-E('cc','ligand_ghost_ion','energy_hartree'))*HA
                exchange={k:interaction[k]-dec(water['per_metal'][metal][k])for k in methods}
                pair[metal]={'interaction':interaction,'water_to_donor_exchange':exchange,'D4_interaction_kcal_mol':disp}
            paired={k:pair['Dy']['water_to_donor_exchange'][k]-pair['La']['water_to_donor_exchange'][k]for k in methods}
            # Independent algebraic arrangement: donor gap minus water gap.
            for k in methods:
                alternative=(pair['Dy']['interaction'][k]-pair['La']['interaction'][k])-(dec(water['per_metal']['Dy'][k])-dec(water['per_metal']['La'][k]))
                check(paired[k]==alternative,'Balanced gap algebra')
            donors[donor]={'per_metal':pair,'paired_exchange':paired}
            if a.mode=='full':
                errors={k:paired[k]-paired['CCSD_T_CP_kcal_mol']for k in methods if k!='CCSD_T_CP_kcal_mol'}
                watererror=(dec(water['per_metal']['Dy']['PBE0_D4_CP_kcal_mol'])-dec(water['per_metal']['La']['PBE0_D4_CP_kcal_mol']))-(dec(water['per_metal']['Dy']['CCSD_T_CP_kcal_mol'])-dec(water['per_metal']['La']['CCSD_T_CP_kcal_mol']))
                donors[donor].update(method_minus_CC_paired_exchange=errors,PBE0_D4_within_declared_resolution=abs(errors['PBE0_D4_CP_kcal_mol'])<=Decimal('.25'),
                    water_only_signed_method_gap=watererror,paired_exchange_minus_water_only_signed_method_gap=errors['PBE0_D4_CP_kcal_mol']-watererror,
                    paired_gap_absolute_magnitude_reduced_from_water_only=abs(errors['PBE0_D4_CP_kcal_mol'])<abs(watererror))
        # Only the already-complete DFT subset is read in subset mode. Full mode
        # exercises the frozen collector in memory after all strict admissions.
        if a.mode=='full':existing=ref.collect(a.panel,a.results)
        else:
            prior=Path(a.panel).parent.parent/'dft_subset_collection_v1/RESULT.json';existing=json.loads(prior.read_text());pins[str(prior.resolve())]=pin(prior)
            check(existing['status']=='DFT_subset_complete_correlated_comparison_pending'and existing['panel']==result['panel'],'Original subset provenance')
        for donor in ref.DONORS:
            for metal in ref.METALS:
                for group in ('interaction','water_to_donor_exchange'):
                    for k in methods:
                        error=abs(donors[donor]['per_metal'][metal][group][k]-dec(existing['donors'][donor]['per_metal'][metal][group][k]));max_arithmetic_error=max(max_arithmetic_error,error);check(error<Decimal('1e-7'),'Independent decimal '+group)
            for k in methods:
                error=abs(donors[donor]['paired_exchange'][k]-dec(existing['donors'][donor]['paired_exchange'][k]));max_arithmetic_error=max(max_arithmetic_error,error);check(error<Decimal('1e-7'),'Independent paired arithmetic')
            if a.mode=='full':check(donors[donor]['PBE0_D4_within_declared_resolution']==existing['donors'][donor]['PBE0_D4_within_declared_resolution'],'Frozen decision predicate')
        # No scientific source is edited; verify observed evidence remained fixed.
        for v in list(pins.values()):verify(v)
        result.update(status='full22_strict_review_passed'if a.mode=='full'else 'DFT10_strict_review_passed_CC_unavailable',
            full22_review_passed=a.mode=='full',reviewed_rows=len(raw_rows),missing_full22_cells=sorted(set(cells)-set(raw_rows)),
            raw_replay_reviews=reviews,independent_decimal_cycles=donors,max_independent_arithmetic_residual_kcal_mol=str(max_arithmetic_error),
            full_collector_admission_gap_closed_for_these_rows=True)
    except Exception as exc:result['error']=repr(exc)
    result.update(checks=checks,verified_pins=list(pins.values()),elapsed_seconds=time.monotonic()-start)
    def convert(x):
        if isinstance(x,Decimal):return str(x)
        raise TypeError(type(x).__name__)
    with (out/'RESULT.json').open('x')as f:json.dump(result,f,indent=2,default=convert,allow_nan=False);f.write('\n')
    print(json.dumps({'status':result['status'],'checks':checks,'result':pin(out/'RESULT.json'),'error':result.get('error')}))
    return 0 if result['status']!='review_failed'else 1


if __name__=='__main__':raise SystemExit(main())
