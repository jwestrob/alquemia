#!/usr/bin/env python3
"""Saved-only exact-state and numerical operator replay; never launches native."""
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import sys
import time
import numpy as np


def pin(path):
    path=Path(path).resolve()
    return dict(path=str(path),sha256=hashlib.sha256(path.read_bytes()).hexdigest())


def load_reader(item):
    if pin(item['path'])!=item: raise ValueError('Snapshot reader pin changed')
    spec=importlib.util.spec_from_file_location('pinned_snapshot',item['path'])
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module


def tokens(path):
    return [line.split() for line in Path(path).read_text().splitlines()]


def replay(manifest_path,off,on):
    start,cpu=time.monotonic(),time.process_time()
    m=json.loads(Path(manifest_path).read_text());reader=load_reader(m['snapshot_reader'])
    off,on=Path(off),Path(on)
    checks={}
    def require(label,condition):
        entry=checks.setdefault(label,dict(count=0,passed=True))
        entry['count']+=1;entry['passed'] &= bool(condition)
        if not condition: raise ValueError(label)
    for mode,folder in [('OFF',off),('ON',on)]:
        for key,name in [('xyz','system.xyz'),('key','system.key'),('parameters','amoeba09.prm')]:
            require('exact prepared '+mode+' '+key,pin(folder/name)['sha256']==m['cells'][mode][key]['sha256'])
    require('literal off/on energy and all force tokens', (off/'harness_ef.dat').read_bytes()==(on/'harness_ef.dat').read_bytes())
    ef=tokens(on/'harness_ef.dat')
    require('complete EF format',len(ef)==9 and ef[0]==['PME_HARNESS_EF_V1'] and ef[-1]==['END','6'] and ef[1][0]=='ENERGY' and len(ef[1])==2)
    require('finite energy',math.isfinite(float(ef[1][1])))
    for i,row in enumerate(ef[2:8],1):
        require(f'finite gradient {i}',len(row)==5 and row[:2]==['GRADIENT',str(i)] and all(math.isfinite(float(v)) for v in row[2:]))
    require('OFF has no export/state/application',all(not(off/name).exists() for name in ('pme_norm_export_001.dat','harness_state.dat','harness_apply.dat')))
    with (on/'pme_norm_export_001.dat').open() as stream:s=reader.read(stream)
    require('exact declared initialized state',s['shape']==(32,32,32) and s['p']==5 and s['fft']=='FFTW' and s['ewald']==.3 and s['volume']==44.**3 and len(s['sites'])==6)
    norm=reader.replay(s)
    reconstructed=[norm[k] for k in ('c1','cinf','cF','U_R_nominal','evenness_perturbation_nominal')]
    qsummary=[norm[k] for k in ('q_min','q_max','q_abs_max','q_zero','q_conjugate_difference_max')]
    for i,(actual,expected) in enumerate(zip(s['native_norms']+s['q_stats'],reconstructed+qsummary)):
        require(f'native summary replay {i}',abs(actual-expected)<=1e-12*max(1.,abs(expected)))
    state=tokens(on/'harness_state.dat');position=0
    def row(tag,count):
        nonlocal position
        if position>=len(state):raise ValueError('Incomplete state')
        value=state[position];position+=1
        require(f'state {tag} row format',len(value)==count+1 and value[0]==tag)
        return value[1:]
    row('PME_HARNESS_STATE_V1',0)
    rr=row('STATE',6)
    require('current state grid/alpha/volume exact',list(map(int,rr[:4]))==[32,32,32,5] and list(map(float,rr[4:]))==[s['ewald'],s['volume']])
    for axis in range(3):
        rr=row('TRANS',4)
        require('current transform exact',int(rr[0])==axis+1 and list(map(float,rr[1:]))==s['trans'][axis].tolist())
    for site in s['sites']:
        rr=row('SITE',9)
        require('current site/alpha/coordinate exact',list(map(int,rr[:5]))==[site['atom'],site['type'],*site['igrid']] and list(map(float,rr[5:]))==[site['alpha'],*site['xyz']])
        for slot in range(5):
            for axis in range(3):
                rr=row('THETA',5)
                require('current spline exact',list(map(int,rr[:3]))==[site['atom'],axis+1,slot+1] and list(map(float,rr[3:]))==site['theta'][axis,slot].tolist())
    for z in range(32):
        for y in range(32):
            for x in range(32):
                rr=row('QFAC',4)
                require('current qfac exact',list(map(int,rr[:3]))==[x+1,y+1,z+1] and float(rr[3])==s['qfac'][x,y,z])
    require('complete state tail',row('END',3)==['6','90','32768'] and position==len(state))
    # Only this six-site fixture: no full-box/dense global-response construction.
    B=np.zeros(s['shape']+(18,))
    for grid,i,b in reader.entries(s):B[grid+(slice(3*i,3*i+3),)]+=b
    transformed=np.fft.fftn(B,axes=(0,1,2))*s['qfac'][...,None]
    back=np.fft.ifftn(transformed,axes=(0,1,2))*math.prod(s['shape'])
    R=B.reshape(-1,18).T@back.reshape(-1,18)
    rows=tokens(on/'harness_apply.dat');pos=0
    require('application header',rows[0]==['PME_HARNESS_APPLY_V1']);pos+=1
    maximum=0.;scale=max(1.,float(np.max(np.abs(R))))
    leakage=0.
    for mode in range(1,4):
        for col in range(18):
            expected=R[:,col] if mode==1 else 1j*R[:,col] if mode==2 else R[:,col]+1j*R[:,(col+1)%18]
            scale=max(scale,float(np.max(np.abs(expected))))
            observed=np.empty(18,dtype=complex)
            for atom in range(1,7):
                rr=rows[pos];pos+=1
                require('complete application row',len(rr)==10 and rr[:4]==['APPLY',str(mode),str(col+1),str(atom)])
                values=np.array(list(map(float,rr[4:])))
                require('finite application values',bool(np.isfinite(values).all()))
                observed[(atom-1)*3:atom*3]=values[:3]+1j*values[3:]
            maximum=max(maximum,float(np.max(np.maximum(np.abs(observed.real-expected.real),np.abs(observed.imag-expected.imag)))))
            if mode==1:leakage=max(leakage,float(np.max(np.abs(observed.imag))))
            if mode==2:leakage=max(leakage,float(np.max(np.abs(observed.real))))
    require('exact application count/end',rows[pos:]==[['END','54']])
    require('actual native packed operator replay',maximum<=1e-11*scale)
    return dict(status='passed_initialized_PME_export_numerical_qualification',manifest=pin(manifest_path),implementation=pin(__file__),
                checks_passed=sum(x['count'] for x in checks.values()),checks=checks,prepared_native_programs=2,
                complete_gradient_records=2,observed_internal_reciprocal_applications=54,execution_receipts_admitted=False,
                raw_inputs={mode:{k:pin(folder/k) for k in ('system.xyz','system.key','amoeba09.prm')} for mode,folder in [('OFF',off),('ON',on)]},
                raw_artifacts={str(p):pin(p) for p in [off/'harness_ef.dat',on/'harness_ef.dat',on/'pme_norm_export_001.dat',on/'harness_state.dat',on/'harness_apply.dat']},
                operator_max_abs_difference=maximum,operator_scale=scale,operator_tolerance=1e-11*scale,channel_leakage_max=leakage,
                norm_replay=norm,full_operator_stability_admitted=False,physical_accuracy_admitted=False,
                cost=dict(wall_seconds=time.monotonic()-start,CPU_seconds=time.process_time()-cpu,threads=1),new_native_calls=0)


if __name__=='__main__':
    manifest,off,on,out=sys.argv[1:]
    destination=Path(out);destination.mkdir(parents=True,exist_ok=False)
    try:result=replay(manifest,off,on)
    except Exception as error:result=dict(status='failed_saved_PME_export_replay',error=repr(error),new_native_calls=0)
    (destination/'RESULT.json').write_text(json.dumps(result,indent=2)+'\n')
    print(result['status'])
    sys.exit(0 if result['status'].startswith('passed_') else 1)
