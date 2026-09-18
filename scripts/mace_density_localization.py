"""Localize actual frozen-density hybrid components, without changing any score."""
from __future__ import annotations
import argparse
import math
from pathlib import Path
import resource
import shutil
import time
import numpy as np
from affordable_common import BOHR_TO_A,HA_TO_KCAL,InvalidArtifact,read_json,record,verify,write_new
from mace_density_gk_hybrid import accepted

PROTOCOL='saved_density_direct_induction_localization_v1'
BOUNDS=(0.,6.,12.,18.,24.,30.,36.,float('inf'))
TOL=1e-7


def reconstruct(collection,case,variant):
    r=read_json(collection);m=read_json(verify(r['manifest']));root=Path(r['manifest']['path']).parent
    if not r['complete'] or not r['numerical_pass']:raise InvalidArtifact('complete qualified actual collection required')
    reported=r['variants'][variant]['cases'][case]
    def task(state):
        t=next(t for t in m['response_tasks'] if t['case_id']==case and t['state']==state and t['variant']==variant and t['poleps']==1e-9)
        a=accepted(t,root/'tasks'/t['task_id'])
        if a is None:raise InvalidArtifact('actual solved response unavailable')
        return t,a[0]
    def induced(s):
        fields=np.loadtxt(verify(s['supplied_fields']),skiprows=1)[:,1:]
        mu=np.array([s['response'][str(i)] for i in range(1,len(fields)+1)])
        return -.5*s['electric']/s['dielec']*np.sum(mu[:,7:10]*fields[:,10:13],axis=1)
    te,se=task('environment');base=induced(se);boundary=read_json(verify(te['boundary']))
    ids=boundary['physical_ids'];index={i:j for j,i in enumerate(ids)};n=len(ids)
    coords=np.array(boundary['positions_A']);distance=np.linalg.norm(coords-coords[index['metal']],axis=1)
    values={};checks=[]
    def check(name,error):
        if not np.isfinite(error):raise InvalidArtifact('nonfinite component reconstruction')
        checks.append(dict(name=name,error_kcal=float(error),pass_=abs(error)<=TOL))
    for metal in ('Ca','La'):
        t,s=task(metal)
        if read_json(verify(t['boundary']))!=boundary:raise InvalidArtifact('endpoint physical boundary differs')
        probes=read_json(verify(t['density_task']['probes']));ix=np.array([index[i] for i in probes['physical_ids']])
        if len(ix)!=len(set(ix)) or set(ix)!=(set(range(n))-{i-1 for i in t['frozen_indices']}):
            raise InvalidArtifact('density field coverage differs from physical exterior')
        moments=np.array([a['global_'] for a in s['moments']])[ix]
        with np.load(verify(t['density_arrays'])) as a:phi=a['phi'];H=a['H']
        with np.load(verify(t['density_task']['field_arrays'])) as a:E=a['exact']
        mat=np.array(t['rotation']);E=E@mat.T;H=np.einsum('ij,njk,lk->nil',mat,H,mat)
        components=dict(charge=moments[:,0]*phi,
            dipole=-np.sum(moments[:,1:4]/BOHR_TO_A*E,axis=1),
            quadrupole=np.sum(moments[:,4:].reshape(-1,3,3)/BOHR_TO_A**2*H,axis=(1,2)))
        arrays={}
        for name,a in components.items():
            arrays[name]=np.zeros(n);arrays[name][ix]=HA_TO_KCAL*a
            check(metal+'_'+name,math.fsum(arrays[name])-reported['endpoints'][metal]['direct_components_kcal'][name])
        arrays['direct']=arrays['charge']+arrays['dipole']+arrays['quadrupole']
        arrays['induction']=induced(s)-base
        for name,component in [('direct','direct_density_kcal'),('induction','environment_induction_transfer_kcal')]:
            check(metal+'_'+name,math.fsum(arrays[name])-reported['endpoints'][metal]['components'][component])
        values[metal]=arrays
    paired={k:values['Ca'][k]-values['La'][k] for k in values['Ca']}
    residues={}
    support=set(boundary['source_support_ids'])
    for i,pid in enumerate(ids):
        resid=pid.rsplit('/',1)[0]
        if resid not in residues:residues[resid]=dict(atom_indices=[],source_support_count=0,modified_exterior_count=0,charge_redistribution_e=0.)
        row=residues[resid];row['atom_indices'].append(i)
        row['source_support_count']+=int(pid in support)
        delta=boundary['environment_charge_by_id'][pid]-boundary['original_charge_by_id'].get(pid,0.)
        if pid not in support:
            row['charge_redistribution_e']+=delta;row['modified_exterior_count']+=int(abs(delta)>1e-12)
    for resid,row in residues.items():
        ix=row.pop('atom_indices');row.update(atom_count=len(ix),min_distance_A=float(min(distance[ix])),max_distance_A=float(max(distance[ix])),
            R_components_kcal={k:math.fsum(a[ix]) for k,a in paired.items()})
    bins=[]
    for lo,hi in zip(BOUNDS[:-1],BOUNDS[1:]):
        ix=(distance>=lo)&(distance<hi)
        bins.append(dict(from_A=lo,to_A=hi if math.isfinite(hi) else None,atom_count=int(ix.sum()),
            R_components_kcal={k:math.fsum(a[ix]) for k,a in paired.items()}))
    for k,a in paired.items():
        check(k+'_pair',math.fsum(a)-(math.fsum(values['Ca'][k])-math.fsum(values['La'][k])))
        check(k+'_residues',math.fsum(v['R_components_kcal'][k] for v in residues.values())-math.fsum(a))
        check(k+'_bins',math.fsum(v['R_components_kcal'][k] for v in bins)-math.fsum(a))
    return dict(case_id=case,variant=variant,collection=record(collection),boundary=te['boundary'],
        physical_ids=ids,distance_A=distance.tolist(),per_atom_R_kcal={k:v.tolist() for k,v in paired.items()},
        residues=residues,radial_bins=bins,checks=checks,checks_pass=all(v['pass_'] for v in checks),
        full_R_kcal=reported['R_kcal'],reported_components=reported['components_R_kcal'] if 'components_R_kcal' in reported else reported)


def audit(config,output):
    start=time.monotonic();cpu=time.process_time();cfg=read_json(config);verify(cfg['plan'])
    if cfg['protocol']!=PROTOCOL:raise InvalidArtifact('unsupported localization protocol')
    cases={};checks=[]
    for name,pin in cfg['collections'].items():
        p=verify(pin);a=reconstruct(p,name,'primary');b=reconstruct(p,name,'rigid')
        checks.extend(dict(case_id=name,**c) for c in a['checks']+b['checks'])
        for resid,row in a['residues'].items():
            error=max(abs(v-b['residues'][resid]['R_components_kcal'][k]) for k,v in row['R_components_kcal'].items())
            checks.append(dict(name=name+'_'+resid+'_rigid',error_kcal=error,pass_=error<=TOL))
        cases[name]=a
    contrasts={}
    for left,right in cfg['contrasts']:
        a,b=cases[left],cases[right];all_res=sorted(set(a['residues'])|set(b['residues']))
        rows={}
        for resid in all_res:
            aa=a['residues'].get(resid);bb=b['residues'].get(resid)
            rows[resid]=dict(left_present=aa is not None,right_present=bb is not None,
                components_kcal={k:(aa['R_components_kcal'][k] if aa else 0.)-(bb['R_components_kcal'][k] if bb else 0.) for k in ('charge','dipole','quadrupole','direct','induction')},
                left=aa,right=bb)
        # An absent residue contributes no atoms to that archived structure;
        # explicit presence flags distinguish this inventory fact from missing energy.
        contrasts[left+'__minus__'+right]=dict(full_R_difference_kcal=a['full_R_kcal']-b['full_R_kcal'],residues=rows,
            unmatched_residue_ids=[k for k,v in rows.items() if not(v['left_present'] and v['right_present'])])
    root=Path(output).resolve();root.mkdir(parents=True,exist_ok=False);shutil.copyfile(__file__,root/Path(__file__).name)
    result=dict(protocol=PROTOCOL,config=record(config),implementation=record(root/Path(__file__).name),cases=cases,
        contrasts=contrasts,checks=checks,checks_pass=all(c['pass_'] for c in checks),tolerance_kcal=TOL,
        new_DFT_calls=0,new_MACE_calls=0,new_solver_calls=0,new_score=None,calibrated_class=None,
        wall_seconds=time.monotonic()-start,CPU_seconds=time.process_time()-cpu,peak_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    write_new(root/'result.json',result)
    if not result['checks_pass']:raise InvalidArtifact('component accounting does not close')
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--config',required=True);p.add_argument('--output',required=True)
    a=p.parse_args();r=audit(a.config,a.output);print({k:r[k] for k in ('checks_pass','wall_seconds','CPU_seconds')})
