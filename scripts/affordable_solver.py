"""Bounded real ESP/APBS validation; never substitutes a baseline score on failure."""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
import copy
import json
import math
import os
from pathlib import Path
import subprocess
import time

import numpy as np

from affordable_common import InvalidArtifact, read_json, write_new, record, verify, xyz, BOHR_TO_A, cache_key, snapshot_implementation
from affordable_environment import mbis_charges, prepare, collect


def real_esp_points(atoms, radii):
    """Two fixed Fibonacci shells per real atom; exclude every atomic interior."""
    coords=np.array([a[1:] for a in atoms]); points=[]
    n=32; indices=np.arange(n); z=1-2*(indices+.5)/n
    theta=indices*math.pi*(3-math.sqrt(5)); radial=np.sqrt(1-z*z)
    directions=np.column_stack((radial*np.cos(theta),radial*np.sin(theta),z))
    atom_radii=np.array([radii[a[0]] for a in atoms])
    for origin,r in zip(coords,atom_radii):
        for scale in (1.4,1.8):
            for point in origin+scale*r*directions:
                if np.all(np.linalg.norm(coords-point,axis=1)>=1.4*atom_radii-1e-8): points.append(point)
    if len(points)<32: raise InvalidArtifact('insufficient exterior ESP points')
    return np.array(points)/BOHR_TO_A


def run_command(command,cwd,log,cost):
    start=time.monotonic()
    with Path(log).open('x') as f:
        r=subprocess.run(['/usr/bin/time','-v','-o',str(Path(cost).resolve()),*command],cwd=cwd,stdout=f,stderr=subprocess.STDOUT,
                         env={**os.environ,'OMP_NUM_THREADS':'1','OPENBLAS_NUM_THREADS':'1','MKL_NUM_THREADS':'1'},check=False)
    result={'command':command,'returncode':r.returncode,'wall_seconds':time.monotonic()-start,
            'log':record(log),'resource_usage':record(cost),'slurm_job_id':os.environ.get('SLURM_JOB_ID')}
    return result


def esp_check(task,orca_dir,radii):
    endpoint=Path(task['output_path']); d=endpoint.parent/'esp_check'
    if (d/'quality.json').exists():
        raise InvalidArtifact('existing ESP attempt requires explicit collection, not overwrite')
    d.mkdir(parents=True,exist_ok=True)
    charges=mbis_charges(endpoint,verify(task['xyz']),task['charge'])
    atoms=xyz(task['xyz']['path']); points=real_esp_points(atoms,radii)
    pp=d/'points_bohr.xyz'; pp.write_text(str(len(points))+'\n'+'\n'.join(' '.join(f'{v:.12f}' for v in row) for row in points)+'\n')
    gbw=endpoint.parent/'endpoint.runtime.gbw'; densities=endpoint.parent/'endpoint.runtime.densities'
    if not gbw.exists() or not densities.exists(): raise InvalidArtifact('wavefunction/density unavailable')
    potential=d/'potential.out'; utility=Path(orca_dir)/'orca_vpot'
    receipt=run_command([str(utility),str(gbw),'endpoint.runtime.scfp',str(pp.resolve()),str(potential.resolve())],
                        endpoint.parent,d/'utility.log',d/'resources.txt')
    write_new(d/'execution.json',receipt)
    if receipt['returncode']!=0: raise InvalidArtifact('orca_vpot failed; see actual utility output')
    values=np.loadtxt(potential,skiprows=1)
    if values.shape!=(len(points),4): raise InvalidArtifact('unexpected potential row count')
    if np.allclose(values[:,:3],points,atol=1e-5,rtol=0): qm=values[:,3]
    elif np.allclose(values[:,1:],points,atol=1e-5,rtol=0): qm=values[:,0]
    else: raise InvalidArtifact('potential coordinate/order mismatch')
    coords=np.array([a[1:] for a in atoms])/BOHR_TO_A
    approximation=np.sum(np.array(charges['charge_e'])[None,:]/np.linalg.norm(points[:,None,:]-coords[None,:,:],axis=2),axis=1)
    rms=float(np.sqrt(np.mean((qm-approximation)**2))); norm=float(np.sqrt(np.mean(qm**2)))
    passed=rms<=.005 or (norm>0 and rms/norm<=.10)
    quality={'status':'passed' if passed else 'failed','RMS_potential_error_au':rms,'relative_RMS':rms/norm if norm else None,
             'point_count':len(points),'actual_quantum_potential':record(potential),'points':record(pp),
             'gbw':record(gbw),'density':record(densities),'utility':record(utility),'charges':charges,
             'execution_receipt':record(d/'execution.json'),'units':'atomic_units_e_per_bohr',
             'rule':'relative_RMS<=0.10 OR absolute_RMS<=0.005_au'}
    write_new(d/'quality.json',quality)
    return quality,record(d/'quality.json')


def settings(physical,spacing=.5,padding=20.,extent=None,center=None):
    coords=np.array([a['xyz_A'] for a in physical])
    if center is None: center=(coords.min(0)+coords.max(0))/2
    if extent is None:
        # Multiples of 64 A keep both 0.5 and 0.4 A multigrid lattices at the
        # exact same extent, independently of endpoint and partition.
        extent=np.ceil((coords.max(0)-coords.min(0)+2*padding)/64)*64
    dims=np.rint(np.asarray(extent)/spacing).astype(int)+1
    return {'solute_dielectric':1.,'solvent_dielectric':78.54,'salt_molar':0.,'temperature_K':298.15,
            'grid_spacing_A':spacing,'grid_dimensions':dims.tolist(),'grid_center_A':np.asarray(center).tolist(),
            'surface':'mol','probe_radius_A':1.4,'charge_discretization':'spl2','boundary_condition':'mdh',
            'radii_policy':'Bondi_CHNOS_common_1p8A_metal_v1'}


def run_state(state,root,label,apbs):
    d=root/label; d.mkdir(parents=True,exist_ok=True)
    sp=d/'state.json'; write_new(sp,state)
    prepare(sp,d/'calculation')
    c=d/'calculation';receipt=run_command([str(apbs),'transfer.in'],c,c/'apbs.out',c/'resources.txt')
    write_new(c/'execution.json',receipt)
    if receipt['returncode']!=0: raise InvalidArtifact('APBS failed; raw output retained')
    result=collect(c/'apbs_manifest.json',c/'apbs.out');result['execution_receipt']=record(c/'execution.json')
    result['execution_cache_key']=cache_key({'preparation':record(c/'apbs_manifest.json'),'solver':record(apbs)})
    write_new(d/'result.json',result)
    return result


def transform(state,rotation=None,translation=None):
    s=copy.deepcopy(state); matrix=np.eye(3) if rotation is None else np.asarray(rotation)
    shift=np.zeros(3) if translation is None else np.asarray(translation)
    center=np.asarray(s['settings']['grid_center_A'])
    for field in ('physical_atoms','core_atoms','environment_atoms'):
        for a in s[field]: a['xyz_A']=(matrix@(np.asarray(a['xyz_A'])-center)+center+shift).tolist()
    # Keep the box fixed under sub-grid translation to measure grid placement
    # error, rather than making the test trivially identical by translating it.
    return s


def execute_checks(pilot_path,root,apbs):
    if not os.environ.get('SLURM_JOB_ID'): raise InvalidArtifact('solver execution requires SLURM')
    m=read_json(pilot_path);root=Path(root).resolve();root.mkdir(parents=True,exist_ok=True)
    cpus=int(os.environ['SLURM_CPUS_ON_NODE']); cap=m['budget']['solver_allocated_core_seconds']
    started=time.monotonic(); outcomes=[]; states={}; numerical=[]
    implementation=snapshot_implementation(root/'implementation')
    from affordable_state import validate_skeleton_pair
    paired_states=[]
    for case in m['cases']:
        paths=[Path(pilot_path).parent.parent/'environment'/f"{case['case']}_{metal}"/'skeleton.json' for metal in ('La','Ca')]
        if all(p.exists() for p in paths):
            paired_states.append({'case':case['case'],**validate_skeleton_pair(*[read_json(p) for p in paths])})
        elif any(p.exists() for p in paths):
            raise InvalidArtifact('unpaired prepared environment')
    write_new(root/'paired_state_checks.json',paired_states)
    def spent(): return (time.monotonic()-started)*cpus
    def allowed(estimated_wall): return spent()+estimated_wall*cpus<=cap
    # The first low-level batch is admitted conservatively using a 30 s estimate;
    # subsequent solver batches use measured APBS time and recorded headroom.
    if not allowed(30): raise InvalidArtifact('solver budget unavailable')
    orca_dir=verify(m['orca']).parent
    def quality_one(t):
        try:
            q,r=esp_check(t,orca_dir,m['environment_model']['radii_A'])
            return t,q,r,None
        except Exception as exc: return t,None,None,str(exc)
    with ThreadPoolExecutor(max_workers=min(len(m['tasks']),cpus)) as pool:
        quality_results=list(pool.map(quality_one,m['tasks']))
    for task,quality,receipt,error in quality_results:
        row={'task_id':task['task_id'],'ESP_status':quality['status'] if quality else 'unavailable','reason':error}
        skeleton=Path(pilot_path).parent.parent/'environment'/task['task_id']/'skeleton.json'
        if quality and quality['status']=='passed' and skeleton.exists():
            s=read_json(skeleton)
            for a,q in zip(s['core_atoms'],quality['charges']['charge_e']): a['charge_e']=q
            s['charge_quality']={'status':'passed','receipt':receipt}
            s['settings']=settings(s['physical_atoms'])
            states[task['task_id']]=s;row['environment_status']='ready_for_physical_checks'
        else:
            row['environment_status']='unavailable';row['reason']=error or ('charge_quality_failed' if quality and quality['status']!='passed' else 'source_environment_preparation_unsupported')
        outcomes.append(row)
    write_new(root/'charge_quality_and_preparation.json',outcomes)
    # Freeze all variants before looking at any APBS value. Fixed-core may be
    # structurally unsupported, which remains a recorded denominator failure.
    schedule=[]
    for name in states: schedule.append((name+'/primary',states[name]))
    for name,s in states.items():
        refined=copy.deepcopy(s);refined['settings']=settings(s['physical_atoms'],spacing=.4)
        extended=copy.deepcopy(s);extended['settings']=settings(s['physical_atoms'],padding=30.)
        schedule.extend([(name+'/refined',refined),(name+'/extended',extended)])
        if name.startswith('1h4i_qm33'):
            theta=.37; rot=[[math.cos(theta),-math.sin(theta),0],[math.sin(theta),math.cos(theta),0],[0,0,1]]
            schedule.extend([(name+'/translated',transform(s,translation=[.173,.117,.231])),
                             (name+'/rotated',transform(s,rotation=rot))])
    anchor=states.get('1h4i_qm33_La')
    if anchor:
        identity=copy.deepcopy(anchor);identity['physical_atoms']=[dict(a,charge_e=0.) for a in identity['core_atoms']]
        identity['environment_atoms']=[];identity['expected_environment_charge_e']=0.
        schedule.extend([('identity/1h4i_qm33_La',identity),('repeat/1h4i_qm33_La',copy.deepcopy(anchor))])
    write_new(root/'numerical_schedule.json',{'entries':[{'label':n,'state_cache_key':cache_key(s),'grid':s['settings']} for n,s in schedule],
                                             'budget_core_seconds':cap,'state_count':len(schedule),'charging_solves_per_state':6})
    estimate=60.0  # first real six-solve state admission estimate, not a timing result
    stop_after_failure=False
    for name,s in schedule:
        if stop_after_failure:
            numerical.append({'label':name,'status':'not_run_after_failure','result':None});continue
        if not allowed(estimate):
            numerical.append({'label':name,'status':'budget_exhausted','result':None});continue
        before=time.monotonic()
        try:
            result=run_state(s,root,name,apbs)
            numerical.append({'label':name,'status':'computed','result':result})
        except Exception as exc:
            numerical.append({'label':name,'status':'failed','reason':str(exc),'result':None})
        estimate=max(estimate,(time.monotonic()-before)*1.5)
        # Infrastructure/syntax failure is not repeated across the entire panel.
        if numerical[-1]['status']=='failed': stop_after_failure=True
    summary={'schema_version':'alquemia.apbs_physical_pilot.v1','pilot':record(pilot_path),'apbs':record(apbs),
             'actual_implementation':implementation,'preparation':outcomes,'numerical_checks':numerical,'elapsed_seconds':time.monotonic()-started,
             'allocated_core_seconds':spent(),'allocated_cpus':cpus,'slurm_job_id':os.environ['SLURM_JOB_ID'],
             'budget_exceeded':spent()>cap,'predictive_claim':'none_development_only'}
    write_new(root/'solver_result.json',summary)
    return summary


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--pilot-manifest',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--apbs',type=Path,required=True)
    a=p.parse_args();r=execute_checks(a.pilot_manifest,a.output,a.apbs.resolve())
    from affordable_compare import compare,report
    repo=Path(__file__).resolve().parents[1]
    data=compare(repo)
    write_new(a.output/'comparison.json',data)
    with (a.output/'REPORT.md').open('x') as f: f.write(report(data))
    print(json.dumps({'allocated_core_seconds':r['allocated_core_seconds'],'states':len(r['numerical_checks'])}))


if __name__=='__main__':main()
