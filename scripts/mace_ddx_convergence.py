"""Logged replay of the unchanged real state that failed native iteration."""
import argparse
import json
import math
import os
from pathlib import Path
import resource
import shutil
import sys
import time
import numpy as np
from affordable_common import BOHR_TO_A,HA_TO_KCAL,InvalidArtifact,cache_key,read_json,record,verify,write_new
from mace_ddx_source_self import validate,direct_potential,potential_error,DIRECT_PROTOCOL,TOL

PROTOCOL='unchanged_GGR_2FW0_La_ddPCM_logged_convergence_v1'

def prepare(config,output):
    c=read_json(config);verify(c['plan']);parent=validate(verify(c['parent']))
    if c['protocol']!=PROTOCOL or parent['protocol']!=DIRECT_PROTOCOL:raise InvalidArtifact('convergence scope differs')
    group=next(g for g in parent['groups'] if g['group_id']=='GGR_2FW0_primary_coarse')
    root=Path(output).resolve();root.mkdir(parents=True,exist_ok=False);impl=root/'implementation';impl.mkdir()
    for n,pin in parent['implementation'].items():shutil.copyfile(verify(pin),impl/n)
    shutil.copyfile(__file__,impl/Path(__file__).name)
    m=dict(protocol=PROTOCOL,config=record(config),parent=c['parent'],plan=c['plan'],source=group['source'],state='La',
           model={**parent['model'],'maxiter':1200},grid=parent['grids']['coarse'],tolerances=parent['tolerances'],python=parent['python'],module=parent['module'],
           implementation={p.name:record(p) for p in impl.glob('*.py')},requested_forward_solves=1,new_DFT_calls=0,new_MACE_calls=0)
    m['cache_key']=cache_key(m);write_new(root/'manifest.json',m);return preflight(root/'manifest.json')

def preflight(path):
    m=read_json(path);parent=validate(verify(m['parent']));payload={k:v for k,v in m.items() if k!='cache_key'}
    group=next(g for g in parent['groups'] if g['group_id']=='GGR_2FW0_primary_coarse')
    if m['protocol']!=PROTOCOL or m['cache_key']!=cache_key(payload) or m['model']!={**parent['model'],'maxiter':1200} or m['grid']!=parent['grids']['coarse'] or m['source']!=group['source'] or m['state']!='La' or m['tolerances']!=parent['tolerances']:raise InvalidArtifact('changed physical state or convergence scope')
    for pin in [m['config'],m['plan'],m['source'],m['python'],m['module'],*m['implementation'].values()]:verify(pin)
    return m

def execute(path,output):
    import pyddx
    m=preflight(path);root=Path(path).resolve().parent;attempt=root/'attempt_0001'
    if not os.environ.get('SLURM_JOB_ID') or int(os.environ.get('SLURM_CPUS_ON_NODE','0'))<64 or record(pyddx.__file__)!=m['module']:raise InvalidArtifact('unqualified runtime')
    attempt.mkdir(exist_ok=False);start=time.monotonic();cpu=time.process_time();data=read_json(verify(m['source']));arrays={}
    r=dict(protocol=PROTOCOL,manifest=record(path),status='failed',forward_solve_started=False,energy_hartree=None,energy_kcal_mol=None,new_score=None,baseline_changed=False)
    try:
        centres=np.asfortranarray(np.array(data['coordinates_A']).T/BOHR_TO_A);radii=np.array(data['radii_A'])/BOHR_TO_A;q=np.array(data['endpoints']['La']['charge_e'])
        model=pyddx.Model(sphere_centres=centres,sphere_radii=radii,**m['model'],**m['grid'],logfile=str(attempt/'solver.log'))
        r['setup_wall_seconds']=time.monotonic()-start;cavity=np.array(model.cavity);phi=direct_potential(cavity,centres,q);psi=np.array(model.multipole_psi(np.asfortranarray(q.reshape(1,-1)/np.sqrt(4*np.pi))))
        arrays.update(cavity_bohr=cavity,phi=phi,psi=psi,source_charge_e=q)
        expected=np.zeros_like(psi);expected[0]=np.sqrt(4*np.pi)*q;pe=potential_error(cavity,centres,q,phi);se=float(np.max(abs(psi-expected)))
        r.update(potential_max_error_au=pe,psi_max_error=se,charge_sum_e=math.fsum(q))
        if pe>TOL['potential_au'] or se>TOL['psi']:raise InvalidArtifact('unchanged source check failed')
        state=pyddx.State(model,np.asfortranarray(psi),phi);r['forward_solve_started']=True;solve=time.monotonic();state.solve(tol=m['tolerances']['solver'])
        x=np.array(state.x);energy=float(state.energy());arrays['x']=x
        r.update(status='computed',energy_hartree=energy,energy_kcal_mol=energy*HA_TO_KCAL,solve_wall_seconds=time.monotonic()-solve,
                 final_single_layer_iterations=state.x_n_iter,contraction_error_kcal=abs(.5*math.fsum((psi*x).ravel())-energy)*HA_TO_KCAL)
    except Exception as exc:r['failure_reason']=str(exc)
    np.savez_compressed(attempt/'arrays.npz',**arrays);r['arrays']=record(attempt/'arrays.npz')
    r.update(wall_seconds=time.monotonic()-start,CPU_seconds=time.process_time()-cpu,peak_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
             solver_log=record(attempt/'solver.log') if (attempt/'solver.log').exists() else None,slurm_job_id=os.environ['SLURM_JOB_ID'],
             convergence_criterion='native relative iterate change; not independently evaluated residual')
    write_new(output,r);return r

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest='command',required=True)
    q=sub.add_parser('prepare');q.add_argument('--config',required=True);q.add_argument('--output',required=True)
    for op in ('dry-run','execute'):
        q=sub.add_parser(op);q.add_argument('--manifest',required=True)
        if op=='execute':q.add_argument('--output',required=True)
    a=p.parse_args()
    if a.command=='prepare':m=prepare(a.config,a.output);r={'status':'prepared','solves':1}
    elif a.command=='dry-run':m=preflight(a.manifest);r={'status':'pass','solves':1}
    else:r=execute(a.manifest,a.output)
    print(json.dumps(r))
