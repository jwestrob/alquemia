"""Qualify real frozen AMOEBA multipoles in the native ddX source interface."""
from pathlib import Path
import argparse
import gc
import json
import math
import os
import resource
import shutil
import time
import numpy as np
from affordable_common import BOHR_TO_A,InvalidArtifact,cache_key,read_json,record,verify,write_new
from mace_frozen_response import source_arrays

PROTOCOL='native_ddX_frozen_AMOEBA_multipole_source_bridge_v1'
TOL=dict(potential_au=1e-10,psi=1e-12,rigid_potential_au=1e-8)
GRID=dict(lmax=6,n_lebedev=194)


def distributions(static,response):
    _,md,mp,_,_=source_arrays(static,response)
    p=np.array([a['global_'] for a in static['moments']]);average=p.copy();difference=np.zeros_like(p)
    average[:,1:4]+=.5*(md+mp);difference[:,1:4]=.5*(md-mp)
    for a in (average,difference):
        a[:,1:4]/=BOHR_TO_A;a[:,4:]/=BOHR_TO_A**2
    return dict(average=average,difference=difference)


def coefficients(p):
    q=p[:,0];d=p[:,1:4];Q=p[:,4:].reshape(-1,3,3)
    out=np.empty((9,len(p)));out[0]=q/np.sqrt(4*np.pi)
    out[1:4]=np.sqrt(3/(4*np.pi))*d[:,[1,2,0]].T
    out[4:]=np.sqrt(15/(4*np.pi))*np.array([2*Q[:,0,1],2*Q[:,1,2],
                (2*Q[:,2,2]-Q[:,0,0]-Q[:,1,1])/np.sqrt(3),2*Q[:,0,2],Q[:,0,0]-Q[:,1,1]])
    return np.asfortranarray(out)


def cartesian_potential(points,centers,p):
    q=p[:,0];d=p[:,1:4];Q=p[:,4:].reshape(-1,3,3);trace=np.trace(Q,axis1=1,axis2=2)
    out=np.empty(len(points))
    for start in range(0,len(points),16):
        r=points[start:start+16,None,:]-centers[None,:,:];norm=np.linalg.norm(r,axis=2)
        if np.min(norm)<1e-8:raise InvalidArtifact('probe overlaps point multipole')
        dr=np.einsum('nij,ij->ni',r,d);rQr=np.einsum('nij,ijk,nik->ni',r,Q,r)
        v=q/norm+dr/norm**3+3*rQr/norm**5-trace/norm**3
        out[start:start+len(v)]=[math.fsum(row) for row in v]
    return out


def prepare(parent,functional,collection,plan,output):
    p=read_json(parent);f=read_json(functional);co=read_json(collection);pm=read_json(verify(co['manifest']))
    if not f['checks_pass'] or not co['complete'] or not co['numerical_pass']:
        raise InvalidArtifact('qualified native functional and archived states required')
    if p['protocol']!='fixed_source_full_protein_ddCPCM_component_refinement_v2':
        raise InvalidArtifact('pinned conductor parent required')
    root=Path(output).resolve();root.mkdir(parents=True,exist_ok=False);impl=root/'implementation';impl.mkdir()
    for name,pin in pm['implementation'].items():shutil.copyfile(verify(pin),impl/name)
    for name in ('mace_ddx_multipole_bridge.py','mace_frozen_response.py'):
        shutil.copyfile(Path(__file__).with_name(name),impl/name)
    m=dict(protocol=PROTOCOL,parent=record(parent),functional=record(functional),collection=record(collection),plan=record(plan),
           python=p['python'],module=p['module'],model={**p['model'],'enable_fmm':False},grid=GRID,tolerances=TOL,
           implementation={x.name:record(x) for x in impl.glob('*.py')},groups=[],requested_setups=4,
           requested_potentials=16,requested_integrals=16,new_solver_calls=0,new_DFT_calls=0,new_MACE_calls=0)
    for case in ('GGR_2FW0','GGR_2FVY'):
        for variant in ('primary','rigid'):
            g=dict(group_id=case+'_'+variant,case_id=case,variant=variant,endpoints={})
            for metal in ('Ca','La'):
                st=next(t for t in pm['static_tasks'] if t['case_id']==case and t['state']==metal and t['variant']==variant)
                rt=next(t for t in pm['response_tasks'] if t['case_id']==case and t['state']==metal and t['variant']==variant and not t['task_id'].endswith('_standard'))
                sp=co['tasks'][st['task_id']]['result'];rp=co['tasks'][rt['task_id']]['result']
                s=read_json(verify(sp));r=read_json(verify(rp));b=read_json(verify(st['boundary']));_,md,mp,_,_=source_arrays(s,r)
                if np.any(md[np.array(st['frozen_indices'])-1]) or np.any(mp[np.array(st['frozen_indices'])-1]):
                    raise InvalidArtifact('quantum source has induced response')
                coords=[a['xyz_A'] for a in s['parameters']['atoms']];radii=[a['radius_A'] for a in s['parameters']['atoms']]
                physical=dict(coordinates_A=coords,radii_A=radii,physical_ids=b['physical_ids'],frozen_indices=st['frozen_indices'],
                              rotation=st['rotation'],translation_A=st['translation_A'])
                if 'physical' in g and g['physical']!=physical:raise InvalidArtifact('paired physical systems differ')
                g['physical']=physical;g['endpoints'][metal]=dict(static=sp,response=rp)
            g['cache_key']=cache_key(dict(group=g,model=m['model'],grid=GRID,implementation=m['implementation'],plan=m['plan']))
            m['groups'].append(g)
    m['cache_key']=cache_key(m);write_new(root/'manifest.json',m);return preflight(root/'manifest.json')


def preflight(path):
    m=read_json(path);p=read_json(verify(m['parent']));f=read_json(verify(m['functional']))
    if m['protocol']!=PROTOCOL or m['grid']!=GRID or m['tolerances']!=TOL or m['model']!={**p['model'],'enable_fmm':False} or not f['checks_pass']:
        raise InvalidArtifact('source representation scope differs')
    if m['cache_key']!=cache_key({k:v for k,v in m.items() if k!='cache_key'}) or len(m['groups'])!=4:
        raise InvalidArtifact('source representation manifest changed')
    if m['python']!=p['python'] or m['module']!=p['module'] or {g['group_id'] for g in m['groups']}!={c+'_'+v for c in ('GGR_2FW0','GGR_2FVY') for v in ('primary','rigid')}:
        raise InvalidArtifact('source geometry inventory or executable changed')
    for pin in (m['collection'],m['plan'],m['python'],m['module'],*m['implementation'].values()):verify(pin)
    for g in m['groups']:
        if g['cache_key']!=cache_key(dict(group={k:v for k,v in g.items() if k!='cache_key'},model=m['model'],grid=GRID,implementation=m['implementation'],plan=m['plan'])):
            raise InvalidArtifact('source group cache changed')
        for e in g['endpoints'].values():
            s=read_json(verify(e['static']));r=read_json(verify(e['response']));_,md,mp,_,_=source_arrays(s,r)
            if np.any(md[np.array(g['physical']['frozen_indices'])-1]) or np.any(mp[np.array(g['physical']['frozen_indices'])-1]):
                raise InvalidArtifact('frozen quantum source dipoles changed')
            for state in (s,r):verify(state['receipt']['log'])
            if [a['xyz_A'] for a in s['parameters']['atoms']]!=g['physical']['coordinates_A'] or [a['radius_A'] for a in s['parameters']['atoms']]!=g['physical']['radii_A']:
                raise InvalidArtifact('source cavity differs from native archive')
    return m


def execute(path,output):
    import pyddx
    m=preflight(path);root=Path(path).resolve().parent;rows={};start=time.monotonic();cpu=time.process_time()
    if not os.environ.get('SLURM_JOB_ID') or int(os.environ.get('SLURM_CPUS_ON_NODE','0'))<64 or record(pyddx.__file__)!=m['module'] or pyddx.__version__!='0.9.0' or np.__version__!='1.26.4':
        raise InvalidArtifact('qualified64CPU/native runtime required')
    for g in m['groups']:
        d=root/'groups'/g['group_id']/'attempt_0001';d.mkdir(parents=True,exist_ok=False);tic=time.monotonic()
        centers=np.array(g['physical']['coordinates_A'])/BOHR_TO_A;radii=np.array(g['physical']['radii_A'])/BOHR_TO_A
        model=pyddx.Model(sphere_centres=np.asfortranarray(centers.T),sphere_radii=radii,**m['model'],**GRID)
        cavity=np.array(model.cavity);indices=np.unique(np.linspace(0,cavity.shape[1]-1,min(256,cavity.shape[1]),dtype=int))
        arrays=dict(cavity_bohr=cavity,probe_indices=indices);row=dict(group=g['group_id'],cache_key=g['cache_key'],states={},
             actual_parameters={k:v for k,v in model.input_parameters.items() if k not in ('sphere_centres','sphere_radii')},setup_wall_seconds=time.monotonic()-tic)
        for metal,e in g['endpoints'].items():
            s=read_json(verify(e['static']));r=read_json(verify(e['response']))
            for kind,poles in distributions(s,r).items():
                key=metal+'_'+kind;before=time.monotonic();result=dict(status='failed')
                try:
                    coeff=coefficients(poles);native=model.multipole_electrostatics(coeff,derivative_order=0)
                    phi=np.array(native['phi']);psi=np.array(model.multipole_psi(coeff))
                    expected=np.zeros_like(psi)
                    for l in range(3):expected[l*l:(l+1)**2]=4*np.pi*coeff[l*l:(l+1)**2]/((2*l+1)*radii**l)
                    direct=cartesian_potential(cavity[:,indices].T,centers,poles)
                    pe=float(np.max(abs(phi[indices]-direct)));se=float(np.max(abs(psi-expected)))
                    if not all(np.isfinite(a).all() for a in (coeff,phi,psi,direct)):raise InvalidArtifact('nonfinite source array')
                    result.update(status='computed',potential_error_au=pe,psi_error=se,charge_sum_e=math.fsum(poles[:,0]),
                         maximum_quadrupole_trace_e_bohr2=float(np.max(abs(np.trace(poles[:,4:].reshape(-1,3,3),axis1=1,axis2=2)))),
                         checks_pass=pe<=TOL['potential_au'] and se<=TOL['psi'])
                    arrays.update({key+'_phi':phi,key+'_psi':psi,key+'_coefficients':coeff,key+'_poles_au':poles})
                except Exception as exc:result.update(failure_reason=str(exc))
                result['wall_seconds']=time.monotonic()-before;row['states'][key]=result
                write_new(d/(key+'.json'),result);print(g['group_id'],key,result,flush=True)
        np.savez_compressed(d/'arrays.npz',**arrays);row['arrays']=record(d/'arrays.npz');row['wall_seconds']=time.monotonic()-tic
        write_new(d/'result.json',row);rows[g['group_id']]=row;model=None;arrays.clear();gc.collect()
    write_new(root/'execution_receipt.json',dict(wall_seconds=time.monotonic()-start,CPU_seconds=time.process_time()-cpu,
               peak_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,slurm_job_id=os.environ['SLURM_JOB_ID'],
               python=os.sys.version,numpy=np.__version__,pyddx=pyddx.__version__,new_solver_calls=0))
    return collect(path,output)


def collect(path,output):
    m=preflight(path);root=Path(path).resolve().parent;groups={};checks=[];transform_errors={}
    for g in m['groups']:
        pin=record(root/'groups'/g['group_id']/'attempt_0001/result.json');r=read_json(verify(pin));verify(r['arrays'])
        if r['cache_key']!=g['cache_key']:raise InvalidArtifact('executed source group changed')
        groups[g['group_id']]=dict(receipt=pin,**r)
        params={**m['model'],**GRID}
        # ddx_parameters.f90 normalizes both ignored FMM orders to -2 when
        # FMM is disabled. This is native input canonicalization, not a mismatch.
        if not params['enable_fmm']:
            params.update(fmm_local_lmax=-2,fmm_multipole_lmax=-2)
        ok=all(r['actual_parameters'].get(k)==v for k,v in params.items() if k not in ('model','enable_force'))
        checks.append(dict(name=g['group_id']+'_parameters',pass_=ok))
        for label,state in r['states'].items():checks.append(dict(name=g['group_id']+'_'+label,pass_=state.get('checks_pass',False)))
    for case in ('GGR_2FW0','GGR_2FVY'):
        g=next(g for g in m['groups'] if g['group_id']==case+'_rigid');physical=g['physical'];primary=groups[case+'_primary'];rigid=groups[case+'_rigid']
        with np.load(verify(primary['arrays']),allow_pickle=False) as a,np.load(verify(rigid['arrays']),allow_pickle=False) as b:
            ix=a['probe_indices'];points=a['cavity_bohr'][:,ix].T@np.array(physical['rotation']).T+np.array(physical['translation_A'])/BOHR_TO_A
            centers=np.array(physical['coordinates_A'])/BOHR_TO_A
            for label in primary['states']:
                if primary['states'][label]['status']!='computed' or rigid['states'][label]['status']!='computed':
                    checks.append(dict(name=case+'_'+label+'_rigid',pass_=False,error=None));continue
                potential=cartesian_potential(points,centers,b[label+'_poles_au']);err=float(np.max(abs(potential-a[label+'_phi'][ix])))
                original=a[label+'_poles_au'];rotated=b[label+'_poles_au'];R=np.array(physical['rotation'])
                expected_Q=np.einsum('ij,njk,lk->nil',R,original[:,4:].reshape(-1,3,3),R)
                transform_errors[case+'_'+label]=dict(charge_e=float(np.max(abs(original[:,0]-rotated[:,0]))),
                    dipole_e_bohr=float(np.max(abs(original[:,1:4]@R.T-rotated[:,1:4]))),
                    quadrupole_e_bohr2=float(np.max(abs(expected_Q-rotated[:,4:].reshape(-1,3,3)))))
                checks.append(dict(name=case+'_'+label+'_rigid',error=err,tolerance=TOL['rigid_potential_au'],pass_=err<=TOL['rigid_potential_au']))
    result=dict(protocol=PROTOCOL,manifest=record(path),reporter=record(__file__),groups=groups,checks=checks,transform_errors=transform_errors,checks_pass=all(c['pass_'] for c in checks),
                complete=all(s['status']=='computed' for g in groups.values() for s in g['states'].values()),
                new_solver_calls=0,new_DFT_calls=0,new_MACE_calls=0,new_score=None,baseline_changed=False)
    write_new(output,result);return result


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);sub=parser.add_subparsers(dest='command',required=True)
    p=sub.add_parser('prepare')
    for k in ('parent','functional','collection','plan','output'):p.add_argument('--'+k,required=True)
    for op in ('dry-run','execute','collect'):
        p=sub.add_parser(op);p.add_argument('--manifest',required=True)
        if op!='dry-run':p.add_argument('--output',required=True)
    a=parser.parse_args()
    if a.command=='prepare':r=prepare(a.parent,a.functional,a.collection,a.plan,a.output)
    elif a.command=='dry-run':r=preflight(a.manifest)
    elif a.command=='execute':r=execute(a.manifest,a.output)
    else:r=collect(a.manifest,a.output)
    print(json.dumps({k:v for k,v in r.items() if k not in ('groups','implementation')}))
