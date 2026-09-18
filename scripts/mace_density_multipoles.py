"""Actual saved-density coupling to native permanent protein multipoles."""
from __future__ import annotations
import argparse
from concurrent.futures import ThreadPoolExecutor
import fcntl
import json
import math
import os
from pathlib import Path
import shutil
import socket
import time
import numpy as np
from affordable_common import BOHR_TO_A,HA_TO_KCAL,InvalidArtifact,cache_key,read_json,record,verify,write_new
from affordable_solver import run_command
from density_embedding import parse_potential
from mace_hybrid import rotation
from mace_qm_field import validate as validate_fields,accepted as accepted_utility,NULLS
from mace_tinker_framework_solver import validate as validate_frameworks,native_call,parse_parameters
from mace_native_field_accounting import build,check_runtime_parameters

PROTOCOL='vacuum_density_native_AMOEBA_permanent_multipole_coupling_diagnostic_v1'
STEPS=(.01,.005)
TOL=dict(center_potential_au=1e-8,moment_rotation_native=1e-10,scalar_rigid_kcal=1e-8,
         endpoint_quadrupole_refinement_kcal=.02,pair_quadrupole_refinement_kcal=.01)
PAIRS=((0,1),(0,2),(1,2))


def stencil():
    offsets=[]
    for axis in np.eye(3):offsets.extend([axis,-axis])
    for i,j in PAIRS:
        for si,sj in ((1,1),(1,-1),(-1,1),(-1,-1)):
            v=np.zeros(3);v[i]=si;v[j]=sj;offsets.append(v)
    return np.array(offsets)


def points(centers):
    p=np.asarray(centers)
    return np.concatenate([p]+[(p[:,None,:]+h*stencil()[None,:,:]).reshape(-1,3) for h in STEPS])


def hessian(values,n):
    v=np.asarray(values);center=v[:n];offsets=v[n:].reshape(2,n,18);out=np.empty((2,n,3,3))
    for s,h in enumerate(STEPS):
        for i in range(3):out[s,:,i,i]=(offsets[s,:,2*i]-2*center+offsets[s,:,2*i+1])/h**2
        for k,(i,j) in enumerate(PAIRS):
            a,b,c,d=np.moveaxis(offsets[s,:,6+4*k:10+4*k],1,0)
            out[s,:,i,j]=out[s,:,j,i]=(a-b-c+d)/(4*h*h)
    return center,out


def monopole_hessian(q,source,probes):
    d=np.asarray(probes)[:,None,:]-np.asarray(source)[None,:,:];r=np.linalg.norm(d,axis=2)
    if np.min(r)<1e-8:raise InvalidArtifact('source/probe overlap')
    return np.sum(np.asarray(q)[None,:,None,None]*(3*d[:,:,:,None]*d[:,:,None,:]/r[:,:,None,None]**5
        -np.eye(3)[None,None,:,:]/r[:,:,None,None]**3),axis=1)


def coupling(q,mu,quad,phi,field,H):
    return dict(charge=HA_TO_KCAL*math.fsum(q*phi),
        dipole=-HA_TO_KCAL*math.fsum((mu*field).ravel()),
        quadrupole=HA_TO_KCAL*math.fsum((quad*H).ravel()),
        isotropic_trace_part_of_quadrupole=HA_TO_KCAL*math.fsum(np.trace(quad,axis1=1,axis2=2)*np.trace(H,axis1=1,axis2=2)/3))


def native_effective_local(pole,axis,planar):
    p=np.array(pole,copy=True)
    if planar:
        if axis=='Z-Bisect':
            p[1]=p[6]=p[10]=0.;p[4]=p[8]=.5*(p[4]+p[8])
        elif axis=='3-Fold':p[1:]=0.
        else:raise InvalidArtifact('unexpected planar native axis')
    return p


def parse_moments(text,task):
    if 'ALQUEMIA_MOMENTS_COMPLETE' not in text or 'ENERGIES' in text or 'RMS Change' in text:
        raise InvalidArtifact('incomplete/unexpected native moment operation')
    params=parse_parameters(text)
    if task is not None:check_runtime_parameters(params,task,1)
    kinds=('RAW_POLE','LOCAL_POLE','GLOBAL_POLE','FRAME','AXIS');data={k:{} for k in kinds}
    for line in text.splitlines():
        f=line.split()
        if f and f[0] in data:
            i=int(f[1])
            if i in data[f[0]]:raise InvalidArtifact('duplicate native moment row')
            data[f[0]][i]=f[2:]
    n=params['inventory'][0]
    if any(set(v)!=set(range(1,n+1)) for v in data.values()):raise InvalidArtifact('missing native moment atoms')
    errors=dict(rotation=0.,trace_preservation=0.,orthogonality=0.,raw_local_chiral=0.)
    native=[];max_trace=0.;planar_count=0
    for i in range(1,n+1):
        raw=np.array(data['RAW_POLE'][i],dtype=float);loc=np.array(data['LOCAL_POLE'][i],dtype=float)
        glob=np.array(data['GLOBAL_POLE'][i],dtype=float);frame=data['FRAME'][i]
        if len(raw)!=13 or len(loc)!=13 or len(glob)!=13 or len(frame)!=10:raise InvalidArtifact('native moment dimensions')
        axis=data['AXIS'][i];a=np.array(frame[1:],dtype=float).reshape(3,3);planar=frame[0]=='T'
        p=native_effective_local(loc,axis[0],planar);planar_count+=planar
        rebuilt=np.r_[p[0],a@p[1:4],(a@p[4:].reshape(3,3)@a.T).ravel()]
        errors['rotation']=max(errors['rotation'],float(np.max(abs(rebuilt-glob))))
        errors['orthogonality']=max(errors['orthogonality'],float(np.max(abs(a@a.T-np.eye(3)))))
        errors['trace_preservation']=max(errors['trace_preservation'],abs(float(np.trace(glob[4:].reshape(3,3))-np.trace(p[4:].reshape(3,3)))))
        # chkpole may only flip y-bearing components; retain both actual records.
        chiral=np.array(raw,copy=True);chiral[[2,5,7,9,11]]*=-1
        errors['raw_local_chiral']=max(errors['raw_local_chiral'],min(float(np.max(abs(loc-raw))),float(np.max(abs(loc-chiral)))))
        max_trace=max(max_trace,abs(float(np.trace(glob[4:].reshape(3,3)))))
        native.append(dict(index=i,raw_local=raw.tolist(),local_after_chkpole=loc.tolist(),global_=glob.tolist(),
            effective_local=p.tolist(),rotation=a.tolist(),planar=planar,axis=axis[0],axis_indices=list(map(int,axis[1:]))))
    result=dict(status='native_global_moments_exported',parameters=params,atoms=native,max_errors=errors,
        max_absolute_native_quadrupole_trace_eA2=max_trace,planar_sites=planar_count,
        pass_=max(errors.values())<=TOL['moment_rotation_native'],new_energy_calls=0,new_response_solves=0,**NULLS)
    cache_key(result);return result


def key(task,m):
    return cache_key(dict(task={k:v for k,v in task.items() if k!='cache_key'},
        **{k:m[k] for k in ('protocol','steps_bohr','tolerances','sources','utility','implementation','plan','moment_exports')}))


def prepare(fields,centers,frameworks,moment_software,plan,output):
    start=time.monotonic();cpu=time.process_time()
    fr=read_json(fields);fm=validate_fields(verify(fr['manifest']));cm=read_json(centers)
    tm=validate_frameworks(frameworks);sw=read_json(moment_software)
    if not fr['complete'] or not fr['numerical_pass']:raise InvalidArtifact('qualified actual field source required')
    if fr['manifest']['sha256']!='0b639b8502fdade3b8dd7727563908cbdf1d2237ae4008b8cfe329ce872d0673':
        raise InvalidArtifact('declared field experiment differs')
    if record(centers)['sha256']!='90d553d2480a1efb3911fdc9169774610c9b124d40db7e17e6b4ad45af32d439':
        raise InvalidArtifact('declared potential centers differ')
    if sw['returncode'] or sw['library']!=read_json(verify(tm['software']))['library']:
        raise InvalidArtifact('wrong/unavailable native moment software')
    root=Path(output).resolve();root.mkdir(parents=True,exist_ok=False);impl=root/'implementation';impl.mkdir()
    for name,pin in fm['implementation'].items():shutil.copyfile(verify(pin),impl/name)
    for name in ('mace_hybrid.py','mace_native_field_accounting.py','mace_tinker_framework_solver.py',
                 'mace_tinker_capability.py','mace_amoeba_capability.py','mace_density_multipoles.py'):
        shutil.copyfile(Path(__file__).with_name(name),impl/name)
    exports={}
    for old in tm['tasks']:
        if old['variant']!='all_standard':continue
        d=root/'moments'/old['case_id'];d.mkdir(parents=True)
        task=dict(old)
        for k in ('xyz','mask','key'):
            src=verify(task[k]);shutil.copyfile(src,d/src.name);task[k]=record(d/src.name)
        receipt=native_call(sw['executable'],task,'moments',d/'native.log',1)
        result=dict(task=task,receipt=receipt,status='failed')
        try:
            if receipt['returncode']:raise InvalidArtifact('native moment export failed')
            result.update(parse_moments(verify(receipt['log']).read_text(),task))
        except Exception as e:result['failure_reason']=str(e)
        write_new(d/'result.json',result)
        if result['status']!='native_global_moments_exported' or not result['pass_']:
            write_new(root/'preparation_failure.json',dict(reason='native moment qualification failed',
                result=record(d/'result.json'),wall_seconds=time.monotonic()-start,process_CPU_seconds=time.process_time()-cpu))
            raise InvalidArtifact('native moment qualification failed')
        exports[old['case_id']]=record(d/'result.json')
    tasks=[]
    for old in fm['tasks']:
        tid=old['task_id'];ct=next(t for t in cm['tasks'] if t['task_id']==tid);row=fr['rows'][tid]
        if any(old[k]!=ct[k] for k in ('case_id','metal','state','xyz','files','source_receipt')):
            raise InvalidArtifact('field/potential density or preparation differs')
        if fm['utility']!=cm['utilities']['orca_vpot']:raise InvalidArtifact('native potential utility differs')
        accepted_centers=[]
        for d in sorted((Path(centers).parent/'potential_execution'/tid).glob('attempt_*')):
            rp=d/'execution.json'
            if not rp.exists():continue
            r=read_json(rp)
            if r['status']=='complete' and r['task']==ct and r['manifest']==record(centers):
                u=r['vpot']
                if u['returncode'] or u['slurm_job_id']!=r['slurm_job_id'] or not r['slurm_job_id']:continue
                for pin in (r['potential'],u['log'],u['resource_usage']):verify(pin)
                accepted_centers.append((rp,r))
        if not accepted_centers:raise InvalidArtifact('missing real center potential')
        rp,r=accepted_centers[-1];data=read_json(verify(old['probes']));weights=read_json(verify(ct['weights']))
        pos=np.loadtxt(verify(ct['points']),skiprows=1)
        if weights['physical_ids']!=data['physical_ids'] or not np.allclose(pos,data['positions_bohr'],atol=1e-11,rtol=0):
            raise InvalidArtifact('center atom order/geometry differs')
        gid='GGR_1GLG' if old['case_id'].startswith('GGR_') else old['case_id']
        native=read_json(verify(exports[gid]));mapping=read_json(verify(native['task']['mapping']))
        amap=dict(zip(mapping['system_atom_ids'],native['atoms']));pmap={a['id']:a for a in mapping['physical_atoms']}
        glob=np.array([amap[i]['global_'] for i in data['physical_ids']])
        if any(not np.allclose(np.array(pmap[pid]['xyz_A'])/BOHR_TO_A,xyz,atol=1e-12,rtol=0)
               for pid,xyz in zip(data['physical_ids'],data['positions_bohr'])):
            raise InvalidArtifact('native multipole observation coordinates differ')
        inp=root/'inputs'/tid;inp.mkdir(parents=True)
        ps=points(data['positions_bohr']);pp=inp/'points_bohr.xyz'
        pp.write_text(str(len(ps))+'\n'+''.join(' '.join(format(float(v),'.14f') for v in p)+'\n' for p in ps))
        with np.load(verify(row['arrays'])) as a:ef=a['exact'].copy()
        np.savez_compressed(inp/'inputs.npz',center_phi=parse_potential(verify(r['potential']),pos),exact_field=ef,
            charge=glob[:,0],dipole=glob[:,1:4]/BOHR_TO_A,quadrupole=glob[:,4:].reshape(-1,3,3)/BOHR_TO_A**2)
        tasks.append({**old,'points':record(pp),'multipole_inputs':record(inp/'inputs.npz'),
            'native_moments':exports[gid],'center_execution':record(rp),'field_arrays':row['arrays'],
            'field_execution':row['execution_receipt']})
        tasks[-1].pop('cache_key')
    m=dict(protocol=PROTOCOL,steps_bohr=list(STEPS),tolerances=TOL,plan=record(plan),tasks=tasks,moment_exports=exports,
        sources=dict(fields=record(fields),field_manifest=fr['manifest'],centers=record(centers),frameworks=record(frameworks),moment_software=record(moment_software)),
        implementation={p.name:record(p) for p in sorted(impl.glob('*.py'))},utility=fm['utility'],utility_calls=8,
        moment_export_calls=3,new_DFT=0,new_MACE=0,new_FF_energy=0,new_nuclear_Hessian=0,
        preparation_receipt=dict(wall_seconds=time.monotonic()-start,process_CPU_seconds=time.process_time()-cpu),**NULLS)
    for task in tasks:task['cache_key']=key(task,m)
    write_new(root/'manifest.json',m);return validate(root/'manifest.json')


def validate(path):
    m=read_json(path)
    if m['protocol']!=PROTOCOL or m['steps_bohr']!=list(STEPS) or m['tolerances']!=TOL or len(m['tasks'])!=8 or len(m['moment_exports'])!=3:
        raise InvalidArtifact('changed density multipole configuration')
    for pin in [m['plan'],m['utility'],*m['sources'].values(),*m['implementation'].values(),*m['moment_exports'].values()]:verify(pin)
    sw=read_json(verify(m['sources']['moment_software']))
    for k in ('executable','source','library','log','parent_software'):verify(sw[k])
    for pin in m['moment_exports'].values():
        r=read_json(verify(pin));verify(r['receipt']['log'])
        if r['status']!='native_global_moments_exported' or not r['pass_']:raise InvalidArtifact('unqualified native moments')
    pairs={}
    for t in m['tasks']:
        for k in ('probes','points','state','framework_mapping','xyz','charge_execution','source_receipt','projection',
                  'multipole_inputs','native_moments','center_execution','field_arrays','field_execution'):verify(t[k])
        for pin in t['files'].values():verify(pin)
        if t['cache_key']!=key(t,m):raise InvalidArtifact('density multipole cache key changed')
        data=read_json(verify(t['probes']));p=np.loadtxt(verify(t['points']),skiprows=1);expected=points(data['positions_bohr'])
        if p.shape!=expected.shape or not np.allclose(p,expected,atol=2e-13,rtol=0):raise InvalidArtifact('observation stencil differs')
        pairs.setdefault(t['case_id'],[]).append({k:data[k] for k in ('physical_ids','positions_bohr')})
    if len(pairs)!=4 or any(len(v)!=2 or v[0]!=v[1] for v in pairs.values()):raise InvalidArtifact('paired observation identity differs')
    return m


def execute(path,retry=False):
    m=validate(path);mp=Path(path).resolve();root=mp.parent/'execution';root.mkdir(exist_ok=True)
    if not os.environ.get('SLURM_JOB_ID') or int(os.environ.get('SLURM_NTASKS','0'))!=8:
        raise InvalidArtifact('eight-worker CPU allocation required')
    def one(t):
        directory=root/t['task_id'];directory.mkdir(exist_ok=True);attempts=sorted(directory.glob('attempt_*'))
        for a in reversed(attempts):
            if accepted_utility(a,t,mp):return dict(task_id=t['task_id'],status='reused')
        if attempts and not retry:return dict(task_id=t['task_id'],status='failed_attempt_requires_explicit_retry')
        d=directory/f'attempt_{len(attempts)+1:04d}';d.mkdir();start=time.monotonic()
        r=dict(task=t,manifest=record(mp),job_id=os.environ['SLURM_JOB_ID'],host=socket.gethostname(),status='failed')
        try:
            for pin in t['files'].values():shutil.copyfile(verify(pin),d/Path(pin['path']).name)
            r['utility']=run_command([str(verify(m['utility'])),str(d/'endpoint.runtime.gbw'),'endpoint.runtime.scfp',
                str(verify(t['points'])),str(d/'potential.out')],d,d/'utility.log',d/'resources.txt')
            if r['utility']['returncode']:raise InvalidArtifact('native potential utility failed')
            parse_potential(d/'potential.out',np.loadtxt(verify(t['points']),skiprows=1))
            for pin in t['files'].values():
                if record(d/Path(pin['path']).name)['sha256']!=pin['sha256']:raise InvalidArtifact('saved density mutated')
            r.update(status='complete',potential=record(d/'potential.out'))
        except Exception as e:r['failure_reason']=str(e)
        r['wall_seconds']=time.monotonic()-start;write_new(d/'execution.json',r)
        return dict(task_id=t['task_id'],status=r['status'],receipt=record(d/'execution.json'))
    with (root/'execute.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        with ThreadPoolExecutor(max_workers=8) as pool:rows=list(pool.map(one,m['tasks']))
    return dict(rows=rows,complete=all(r['status'] in ('complete','reused') for r in rows))


def collect(path,output):
    m=validate(path);mp=Path(path).resolve();out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False)
    rows={};pairs={};attempts=[]
    for t in m['tasks']:
        good=[];tid=t['task_id']
        for d in sorted((mp.parent/'execution'/tid).glob('attempt_*')):
            r=accepted_utility(d,t,mp);attempts.append(dict(task_id=tid,attempt=str(d),accepted=r is not None,
                receipt=record(d/'execution.json') if (d/'execution.json').exists() else None))
            if r:good.append((d,r))
        if not good:rows[tid]=dict(status='unavailable');continue
        d,r=good[-1];data=read_json(verify(t['probes']));n=len(data['physical_ids'])
        v=parse_potential(verify(r['potential']),np.loadtxt(verify(t['points']),skiprows=1));phi,H=hessian(v,n)
        with np.load(verify(t['multipole_inputs'])) as a:inputs={k:a[k].copy() for k in a.files}
        q,mu,Q,E=[inputs[k] for k in ('charge','dipole','quadrupole','exact_field')]
        components=[coupling(q,mu,Q,phi,E,h) for h in H]
        rot=rotation();Qr=np.einsum('ij,njk,lk->nil',rot,Q,rot);Hr=np.einsum('ij,njk,lk->nil',rot,H[1],rot)
        rotated=coupling(q,mu@rot.T,Qr,phi,E@rot.T,Hr)
        rigid={k:rotated[k]-components[1][k] for k in rotated}
        center_error=float(np.max(abs(phi-inputs['center_phi'])))
        delta=components[1]['quadrupole']-components[0]['quadrupole']
        analytic=monopole_hessian(data['projected_charge_e'],data['projected_positions_bohr'],data['positions_bohr'])
        hdiff=H[1]-H[0];hdiff_proj=analytic-H[1]
        np.savez_compressed(out/(tid+'.npz'),phi=phi,H_coarse=H[0],H=H[1],projected_monopole_H=analytic)
        rows[tid]=dict(status='computed_component_diagnostic',center_potential_error_au=center_error,
            center_pass=center_error<=TOL['center_potential_au'],quadrupole_refinement_kcal_mol=delta,
            quadrupole_refinement_pass=abs(delta)<=TOL['endpoint_quadrupole_refinement_kcal'],
            max_H_refinement_au=float(np.max(abs(hdiff))),RMS_H_refinement_au=float(np.sqrt(np.mean(hdiff**2))),
            max_projected_monopole_H_error_au=float(np.max(abs(hdiff_proj))),
            RMS_projected_monopole_H_error_au=float(np.sqrt(np.mean(hdiff_proj**2))),
            rigid_component_errors_kcal_mol=rigid,rigid_pass=max(abs(v) for v in rigid.values())<=TOL['scalar_rigid_kcal'],
            components_coarse_kcal_mol=components[0],components_fine_kcal_mol=components[1],
            total_component_diagnostic_kcal_mol=sum(components[1][k] for k in ('charge','dipole','quadrupole')),
            projected_monopole_quadrupole_component_kcal_mol=coupling(q,mu,Q,phi,E,analytic)['quadrupole'],
            arrays=record(out/(tid+'.npz')),execution_receipt=record(d/'execution.json'),potential=r['potential'])
    for case in sorted({t['case_id'] for t in m['tasks']}):
        ca,la=[rows[case+'_'+metal] for metal in ('Ca','La')]
        if ca['status']=='unavailable' or la['status']=='unavailable':pairs[case]=dict(status='unavailable');continue
        delta=ca['quadrupole_refinement_kcal_mol']-la['quadrupole_refinement_kcal_mol']
        pairs[case]=dict(status='computed_component_diagnostic',quadrupole_refinement_kcal_mol=delta,
            quadrupole_refinement_pass=abs(delta)<=TOL['pair_quadrupole_refinement_kcal'],
            Ca_minus_La_component_diagnostic_kcal_mol={k:ca['components_fine_kcal_mol'][k]-la['components_fine_kcal_mol'][k] for k in ca['components_fine_kcal_mol']},
            projected_monopole_quadrupole_error_kcal_mol=(ca['projected_monopole_quadrupole_component_kcal_mol']-la['projected_monopole_quadrupole_component_kcal_mol'])-
                (ca['components_fine_kcal_mol']['quadrupole']-la['components_fine_kcal_mol']['quadrupole']))
    complete=all(r['status']!='unavailable' for r in rows.values())
    result=dict(protocol=PROTOCOL,manifest=record(mp),rows=rows,pairs=pairs,attempts=attempts,complete=complete,
        numerical_pass=complete and all(r['center_pass'] and r['quadrupole_refinement_pass'] and r['rigid_pass'] for r in rows.values()) and
                       all(p['quadrupole_refinement_pass'] for p in pairs.values()),
        units=dict(phi='Hartree/e',field='Hartree/(e bohr)',H='Hartree/(e bohr^2)',components='kcal/mol'),
        scope='Unchanged exterior native multipole diagnostic; not charge-closed hybrid/environment correction',**NULLS)
    write_new(out/'result.json',result);return result


def main():
    p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest='command',required=True)
    s=sub.add_parser('build-moments')
    for name in ('parent-software','source','output'):s.add_argument('--'+name,required=True)
    s=sub.add_parser('prepare')
    for name in ('fields','centers','frameworks','moment-software','plan','output'):s.add_argument('--'+name,required=True)
    for name in ('dry-run','execute','collect'):
        s=sub.add_parser(name);s.add_argument('--manifest',required=True)
        if name!='dry-run':s.add_argument('--output',required=True)
        if name=='execute':s.add_argument('--retry-failed',action='store_true')
    a=p.parse_args()
    if a.command=='build-moments':
        r=build(a.parent_software,a.source,a.output);print(json.dumps({'returncode':r['returncode']}));raise SystemExit(r['returncode'])
    elif a.command=='prepare':print(json.dumps({'tasks':len(prepare(a.fields,a.centers,a.frameworks,a.moment_software,a.plan,a.output)['tasks'])}))
    elif a.command=='dry-run':print(json.dumps({'tasks':len(validate(a.manifest)['tasks']),'new_native_calls':0}))
    elif a.command=='execute':r=execute(a.manifest,a.retry_failed);write_new(a.output,r);print(json.dumps(r))
    else:
        r=collect(a.manifest,a.output);print(json.dumps({k:r[k] for k in ('complete','numerical_pass')}))


if __name__=='__main__':main()
