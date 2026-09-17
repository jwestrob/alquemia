"""Constrained charge/dipole source fit; independent spatial validation, no score."""
from __future__ import annotations
import argparse
from concurrent.futures import ThreadPoolExecutor
import fcntl
import json
import os
from pathlib import Path
import resource
import shutil
import socket
import time
import numpy as np
import scipy
from scipy.linalg import null_space,svd
from affordable_common import BOHR_TO_A,InvalidArtifact,cache_key,read_json,record,verify,write_new
from affordable_solver import run_command
from density_embedding import parse_potential
from mace_hybrid import rotation
from mace_qm_field import (STEPS,TOL as FIELD_TOL,NULLS,offset_points,derivative,metric,
    point_field,diagonal_response,accepted as accepted_utility,validate as validate_fields)

PROTOCOL='vacuum_density_physical_charge_dipole_source_v1'
SAMPLING_PROTOCOL='vacuum_density_physical_charge_dipole_two_center_source_v2'
ELL=1/BOHR_TO_A
SHIFT_A=np.array([.17,-.11,.13])
TOL={**FIELD_TOL,'potential_RMS_au':.005,'potential_relative_RMS':.10,'charge_e':1e-9,'replay_au':1e-8}
MODEL=dict(ell_A=1.,regularization_fraction=1e-4,field_row_weights='unweighted',
    source_sites='unchanged_physical_projection_support',parameters='charge_and_Cartesian_dipole',
    total_charge_constraint='formal_QM_charge',prior='projected_CHELPG_plus_zero_dipoles')


def model_for(protocol):
    if protocol==PROTOCOL:return MODEL
    if protocol==SAMPLING_PROTOCOL:return {**MODEL,'training_centers':'original_and_positive_offset',
        'validation_shift_A':(-SHIFT_A).tolist()}
    raise InvalidArtifact('unknown distributed-source protocol')


def matrices(source,probes):
    """x=[q(e),mu/(e*ell)]; return phi and Cartesian field linear maps."""
    d=np.asarray(probes)[:,None,:]-np.asarray(source)[None,:,:];r=np.linalg.norm(d,axis=2)
    if np.min(r)<1e-8:raise InvalidArtifact('source/probe overlap')
    n,k=r.shape
    p=np.concatenate([1/r,(ELL*d/r[:,:,None]**3).reshape(n,3*k)],axis=1)
    tensor=ELL*(3*d[:,:,:,None]*d[:,:,None,:]/r[:,:,None,None]**5-np.eye(3)[None,None,:,:]/r[:,:,None,None]**3)
    e=np.concatenate([(d/r[:,:,None]**3).transpose(0,2,1),tensor.transpose(0,2,1,3).reshape(n,3,3*k)],axis=2)
    return p,e


def evaluate(source,probes,charge,dipole_e_bohr):
    p,e=matrices(source,probes);x=np.r_[charge,np.asarray(dipole_e_bohr).ravel()/ELL]
    return p@x,np.einsum('ijk,k->ij',e,x)


def fit_source(data,potential,field,formal_charge):
    start=time.monotonic();cpu=time.thread_time()
    src=np.array(data['projected_positions_bohr']);k=len(src)
    p,e=matrices(src,data['positions_bohr']);a=np.vstack([p,ELL*e.reshape(-1,4*k)])
    y=np.r_[potential,ELL*np.asarray(field).ravel()];prior=np.r_[data['projected_charge_e'],np.zeros(3*k)]
    c=np.r_[np.ones(k),np.zeros(3*k)];z=null_space(c[None,:])
    constrained=prior+(formal_charge-c@prior)*c/(c@c)
    an=a@z;u,s,vh=svd(an,full_matrices=False,check_finite=True)
    lam=MODEL['regularization_fraction']*s[0]
    x=constrained+z@(vh.T@((s/(s*s+lam*lam))*(u.T@(y-a@constrained))))
    q=x[:k];mu=x[k:].reshape(k,3)*ELL;closure=float(q.sum()-formal_charge)
    if abs(closure)>TOL['charge_e'] or not np.all(np.isfinite(x)):raise InvalidArtifact('fit charge/nonfinite failure')
    pred=p@x;ef=np.einsum('ijk,k->ij',e,x)
    # A transformed evaluation of these same moments is not another fit.
    rot=rotation();shift=np.array([17.,-23.,9.])/BOHR_TO_A
    phi_r,field_r=evaluate(src@rot.T+shift,np.array(data['positions_bohr'])@rot.T+shift,q,mu@rot.T)
    rigid=dict(potential_max_au=float(np.max(np.abs(phi_r-pred))),
        field_max_vector_au=float(np.max(np.linalg.norm(field_r-ef@rot.T,axis=1))))
    if max(rigid.values())>TOL['replay_au']:raise InvalidArtifact('fitted representation rigid replay failed')
    r=dict(status='fitted_density_representation',charge_e=q.tolist(),dipole_e_bohr=mu.tolist(),
        formal_charge=formal_charge,charge_residual_e=closure,lambda_=float(lam),singular_values=s.tolist(),
        effective_rank=int(sum(s>lam)),effective_degrees_of_freedom=float(sum(s*s/(s*s+lam*lam))),
        max_absolute_charge_e=float(max(abs(q))),max_dipole_eA=float(max(np.linalg.norm(mu,axis=1))*BOHR_TO_A),
        prior_displacement_norm_e=float(np.linalg.norm(x-prior)),training_scaled_residual_norm_au=float(np.linalg.norm(a@x-y)),
        rigid_replay=rigid,wall_seconds=time.monotonic()-start,thread_CPU_seconds=time.thread_time()-cpu,
        shared_process_RSS_highwater_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,**NULLS)
    return r,pred,ef


def potential_metric(pred,exact):
    rms=float(np.sqrt(np.mean((pred-exact)**2)));norm=float(np.sqrt(np.mean(exact**2)))
    return dict(RMS_au=rms,relative_RMS=rms/norm if norm else None,
        pass_=rms<=TOL['potential_RMS_au'] or norm>0 and rms/norm<=TOL['potential_relative_RMS'])


def key(task,m):
    return cache_key(dict(task={k:v for k,v in task.items() if k!='cache_key'},
        **{k:m[k] for k in ('protocol','model','tolerances','plan','sources','implementation','utility','versions')}))


def prepare(fields,center_manifest,plan,output,additional_training_report=None):
    fr=read_json(fields);fm=validate_fields(verify(fr['manifest']));cm=read_json(center_manifest)
    if not fr['complete'] or not fr['numerical_pass']:raise InvalidArtifact('qualified actual field source required')
    if fr['manifest']['sha256']!='0b639b8502fdade3b8dd7727563908cbdf1d2237ae4008b8cfe329ce872d0673':
        raise InvalidArtifact('declared field experiment differs')
    if record(center_manifest)['sha256']!='90d553d2480a1efb3911fdc9169774610c9b124d40db7e17e6b4ad45af32d439':
        raise InvalidArtifact('declared center potentials differ')
    additional=None;previous=None;previous_path=None
    protocol=SAMPLING_PROTOCOL if additional_training_report else PROTOCOL
    if additional_training_report:
        additional=read_json(additional_training_report);previous_path=verify(additional['manifest']);previous=validate(previous_path)
        if (previous['protocol']!=PROTOCOL or not additional['complete'] or not additional['numerical_pass']
                or additional['manifest']['sha256']!='98820483e24705f632a993400f01aaabbd55ffe6551005c8a7db998040e67282'
                or previous['sources']['fields']!=record(fields) or previous['sources']['center_manifest']!=record(center_manifest)):
            raise InvalidArtifact('declared V1 native spatial observations required')
    root=Path(output).resolve();root.mkdir(parents=True,exist_ok=False);impl=root/'implementation';impl.mkdir()
    for name,pin in fm['implementation'].items():shutil.copyfile(verify(pin),impl/name)
    shutil.copyfile(__file__,impl/Path(__file__).name)
    # The prior field implementation obtains this helper through its parent tree.
    if not (impl/'mace_hybrid.py').exists():shutil.copyfile(Path(__file__).with_name('mace_hybrid.py'),impl/'mace_hybrid.py')
    tasks=[];start=time.monotonic();cpu=time.process_time()
    for old in fm['tasks']:
        tid=old['task_id'];ct=next(t for t in cm['tasks'] if t['task_id']==tid);row=fr['rows'][tid]
        if any(old[k]!=ct[k] for k in ('case_id','metal','state','xyz','files','source_receipt')):
            raise InvalidArtifact('field/potential density or preparation differs')
        if fm['utility']!=cm['utilities']['orca_vpot']:raise InvalidArtifact('native potential executable differs')
        for pin in old['files'].values():verify(pin)
        center_attempts=[]
        for d in sorted((Path(center_manifest).parent/'potential_execution'/tid).glob('attempt_*')):
            rp=d/'execution.json'
            if not rp.exists():continue
            r=read_json(rp)
            if r['status']=='complete' and r['task']==ct and r['manifest']==record(center_manifest):
                u=r['vpot']
                if u['returncode'] or u['slurm_job_id']!=r['slurm_job_id'] or not r['slurm_job_id']:continue
                for pin in (r['potential'],u['log'],u['resource_usage']):verify(pin)
                center_attempts.append((rp,r))
        if not center_attempts:raise InvalidArtifact('missing accepted center potential')
        rp,r=center_attempts[-1];data=read_json(verify(old['probes']));weights=read_json(verify(ct['weights']))
        centers=np.loadtxt(verify(ct['points']),skiprows=1)
        if weights['physical_ids']!=data['physical_ids'] or not np.allclose(centers,data['positions_bohr'],atol=1e-11,rtol=0):
            raise InvalidArtifact('center potential physical identity/order differs')
        potential=parse_potential(verify(r['potential']),centers)
        with np.load(verify(row['arrays'])) as arr:field=arr['exact'].copy()
        projection=read_json(verify(old['projection']))
        if len(projection['physical_ids'])!=len(data['projected_charge_e']):raise InvalidArtifact('source identity coverage differs')
        physical={a['id']:a for a in read_json(verify(old['state']))['physical_atoms']}
        for i,pid in enumerate(projection['physical_ids']):
            atom=physical[pid]
            if not np.allclose(np.array(atom['xyz_A'])/BOHR_TO_A,data['projected_positions_bohr'][i],atol=1e-12,rtol=0):
                raise InvalidArtifact('source physical coordinates differ')
        training_positions=np.array(data['positions_bohr']);extra_pins={}
        if additional is not None:
            previous_task=next(t for t in previous['tasks'] if t['task_id']==tid)
            extra_rp=verify(additional['rows'][tid]['execution_receipt'])
            extra=accepted(extra_rp.parent,previous_task,previous_path)
            if extra is None or previous_task['probes']!=old['probes'] or previous_task['files']!=old['files'] or previous_task['formal_charge']!=ct['charge']:
                raise InvalidArtifact('additional native training state differs')
            extra_points=np.loadtxt(verify(previous_task['points']),skiprows=1);n=len(training_positions)
            if not np.allclose(extra_points[:n],training_positions+SHIFT_A/BOHR_TO_A,atol=2e-13,rtol=0):
                raise InvalidArtifact('additional native spatial centers differ')
            native=parse_potential(verify(extra['potential']),extra_points)
            training_positions=np.vstack([training_positions,extra_points[:n]])
            potential=np.r_[potential,native[:n]];field=np.vstack([field,derivative(native[n:],n)[1]])
            extra_pins=dict(additional_training_receipt=record(extra_rp),additional_training_potential=extra['potential'],
                additional_training_points=previous_task['points'])
        d=root/'inputs'/tid;d.mkdir(parents=True)
        np.savez_compressed(d/'training.npz',potential=potential,field=field,positions_bohr=training_positions)
        shift=SHIFT_A if protocol==PROTOCOL else -SHIFT_A
        probe=np.array(data['positions_bohr'])+shift/BOHR_TO_A
        points=np.vstack([probe,offset_points(probe)])
        (d/'points_bohr.xyz').write_text(str(len(points))+'\n'+''.join(' '.join(format(float(v),'.14f') for v in p)+'\n' for p in points))
        tasks.append(dict(task_id=tid,case_id=old['case_id'],metal=old['metal'],formal_charge=ct['charge'],
            source_ids=projection['physical_ids'],probes=old['probes'],state=old['state'],source_receipt=old['source_receipt'],
            files=old['files'],training=record(d/'training.npz'),points=record(d/'points_bohr.xyz'),
            center_receipt=record(rp),field_receipt=row['execution_receipt'],field_arrays=row['arrays'],**extra_pins))
    m=dict(protocol=protocol,model=model_for(protocol),tolerances=TOL,plan=record(plan),
        sources=dict(fields=record(fields),center_manifest=record(center_manifest),field_manifest=fr['manifest']),
        implementation={p.name:record(p) for p in sorted(impl.glob('*.py'))},utility=fm['utility'],
        versions=dict(numpy=np.__version__,scipy=scipy.__version__),tasks=tasks,
        new_fits=8,new_potential_calls=8,new_DFT=0,new_MACE=0,new_FF_energy=0,
        preparation=dict(wall_seconds=time.monotonic()-start,process_CPU_seconds=time.process_time()-cpu),**NULLS)
    if additional_training_report:m['sources']['development_spatial_report']=record(additional_training_report)
    for t in tasks:t['cache_key']=key(t,m)
    write_new(root/'manifest.json',m);return validate(root/'manifest.json')


def validate(path):
    m=read_json(path)
    if m['model']!=model_for(m['protocol']) or m['tolerances']!=TOL or len(m['tasks'])!=8:
        raise InvalidArtifact('distributed-source configuration differs')
    if m['versions']!=dict(numpy=np.__version__,scipy=scipy.__version__):raise InvalidArtifact('numerical library version differs')
    for pin in [m['plan'],m['utility'],*m['sources'].values(),*m['implementation'].values()]:verify(pin)
    pairs={}
    for t in m['tasks']:
        for k in ('probes','state','source_receipt','training','points','center_receipt','field_receipt','field_arrays'):verify(t[k])
        if m['protocol']==SAMPLING_PROTOCOL:
            for k in ('additional_training_receipt','additional_training_potential','additional_training_points'):verify(t[k])
        for pin in t['files'].values():verify(pin)
        shift=SHIFT_A if m['protocol']==PROTOCOL else -SHIFT_A
        data=read_json(verify(t['probes']));probe=np.array(data['positions_bohr'])+shift/BOHR_TO_A
        expected=np.vstack([probe,offset_points(probe)]);p=np.loadtxt(verify(t['points']),skiprows=1)
        if p.shape!=expected.shape or not np.allclose(p,expected,atol=2e-13,rtol=0):raise InvalidArtifact('validation probe coordinates differ')
        if t['cache_key']!=key(t,m):raise InvalidArtifact('source scientific cache key differs')
        pairs.setdefault(t['case_id'],[]).append((p,t['source_ids'],t['formal_charge']))
    for v in pairs.values():
        if len(v)!=2 or not np.array_equal(v[0][0],v[1][0]) or v[0][1]!=v[1][1] or abs(v[0][2]-v[1][2])!=1:
            raise InvalidArtifact('paired source geometry/charge differs')
    return m


def accepted(directory,task,path):
    r=accepted_utility(directory,task,path)
    if r:
        for k in ('fit','fit_training_arrays'):verify(r[k])
        fit=read_json(verify(r['fit']))
        if fit['cache_key']!=task['cache_key'] or abs(fit['charge_residual_e'])>TOL['charge_e']:return None
    return r


def execute(path,retry=False):
    m=validate(path);mp=Path(path).resolve();root=mp.parent/'execution';root.mkdir(exist_ok=True)
    if not os.environ.get('SLURM_JOB_ID') or int(os.environ.get('SLURM_NTASKS','0'))!=8:
        raise InvalidArtifact('eight-worker CPU allocation required')
    if any(os.environ.get(k)!='1' for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS')):
        raise InvalidArtifact('one thread per fit/utility required')
    def one(task):
        td=root/task['task_id'];td.mkdir(exist_ok=True);attempts=sorted(td.glob('attempt_*'))
        for d in reversed(attempts):
            if accepted(d,task,mp):return dict(task_id=task['task_id'],status='reused')
        if attempts and not retry:return dict(task_id=task['task_id'],status='failed_attempt_requires_explicit_retry')
        d=td/f'attempt_{len(attempts)+1:04d}';d.mkdir();start=time.monotonic()
        r=dict(task=task,manifest=record(mp),job_id=os.environ['SLURM_JOB_ID'],host=socket.gethostname(),status='failed')
        try:
            data=read_json(verify(task['probes']))
            with np.load(verify(task['training'])) as arr:
                phi=arr['potential'];field=arr['field']
                fit_data={**data,'positions_bohr':arr['positions_bohr']} if 'positions_bohr' in arr else data
            r['fit_attempted']=True
            fit,p,e=fit_source(fit_data,phi,field,task['formal_charge']);fit['cache_key']=task['cache_key']
            write_new(d/'fit.json',fit);np.savez_compressed(d/'fit_training.npz',potential=p,field=e)
            r.update(fit=record(d/'fit.json'),fit_training_arrays=record(d/'fit_training.npz'))
            # Freeze fitted moments before starting the separate validation call.
            for pin in task['files'].values():shutil.copyfile(verify(pin),d/Path(pin['path']).name)
            r['utility']=run_command([str(verify(m['utility'])),str(d/'endpoint.runtime.gbw'),'endpoint.runtime.scfp',
                str(verify(task['points'])),str(d/'potential.out')],d,d/'utility.log',d/'resources.txt')
            if r['utility']['returncode']:raise InvalidArtifact('native validation potential failed')
            parse_potential(d/'potential.out',np.loadtxt(verify(task['points']),skiprows=1))
            for pin in task['files'].values():
                if record(d/Path(pin['path']).name)['sha256']!=pin['sha256']:raise InvalidArtifact('saved density mutated')
            r.update(status='complete',potential=record(d/'potential.out'))
        except Exception as exc:r['failure_reason']=str(exc)
        r['wall_seconds']=time.monotonic()-start;write_new(d/'execution.json',r)
        return dict(task_id=task['task_id'],status=r['status'],receipt=record(d/'execution.json'))
    with (root/'execute.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        with ThreadPoolExecutor(max_workers=8) as pool:rows=list(pool.map(one,m['tasks']))
    return dict(rows=rows,complete=all(r['status'] in ('complete','reused') for r in rows))


def compare_arrays(arr,alpha):
    return dict(field_quality={k:metric(arr[k+'_field'],arr['exact_field'],alpha) for k in ('distributed','fitted','projected')},
        potential_quality={k:potential_metric(arr[k+'_potential'],arr['exact_potential']) for k in ('distributed','fitted','projected')},
        U0_kcal_mol={k:diagonal_response(arr[k+'_field'],alpha) for k in ('exact','coarse','distributed','fitted','projected')})


def collect(path,output):
    m=validate(path);mp=Path(path).resolve();out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False)
    rows={};vectors={};attempts=[]
    for t in m['tasks']:
        tid=t['task_id'];good=[]
        for d in sorted((mp.parent/'execution'/tid).glob('attempt_*')):
            r=accepted(d,t,mp);attempts.append(dict(task_id=tid,directory=str(d),accepted=r is not None,
                receipt=record(d/'execution.json') if (d/'execution.json').exists() else None))
            if r:good.append((d,r))
        if not good:rows[tid]=dict(status='unavailable');continue
        d,r=good[-1];fit=read_json(verify(r['fit']));data=read_json(verify(t['probes']));alpha=np.array(data['alpha_bohr3'])
        n=len(alpha);points=np.loadtxt(verify(t['points']),skiprows=1);pot=parse_potential(verify(r['potential']),points)
        fields=derivative(pot[n:],n);probe=points[:n]
        distributed_phi,distributed_field=evaluate(data['projected_positions_bohr'],probe,fit['charge_e'],fit['dipole_e_bohr'])
        arr=dict(exact_potential=pot[:n],exact_field=fields[1],coarse_field=fields[0],
            distributed_potential=distributed_phi,distributed_field=distributed_field)
        for label,q,src in [('fitted',data['charge_e'],data['QM_positions_bohr']),('projected',data['projected_charge_e'],data['projected_positions_bohr'])]:
            arr[label+'_potential'],arr[label+'_field']=evaluate(src,probe,q,np.zeros((len(q),3)))
        np.savez_compressed(out/(tid+'.npz'),**arr);vectors[tid]=(arr,alpha)
        with np.load(verify(t['training'])) as src,np.load(verify(r['fit_training_arrays'])) as pred:
            factor=1 if m['protocol']==PROTOCOL else 2
            if pred['field'].shape!=(factor*n,3):raise InvalidArtifact('training observation count differs')
            training=dict(field_quality=metric(pred['field'],src['field'],np.tile(alpha,factor)),potential_quality=potential_metric(pred['potential'],src['potential']))
        numerical=float(np.max(np.linalg.norm(fields[1]-fields[0],axis=1)))
        rows[tid]=dict(status='computed_representation',validation=compare_arrays(arr,alpha),training=training,
            numerical_field_max_change_au=numerical,numerical_pass=numerical<=TOL['numerical_field_max_au'],
            fit=r['fit'],execution_receipt=record(d/'execution.json'),arrays=record(out/(tid+'.npz')),
            max_absolute_charge_e=fit['max_absolute_charge_e'],max_dipole_eA=fit['max_dipole_eA'],
            charge_residual_e=fit['charge_residual_e'],rigid_replay=fit['rigid_replay'])
    pairs={}
    for case in sorted(set(t['case_id'] for t in m['tasks'])):
        ca,la=case+'_Ca',case+'_La'
        if ca not in vectors or la not in vectors:pairs[case]=dict(status='unavailable');continue
        ac,alpha=vectors[ca];al,_=vectors[la];delta={k:ac[k]-al[k] for k in ac}
        checks=compare_arrays(delta,alpha)
        # The difference of response energies is not the response to the difference field.
        u={k:rows[ca]['validation']['U0_kcal_mol'][k]-rows[la]['validation']['U0_kcal_mol'][k] for k in ('exact','coarse','distributed','fitted','projected')}
        checks.pop('U0_kcal_mol')
        errors={k:u[k]-u['exact'] for k in ('distributed','fitted','projected')}
        pairs[case]=dict(status='computed_representation',**checks,Ca_minus_La_U0_kcal_mol=u,
            paired_U0_errors_kcal_mol=errors,paired_U0_flag_pass=abs(errors['distributed'])<=TOL['paired_U0_error_kcal'],
            numerical_pair_U0_change_kcal_mol=u['exact']-u['coarse'],
            numerical_pass=abs(u['exact']-u['coarse'])<=TOL['numerical_pair_U0_kcal'])
    complete=len(vectors)==8
    result=dict(protocol=m['protocol'],manifest=record(mp),rows=rows,pairs=pairs,attempts=attempts,complete=complete,
        numerical_pass=complete and all(x['numerical_pass'] for x in [*rows.values(),*pairs.values()]),
        heldout_field_screen_pass=complete and all(x['validation']['field_quality']['distributed']['pass_'] for x in rows.values()) and all(x['field_quality']['distributed']['pass_'] for x in pairs.values()),
        heldout_potential_screen_pass=complete and all(x['validation']['potential_quality']['distributed']['pass_'] for x in rows.values()) and all(x['potential_quality']['distributed']['pass_'] for x in pairs.values()),
        paired_U0_flags_pass=complete and all(x['paired_U0_flag_pass'] for x in pairs.values()),
        U0_scope='undamped diagonal diagnostic; not environmental energy',**NULLS)
    write_new(out/'result.json',result);return result


def main():
    p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest='op',required=True)
    s=sub.add_parser('prepare')
    for name in ('fields','center-manifest','plan','output'):s.add_argument('--'+name,required=True)
    s.add_argument('--additional-training-report')
    for name in ('dry-run','execute','collect'):
        s=sub.add_parser(name);s.add_argument('--manifest',required=True)
        if name!='dry-run':s.add_argument('--output',required=True)
        if name=='execute':s.add_argument('--retry-failed',action='store_true')
    a=p.parse_args()
    if a.op=='prepare':print(json.dumps({'tasks':len(prepare(a.fields,a.center_manifest,a.plan,a.output,a.additional_training_report)['tasks'])}))
    elif a.op=='dry-run':print(json.dumps({'tasks':len(validate(a.manifest)['tasks']),'new_fits':0}))
    elif a.op=='execute':r=execute(a.manifest,a.retry_failed);write_new(a.output,r);print(json.dumps(r))
    else:
        r=collect(a.manifest,a.output);print(json.dumps({k:r[k] for k in ('complete','numerical_pass','heldout_field_screen_pass','heldout_potential_screen_pass','paired_U0_flags_pass')}))


if __name__=='__main__':main()
