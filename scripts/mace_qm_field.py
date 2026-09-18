"""Real saved-density electric-field diagnostic; never a scoring correction."""
from __future__ import annotations
import argparse
from concurrent.futures import ThreadPoolExecutor
import fcntl
import json
import os
from pathlib import Path
import shutil
import socket
import time
import numpy as np
from affordable_common import BOHR_TO_A,HA_TO_KCAL,InvalidArtifact,cache_key,read_json,record,verify,write_new,xyz
from affordable_solver import run_command
from density_embedding import parse_potential
from mace_omol_charges import validate as validate_charges

PROTOCOL='normalized_vacuum_QM_electric_field_diagnostic_v1'
CHARGES_SHA='2689bcdfa9b3f2d0358c5e21fa8b5b38533b36502c7e5c3f9826ee7487041d10'
STEPS=(.001,.0005)
TOL=dict(numerical_field_max_au=1e-6,numerical_pair_U0_kcal=.01,field_RMS_au=1e-4,
         field_relative_RMS=.10,paired_U0_error_kcal=1.,far_distance_A=3.)
NULLS=dict(numerical_score=None,environment_correction=None,classification=None,baseline_changed=False)


def offset_points(positions_bohr):
    p=np.asarray(positions_bohr);axis=np.eye(3)
    return np.array([p[:,None,None,:]+h*axis[None,:,None,:]*np.array([1.,-1.])[None,None,:,None]
                     for h in STEPS]).reshape(-1,3)


def derivative(values,n):
    v=np.asarray(values).reshape(2,n,3,2)
    return -(v[:,:,:,0]-v[:,:,:,1])/(2*np.array(STEPS)[:,None,None])


def point_field(charges,source_bohr,probes_bohr):
    delta=np.asarray(probes_bohr)[:,None,:]-np.asarray(source_bohr)[None,:,:]
    r=np.linalg.norm(delta,axis=2)
    if np.min(r)<1e-8:raise InvalidArtifact('source/probe overlap')
    return np.sum(np.asarray(charges)[None,:,None]*delta/r[:,:,None]**3,axis=1)


def metric(candidate,exact,alpha):
    error=np.linalg.norm(candidate-exact,axis=1);norm=np.linalg.norm(exact,axis=1)
    wrms=float(np.sqrt(np.dot(alpha,error**2)/sum(alpha)))
    wnorm=float(np.sqrt(np.dot(alpha,norm**2)/sum(alpha)))
    return dict(atom_count=len(alpha),RMS_vector_error_au=float(np.sqrt(np.mean(error**2))),
        max_vector_error_au=float(max(error)),weighted_RMS_error_au=wrms,
        relative_weighted_RMS=None if wnorm==0 else wrms/wnorm,
        pass_=wrms<=TOL['field_RMS_au'] or wnorm>0 and wrms/wnorm<=TOL['field_relative_RMS'])


def diagonal_response(field,alpha):
    return -.5*float(np.einsum('i,ij,ij->',alpha,field,field))*HA_TO_KCAL


def task_key(task,m):
    return cache_key(dict(task={k:v for k,v in task.items() if k!='cache_key'},
        **{k:m[k] for k in ('protocol','steps_bohr','tolerances','sources','utility','implementation','plan')}))


def prepare(charges,frameworks,states,plan,output):
    if record(charges)['sha256']!=CHARGES_SHA:raise InvalidArtifact('declared charge report required')
    c=read_json(charges);validate_charges(verify(c['manifest']));cm=read_json(verify(c['manifest']))
    if c['status']!='complete' or not c['projection_gate_pass']:raise InvalidArtifact('charge source incomplete')
    fm=read_json(frameworks);by_id={r['case_id']:r for r in fm['cases']}
    root=Path(output).resolve();root.mkdir(parents=True,exist_ok=False);impl=root/'implementation';impl.mkdir()
    for name,pin in cm['implementation'].items():shutil.copyfile(verify(pin),impl/name)
    shutil.copyfile(__file__,impl/Path(__file__).name)
    tasks=[];start=time.monotonic()
    for old in cm['tasks']:
        name,metal=old['case_id'],old['metal'];gid='GGR_1GLG' if name.startswith('GGR_') else name
        fr=by_id[gid];mapping=read_json(verify(fr['mapping']));sp=Path(states)/name/'state.json';s=read_json(sp)
        if mapping['physical_atoms']!=s['physical_atoms']:raise InvalidArtifact('physical framework differs')
        support=set(s['projection_support_ids']);par={p['id']:p for p in mapping['parameters']}
        atoms=[a for a in s['physical_atoms'] if a['id'] not in support]
        alpha=[par[a['id']]['multipole_md_units'][9]*(10/BOHR_TO_A)**3 for a in atoms]
        if not atoms or min(alpha)<=0:raise InvalidArtifact('missing/nonpositive real polarizabilities')
        coords=np.array([a['xyz_A'] for a in atoms])/BOHR_TO_A;qm=xyz(verify(old['xyz']))
        near=np.min(np.linalg.norm(coords[:,None,:]*BOHR_TO_A-np.array([a[1:] for a in qm])[None,:,:],axis=2),axis=1)
        proj=read_json(verify(old['projection']));row=c['rows'][old['task_id']]
        if set(proj['physical_ids'])!=support:raise InvalidArtifact('projected support mismatch')
        if not np.allclose(np.array(proj['weights'])@row['charge_e'],row['projected_charge_e'],atol=1e-14,rtol=0):
            raise InvalidArtifact('charge projection differs')
        data=dict(physical_ids=[a['id'] for a in atoms],positions_bohr=coords.tolist(),alpha_bohr3=alpha,
            nearest_QM_or_cap_A=near.tolist(),charge_e=row['charge_e'],projected_charge_e=row['projected_charge_e'],
            QM_positions_bohr=(np.array([a[1:] for a in qm])/BOHR_TO_A).tolist(),
            projected_positions_bohr=(np.array(proj['coordinates_A'])/BOHR_TO_A).tolist())
        d=root/'inputs'/old['task_id'];d.mkdir(parents=True);write_new(d/'probes.json',data)
        points=offset_points(coords);pp=d/'points_bohr.xyz'
        pp.write_text(str(len(points))+'\n'+''.join(' '.join(format(float(v),'.14f') for v in p)+'\n' for p in points))
        tasks.append(dict(task_id=old['task_id'],case_id=name,metal=metal,probes=record(d/'probes.json'),
            points=record(pp),state=record(sp),framework_mapping=fr['mapping'],xyz=old['xyz'],files=old['files'],
            charge_execution=row['execution_receipt'],source_receipt=old['source_receipt'],projection=old['projection']))
    m=dict(protocol=PROTOCOL,steps_bohr=list(STEPS),tolerances=TOL,plan=record(plan),tasks=tasks,
        sources=dict(charges=record(charges),frameworks=record(frameworks),charge_manifest=c['manifest']),
        implementation={p.name:record(p) for p in sorted(impl.glob('*.py'))},utility=cm['utilities']['orca_vpot'],
        utility_calls=8,new_DFT=0,new_MACE=0,new_charge_fit=0,new_nuclear_gradients=0,
        preparation_wall_seconds=time.monotonic()-start,**NULLS)
    for task in tasks:task['cache_key']=task_key(task,m)
    write_new(root/'manifest.json',m);return validate(root/'manifest.json')


def validate(path):
    m=read_json(path)
    if m['protocol']!=PROTOCOL or m['steps_bohr']!=list(STEPS) or m['tolerances']!=TOL or len(m['tasks'])!=8:
        raise InvalidArtifact('undeclared field configuration')
    if m['sources']['charges']['sha256']!=CHARGES_SHA:raise InvalidArtifact('charge source changed')
    for pin in [m['plan'],m['utility'],*m['sources'].values(),*m['implementation'].values()]:verify(pin)
    pairs={}
    for task in m['tasks']:
        for k in ('probes','points','state','framework_mapping','xyz','charge_execution','source_receipt','projection'):verify(task[k])
        for pin in task['files'].values():verify(pin)
        data=read_json(verify(task['probes']));p=np.loadtxt(verify(task['points']),skiprows=1)
        expected=offset_points(data['positions_bohr'])
        if p.shape!=expected.shape or not np.allclose(p,expected,atol=2e-13,rtol=0):raise InvalidArtifact('field probes changed')
        if task['cache_key']!=task_key(task,m):raise InvalidArtifact('field cache key changed')
        pairs.setdefault(task['case_id'],[]).append({k:data[k] for k in ('physical_ids','positions_bohr','alpha_bohr3','nearest_QM_or_cap_A')})
    if len(pairs)!=4 or any(len(v)!=2 or v[0]!=v[1] for v in pairs.values()):raise InvalidArtifact('paired field geometry differs')
    return m


def accepted(directory,task,path):
    p=directory/'execution.json'
    if not p.exists():return None
    r=read_json(p)
    if r['status']!='complete' or r['task']!=task or r['manifest']!=record(path) or not r['job_id']:return None
    if r['utility']['returncode'] or r['utility']['slurm_job_id']!=r['job_id']:return None
    for pin in (r['potential'],r['utility']['log'],r['utility']['resource_usage']):verify(pin)
    return r


def execute(path,retry=False):
    return execute_prepared(path,validate(path),retry)


def execute_prepared(path,m,retry=False):
    """Shared utility executor; caller supplies its protocol-validated manifest."""
    if read_json(path)!=m:raise InvalidArtifact('validated utility manifest changed')
    mp=Path(path).resolve();root=mp.parent/'execution';root.mkdir(exist_ok=True)
    if not os.environ.get('SLURM_JOB_ID') or int(os.environ.get('SLURM_NTASKS','0'))!=8:
        raise InvalidArtifact('eight-worker CPU allocation required')
    def one(task):
        directory=root/task['task_id'];directory.mkdir(exist_ok=True);attempts=sorted(directory.glob('attempt_*'))
        for a in reversed(attempts):
            if accepted(a,task,mp):return dict(task_id=task['task_id'],status='reused')
        if attempts and not retry:return dict(task_id=task['task_id'],status='failed_attempt_requires_explicit_retry')
        d=directory/f'attempt_{len(attempts)+1:04d}';d.mkdir();start=time.monotonic()
        r=dict(task=task,manifest=record(mp),job_id=os.environ['SLURM_JOB_ID'],host=socket.gethostname(),status='failed')
        try:
            for pin in task['files'].values():shutil.copyfile(verify(pin),d/Path(pin['path']).name)
            r['utility']=run_command([str(verify(m['utility'])),str(d/'endpoint.runtime.gbw'),'endpoint.runtime.scfp',
                str(verify(task['points'])),str(d/'potential.out')],d,d/'utility.log',d/'resources.txt')
            if r['utility']['returncode']:raise InvalidArtifact('native potential utility failed')
            parse_potential(d/'potential.out',np.loadtxt(verify(task['points']),skiprows=1))
            for pin in task['files'].values():
                if record(d/Path(pin['path']).name)['sha256']!=pin['sha256']:raise InvalidArtifact('utility mutated saved density')
            r.update(status='complete',potential=record(d/'potential.out'))
        except Exception as exc:r['failure_reason']=str(exc)
        r['wall_seconds']=time.monotonic()-start;write_new(d/'execution.json',r)
        return dict(task_id=task['task_id'],status=r['status'],receipt=record(d/'execution.json'))
    with (root/'execute.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        with ThreadPoolExecutor(max_workers=8) as pool:rows=list(pool.map(one,m['tasks']))
    return dict(rows=rows,complete=all(r['status'] in ('complete','reused') for r in rows))


def collect(path,output):
    m=validate(path);mp=Path(path).resolve();out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False)
    rows={};arrays={};pairs={};attempts=[]
    for task in m['tasks']:
        tid=task['task_id'];good=[]
        for d in sorted((mp.parent/'execution'/tid).glob('attempt_*')):
            r=accepted(d,task,mp);attempts.append(dict(task_id=tid,attempt=str(d),accepted=r is not None,
                receipt=record(d/'execution.json') if (d/'execution.json').exists() else None))
            if r:good.append((d,r))
        if not good:rows[tid]=dict(status='unavailable');continue
        accepted_dir,r=good[-1];data=read_json(verify(task['probes']));p=np.loadtxt(verify(task['points']),skiprows=1)
        v=parse_potential(verify(r['potential']),p);n=len(data['physical_ids']);f=derivative(v,n)
        alpha=np.array(data['alpha_bohr3']);a=dict(exact_coarse=f[0],exact=f[1],
            fitted=point_field(data['charge_e'],data['QM_positions_bohr'],data['positions_bohr']),
            projected=point_field(data['projected_charge_e'],data['projected_positions_bohr'],data['positions_bohr']))
        if not all(np.all(np.isfinite(x)) for x in a.values()):raise InvalidArtifact('nonfinite field')
        arrays[tid]=(a,data);np.savez_compressed(out/(tid+'.npz'),**a)
        delta=float(np.max(np.linalg.norm(f[1]-f[0],axis=1)));strata={}
        for label,mask in [('all',np.ones(n,dtype=bool)),('far',np.array(data['nearest_QM_or_cap_A'])>=TOL['far_distance_A'])]:
            if not np.any(mask):raise InvalidArtifact('empty declared field stratum')
            strata[label]=dict(atom_count=int(mask.sum()),
                quality={k:metric(a[k][mask],a['exact'][mask],alpha[mask]) for k in ('fitted','projected')},
                U0_kcal_mol={k:diagonal_response(val[mask],alpha[mask]) for k,val in a.items()})
        rows[tid]=dict(status='computed_diagnostic',numerical_max_vector_change_au=delta,
            numerical_pass=delta<=TOL['numerical_field_max_au'],strata=strata,arrays=record(out/(tid+'.npz')),
            execution_receipt=record(accepted_dir/'execution.json'),
            potential=r['potential'],probes=task['probes'],utility_wall_seconds=r['utility']['wall_seconds'])
    for case in sorted(set(t['case_id'] for t in m['tasks'])):
        ca,la=case+'_Ca',case+'_La'
        if ca not in arrays or la not in arrays:pairs[case]=dict(status='unavailable');continue
        ac,data=arrays[ca];al,_=arrays[la];alpha=np.array(data['alpha_bohr3']);pr={}
        for label,mask in [('all',np.ones(len(alpha),dtype=bool)),('far',np.array(data['nearest_QM_or_cap_A'])>=TOL['far_distance_A'])]:
            exact=ac['exact'][mask]-al['exact'][mask]
            u={k:rows[ca]['strata'][label]['U0_kcal_mol'][k]-rows[la]['strata'][label]['U0_kcal_mol'][k] for k in ac}
            errors={k:u[k]-u['exact'] for k in ('fitted','projected')}
            pr[label]=dict(quality={k:metric(ac[k][mask]-al[k][mask],exact,alpha[mask]) for k in ('fitted','projected')},
                Ca_minus_La_U0_kcal_mol=u,paired_U0_error_kcal_mol=errors,
                paired_projected_U0_flag_pass=abs(errors['projected'])<=TOL['paired_U0_error_kcal'],
                numerical_U0_change_kcal_mol=u['exact']-u['exact_coarse'],
                numerical_U0_pass=abs(u['exact']-u['exact_coarse'])<=TOL['numerical_pair_U0_kcal'])
        pairs[case]=dict(status='computed_diagnostic',strata=pr)
    complete=len(arrays)==8
    result=dict(protocol=PROTOCOL,manifest=record(mp),rows=rows,pairs=pairs,attempts=attempts,complete=complete,
        numerical_pass=complete and all(r['numerical_pass'] for r in rows.values()) and all(s['numerical_U0_pass'] for p in pairs.values() for s in p['strata'].values()),
        fitted_field_screen_pass=complete and all(r['strata']['all']['quality']['fitted']['pass_'] for r in [*rows.values(),*pairs.values()]),
        projected_field_screen_pass=complete and all(r['strata']['all']['quality']['projected']['pass_'] for r in [*rows.values(),*pairs.values()]),
        paired_U0_flags_pass=complete and all(p['strata']['all']['paired_projected_U0_flag_pass'] for p in pairs.values()),
        units=dict(field='Hartree/(e bohr)',alpha='bohr^3',U0='kcal/mol'),
        U0_scope='diagonal undamped isolated-source response diagnostic; not environmental energy',**NULLS)
    write_new(out/'result.json',result);return result


def main():
    p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest='op',required=True)
    s=sub.add_parser('prepare')
    for name in ('charges','frameworks','states','plan','output'):s.add_argument('--'+name,required=True)
    for name in ('dry-run','execute','collect'):
        s=sub.add_parser(name);s.add_argument('--manifest',required=True)
        if name!='dry-run':s.add_argument('--output',required=True)
        if name=='execute':s.add_argument('--retry-failed',action='store_true')
    a=p.parse_args()
    if a.op=='prepare':r=prepare(a.charges,a.frameworks,a.states,a.plan,a.output);print(json.dumps({'tasks':len(r['tasks'])}))
    elif a.op=='dry-run':print(json.dumps({'tasks':len(validate(a.manifest)['tasks']),'new_utility_calls':0}))
    elif a.op=='execute':r=execute(a.manifest,a.retry_failed);write_new(a.output,r);print(json.dumps(r))
    else:
        r=collect(a.manifest,a.output);print(json.dumps({k:r[k] for k in ('complete','numerical_pass','fitted_field_screen_pass','projected_field_screen_pass','paired_U0_flags_pass')}))


if __name__=='__main__':main()
