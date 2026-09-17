"""Saved-wavefunction charge and cap-projection checks; no solvent/affinity scalar."""
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
from affordable_common import BOHR_TO_A, InvalidArtifact, cache_key, read_json, record, verify, write_new, xyz
from affordable_solver import real_esp_points, run_command
from affordable_response import source_key
from density_embedding import parse_chelpg, parse_potential

PROTOCOL='normalized_vacuum_CHELPG_source_projection_diagnostic_v1'
QUANTUM_SHA='d63801217ed7b6092cf67f4224b34b28574199aefc1ca291abde34c07fee6275'
RADII={'H':1.2,'C':1.7,'N':1.55,'O':1.52,'S':1.8,'Ca':1.8,'La':1.8}
TOL={'charge_e':5e-5,'ESP_RMS_au':.005,'ESP_relative_RMS':.10,'projection_charge_e':1e-9,'projection_dipole_eA':1e-8}
MODEL={'scheme':'native_ORCA6.1.1_CHELPG','grid_A':.3,'extent_A':2.8,'radii':'COSMO','dipole_constraint':False,
       'probes':'existing_32_Fibonacci_1.4_1.8_shells','probe_radii_A':RADII,
       'projection':'source_identity_sigma_cap_geometric_barycentric_v1','fit_charge_renormalization':False}
NULLS={'solvent_correction':None,'affinity_score':None,'calibrated_class':None,'full_charge_model':None,'baseline_changed':False}


def projection(case,atoms):
    """Map only the QM distribution, preserving its charge and dipole."""
    prep=read_json(verify(case['normalized_global_preparation']))
    old=read_json(verify(case['physical_atoms']))
    source_to_id={a['source_key']:a['id'] for a in old}
    physical={a['id']:a for a in prep['physical_atoms']}
    if len(physical)!=len(prep['physical_atoms']):raise InvalidArtifact('duplicate physical IDs')
    weights=[{'metal':1.}];caps=[]
    maps={a['qm_index']:a for a in case['mapping']}
    if set(maps)!=set(range(1,len(atoms))):raise InvalidArtifact('source mapping atom coverage differs')
    for i in range(1,len(atoms)):
        a=maps[i]
        if a['kind']=='source':
            pid=a['physical_id'];p=physical[pid]
            if p['element']!=atoms[i][0] or not np.allclose(p['xyz_A'],atoms[i][1:],atol=1e-9,rtol=0):
                raise InvalidArtifact('mapped source coordinates differ')
            weights.append({pid:1.})
        elif a['kind']=='sigma_link_H':
            left,right=(source_to_id[source_key(a[k])] for k in ('retained','omitted'))
            x,y=(np.array(physical[k]['xyz_A']) for k in (left,right))
            lam=a['length_A']/float(np.linalg.norm(y-x))
            if not 0<lam<1 or atoms[i][0]!='H' or not np.allclose((1-lam)*x+lam*y,atoms[i][1:],atol=1e-9,rtol=0):
                raise InvalidArtifact('cap anchors/geometry invalid')
            weights.append({left:1-lam,right:lam})
            caps.append({'qm_index':i,'retained':left,'omitted':right,'lambda':lam})
        else:raise InvalidArtifact('unsupported charge mapping')
    ids=sorted({k for w in weights for k in w});lookup={k:i for i,k in enumerate(ids)}
    matrix=np.zeros((len(ids),len(atoms)))
    for j,w in enumerate(weights):
        for k,v in w.items():matrix[lookup[k],j]+=v
    coords=np.array([physical[k]['xyz_A'] for k in ids]);qm=np.array([a[1:] for a in atoms])
    if np.max(np.abs(matrix.sum(0)-1))>1e-12 or np.max(np.abs(coords.T@matrix-qm.T))>1e-9:
        raise InvalidArtifact('projection does not conserve monopole/dipole geometry')
    return {'physical_ids':ids,'coordinates_A':coords.tolist(),'weights':matrix.tolist(),'caps':caps,
            'scope':'projected_QM_distribution_only; no_classical_environment_charges'}


def prepare(quantum,plan,output):
    from mace_omol_matched_h import collect_quantum
    if record(quantum)['sha256']!=QUANTUM_SHA:raise InvalidArtifact('declared quantum source differs')
    q=collect_quantum(quantum)
    if q['status']!='complete':raise InvalidArtifact('eight actual quantum endpoints required')
    parent=read_json(quantum);p=read_json(verify(parent['preparation']))
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False);impl=out/'implementation';impl.mkdir()
    # Preserve the executed parent implementations and add only this utility interface.
    names={k:verify(v) for k,v in parent['implementation'].items()}
    for name in ('mace_omol_charges.py','density_embedding.py','global_electrostatic.py'):
        names[name]=Path(__file__).with_name(name)
    pins={}
    for name,source in names.items():shutil.copyfile(source,impl/name);pins[name]=record(impl/name)
    od=verify(parent['orca']).parent
    utilities={k:record(od/k) for k in ('orca_chelpg','orca_vpot')};tasks=[]
    for task in parent['tasks']:
        tid=task['task_id'];case=read_json(verify(p['cases'][task['case_id']]))
        atoms=xyz(verify(task['xyz']));d=out/'inputs'/tid;d.mkdir(parents=True)
        pp=d/'points_bohr.xyz';points=real_esp_points(atoms,RADII)
        pp.write_text(str(len(points))+'\n'+''.join(' '.join(format(v,'.12f') for v in x)+'\n' for x in points))
        proj=d/'projection.json';write_new(proj,projection(case,atoms));files={};sources={}
        source_dir=Path(task['output_path']).parent
        for key,suffix in [('gbw','gbw'),('density','densities'),('density_info','densitiesinfo')]:
            source=source_dir/('endpoint.runtime.'+suffix);destination=d/source.name
            sources[key]=record(source);shutil.copyfile(source,destination);files[key]=record(destination)
        t={'task_id':tid,'case_id':task['case_id'],'metal':task['metal'],'charge':task['charge'],
           'source_receipt':q['rows'][tid]['receipt'],'source_output':q['rows'][tid]['output'],
           'source_mapping':p['cases'][task['case_id']],'xyz':task['xyz'],'points':record(pp),
           'projection':record(proj),'files':files,'source_wavefunctions':sources}
        t['cache_key']=cache_key({'task':t,'model':MODEL,'tolerances':TOL,'utilities':utilities,'implementation':pins});tasks.append(t)
    m={'protocol_id':PROTOCOL,'quantum':record(quantum),'preparation':parent['preparation'],'plan':record(plan),'tasks':tasks,
       'model':MODEL,'tolerances':TOL,'utilities':utilities,'implementation':pins,'python':record(os.path.realpath(os.sys.executable)),
       'new_DFT_calls':0,'new_MACE_calls':0,'charge_utility_calls':8,'potential_utility_calls':8,
       'compute_budget':None,'wall_time_limit':None,**NULLS}
    write_new(out/'manifest.json',m);return validate(out/'manifest.json')


def validate(manifest):
    if read_json(manifest).get('protocol_id')=='embedded_CHELPG_source_projection_diagnostic_v1':
        from mace_responsive_charges import validate
        return validate(manifest)
    m=read_json(manifest)
    if m['protocol_id']!=PROTOCOL or m['model']!=MODEL or m['tolerances']!=TOL or len(m['tasks'])!=8 or m['quantum']['sha256']!=QUANTUM_SHA:
        raise InvalidArtifact('charge protocol/source/settings differ')
    for pin in [m['quantum'],m['preparation'],m['plan'],m['python'],*m['utilities'].values(),*m['implementation'].values()]:verify(pin)
    p=read_json(verify(m['preparation']));parent=read_json(verify(m['quantum']));expected={t['task_id']:t for t in parent['tasks']};pairs={}
    if set(expected)!={t['task_id'] for t in m['tasks']}:raise InvalidArtifact('endpoint inventory differs')
    for t in m['tasks']:
        for pin in [t['source_receipt'],t['source_output'],t['source_mapping'],t['xyz'],t['points'],t['projection'],*t['files'].values(),*t['source_wavefunctions'].values()]:verify(pin)
        old=expected[t['task_id']];receipt=read_json(verify(t['source_receipt']))
        if (t['source_mapping']!=p['cases'][t['case_id']] or t['xyz']!=old['xyz'] or t['charge']!=old['charge']
            or t['case_id']!=old['case_id'] or t['metal']!=old['metal'] or receipt['task_id']!=t['task_id']
            or receipt['manifest']!=m['quantum'] or receipt['artifacts']['output']!=t['source_output']
            or receipt['returncode']!=0 or not receipt['normal_termination']):raise InvalidArtifact('source quantum receipt/state differs')
        if any(t['files'][k]['sha256']!=pin['sha256'] for k,pin in t['source_wavefunctions'].items()):raise InvalidArtifact('copied density changed')
        atoms=xyz(verify(t['xyz']));points=np.loadtxt(verify(t['points']),skiprows=1)
        if not np.allclose(points,real_esp_points(atoms,RADII),atol=1e-11,rtol=0):raise InvalidArtifact('ESP probes changed')
        saved=read_json(verify(t['projection']));actual=projection(read_json(verify(t['source_mapping'])),atoms)
        if saved['physical_ids']!=actual['physical_ids'] or saved['scope']!=actual['scope'] or not np.allclose(saved['weights'],actual['weights'],atol=1e-12,rtol=0) or saved['coordinates_A']!=actual['coordinates_A']:
            raise InvalidArtifact('source projection changed')
        payload={k:v for k,v in t.items() if k!='cache_key'}
        if t['cache_key']!=cache_key({'task':payload,'model':MODEL,'tolerances':TOL,'utilities':m['utilities'],'implementation':m['implementation']}):raise InvalidArtifact('charge cache differs')
        pairs.setdefault(t['case_id'],[]).append(points)
    for a,b in pairs.values():
        if not np.array_equal(a,b):raise InvalidArtifact('paired ESP probes differ')
    return {'status':'pass','tasks':8,'utility_calls':16,'manifest':record(manifest)}


def accepted(attempt,task,manifest):
    path=attempt/'execution.json'
    if not path.exists():return None
    r=read_json(path)
    if r['task']!=task or r['manifest']!=record(manifest) or r['status']!='complete' or not r['slurm_job_id']:return None
    for step in ('chelpg','vpot'):
        u=r[step]
        if u['returncode']!=0 or u['slurm_job_id']!=r['slurm_job_id']:return None
        verify(u['log']);verify(u['resource_usage'])
    verify(r['potential']);return r


def execute(manifest,retry_failed=False):
    validate(manifest);m=read_json(manifest);mp=Path(manifest).resolve();root=mp.parent/'execution';root.mkdir(exist_ok=True)
    if not os.environ.get('SLURM_JOB_ID'):raise InvalidArtifact('utilities require an allocation')
    def one(t):
        td=root/t['task_id'];td.mkdir(exist_ok=True);attempts=sorted(td.glob('attempt_*'))
        for a in reversed(attempts):
            if accepted(a,t,mp):return {'task_id':t['task_id'],'status':'reused','receipt':record(a/'execution.json')}
        if attempts and not retry_failed:return {'task_id':t['task_id'],'status':'failed_attempt_requires_explicit_retry'}
        d=td/f'attempt_{len(attempts)+1:04d}';d.mkdir();start=time.monotonic();r={'task':t,'manifest':record(mp),'slurm_job_id':os.environ['SLURM_JOB_ID'],'host':socket.gethostname(),'status':'failed'}
        try:
            for pin in t['files'].values():shutil.copyfile(verify(pin),d/Path(pin['path']).name)
            gbw=d/'endpoint.runtime.gbw'
            r['chelpg']=run_command([str(verify(m['utilities']['orca_chelpg'])),str(gbw)],d,d/'chelpg.log',d/'chelpg.resources.txt')
            if r['chelpg']['returncode']:raise InvalidArtifact('CHELPG utility failed')
            r['vpot']=run_command([str(verify(m['utilities']['orca_vpot'])),str(gbw),'endpoint.runtime.scfp',str(verify(t['points'])),str(d/'potential.out')],d,d/'vpot.log',d/'vpot.resources.txt')
            if r['vpot']['returncode']:raise InvalidArtifact('potential utility failed')
            for pin in t['files'].values():
                if record(d/Path(pin['path']).name)['sha256']!=pin['sha256']:raise InvalidArtifact('utility mutated saved wavefunction')
            r.update(status='complete',potential=record(d/'potential.out'))
        except Exception as exc:r['reason']=str(exc)
        r['wall_seconds']=time.monotonic()-start;write_new(d/'execution.json',r)
        return {'task_id':t['task_id'],'status':r['status'],'receipt':record(d/'execution.json')}
    with (mp.parent/'execute.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        workers=min(8,int(os.environ['SLURM_CPUS_ON_NODE']))
        with ThreadPoolExecutor(max_workers=workers) as pool:rows=list(pool.map(one,m['tasks']))
    return {'status':'complete' if all(r['status'] in ('complete','reused') for r in rows) else 'incomplete','rows':rows,'manifest':record(mp)}


def potential(charges,coords_A,points_bohr):
    distance=np.linalg.norm(points_bohr[:,None,:]-np.array(coords_A)[None,:,:]/BOHR_TO_A,axis=2)
    if np.min(distance)<1e-8:raise InvalidArtifact('potential probe/nucleus overlap')
    return np.sum(np.array(charges)[None,:]/distance,axis=1)


def quality(predicted,exact):
    rms=float(np.sqrt(np.mean((predicted-exact)**2)));norm=float(np.sqrt(np.mean(exact**2)))
    return {'RMS_au':rms,'relative_RMS':rms/norm if norm else None,
            'pass':rms<=TOL['ESP_RMS_au'] or norm>0 and rms/norm<=TOL['ESP_relative_RMS']}


def report(manifest,output):
    validate(manifest);mp=Path(manifest).resolve();m=read_json(mp);rows={};arrays={};attempts=[]
    for t in m['tasks']:
        tid=t['task_id'];row={'status':'unavailable'};good=[]
        for a in sorted((mp.parent/'execution'/tid).glob('attempt_*')):
            r=accepted(a,t,mp);attempts.append({'task_id':tid,'directory':str(a),'accepted':r is not None,
                'receipt':record(a/'execution.json') if (a/'execution.json').exists() else None})
            if r:good.append((a,r))
        try:
            if not good:raise InvalidArtifact('no complete actual utility receipt')
            a,r=good[-1];atoms=xyz(verify(t['xyz']));points=np.loadtxt(verify(t['points']),skiprows=1)
            q=parse_chelpg(verify(r['chelpg']['log']).read_text(),[v[0] for v in atoms],t['charge'])
            exact=parse_potential(verify(r['potential']),points);fit=potential(q,[a[1:] for a in atoms],points)
            proj=read_json(verify(t['projection']));matrix=np.array(proj['weights']);pq=matrix@q
            mapped=potential(pq,proj['coordinates_A'],points)
            cq=abs(float(pq.sum()-q.sum()));dq=float(np.max(np.abs(np.array(proj['coordinates_A']).T@pq-np.array([v[1:] for v in atoms]).T@q)))
            row={'status':'computed','charge_e':q.tolist(),'charge_sum_e':float(q.sum()),'projected_charge_e':pq.tolist(),
                 'fit_quality':quality(fit,exact),'projected_quality':quality(mapped,exact),'projection_only_error':quality(mapped,fit),
                 'projection_charge_error_e':cq,'projection_dipole_error_eA':dq,
                 'projection_conservation_pass':cq<=TOL['projection_charge_e'] and dq<=TOL['projection_dipole_eA'],
                 'execution_receipt':record(a/'execution.json'),'potential':r['potential'],'projection':t['projection'],'points':t['points'],
                 'utility_wall_seconds':r['chelpg']['wall_seconds']+r['vpot']['wall_seconds']}
            arrays[tid]={'exact':exact.tolist(),'fitted':fit.tolist(),'projected':mapped.tolist()}
        except (ValueError,OSError,KeyError) as exc:row['reason']=str(exc)
        rows[tid]=row
    pairs={}
    for name in sorted({t['case_id'] for t in m['tasks']}):
        ca,la=(arrays.get(name+'_'+metal) for metal in ('Ca','La'))
        if ca is not None and la is not None:
            delta={k:np.array(ca[k])-la[k] for k in ca}
            pairs[name]={'fit_quality':quality(delta['fitted'],delta['exact']),'projected_quality':quality(delta['projected'],delta['exact']),
                         'projection_only_error':quality(delta['projected'],delta['fitted'])}
    complete=len(arrays)==8
    result={'protocol_id':PROTOCOL,'manifest':record(mp),'rows':rows,'paired_Ca_minus_La':pairs,'attempts':attempts,
            'status':'complete' if complete else 'incomplete','tolerances':TOL,
            'fit_quality_gate_pass':complete and all(r['fit_quality']['pass'] for r in rows.values()) and all(r['fit_quality']['pass'] for r in pairs.values()),
            'projection_gate_pass':complete and all(r['projected_quality']['pass'] and r['projection_conservation_pass'] for r in rows.values()) and all(r['projected_quality']['pass'] for r in pairs.values()),
            'implementation':record(__file__),**NULLS}
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False);write_new(out/'potentials.json',arrays);result['potential_arrays']=record(out/'potentials.json');write_new(out/'result.json',result)
    return {k:result[k] for k in ('status','fit_quality_gate_pass','projection_gate_pass','paired_Ca_minus_La')}


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);s=p.add_subparsers(dest='op',required=True)
    a=s.add_parser('prepare')
    for k in ('quantum','plan','output'):a.add_argument('--'+k,required=True)
    for op in ('dry-run','execute','report'):
        a=s.add_parser(op);a.add_argument('--manifest',required=True);a.add_argument('--output',required=op!='dry-run')
        if op=='execute':a.add_argument('--retry-failed',action='store_true')
    a=p.parse_args()
    if a.op=='prepare':r=prepare(a.quantum,a.plan,a.output)
    elif a.op=='dry-run':r=validate(a.manifest)
    elif a.op=='execute':r=execute(a.manifest,a.retry_failed);write_new(a.output,r)
    else:r=report(a.manifest,a.output)
    print(json.dumps(r,indent=2))
