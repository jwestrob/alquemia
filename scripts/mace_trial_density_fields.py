"""Actual potential, electric field and spatial Hessian of saved trial densities."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import shutil
import time
import numpy as np
from affordable_common import BOHR_TO_A,InvalidArtifact,cache_key,read_json,record,verify,write_new
from density_embedding import parse_potential
from mace_hybrid import rotation
from mace_qm_field import offset_points,derivative,diagonal_response,execute_prepared,accepted
from mace_density_multipoles import points as quadrupole_points,hessian,coupling
from mace_trial_density import PROTOCOL as SOURCE_PROTOCOL

PROTOCOL='saved_responsive_trial_density_field_multipole_observations_v1'
TOL=dict(center_potential_au=1e-8,field_refinement_max_au=1e-6,pair_U0_refinement_kcal=.01,
         endpoint_quadrupole_refinement_kcal=.02,pair_quadrupole_refinement_kcal=.01,rigid_kcal=1e-8)
NULLS=dict(numerical_score=None,classification=None,environment_correction=None,baseline_changed=False)


def points(centers):
    return np.concatenate([quadrupole_points(centers),offset_points(centers)])


def key(t,m):
    return cache_key(dict(task={k:v for k,v in t.items() if k!='cache_key'},
        **{k:m[k] for k in ('protocol','source_audit','plan','utility','implementation','tolerances','steps_bohr')}))


def prepare(source_audit,plan,output):
    start=time.monotonic();cpu=time.process_time();a=read_json(source_audit)
    if a['protocol']!=SOURCE_PROTOCOL or not a['complete'] or not a['gates_pass'] or a['plan']!=record(plan):
        raise InvalidArtifact('qualified declared trial density audit required')
    verify(a['implementation']);pm=read_json(verify(a['parent_manifest']))
    dm=read_json(verify(pm['sources']['density_manifest']));br=read_json(verify(pm['sources']['boundary']))
    bm=read_json(verify(br['manifest']));cm=read_json(verify(read_json(verify(a['charges']))['manifest']))
    root=Path(output).resolve();root.mkdir(parents=True,exist_ok=False);impl=root/'implementation';impl.mkdir()
    for name,pin in cm['implementation'].items():shutil.copyfile(verify(pin),impl/name)
    for name,pin in pm['implementation'].items():shutil.copyfile(verify(pin),impl/name)
    for name in ('mace_qm_field.py','mace_trial_density.py','mace_trial_density_fields.py'):
        shutil.copyfile(Path(__file__).with_name(name),impl/name)
    m=dict(protocol=PROTOCOL,source_audit=record(source_audit),plan=record(plan),utility=cm['utilities']['orca_vpot'],
        implementation={p.name:record(p) for p in sorted(impl.glob('*.py'))},tolerances=TOL,
        steps_bohr=dict(field=[.001,.0005],H=[.01,.005]),observations_per_site=49,tasks=[],
        utility_calls=8,new_DFT_calls=0,new_charge_fits=0,new_MACE_calls=0,**NULLS)
    for tid,row in a['rows'].items():
        ct=row['source_task'];old=next(t for t in dm['tasks'] if t['task_id']==tid)
        data=read_json(verify(old['probes']));meta=read_json(verify(row['boundary']))
        bt=next(t for t in bm['tasks'] if t['case_id']==row['case_id'] and t['metal']==row['metal'] and t['mode']=='source')
        native=read_json(verify(br['tasks'][bt['task_id']]['result']))
        lookup={pid:i for i,pid in enumerate(meta['physical_ids'])};indices=[lookup[pid] for pid in data['physical_ids']]
        glob=np.array([atom['global_'] for atom in native['atoms']])[indices]
        expected=set(meta['physical_ids'])-set(meta['source_support_ids'])
        if set(data['physical_ids'])!=expected:raise InvalidArtifact('exterior source inventory differs')
        positions=np.array(meta['positions_A'])[indices]/BOHR_TO_A
        if not np.allclose(positions,data['positions_bohr'],atol=1e-12,rtol=0):raise InvalidArtifact('physical observation positions differ')
        groups=read_json(verify(ct['probe_groups']));weights=read_json(verify(groups['environment_weights']))
        oldpoints=np.loadtxt(verify(ct['points']),skiprows=1);ix=np.array(groups['environment_indices'])
        if weights['physical_ids']!=data['physical_ids'] or not np.allclose(oldpoints[ix],positions,atol=1e-11,rtol=0):
            raise InvalidArtifact('archived responsive center order/geometry differs')
        phi=parse_potential(verify(row['existing_center_potential']),oldpoints)[ix]
        d=root/'inputs'/tid;d.mkdir(parents=True)
        p=d/'points_bohr.xyz';ps=points(positions)
        p.write_text(str(len(ps))+'\n'+''.join(' '.join(format(float(v),'.14f') for v in x)+'\n' for x in ps))
        probes=dict(physical_ids=data['physical_ids'],positions_bohr=positions.tolist(),alpha_bohr3=data['alpha_bohr3'],
            physical_boundary=row['boundary'],source_support_ids=meta['source_support_ids'])
        write_new(d/'probes.json',probes)
        np.savez_compressed(d/'multipole_inputs.npz',center_phi=phi,charge=glob[:,0],dipole=glob[:,1:4]/BOHR_TO_A,
            quadrupole=glob[:,4:].reshape(-1,3,3)/BOHR_TO_A**2)
        t=dict(task_id=tid,case_id=row['case_id'],metal=row['metal'],points=record(p),probes=record(d/'probes.json'),
            multipole_inputs=record(d/'multipole_inputs.npz'),files=ct['files'],source_receipt=ct['source_receipt'],
            source_output=ct['source_output'],xyz=ct['xyz'],state=ct['state'],projection=ct['projection'],
            source_charge_execution=row['source_charge_receipt'],native_moments=br['tasks'][bt['task_id']]['result'],
            old_center_potential=row['existing_center_potential'],old_center_points=ct['points'])
        t['cache_key']=key(t,m);m['tasks'].append(t)
    m['preparation_receipt']=dict(wall_seconds=time.monotonic()-start,process_CPU_seconds=time.process_time()-cpu)
    write_new(root/'manifest.json',m);return validate(root/'manifest.json')


def validate(path):
    m=read_json(path)
    if m['protocol']!=PROTOCOL or m['tolerances']!=TOL or len(m['tasks'])!=8 or m['observations_per_site']!=49:
        raise InvalidArtifact('trial density query configuration changed')
    if m['steps_bohr']!=dict(field=[.001,.0005],H=[.01,.005]):raise InvalidArtifact('trial density steps changed')
    for pin in [m['source_audit'],m['plan'],m['utility'],*m['implementation'].values()]:verify(pin)
    a=read_json(verify(m['source_audit']));paired={}
    for t in m['tasks']:
        if t['cache_key']!=key(t,m):raise InvalidArtifact('trial density scientific cache key changed')
        for k in ('points','probes','multipole_inputs','source_receipt','source_output','xyz','state','projection',
                  'source_charge_execution','native_moments','old_center_potential','old_center_points'):verify(t[k])
        ct=a['rows'][t['task_id']]['source_task']
        if any(t[k]!=ct[k] for k in ('files','source_receipt','source_output','xyz','state','projection')):
            raise InvalidArtifact('trial density source substitution')
        for pin in t['files'].values():verify(pin)
        data=read_json(verify(t['probes']));p=np.loadtxt(verify(t['points']),skiprows=1);expected=points(data['positions_bohr'])
        if p.shape!=expected.shape or not np.allclose(p,expected,atol=2e-13,rtol=0):raise InvalidArtifact('trial density stencil changed')
        paired.setdefault(t['case_id'],[]).append(data)
    if len(paired)!=4 or any(len(v)!=2 or v[0]!=v[1] for v in paired.values()):raise InvalidArtifact('paired physical observation mismatch')
    return m


def collect(path,output):
    m=validate(path);mp=Path(path).resolve();root=mp.parent;out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False)
    rows={};attempts=[];pairs={}
    for t in m['tasks']:
        tid=t['task_id'];good=[]
        for d in sorted((root/'execution'/tid).glob('attempt_*')):
            r=accepted(d,t,mp);attempts.append(dict(task_id=tid,accepted=r is not None,
                receipt=record(d/'execution.json') if (d/'execution.json').exists() else None))
            if r:good.append((d,r))
        if not good:rows[tid]=dict(status='unavailable');continue
        d,r=good[-1];data=read_json(verify(t['probes']));n=len(data['physical_ids'])
        v=parse_potential(verify(r['potential']),np.loadtxt(verify(t['points']),skiprows=1))
        phi,H=hessian(v[:37*n],n);fields=derivative(v[37*n:],n)
        with np.load(verify(t['multipole_inputs'])) as a:inputs={k:a[k].copy() for k in a.files}
        q,mu,Q=[inputs[k] for k in ('charge','dipole','quadrupole')]
        c=[coupling(q,mu,Q,phi,fields[1],h) for h in H]
        mat=rotation();rot=coupling(q,mu@mat.T,np.einsum('ij,njk,lk->nil',mat,Q,mat),phi,fields[1]@mat.T,
            np.einsum('ij,njk,lk->nil',mat,H[1],mat))
        rigid=max(abs(rot[k]-c[1][k]) for k in rot);center=float(np.max(abs(phi-inputs['center_phi'])))
        field_delta=float(np.max(np.linalg.norm(fields[1]-fields[0],axis=1)))
        quad=c[1]['quadrupole']-c[0]['quadrupole'];u=[diagonal_response(f,np.array(data['alpha_bohr3'])) for f in fields]
        np.savez_compressed(out/(tid+'.npz'),phi=phi,H=H[1],H_coarse=H[0],exact=fields[1],exact_coarse=fields[0])
        rows[tid]=dict(status='computed_trial_density_observations',arrays=record(out/(tid+'.npz')),
            center_error_au=center,field_refinement_max_au=field_delta,quadrupole_refinement_kcal=quad,
            rigid_error_kcal=rigid,U0_coarse_kcal=u[0],U0_fine_kcal=u[1],direct_components_kcal=c[1],
            numerical_pass=center<=TOL['center_potential_au'] and field_delta<=TOL['field_refinement_max_au'] and
                           abs(quad)<=TOL['endpoint_quadrupole_refinement_kcal'] and rigid<=TOL['rigid_kcal'],
            execution_receipt=record(d/'execution.json'),potential=r['potential'],utility_wall_seconds=r['utility']['wall_seconds'])
    for case in sorted({t['case_id'] for t in m['tasks']}):
        ca,la=[rows[case+'_'+metal] for metal in ('Ca','La')]
        if any(x['status']=='unavailable' for x in (ca,la)):pairs[case]=dict(status='unavailable');continue
        quad=ca['quadrupole_refinement_kcal']-la['quadrupole_refinement_kcal']
        du=(ca['U0_fine_kcal']-la['U0_fine_kcal'])-(ca['U0_coarse_kcal']-la['U0_coarse_kcal'])
        pairs[case]=dict(status='computed_trial_density_observations',quadrupole_refinement_kcal=quad,U0_refinement_kcal=du,
            numerical_pass=abs(quad)<=TOL['pair_quadrupole_refinement_kcal'] and abs(du)<=TOL['pair_U0_refinement_kcal'])
    complete=all(x['status']!='unavailable' for x in rows.values())
    result=dict(protocol=PROTOCOL,manifest=record(mp),complete=complete,rows=rows,pairs=pairs,attempts=attempts,
        numerical_pass=complete and all(x['numerical_pass'] for x in [*rows.values(),*pairs.values()]),
        U0_scope='Bare undamped diagonal response numerical diagnostic; not environmental energy',**NULLS)
    write_new(out/'result.json',result);return result


def main():
    p=argparse.ArgumentParser(description=__doc__);s=p.add_subparsers(dest='command',required=True)
    a=s.add_parser('prepare')
    for name in ('source-audit','plan','output'):a.add_argument('--'+name,required=True)
    for name in ('dry-run','execute','collect'):
        a=s.add_parser(name);a.add_argument('--manifest',required=True)
        if name!='dry-run':a.add_argument('--output',required=True)
        if name=='execute':a.add_argument('--retry-failed',action='store_true')
    a=p.parse_args()
    if a.command=='prepare':print(json.dumps({'tasks':len(prepare(a.source_audit,a.plan,a.output)['tasks'])}))
    elif a.command=='dry-run':print(json.dumps({'tasks':len(validate(a.manifest)['tasks']),'new_calls':0}))
    elif a.command=='execute':write_new(a.output,execute_prepared(a.manifest,validate(a.manifest),a.retry_failed))
    else:
        r=collect(a.manifest,a.output);print(json.dumps({k:r[k] for k in ('complete','numerical_pass')}))


if __name__=='__main__':main()
