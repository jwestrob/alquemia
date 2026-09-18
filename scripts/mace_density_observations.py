"""Fresh endpoint CHELPG and exact density observations for a declared panel."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import shutil
import time
import numpy as np
from affordable_common import BOHR_TO_A,HA_TO_KCAL,InvalidArtifact,cache_key,read_json,record,verify,write_new,xyz
from affordable_solver import real_esp_points
from density_embedding import parse_chelpg,parse_potential
from mace_omol_charges import MODEL,RADII,TOL as CHARGE_TOL,accepted,execute_prepared,potential,quality
from mace_trial_density_fields import points,TOL as FIELD_TOL
from mace_density_multipoles import hessian,coupling
from mace_qm_field import derivative,diagonal_response
from mace_hybrid import rotation
from mace_density_boundary_panel import validate as validate_boundary

PROTOCOL='source_graph_responsive_density_CHELPG_field_observations_v1'


def key(t,m):
    return cache_key(dict(task={k:v for k,v in t.items() if k!='cache_key'},
        **{k:m[k] for k in ('protocol','quantum','environment_boundary','query_scope','model','charge_tolerances','field_tolerances','utilities','implementation')}))


def prepare(quantum,environment_boundary,query_scope,output):
    start=time.monotonic();cpu=time.process_time();qr=read_json(quantum)
    if qr['status']!='complete':raise InvalidArtifact('complete actual quantum endpoints required')
    for n in ('collection_implementation','parser_implementation'):verify(qr[n])
    qm=read_json(verify(qr['manifest']));p=read_json(verify(qm['preparation']))
    br=read_json(environment_boundary);bm=validate_boundary(verify(br['manifest']))
    if not br['complete'] or not br['gates_pass'] or br['mode']!='zero_source' or bm['preparation']!=qm['preparation']:
        raise InvalidArtifact('same physical environment-only boundary required')
    root=Path(output).resolve();root.mkdir(parents=True,exist_ok=False);impl=root/'implementation';impl.mkdir()
    paths={k:verify(v) for parent in (qm,bm) for k,v in parent['implementation'].items()}
    for n in ('mace_density_observations.py','mace_omol_charges.py','mace_trial_density_fields.py','mace_density_inputs.py',
              'mace_responsive_charges.py','mace_omol_vacuum.py','mace_density_boundary_panel.py','mace_density_frameworks.py',
              'mace_trial_density.py','mace_density_gk_hybrid.py','mace_native_field_input.py'):
        paths[n]=Path(__file__).with_name(n)
    for n,src in paths.items():shutil.copyfile(src,impl/n)
    od=verify(qm['orca']).parent
    m=dict(protocol=PROTOCOL,quantum=record(quantum),environment_boundary=record(environment_boundary),
        query_scope=record(query_scope),plan=p['agreement'],model=MODEL,charge_tolerances=CHARGE_TOL,field_tolerances=FIELD_TOL,
        utilities={n:record(od/n) for n in ('orca_chelpg','orca_vpot')},implementation={p.name:record(p) for p in impl.glob('*.py')},
        tasks=[],new_DFT_calls=0,new_MACE_calls=0,baseline_changed=False)
    for qt in qm['tasks']:
        tid=qt['task_id'];name,metal=qt['case_id'],qt['metal'];case=read_json(verify(p['cases'][name]))
        meta=read_json(verify(br['cases'][name]));state=read_json(verify(meta['state']))
        bt=next(t for t in bm['tasks'] if t['case_id']==name and t['metal']==metal)
        nr=read_json(verify(br['tasks'][bt['task_id']]['result']));ids=meta['physical_ids'];support=set(meta['source_support_ids'])
        selected=[i for i,pid in enumerate(ids) if pid not in support];names=[ids[i] for i in selected]
        positions=np.array(meta['positions_A'])[selected]/BOHR_TO_A
        glob=np.array([a['global_'] for a in nr['atoms']])[selected]
        alpha=np.array([a['polarizability_A3'] for a in nr['parameters']['atoms']])[selected]/BOHR_TO_A**3
        if not len(selected) or np.min(alpha)<=0:raise InvalidArtifact('nonpositive environment response parameter')
        oldq=dict(zip([a['id'] for a in state['physical_atoms']],state['environment_charges_e']))
        if any(q and pid not in names for pid,q in oldq.items()):raise InvalidArtifact('old generating charge missing from observations')
        atoms=xyz(verify(qt['xyz']));esp=real_esp_points(atoms,RADII);allpoints=np.concatenate([points(positions),esp])
        d=root/'inputs'/tid;d.mkdir(parents=True)
        pp=d/'points_bohr.xyz';pp.write_text(str(len(allpoints))+'\n'+''.join(' '.join(format(float(v),'.14f') for v in row)+'\n' for row in allpoints))
        probes=dict(physical_ids=names,positions_bohr=positions.tolist(),alpha_bohr3=alpha.tolist(),
            old_generating_charges_e=[oldq[pid] for pid in names],exterior_start=49*len(names),exterior_count=len(esp),
            charge=glob[:,0].tolist(),dipole_au=(glob[:,1:4]/BOHR_TO_A).tolist(),quadrupole_au=(glob[:,4:].reshape(-1,3,3)/BOHR_TO_A**2).tolist())
        write_new(d/'probes.json',probes);files={}
        for k,suffix in [('gbw','gbw'),('density','densities'),('density_info','densitiesinfo')]:
            src=Path(qt['output_path']).parent/('endpoint.runtime.'+suffix);dst=d/src.name
            verify(record(src));shutil.copyfile(src,dst);files[k]=record(dst)
        t=dict(task_id=tid,case_id=name,metal=metal,charge=qt['charge'],xyz=qt['xyz'],projection=case['projection'],
            points=record(pp),probes=record(d/'probes.json'),files=files,source_task=qt,
            source_receipt=qr['rows'][tid]['receipt'],source_output=qr['rows'][tid]['output'],
            quantum_energy_hartree=qr['rows'][tid]['energy_hartree'],native_environment=br['tasks'][bt['task_id']]['result'])
        t['cache_key']=key(t,m);m['tasks'].append(t)
    m['preparation_receipt']=dict(wall_seconds=time.monotonic()-start,CPU_seconds=time.process_time()-cpu)
    m['charge_utility_calls']=len(m['tasks']);m['potential_utility_calls']=len(m['tasks'])
    write_new(root/'manifest.json',m);return validate(root/'manifest.json')


def validate(path):
    m=read_json(path)
    if m['protocol']!=PROTOCOL or m['model']!=MODEL or m['charge_tolerances']!=CHARGE_TOL or m['field_tolerances']!=FIELD_TOL:
        raise InvalidArtifact('density observation model/tolerances differ')
    for pin in [m['quantum'],m['environment_boundary'],m['query_scope'],m['plan'],*m['utilities'].values(),*m['implementation'].values()]:verify(pin)
    qr=read_json(verify(m['quantum']));qm=read_json(verify(qr['manifest']));paired={}
    if len(m['tasks'])!=len(qm['tasks']) or {t['task_id'] for t in m['tasks']}!={t['task_id'] for t in qm['tasks']}:
        raise InvalidArtifact('density task inventory differs')
    for t in m['tasks']:
        if t['cache_key']!=key(t,m):raise InvalidArtifact('density observation cache differs')
        for k in ('xyz','projection','points','probes','source_receipt','source_output','native_environment'):verify(t[k])
        for pin in t['files'].values():verify(pin)
        qt=next(x for x in qm['tasks'] if x['task_id']==t['task_id'])
        if t['source_task']!=qt or t['quantum_energy_hartree']!=qr['rows'][t['task_id']]['energy_hartree']:
            raise InvalidArtifact('actual quantum source differs')
        data=read_json(verify(t['probes']));saved=np.loadtxt(verify(t['points']),skiprows=1)
        wanted=np.concatenate([points(data['positions_bohr']),real_esp_points(xyz(verify(t['xyz'])),RADII)])
        if saved.shape!=wanted.shape or not np.allclose(saved,wanted,atol=2e-13,rtol=0):raise InvalidArtifact('density observation stencil changed')
        paired.setdefault(t['case_id'],[]).append(data)
    if any(len(v)!=2 or v[0]!=v[1] for v in paired.values()):raise InvalidArtifact('paired environment observations differ')
    return m


def collect(path,output):
    m=validate(path);mp=Path(path).resolve();out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False)
    rows={};pairs={};attempts=[];esp_arrays={}
    for t in m['tasks']:
        tid=t['task_id'];good=[]
        for d in sorted((mp.parent/'execution'/tid).glob('attempt_*')):
            r=accepted(d,t,mp);attempts.append(dict(task_id=tid,directory=str(d),accepted=r is not None))
            if r:good.append((d,r))
        if not good:rows[tid]=dict(status='unavailable');continue
        d,r=good[-1];data=read_json(verify(t['probes']));n=len(data['physical_ids']);atoms=xyz(verify(t['xyz']))
        coords=np.loadtxt(verify(t['points']),skiprows=1);v=parse_potential(verify(r['potential']),coords)
        q=parse_chelpg(verify(r['chelpg']['log']).read_text(),[a[0] for a in atoms],t['charge'])
        proj=read_json(verify(t['projection']));pq=np.array(proj['weights'])@q
        cq=abs(float(pq.sum()-q.sum()));dq=float(np.max(abs(np.array(proj['coordinates_A']).T@pq-np.array([a[1:] for a in atoms]).T@q)))
        tail=coords[49*n:];exact=v[49*n:];fit=potential(q,[a[1:] for a in atoms],tail);mapped=potential(pq,proj['coordinates_A'],tail)
        esp_arrays[tid]=dict(exact=exact,fitted=fit,projected=mapped)
        phi,H=hessian(v[:37*n],n);fields=derivative(v[37*n:49*n],n)
        eq,mu,Q=[np.array(data[k]) for k in ('charge','dipole_au','quadrupole_au')]
        coupling_values=[coupling(eq,mu,Q,phi,fields[1],h) for h in H];mat=rotation()
        rotated=coupling(eq,mu@mat.T,np.einsum('ij,njk,lk->nil',mat,Q,mat),phi,fields[1]@mat.T,np.einsum('ij,njk,lk->nil',mat,H[1],mat))
        rigid=max(abs(rotated[k]-coupling_values[1][k]) for k in rotated)
        field_delta=float(np.max(np.linalg.norm(fields[1]-fields[0],axis=1)))
        quad=coupling_values[1]['quadrupole']-coupling_values[0]['quadrupole']
        u=[diagonal_response(f,np.array(data['alpha_bohr3'])) for f in fields]
        old=float(np.dot(phi,data['old_generating_charges_e'])*HA_TO_KCAL)
        np.savez_compressed(out/(tid+'.npz'),phi=phi,H=H[1],H_coarse=H[0],exact=fields[1],exact_coarse=fields[0])
        rows[tid]=dict(status='computed',charge_e=q.tolist(),charge_sum_e=float(q.sum()),projected_charge_e=pq.tolist(),
            fit_quality=quality(fit,exact),projected_quality=quality(mapped,exact),projection_charge_error_e=cq,projection_dipole_error_eA=dq,
            projection_conservation_pass=cq<=CHARGE_TOL['projection_charge_e'] and dq<=CHARGE_TOL['projection_dipole_eA'],
            arrays=record(out/(tid+'.npz')),projection=t['projection'],probes=t['probes'],
            field_refinement_max_au=field_delta,quadrupole_refinement_kcal=quad,rigid_error_kcal=rigid,
            U0_coarse_kcal=u[0],U0_fine_kcal=u[1],direct_components_kcal=coupling_values[1],
            old_generating_interaction_kcal=old,intrinsic_trial_core_kcal=t['quantum_energy_hartree']*HA_TO_KCAL-old,
            center_archive_error_au=None,center_archive_status='unavailable_on_new_geometry',
            numerical_pass=field_delta<=FIELD_TOL['field_refinement_max_au'] and abs(quad)<=FIELD_TOL['endpoint_quadrupole_refinement_kcal'] and rigid<=FIELD_TOL['rigid_kcal'],
            execution_receipt=record(d/'execution.json'),potential=r['potential'])
    for name in sorted({t['case_id'] for t in m['tasks']}):
        a,b=[rows[name+'_'+e] for e in ('Ca','La')]
        if any(r['status']!='computed' for r in (a,b)):pairs[name]=dict(status='unavailable');continue
        x,y=[esp_arrays[name+'_'+e] for e in ('Ca','La')];delta={k:x[k]-y[k] for k in x}
        quad=a['quadrupole_refinement_kcal']-b['quadrupole_refinement_kcal']
        du=(a['U0_fine_kcal']-b['U0_fine_kcal'])-(a['U0_coarse_kcal']-b['U0_coarse_kcal'])
        pairs[name]=dict(status='computed',fit_quality=quality(delta['fitted'],delta['exact']),projected_quality=quality(delta['projected'],delta['exact']),
            quadrupole_refinement_kcal=quad,U0_refinement_kcal=du,
            numerical_pass=abs(quad)<=FIELD_TOL['pair_quadrupole_refinement_kcal'] and abs(du)<=FIELD_TOL['pair_U0_refinement_kcal'])
    complete=all(r['status']=='computed' for r in rows.values())
    r=dict(protocol=PROTOCOL,manifest=record(path),status='complete' if complete else 'incomplete',rows=rows,pairs=pairs,attempts=attempts,
        fit_quality_gate_pass=complete and all(r['fit_quality']['pass'] for r in [*rows.values(),*pairs.values()]),
        projection_gate_pass=complete and all(r['projected_quality']['pass'] for r in [*rows.values(),*pairs.values()]) and all(r['projection_conservation_pass'] for r in rows.values()),
        numerical_pass=complete and all(r['numerical_pass'] for r in [*rows.values(),*pairs.values()]),
        U0_scope='bare diagonal numerical diagnostic, not environmental energy',center_archive_check=None,
        baseline_changed=False,environment_correction=None,full_hybrid_score=None,calibrated_class=None)
    write_new(out/'result.json',r);return r


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest='command',required=True)
    s=sub.add_parser('prepare')
    for n in ('quantum','environment-boundary','query-scope','output'):s.add_argument('--'+n,required=True)
    for op in ('dry-run','execute','collect'):
        s=sub.add_parser(op);s.add_argument('--manifest',required=True)
        if op!='dry-run':s.add_argument('--output',required=True)
        if op=='execute':s.add_argument('--retry-failed',action='store_true')
    a=p.parse_args()
    if a.command=='prepare':r=prepare(a.quantum,a.environment_boundary,a.query_scope,a.output);print(json.dumps(dict(tasks=len(r['tasks']))))
    elif a.command=='dry-run':print(json.dumps(dict(tasks=len(validate(a.manifest)['tasks']),status='pass')))
    elif a.command=='execute':r=execute_prepared(a.manifest,validate(a.manifest),a.retry_failed);write_new(a.output,r);print(json.dumps(dict(status=r['status'])))
    else:r=collect(a.manifest,a.output);print(json.dumps({k:r[k] for k in ('status','fit_quality_gate_pass','projection_gate_pass','numerical_pass')}))
