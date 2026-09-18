"""Recover declared real ddPCM states without reusing a poisoned native Model."""
from __future__ import annotations
import argparse
import gc
import json
import math
import os
import platform
from pathlib import Path
import resource
import shutil
import time
import numpy as np
from affordable_common import BOHR_TO_A, HA_TO_KCAL, InvalidArtifact, cache_key, read_json, record, verify, write_new
from mace_ddx_source_self import DIRECT_PROTOCOL, validate, direct_potential, potential_error
from mace_ddx_convergence import preflight as logged_preflight

PROTOCOL=DIRECT_PROTOCOL+'_independent_iteration_recovery_v3'


def prepare(config,output):
    c=read_json(config);verify(c['plan']);old=read_json(verify(c['parent_collection']));parent=validate(verify(old['manifest']))
    replay=read_json(verify(c['logged_replay']));logged=logged_preflight(verify(replay['manifest']))
    if parent['protocol']!=DIRECT_PROTOCOL or logged['parent']!=old['manifest'] or replay['status']!='computed':
        raise InvalidArtifact('matching completed original collection and logged replay required')
    root=Path(output).resolve();root.mkdir(parents=True,exist_ok=False);impl=root/'implementation';impl.mkdir()
    for name in ('affordable_common.py','mace_ddx_source_self.py','mace_ddx_convergence.py','mace_ddx_recovery.py','mace_ddx_source_report.py'):
        shutil.copyfile(Path(__file__).with_name(name),impl/name)
    m=dict(protocol=PROTOCOL,config=record(config),plan=c['plan'],parent_collection=c['parent_collection'],logged_replay=c['logged_replay'],
           parent=old['manifest'],implementation={p.name:record(p) for p in impl.glob('*.py')},
           model={**parent['model'],'maxiter':1200},grids=parent['grids'],tolerances=parent['tolerances'],python=parent['python'],module=parent['module'],groups=[],
           new_DFT_calls=0,new_MACE_calls=0,new_score=None,baseline_changed=False,requested_state_roles=20)
    for g in parent['groups']:
        prior=old['groups'][g['group_id']];reuse={}
        if prior['cache_key']!=g['cache_key'] or prior['source']!=g['source']:raise InvalidArtifact('original group identity differs')
        for label in g['states']:
            if prior['rows'].get(label,{}).get('status')=='computed':
                reuse[label]=dict(kind='original',receipt=c['parent_collection'],group=g['group_id'],state=label,arrays=prior['arrays'],maxiter=300)
            elif g['source']==logged['source'] and parent['grids'][g['grid']]==logged['grid'] and label==logged['state']:
                reuse[label]=dict(kind='logged',receipt=c['logged_replay'],state=label,arrays=replay['arrays'],maxiter=1200)
        group={k:v for k,v in g.items() if k!='cache_key'};group.update(reuse=reuse,new_states=[s for s in g['states'] if s not in reuse])
        group['cache_key']=cache_key(dict(group=group,model=m['model'],config=m['config'],implementation=m['implementation']))
        m['groups'].append(group)
    m['requested_new_forward_solves']=sum(len(g['new_states']) for g in m['groups']);m['reused_state_roles']=20-m['requested_new_forward_solves']
    m['cache_key']=cache_key(m);write_new(root/'manifest.json',m);return preflight(root/'manifest.json')


def preflight(path):
    m=read_json(path);parent=validate(verify(m['parent']));old=read_json(verify(m['parent_collection']));replay=read_json(verify(m['logged_replay']))
    lm=logged_preflight(verify(replay['manifest']))
    if m['protocol']!=PROTOCOL or m['cache_key']!=cache_key({k:v for k,v in m.items() if k!='cache_key'}):raise InvalidArtifact('recovery manifest changed')
    for pin in (m['config'],m['plan'],m['python'],m['module'],*m['implementation'].values()):verify(pin)
    if old['manifest']!=m['parent'] or lm['parent']!=m['parent'] or m['model']!={**parent['model'],'maxiter':1200} or m['grids']!=parent['grids'] or m['tolerances']!=parent['tolerances'] or m['python']!=parent['python'] or m['module']!=parent['module']:
        raise InvalidArtifact('physical model or software changed in numerical recovery')
    if len(m['groups'])!=len(parent['groups']):raise InvalidArtifact('group inventory changed')
    for g,p in zip(m['groups'],parent['groups']):
        if any(g[k]!=v for k,v in p.items() if k!='cache_key'):raise InvalidArtifact('original scientific inventory changed')
        if g['cache_key']!=cache_key(dict(group={k:v for k,v in g.items() if k!='cache_key'},model=m['model'],config=m['config'],implementation=m['implementation'])):raise InvalidArtifact('group cache changed')
        expected={}
        for label in g['states']:
            previous=old['groups'][g['group_id']]
            if previous['rows'].get(label,{}).get('status')=='computed':
                expected[label]=dict(kind='original',receipt=m['parent_collection'],group=g['group_id'],state=label,arrays=previous['arrays'],maxiter=300)
            elif g['source']==lm['source'] and parent['grids'][g['grid']]==lm['grid'] and label==lm['state']:
                expected[label]=dict(kind='logged',receipt=m['logged_replay'],state=label,arrays=replay['arrays'],maxiter=1200)
        if g['reuse']!=expected or g['new_states']!=[s for s in g['states'] if s not in expected]:raise InvalidArtifact('reuse inventory changed')
        for item in g['reuse'].values():verify(item['arrays'])
    if m['requested_new_forward_solves']!=sum(len(g['new_states']) for g in m['groups']) or m['reused_state_roles']+m['requested_new_forward_solves']!=20 or m['new_DFT_calls'] or m['new_MACE_calls']:
        raise InvalidArtifact('execution inventory changed')
    return m


def retained_state(item):
    source=read_json(verify(item['receipt']))
    if item['kind']=='original':
        group=source['groups'][item['group']];original=group['rows'][item['state']];prefix=item['state']+'_';parameters=original.get('actual_parameters',group['actual_parameters'])
    else:
        original=source;prefix='';manifest=logged_preflight(verify(source['manifest']));parameters={**manifest['model'],**manifest['grid']}
    if original['status']!='computed':raise InvalidArtifact('cannot reuse unsuccessful state')
    with np.load(verify(item['arrays']),allow_pickle=False) as a:
        arrays={k:np.array(a[prefix+k]) for k in ('phi','psi','x')};arrays['cavity_bohr']=np.array(a['cavity_bohr'])
    row={k:v for k,v in original.items() if k not in ('arrays','manifest','solver_log')}
    row.update(state=item['state'],execution='reused',reused_from=item,solver_maxiter=item['maxiter'],new_forward_solve_started=False,actual_parameters=parameters,
               parameters_evidence='native_input_parameters' if item['kind']=='original' else 'constructor_parameters_from_pinned_logged_execution')
    return row,arrays


def execute(path,group_id,output):
    return execute_group(preflight(path),path,group_id,output)


def execute_group(m,path,group_id,output):
    """Execute one already validated group; explicit optional conductor prefactor."""
    g=next(g for g in m['groups'] if g['group_id']==group_id);factor=m.get('energy_prefactor',1.)
    pyddx=None
    if g['new_states']:
        import pyddx
        if not os.environ.get('SLURM_JOB_ID') or int(os.environ.get('SLURM_CPUS_ON_NODE','0'))<64 or record(pyddx.__file__)!=m['module'] or pyddx.__version__!='0.9.0' or np.__version__!='1.26.4':raise InvalidArtifact('unqualified allocation or module')
    d=Path(path).resolve().parent/'groups'/group_id/'attempt_0001';d.mkdir(parents=True,exist_ok=False)
    start=time.monotonic();cpu=time.process_time();data=read_json(verify(g['source']));arrays={}
    centres=np.asfortranarray(np.array(data['coordinates_A']).T/BOHR_TO_A);radii=np.array(data['radii_A'])/BOHR_TO_A
    r=dict(group_id=group_id,cache_key=g['cache_key'],source=g['source'],manifest=record(path),rows={},actual_parameters={},reciprocity_error_kcal=None)
    for label in g['states']:
        before=time.monotonic();model=None;state=None;a={};row=dict(state=label,status='failed',execution='new',forward_solve_started=False,new_forward_solve_started=False,solver_maxiter=1200)
        write_new(d/(label+'_admission.json'),dict(manifest=record(path),group=group_id,state=label,reused=label in g['reuse'],time_unix=time.time(),slurm_job_id=os.environ.get('SLURM_JOB_ID')))
        try:
            if label in g['reuse']:
                row,a=retained_state(g['reuse'][label])
            else:
                model=pyddx.Model(sphere_centres=centres,sphere_radii=radii,**m['model'],**m['grids'][g['grid']],logfile=str(d/(label+'_solver.log')))
                row['actual_parameters']={k:v for k,v in model.input_parameters.items() if k not in ('sphere_centres','sphere_radii')};row['parameters_evidence']='native_input_parameters'
                q=np.zeros(len(radii)) if label=='zero' else np.array(data['endpoints'][label.split('_')[0]]['charge_e'])
                cavity=np.array(model.cavity);phi=direct_potential(cavity,centres,q);psi=np.array(model.multipole_psi(np.asfortranarray(q.reshape(1,-1)/np.sqrt(4*np.pi))))
                expected=np.zeros_like(psi);expected[0]=np.sqrt(4*np.pi)*q;pe=potential_error(cavity,centres,q,phi);se=float(np.max(abs(psi-expected)))
                a.update(cavity_bohr=cavity,phi=phi,psi=psi);row.update(potential_max_error_au=pe,psi_max_error=se,charge_sum_e=math.fsum(q),preparation_wall_seconds=time.monotonic()-before)
                if pe>m['tolerances']['potential_au'] or se>m['tolerances']['psi']:raise InvalidArtifact('source potential/integral gate failed')
                state=pyddx.State(model,np.asfortranarray(psi),phi);row['forward_solve_started']=row['new_forward_solve_started']=True;solve=time.monotonic();state.solve(tol=m['tolerances']['solver'])
                x=np.array(state.x);raw_energy=float(state.energy());energy=raw_energy*factor;a['x']=x
                if not state.is_solved or not np.isfinite(energy) or not np.isfinite(x).all():raise InvalidArtifact('native solve incomplete/nonfinite')
                err=abs(factor*.5*math.fsum((psi*x).ravel())-energy)*HA_TO_KCAL
                row.update(status='computed',energy_hartree=energy,energy_kcal_mol=energy*HA_TO_KCAL,contraction_error_kcal=err,
                           raw_native_energy_hartree=raw_energy,energy_prefactor=factor,
                           contraction_pass=err<=m['tolerances']['energy_kcal'],passive_energy_pass=energy*HA_TO_KCAL<=m['tolerances']['energy_kcal'],
                           final_single_layer_iterations=state.x_n_iter,solve_wall_seconds=time.monotonic()-solve)
            if 'cavity_bohr' in a:
                if 'cavity_bohr' in arrays and not np.array_equal(arrays['cavity_bohr'],a['cavity_bohr']):raise InvalidArtifact('matched cavity points changed')
                arrays['cavity_bohr']=a['cavity_bohr']
                for key in ('phi','psi','x'):
                    if key in a:arrays[label+'_'+key]=a[key]
        except Exception as exc:
            row.update(status='failed',failure_reason=str(exc),energy_hartree=None,energy_kcal_mol=None)
        row['current_execution_wall_seconds']=time.monotonic()-before
        if (d/(label+'_solver.log')).exists():row['solver_log']=record(d/(label+'_solver.log'))
        write_new(d/(label+'_result.json'),row);r['rows'][label]=row
        print(group_id,label,row['execution'],row['status'],flush=True)
        state=None;model=None;a.clear();gc.collect()
    if all(r['rows'][label]['status']=='computed' for label in ('Ca','La')):
        r['reciprocity_error_kcal']=factor*.5*abs(math.fsum((arrays['Ca_psi']*arrays['La_x']).ravel())-math.fsum((arrays['La_psi']*arrays['Ca_x']).ravel()))*HA_TO_KCAL
    np.savez_compressed(d/'arrays.npz',**arrays)
    r.update(status='complete' if all(row['status']=='computed' for row in r['rows'].values()) else 'incomplete',arrays=record(d/'arrays.npz'),
             wall_seconds=time.monotonic()-start,CPU_seconds=time.process_time()-cpu,peak_process_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
             actual_new_forward_solve_starts=sum(row['new_forward_solve_started'] for row in r['rows'].values()),slurm_job_id=os.environ.get('SLURM_JOB_ID'))
    r['runtime']=dict(hostname=platform.node(),python=platform.python_version(),numpy=np.__version__,pyddx=pyddx.__version__ if pyddx is not None else None,
                      environment={k:os.environ.get(k) for k in ('SLURM_JOB_ID','SLURMD_NODENAME','SLURM_CPUS_ON_NODE','SLURM_MEM_PER_NODE','OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS')},
                      cpu_models=sorted({line.split(':',1)[1].strip() for line in Path('/proc/cpuinfo').read_text().splitlines() if line.startswith('model name')}))
    write_new(d/'result.json',r);write_new(output,r);return r


def collect(path,output):
    m=preflight(path);root=Path(path).resolve().parent;groups={};receipts={}
    for g in m['groups']:
        p=root/'groups'/g['group_id']/'attempt_0001/result.json'
        if not p.exists():raise InvalidArtifact('group result not yet available: '+g['group_id'])
        r=read_json(p)
        if r['manifest']!=record(path) or r['cache_key']!=g['cache_key'] or set(r['rows'])!=set(g['states']):raise InvalidArtifact('group receipt mismatch')
        verify(r['arrays']);groups[g['group_id']]=r;receipts[g['group_id']]=record(p)
    r=dict(protocol=PROTOCOL,manifest=record(path),groups=groups,receipts=receipts,complete=all(g['status']=='complete' for g in groups.values()),
           requested_state_roles=20,reused_state_roles=m['reused_state_roles'],requested_new_forward_solves=m['requested_new_forward_solves'],
           actual_new_forward_solve_starts=sum(g['actual_new_forward_solve_starts'] for g in groups.values()),
           new_DFT_calls=0,new_MACE_calls=0,new_score=None,baseline_changed=False)
    write_new(output,r);return r


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest='command',required=True)
    q=sub.add_parser('prepare');q.add_argument('--config',required=True);q.add_argument('--output',required=True)
    for op in ('dry-run','execute','collect'):
        q=sub.add_parser(op);q.add_argument('--manifest',required=True)
        if op in ('execute','collect'):q.add_argument('--output',required=True)
        if op=='execute':q.add_argument('--group',required=True)
    a=p.parse_args()
    if a.command=='prepare':r=prepare(a.config,a.output)
    elif a.command=='dry-run':r=preflight(a.manifest)
    elif a.command=='execute':r=execute(a.manifest,a.group,a.output)
    else:r=collect(a.manifest,a.output)
    print(json.dumps({k:v for k,v in r.items() if k not in ('groups','implementation','rows')}))
