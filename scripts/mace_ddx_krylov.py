"""Two real endpoint qualification of GMRES with unchanged native ddX operators."""
from pathlib import Path
import argparse
import gc
import json
import math
import os
import platform
import resource
import shutil
import time
import numpy as np
from affordable_common import BOHR_TO_A,HA_TO_KCAL,InvalidArtifact,cache_key,read_json,record,verify,write_new
from mace_ddx_recovery import preflight as parent_preflight
from mace_ddx_source_self import direct_potential,potential_error

PROTOCOL='fixed_GGR_2FW0_coarse_native_ddPCM_GMRES_qualification_v1'
SETTINGS=dict(restart=40,maxiter=30,rtol=1e-10,atol=0.,callback_type='pr_norm')
TOL=dict(norm_relative=1e-12,residual_relative=1e-9,energy_kcal=1e-6)


def prepare(config,output):
    c=read_json(config);parent=parent_preflight(verify(c['parent']));software=read_json(verify(c['software_receipt']));group=read_json(verify(c['source_group']))
    source=next(g for g in parent['groups'] if g['group_id']=='GGR_2FW0_primary_coarse')
    if software['status']!='built_and_imported' or group['manifest']!=c['parent'] or group['cache_key']!=source['cache_key'] or group['status']!='complete':raise InvalidArtifact('completed native references and isolated adapter required')
    sm=read_json(verify(software['manifest']))
    if sm['protocol']!='isolated_pyddx_native_operator_adapter_build_v1' or sm['plan']!=c['plan']:raise InvalidArtifact('adapter build scope differs')
    root=Path(output).resolve();root.mkdir(parents=True,exist_ok=False);impl=root/'implementation';impl.mkdir()
    for name,pin in parent['implementation'].items():shutil.copyfile(verify(pin),impl/name)
    shutil.copyfile(__file__,impl/Path(__file__).name)
    py=verify(c['software_receipt']).parent/'venv/bin/python';pypin=record(py);pypin['path']=str(py)
    module=list((py.parent.parent/'lib/python3.11/site-packages').glob('pyddx*.so'))
    if len(module)!=1:raise InvalidArtifact('one isolated adapter module required')
    m=dict(protocol=PROTOCOL,config=record(config),plan=c['plan'],parent=c['parent'],source_group=c['source_group'],software_receipt=c['software_receipt'],
           python=pypin,module=record(module[0]),source=source['source'],model=parent['model'],grid=parent['grids']['coarse'],
           solver=SETTINGS,tolerances=TOL,source_tolerances=parent['tolerances'],implementation={p.name:record(p) for p in impl.glob('*.py')},
           states=['Ca','La'],requested_forward_endpoints=2,maximum_linear_system_solves=4,new_DFT_calls=0,new_MACE_calls=0)
    m['cache_key']=cache_key(m);write_new(root/'manifest.json',m);return preflight(root/'manifest.json')


def preflight(path):
    m=read_json(path);parent=parent_preflight(verify(m['parent']));g=next(g for g in parent['groups'] if g['group_id']=='GGR_2FW0_primary_coarse');r=read_json(verify(m['source_group']))
    if m['protocol']!=PROTOCOL or m['cache_key']!=cache_key({k:v for k,v in m.items() if k!='cache_key'}) or m['solver']!=SETTINGS or m['tolerances']!=TOL or m['source_tolerances']!=parent['tolerances'] or m['states']!=['Ca','La'] or m['model']!=parent['model'] or m['grid']!=parent['grids']['coarse'] or m['source']!=g['source'] or r['manifest']!=m['parent'] or r['cache_key']!=g['cache_key'] or r['status']!='complete':raise InvalidArtifact('declared Krylov state or settings changed')
    for pin in (m['config'],m['plan'],m['source'],m['software_receipt'],m['module'],m['python'],*m['implementation'].values(),r['arrays']):verify(pin)
    software=read_json(verify(m['software_receipt']));sm=read_json(verify(software['manifest']))
    if software['status']!='built_and_imported' or sm['protocol']!='isolated_pyddx_native_operator_adapter_build_v1' or sm['plan']!=m['plan']:raise InvalidArtifact('adapter receipt or scope differs')
    for pin in sm['pins']:verify(pin)
    return m


def execute(path,output):
    import pyddx
    import scipy
    from scipy.sparse.linalg import LinearOperator,gmres
    m=preflight(path)
    if not os.environ.get('SLURM_JOB_ID') or int(os.environ.get('SLURM_CPUS_ON_NODE','0'))<64 or record(pyddx.__file__)!=m['module'] or getattr(pyddx,'alquemia_pcm_kernel_version',None)!='native_pcm_operators_v1' or scipy.__version__!='1.17.1' or np.__version__!='1.26.4':raise InvalidArtifact('unqualified allocation or adapter')
    root=Path(path).resolve().parent;attempt=root/'attempt_0001';attempt.mkdir(exist_ok=False);start=time.monotonic();cpu=time.process_time()
    data=read_json(verify(m['source']));old=read_json(verify(m['source_group']))
    with np.load(verify(old['arrays']),allow_pickle=False) as f:saved={k:np.array(f[k]) for k in f.files}
    coords=np.asfortranarray(np.array(data['coordinates_A']).T/BOHR_TO_A);radii=np.array(data['radii_A'])/BOHR_TO_A
    rows={};all_checks=[];total_starts=0
    for metal in m['states']:
        before=time.monotonic();arrays={};model=None;state=None;calls={str(i):0 for i in range(1,6)};checks=[]
        row=dict(state=metal,status='failed',energy_hartree=None,energy_kcal_mol=None,linear_system_starts=0,stages={})
        def check(name,error,tolerance):
            checks.append(dict(name=name,error=float(error) if np.isfinite(error) else None,nonfinite_observation=None if np.isfinite(error) else str(error),tolerance=tolerance,pass_=bool(np.isfinite(error) and abs(error)<=tolerance)))
        try:
            model=pyddx.Model(sphere_centres=coords,sphere_radii=radii,**m['model'],**m['grid']);shape=(model.n_basis,model.n_spheres)
            q=np.array(data['endpoints'][metal]['charge_e']);cavity=np.array(model.cavity)
            if not np.array_equal(cavity,saved['cavity_bohr']):raise InvalidArtifact('native cavity-point identity changed')
            phi=direct_potential(cavity,coords,q);psi=np.array(model.multipole_psi(np.asfortranarray(q.reshape(1,-1)/np.sqrt(4*np.pi))))
            if not np.array_equal(phi,saved[metal+'_phi']) or not np.array_equal(psi,saved[metal+'_psi']):raise InvalidArtifact('source potential/integral identity changed')
            check('source_potential_au',potential_error(cavity,coords,q,phi),m['source_tolerances']['potential_au'])
            state=pyddx.State(model,np.asfortranarray(psi),phi);signed_phi=np.array(state.pcm_source_vector())
            weights=np.repeat([1/math.sqrt((l+1)*model.n_spheres) for l in range(m['grid']['lmax']+1)],[2*l+1 for l in range(m['grid']['lmax']+1)])[:,None]
            def native_norm(a):return float(state.pcm_native_norm(np.asfortranarray(a)))
            def ratio(a,b):
                n=native_norm(a);d=native_norm(b)
                return n/d if d else (0. if n==0. else float('inf'))
            def apply(op,a):
                calls[str(op)]+=1
                return np.array(state.pcm_operator(op,np.asfortranarray(a)))
            def norm_check(name,a):
                actual=native_norm(a);transformed=float(np.linalg.norm((weights*a).ravel()))
                check(name,(transformed-actual)/actual if actual else transformed,TOL['norm_relative'])
            norm_check('source_norm_equivalence',signed_phi);norm_check('reference_norm_equivalence',saved[metal+'_x'])
            rhs=apply(1,signed_phi);old_res=apply(2,apply(3,saved[metal+'_x']))-rhs
            check('saved_native_composed_residual',ratio(old_res,rhs),TOL['residual_relative'])
            arrays.update(phi=phi,psi=psi,cavity_bohr=cavity,native_signed_phi=signed_phi,native_R_infinity_rhs=rhs,saved_native_composed_residual=old_res)
            def solve(label,op,prec,b):
                nonlocal total_starts
                t=time.monotonic();history=[];n=b.size
                def wrapped(kind,v):return (weights*apply(kind,v.reshape(shape,order='F')/weights)).ravel(order='F')
                A=LinearOperator((n,n),matvec=lambda v:wrapped(op,v),dtype=np.float64)
                P=LinearOperator((n,n),matvec=lambda v:wrapped(prec,v),dtype=np.float64)
                log=attempt/(metal+'_'+label+'_iterations.jsonl')
                def callback(value):
                    history.append(float(value))
                    with log.open('a') as f:f.write(json.dumps(dict(iteration=len(history),preconditioned_relative_residual=float(value)))+'\n')
                row['linear_system_starts']+=1;total_starts+=1
                z,info=gmres(A,(weights*b).ravel(order='F'),M=P,callback=callback,**m['solver'])
                x=z.reshape(shape,order='F')/weights;res=apply(op,x)-b;relative=ratio(res,b)
                row['stages'][label]=dict(info=int(info),iterations=len(history),true_relative_native_hnorm_residual=relative,wall_seconds=time.monotonic()-t,
                    log=record(log) if log.exists() else None)
                arrays[label+'_solution']=x;arrays[label+'_residual']=res
                check(label+'_true_residual',relative,TOL['residual_relative']);norm_check(label+'_norm_equivalence',x)
                if info!=0 or relative>TOL['residual_relative']:raise InvalidArtifact('GMRES '+label+' did not satisfy status/residual criteria')
                return x
            phieps=solve('dielectric',2,4,rhs);x=solve('single_layer',3,5,phieps);e=.5*math.fsum((psi*x).ravel())
            arrays['x']=x;check('native_energy_equivalence_kcal',e*HA_TO_KCAL-old['rows'][metal]['energy_kcal_mol'],TOL['energy_kcal'])
            check('passive_energy_kcal',max(0.,e*HA_TO_KCAL),m['source_tolerances']['energy_kcal'])
            row.update(status='computed',energy_hartree=e,energy_kcal_mol=e*HA_TO_KCAL)
        except Exception as exc:row['failure_reason']=str(exc)
        np.savez_compressed(attempt/(metal+'_arrays.npz'),**arrays)
        row.update(arrays=record(attempt/(metal+'_arrays.npz')),operator_and_preconditioner_calls=calls,checks=checks,checks_pass=bool(checks) and all(c['pass_'] for c in checks),wall_seconds=time.monotonic()-before)
        write_new(attempt/(metal+'_result.json'),row);rows[metal]=row;all_checks.extend(dict(c,state=metal) for c in checks)
        print(metal,row['status'],row['checks_pass'],row['stages'],flush=True)
        state=None;model=None;arrays.clear();gc.collect()
    complete=all(r['status']=='computed' for r in rows.values());contrast=None;delta=None
    if complete:
        contrast=rows['Ca']['energy_kcal_mol']-rows['La']['energy_kcal_mol'];native=old['rows']['Ca']['energy_kcal_mol']-old['rows']['La']['energy_kcal_mol'];delta=contrast-native
        all_checks.append(dict(name='native_paired_contrast_equivalence_kcal',error=delta,tolerance=TOL['energy_kcal'],pass_=abs(delta)<=TOL['energy_kcal']))
    result=dict(protocol=PROTOCOL,manifest=record(path),rows=rows,complete=complete,checks=all_checks,checks_pass=complete and all(c['pass_'] for c in all_checks),
        actual_forward_endpoint_starts=sum(r['linear_system_starts']>0 for r in rows.values()),actual_linear_system_starts=total_starts,
        source_self_R_kcal=contrast,native_R_difference_kcal=delta,reference_reciprocity_error_kcal=old['reciprocity_error_kcal'],
        new_DFT_calls=0,new_MACE_calls=0,new_score=None,baseline_changed=False,wall_seconds=time.monotonic()-start,CPU_seconds=time.process_time()-cpu,
        peak_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,slurm_job_id=os.environ['SLURM_JOB_ID'],hostname=platform.node(),scipy=scipy.__version__,
        interpretation='same PCM equations, numerical solver qualification only; coarse reciprocity failure remains')
    write_new(output,result);return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest='command',required=True)
    q=sub.add_parser('prepare');q.add_argument('--config',required=True);q.add_argument('--output',required=True)
    for op in ('dry-run','execute'):
        q=sub.add_parser(op);q.add_argument('--manifest',required=True)
        if op=='execute':q.add_argument('--output',required=True)
    a=p.parse_args()
    if a.command=='prepare':r=prepare(a.config,a.output)
    elif a.command=='dry-run':r=preflight(a.manifest)
    else:r=execute(a.manifest,a.output)
    print(json.dumps({k:v for k,v in r.items() if k not in ('rows','checks','implementation')}))
