"""Bounded nonlinear donor accommodation with warm OMOL and native GFN2 gradients.

Experimental composite energies only. No DFT, entropy, refit or default change.
CPU owns physical mappings/ORCA; the isolated GPU process owns native OMOL.
"""
from __future__ import annotations
import argparse
import fcntl
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import time
import traceback
import numpy as np
from affordable_common import HA_TO_KCAL, InvalidArtifact, cache_key, read_json, record, verify, write_new, xyz
from mace_hybrid import EV_TO_KCAL, check_atoms, write_xyz

PROTOCOL='bounded_composite_donor_accommodation_v1'
CASES=('1H4I','4MAE','PQQSEQ_83440678cbbd658047c9','PQQSEQ_07ab500e3df76b30d71c')
SETTINGS={'method':'L-BFGS-B','bounds_radian':[-.8,.8],'maxiter':24,'maxfun':80,'maxls':10,
          'gtol':.2,'ftol':1e-12,'interior_margin_radian':.02,'maximum_energy_increase_kcal_mol':.1,
          'hessian_steps_radian':[.005,.010],'hessian_absolute_tolerance':.1,'hessian_relative_tolerance':.01,
          'new_overlap_H_A':.55,'new_overlap_heavy_A':1.0}
RESOURCES={'cpus':32,'gpus':1,'host_mem_MiB':200000,'mpi_ranks':8,'concurrent_tasks':2}
ROLE_ORDER=('anchor_glutamate','extra_acidic_ligand_homolog')


def prepare(source,primary_result,agreement,output,cpu_python,gpu_python):
    from mace_site_kinematics import Kinematics
    from accommodation_torsion_profiles import role_mode, geometry_check
    from accommodation_response import recipe
    from affordable_workflow import dry_run
    old=read_json(source)
    primary=read_json(primary_result)
    if read_json(verify(primary['manifest']))['design']!=record(source):raise InvalidArtifact('primary result describes another source design')
    if read_json(verify(primary['MACE_result']))['model']!=old['model']:raise InvalidArtifact('primary native checkpoint differs')
    if tuple(c['case_id'] for c in old['cases'])!=CASES:raise InvalidArtifact('exact original four-case torsion design required')
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False)
    impl=out/'implementation';impl.mkdir();pins={}
    for p in Path(__file__).parent.glob('*.py'):
        dest=impl/p.name;shutil.copyfile(p,dest);pins[p.name]=record(dest)
    tasks=[]
    for c in old['cases']:
        parent=read_json(verify(c['source']['parent']))
        roles=[r for r in ROLE_ORDER if r in c['modes']]
        if roles!=['anchor_glutamate']+([] if c['case_id']=='1H4I' else ['extra_acidic_ligand_homolog']):
            raise InvalidArtifact('original active donor role scope differs')
        expected=[role_mode(parent,r) for r in roles]
        if any(c['modes'][r]!=mid for r,mid in zip(roles,expected)):raise InvalidArtifact('declared mode differs from actual role chemistry')
        for z in ('Ca','La'):
            ep=c['origins'][z];mapping=read_json(verify(c['maps'][z]));kin=Kinematics(mapping['context']);atoms=xyz(verify(ep['xyz']))
            names=[m['id'] for m in kin.modes];active=[names.index(v) for v in expected]
            if any(kin.modes[i]['unit']!='radian' for i in active):raise InvalidArtifact('donor mode is not an angle')
            q=np.zeros(len(names));checks=geometry_check(kin,q,[a[0] for a in atoms]);check_atoms(atoms,ep['charge'])
            if not checks['pass'] or not np.allclose(kin.evaluate(q)[1],[a[1:] for a in atoms],atol=1e-12,rtol=0):
                raise InvalidArtifact('mapping does not reproduce source')
            # Mapping derivative check at a fixed nonzero joint geometry; no energy.
            q[active]=.01;probe=geometry_check(kin,q,[a[0] for a in atoms])
            if not probe['pass']:raise InvalidArtifact('active-coordinate physical mapping check failed')
            tasks.append({'task_id':c['case_id']+'__'+z,'case_id':c['case_id'],'metal':z,'xyz':ep['xyz'],
                'charge':ep['charge'],'multiplicity':ep['multiplicity'],'mapping':c['maps'][z],'source_preparation':c['preparation'],
                'active_roles':roles,'active_mode_ids':expected,'active_indices':active,'mode_count':len(names),
                'origin_geometry_checks':checks,'joint_mapping_check':probe,'label_scope':c['label_scope']})
    manifest={'protocol_id':PROTOCOL,'source':record(source),'primary_result':record(primary_result),'agreement':record(agreement),'settings':SETTINGS,
        'resources':RESOURCES,'tasks':tasks,'model':old['model'],'software':old['software'],'orca':old['orca'],
        'implementation':pins,'cpu_python':str(Path(cpu_python).absolute()),'gpu_python':str(Path(gpu_python).absolute()),
        'cpu_executable':record(cpu_python),'gpu_executable':record(gpu_python),'new_DFT_calls':0,'baseline_changed':False,
        'algorithm_starts':8,'gradient_recipe':{'vacuum':recipe(0,'vacuum',True),'alpb':recipe(0,'alpb',True)}}
    write_new(out/'manifest.json',manifest)
    # Actual origin-only finite manifests allow existing runtime dry-run checks.
    preflights=[]
    for t in tasks:
        low=low_manifest(manifest,t,t['xyz'],out/'preflight'/t['task_id'])
        preflights.append(dry_run(low))
    result=validate(out/'manifest.json');result['origin_manifest_preflights']=preflights
    write_new(out/'PREFLIGHT.json',result);return result


def validate(manifest):
    m=read_json(manifest)
    if m['protocol_id']!=PROTOCOL or m['settings']!=SETTINGS or m['resources']!=RESOURCES:raise InvalidArtifact('frozen nonlinear policy changed')
    if [(t['case_id'],t['metal']) for t in m['tasks']]!=[(c,z) for c in CASES for z in ('Ca','La')]:raise InvalidArtifact('eight-start scope changed')
    for key in ('source','primary_result','agreement','software','orca','cpu_executable','gpu_executable'):verify(m[key])
    for pin in m['implementation'].values():verify(pin)
    for which in ('cpu','gpu'):
        if record(m[which+'_python'])!=m[which+'_executable']:raise InvalidArtifact('interpreter invocation changed')
    from mace_omol import model
    if model(verify(m['software']))!=m['model']:raise InvalidArtifact('native checkpoint changed')
    for t in m['tasks']:
        for key in ('xyz','mapping','source_preparation'):verify(t[key])
        check_atoms(xyz(verify(t['xyz'])),t['charge'])
        if t['multiplicity']!=1:raise InvalidArtifact('unsupported spin')
    from affordable_common import paired
    for case in CASES:
        ca,la=[next(t for t in m['tasks'] if t['case_id']==case and t['metal']==z) for z in ('Ca','La')]
        paired(verify(la['xyz']),verify(ca['xyz']),la['charge'],ca['charge'])
        if (ca['active_mode_ids']!=la['active_mode_ids'] or ca['active_indices']!=la['active_indices']
                or read_json(verify(ca['mapping']))!=read_json(verify(la['mapping']))):
            raise InvalidArtifact('paired physical axes, moving atom groups or nuclear mapping differ')
    return {'status':'prepared_and_validated','manifest':record(manifest),'starts':8,'new_DFT_calls':0}


def low_manifest(m,t,coordinate_pin,directory):
    from accommodation_response import recipe
    directory=Path(directory);tasks=[]
    for medium in ('vacuum','alpb'):
        td=directory/medium;td.mkdir(parents=True,exist_ok=False);xp=td/'core.xyz';shutil.copyfile(verify(coordinate_pin),xp)
        ip=td/'endpoint.inp';ip.write_text(recipe(t['charge'],medium,True))
        tasks.append({'task_id':t['task_id']+'__'+medium,'case':t['case_id'],'case_id':t['case_id'],'metal':t['metal'],
            'medium':medium,'charge':t['charge'],'multiplicity':t['multiplicity'],'xyz':record(xp),'input':record(ip),
            'output_path':str(td/'endpoint.out'),'physical_mapping':t['mapping']})
    mm={'protocol_id':PROTOCOL,'agreement':m['agreement'],'orca':m['orca'],'implementation':m['implementation'],'tasks':tasks,
        'execution_policy':{'task_runner':m['implementation']['run_orca_task_manifest.py'],
                            'runtime_renderer':m['implementation']['render_orca_runtime_input.py']},
        'execution_resources':{'mpi_ranks':RESOURCES['mpi_ranks'],'concurrent_tasks':RESOURCES['concurrent_tasks']}}
    path=directory/'manifest.json';write_new(path,mm);return path


def low_collect(manifest):
    from compact_solvation import completed,diagnostics
    from accommodation_response import parse_gradient
    m=read_json(manifest);rows={}
    for t in m['tasks']:
        row={'status':'unavailable','energy_hartree':None,'task_id':t['task_id']}
        op=Path(t['output_path']);text=op.read_text() if op.exists() else ''
        row.update(ORCA_execution_receipt_present=Path(str(op)+'.execution.json').exists(),
                   native_SCF_started=bool(re.search(r'Iteration\s+Energy \(Eh\)\s+Delta-E',text)))
        pin=completed(manifest,t['task_id'])
        if pin:
            try:
                audit=diagnostics(pin,t);text=verify(pin['output']).read_text()
                tol=re.findall(r'Energy Change\s+TolE\s+\.{4}\s+([-+0-9.eE]+)',text)
                if (audit['charge_sanity_status']!='pass' or 'INFO: Using special xTB SCF mixer' not in text
                        or len(tol)!=1 or float(tol[0])>1e-8):raise InvalidArtifact('native analytic Tight policy not confirmed')
                gradient=parse_gradient(t,pin)
                # Existing parser returns a q0 projection. Reproject its actual
                # Cartesian gradient at the current q below instead.
                gradient.pop('projected_gradient_kcal_per_unit')
                row.update(pin);row.update(audit);row.update(gradient);row['status']='complete'
            except Exception as exc:row.update(reason=str(exc))
        else:
            row['available_artifacts']=[record(p) for p in (Path(t['output_path']),Path(t['output_path']+'.execution.json')) if p.exists()]
        rows[t['medium']]=row
    result={'manifest':record(manifest),'rows':rows,'status':'complete' if all(v['status']=='complete' for v in rows.values()) else 'unavailable'}
    write_new(Path(manifest).parent/'collection.json',result);return result


def gpu_server(manifest):
    """Persistent isolated worker; all requested geometries and results stay on disk."""
    import torch
    from ase import Atoms
    from mace.calculators import mace_omol
    from mace_omol import input_batch
    import resource
    if not os.environ.get('SLURM_JOB_ID') or not torch.cuda.is_available():raise InvalidArtifact('allocated GPU required')
    m=read_json(manifest);torch.set_num_threads(1);torch.set_default_dtype(torch.float64);torch.cuda.reset_peak_memory_stats()
    began=time.monotonic();calc=mace_omol(model=str(verify(m['model']['checkpoint'])),device='cuda',default_dtype='float64');model=calc.models[0]
    if (type(model).__name__!='ScaleShiftMACE' or dict(model.embedding_specs)!=m['model']['embedding_specs'] or list(model.heads)!=['omol'] or calc.energy_units_to_eV!=1.):raise InvalidArtifact('native model differs')
    for p in model.parameters():p.requires_grad_(False)
    load=time.monotonic()-began;calls=0;failed=0
    print('NONLINEAR_JSON '+json.dumps({'status':'ready','load_seconds':load}),flush=True)
    for line in sys.stdin:
        request=json.loads(line);pin=request['request'];rp=verify(pin);r=read_json(rp);result={'request':pin,'status':'failed','energy_eV':None,'model_call_started':False};tick=time.monotonic()
        try:
            if r['manifest']!=record(manifest):raise InvalidArtifact('GPU request manifest mismatch')
            t=next(t for t in m['tasks'] if t['task_id']==r['task_id'])
            if (r['charge'],r['multiplicity'])!=(t['charge'],t['multiplicity']):raise InvalidArtifact('GPU state mismatch')
            rows=xyz(verify(r['xyz']));atoms=Atoms([a[0] for a in rows],positions=[a[1:] for a in rows],pbc=False)
            atoms.info.update(charge=t['charge'],spin=t['multiplicity']);atoms.calc=calc
            state=input_batch(calc,atoms,t['charge'],t['multiplicity']);calls+=1;result['model_call_started']=True
            energy=float(atoms.get_potential_energy());forces=np.asarray(atoms.get_forces(),dtype=float)
            if not np.isfinite(energy) or forces.shape!=(len(rows),3) or not np.isfinite(forces).all():raise InvalidArtifact('native energy/force unavailable')
            fp=rp.parent/'forces_eV_A.npy';np.save(fp,forces)
            result.update(status='complete',energy_eV=energy,forces=record(fp),state=state,force_definition='negative Cartesian energy gradient, eV/angstrom')
        except Exception as exc:
            failed+=1;result.update(reason=str(exc),traceback=traceback.format_exc())
        result.update(wall_seconds=time.monotonic()-tick,model=m['model'],job_id=os.environ['SLURM_JOB_ID'],device=torch.cuda.get_device_name())
        dest=rp.parent/'mace_result.json';write_new(dest,result)
        print('NONLINEAR_JSON '+json.dumps({'result':record(dest)}),flush=True)
    dest=Path(manifest).parent/('gpu_summary_'+os.environ['SLURM_JOB_ID']+'.json')
    write_new(dest,{'manifest':record(manifest),'new_MACE_calls':calls,'failed_requests':failed,'model_load_seconds':load,
        'wall_seconds':time.monotonic()-began,'peak_host_RSS_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        'peak_cuda_allocated_bytes':torch.cuda.max_memory_allocated(),'job_id':os.environ['SLURM_JOB_ID'],'device':torch.cuda.get_device_name()})


class WarmGPU:
    def __init__(self,manifest):
        m=read_json(manifest);self.manifest=manifest
        self.log=(Path(manifest).parent/('gpu_log_'+os.environ['SLURM_JOB_ID']+'.txt')).open('x')
        script=verify(m['implementation']['accommodation_nonlinear.py'])
        self.command=['srun','--exact','--overlap','--ntasks=1','--cpus-per-task=1','--gres=gpu:1',m['gpu_python'],str(script),'gpu-server','--manifest',str(manifest)]
        self.process=subprocess.Popen(self.command,stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,bufsize=1)
        try:
            message=self.receive()
            if message.get('status')!='ready':raise InvalidArtifact('warm worker did not initialize')
        except Exception:
            # Terminate only the child step this instance created; no external
            # job, process, queue or shared runtime is touched.
            self.process.stdin.close()
            if self.process.poll() is None:
                self.process.terminate()
                try:self.process.wait(timeout=5)
                except subprocess.TimeoutExpired:self.process.kill();self.process.wait()
            self.log.close()
            raise
    def receive(self):
        while True:
            line=self.process.stdout.readline()
            if not line:raise InvalidArtifact('GPU worker terminated: '+str(self.process.poll()))
            self.log.write(line);self.log.flush()
            if line.startswith('NONLINEAR_JSON '):return json.loads(line[len('NONLINEAR_JSON '):])
    def evaluate(self,pin):
        self.process.stdin.write(json.dumps({'request':pin})+'\n');self.process.stdin.flush()
        message=self.receive();result=read_json(verify(message['result']))
        if result['request']!=pin:raise InvalidArtifact('GPU response belongs to another geometry')
        return result,message['result']
    def close(self):
        self.process.stdin.close()
        for line in self.process.stdout:self.log.write(line)
        code=self.process.wait();self.log.close()
        if code:raise InvalidArtifact('GPU worker exit code '+str(code))


def relative_components(value,origin):
    mace=(value['MACE_eV']-origin['MACE_eV'])*EV_TO_KCAL
    solvent=((value['GFN2_ALPB_hartree']-origin['GFN2_ALPB_hartree'])-
             (value['GFN2_vacuum_hartree']-origin['GFN2_vacuum_hartree']))*HA_TO_KCAL
    return {'native_MACE_kcal_mol':mace,'solvent_transfer_kcal_mol':solvent,'composite_kcal_mol':mace+solvent}


def project_components(jac,forces_eV_A,vac_gradient_kcal_A,alpb_gradient_kcal_A):
    native=-np.einsum('mij,ij->m',jac,np.asarray(forces_eV_A))*EV_TO_KCAL
    transfer=np.einsum('mij,ij->m',jac,np.asarray(alpb_gradient_kcal_A)-np.asarray(vac_gradient_kcal_A))
    return native,transfer,native+transfer


def candidate_gate(success,q,gradient,work,geometry):
    q=np.asarray(q);gradient=np.asarray(gradient)
    checks={'optimizer_success':bool(success),'absolute_gradient_pass':bool(np.max(np.abs(gradient))<=SETTINGS['gtol']),
            'interior_pass':bool(np.min(SETTINGS['bounds_radian'][1]-np.abs(q))>=SETTINGS['interior_margin_radian']),
            'nonincreasing_within_numerical_scale':bool(work<=SETTINGS['maximum_energy_increase_kcal_mol']),
            'geometry_pass':bool(geometry['pass'])}
    return {**checks,'eligible_for_curvature':all(checks.values()),'max_absolute_gradient_kcal_mol_rad':float(np.max(np.abs(gradient))),
            'boundary_distance_radian':float(np.min(SETTINGS['bounds_radian'][1]-np.abs(q)))}


def curvature_gate(fine,coarse):
    fine=np.asarray(fine);coarse=np.asarray(coarse)
    if fine.shape!=coarse.shape or fine.ndim!=2 or fine.shape[0]!=fine.shape[1] or not np.isfinite([fine,coarse]).all():raise InvalidArtifact('invalid projected curvature')
    scale=float(max(np.max(np.abs(fine)),np.max(np.abs(coarse))))
    tol=max(SETTINGS['hessian_absolute_tolerance'],SETTINGS['hessian_relative_tolerance']*scale)
    symmetry=max(float(np.max(np.abs(fine-fine.T))),float(np.max(np.abs(coarse-coarse.T))))
    error=float(np.max(np.abs(fine-coarse)));spectral=float(np.linalg.norm(fine-coarse,2))
    # Symmetrization is only for a diagnostic eigenvalue; failed asymmetry still
    # disqualifies the matrix. No eigenvalue clipping or soft-mode deletion.
    eig=np.linalg.eigvalsh((fine+fine.T)/2)
    return {'fine_H_kcal_mol_rad2':fine.tolist(),'coarse_H_kcal_mol_rad2':coarse.tolist(),
        'tolerance_kcal_mol_rad2':tol,'symmetry_error':symmetry,'refinement_max_error':error,'refinement_spectral_norm':spectral,
        'eigenvalues_unclipped':eig.tolist(),'status':'pass' if symmetry<=tol and error<=tol and eig[0]>spectral+.1 else 'fail',
        'entropy':None}


class Evaluator:
    def __init__(self,manifest,task,gpu):
        from mace_site_kinematics import Kinematics
        self.path=Path(manifest);self.manifest=read_json(manifest);self.pin=record(manifest);self.task=task;self.gpu=gpu
        self.kin=Kinematics(read_json(verify(task['mapping']))['context']);self.atoms=xyz(verify(task['xyz']))
        self.directory=self.path.parent/'endpoints'/task['task_id'];self.directory.mkdir(parents=True,exist_ok=True)
        self.trace=self.directory/('requests_'+os.environ['SLURM_JOB_ID']+'.jsonl');self.new_calls=0;self.cache_hits=0
    def log(self,row):
        with self.trace.open('a') as f:f.write(json.dumps(row,allow_nan=False)+'\n');f.flush();os.fsync(f.fileno())
    def __call__(self,active_q,purpose='optimizer'):
        from accommodation_torsion_profiles import geometry_check
        from affordable_workflow import execute as existing_execute
        active_q=np.asarray(active_q,dtype=float)
        if active_q.shape!=(len(self.task['active_indices']),) or not np.isfinite(active_q).all() or np.max(np.abs(active_q))>.8+1e-12:
            raise InvalidArtifact('trial outside declared active angular domain')
        q=np.zeros(self.task['mode_count']);q[self.task['active_indices']]=active_q
        _,coords,_,jac=self.kin.evaluate(q)
        key=cache_key({'manifest':self.pin,'task':self.task['task_id'],'active_q':active_q.tolist()})
        directory=self.directory/'evaluations'/key;result_path=directory/'result.json'
        if result_path.exists():
            result=read_json(result_path)
            if result['cache_key']!=key or result['manifest']!=self.pin:raise InvalidArtifact('evaluation cache mismatch')
            for pin in result.get('artifacts',{}).values():verify(pin)
            self.cache_hits+=1;self.log({'purpose':purpose,'q':active_q.tolist(),'cache_reuse':True,'result':record(result_path),'status':result['status']})
            if result['status']!='complete':raise InvalidArtifact('previous failed evaluation retained: '+result.get('reason','unknown'))
            return result
        if directory.exists():raise InvalidArtifact('partial evaluation retained; no implicit retry')
        directory.mkdir(parents=True);start=time.monotonic()
        result={'cache_key':key,'manifest':self.pin,'task_id':self.task['task_id'],'active_q_radian':active_q.tolist(),'full_q':q.tolist(),
                'purpose_first_requested':purpose,'status':'failed','artifacts':{},'new_MACE_calls':0,'new_GFN2_calls':0}
        try:
            checks=geometry_check(self.kin,q,[a[0] for a in self.atoms]);result['geometry_checks']=checks
            if not checks['pass']:raise InvalidArtifact('unsupported physical trial geometry')
            xp=directory/'context.xyz'
            if np.array_equal(active_q,np.zeros_like(active_q)):shutil.copyfile(verify(self.task['xyz']),xp)
            else:write_xyz(xp,[(a[0],*p) for a,p in zip(self.atoms,coords)])
            request={'manifest':self.pin,'task_id':self.task['task_id'],'xyz':record(xp),'charge':self.task['charge'],'multiplicity':self.task['multiplicity']}
            request_path=directory/'mace_request.json';write_new(request_path,request)
            native,pin=self.gpu.evaluate(record(request_path));result['artifacts']['MACE']=pin;result['new_MACE_calls']=1
            if native['status']!='complete':raise InvalidArtifact('native MACE unavailable: '+native.get('reason','unknown'))
            low=low_manifest(self.manifest,self.task,record(xp),directory/'GFN2');result['artifacts']['low_manifest']=record(low)
            result['new_GFN2_calls']=2;error=None
            try:existing_execute(low)
            except Exception as exc:error=str(exc)
            collected=low_collect(low);result['artifacts']['low_collection']=record(low.parent/'collection.json')
            if collected['status']!='complete':raise InvalidArtifact('native analytic solvent unavailable: '+str(error))
            vac,alpb=(collected['rows'][medium] for medium in ('vacuum','alpb'))
            if read_json(verify(vac['parameter_export']))!=read_json(verify(alpb['parameter_export'])):raise InvalidArtifact('matched GFN2 media parameters differ')
            forces=np.load(verify(native['forces']),allow_pickle=False)
            selected_jac=jac[self.task['active_indices']]
            ng,sg,cg=project_components(selected_jac,forces,vac['gradient_kcal_mol_A'],alpb['gradient_kcal_mol_A'])
            components={'MACE_eV':native['energy_eV'],'GFN2_vacuum_hartree':vac['energy_hartree'],'GFN2_ALPB_hartree':alpb['energy_hartree']}
            result.update(status='complete',components=components,coordinate=record(xp),
                native_gradient_kcal_mol_rad=ng.tolist(),solvent_gradient_kcal_mol_rad=sg.tolist(),gradient_kcal_mol_rad=cg.tolist(),
                solvent_transfer_kcal_mol=(alpb['energy_hartree']-vac['energy_hartree'])*HA_TO_KCAL,
                analytic_gradient=True,numerical_DFT_gradients=False)
            self.new_calls+=1
        except Exception as exc:result.update(reason=str(exc),traceback=traceback.format_exc())
        result['wall_seconds']=time.monotonic()-start;write_new(result_path,result)
        self.log({'purpose':purpose,'q':active_q.tolist(),'cache_reuse':False,'result':record(result_path),'status':result['status']})
        if result['status']!='complete':raise InvalidArtifact(result['reason'])
        return result


def optimize_endpoint(manifest,task,gpu):
    from scipy.optimize import minimize
    ev=Evaluator(manifest,task,gpu);path=ev.directory/'result.json'
    if path.exists():
        result=read_json(path)
        if result['manifest']!=record(manifest):raise InvalidArtifact('endpoint completion belongs to a different plan')
        return record(path)
    zero=np.zeros(len(task['active_indices']));origin=None;start=time.monotonic();accepted=[]
    result={'task_id':task['task_id'],'case_id':task['case_id'],'metal':task['metal'],'manifest':record(manifest),'status':'unavailable',
            'candidate':None,'origin':None,'accommodation_work':None,'curvature':{'status':'not_eligible'},'stationary_minimum_qualified':False}
    try:
        origin=ev(zero,'origin');result['origin']=origin
        def objective(q):
            value=ev(q);work=relative_components(value['components'],origin['components'])
            return work['composite_kcal_mol'],np.asarray(value['gradient_kcal_mol_rad'])
        def callback(q):accepted.append({'q':q.tolist(),'relative_energy_kcal_mol':objective(q)[0]})
        opt=minimize(objective,zero,jac=True,method='L-BFGS-B',bounds=[(-.8,.8)]*len(zero),callback=callback,
            options={k:SETTINGS[k] for k in ('maxiter','maxfun','maxls','gtol','ftol')})
        candidate=ev(opt.x,'final_candidate');work=relative_components(candidate['components'],origin['components'])
        gate=candidate_gate(opt.success,opt.x,candidate['gradient_kcal_mol_rad'],work['composite_kcal_mol'],candidate['geometry_checks'])
        result.update(status='candidate_available',candidate=candidate,accommodation_work=work,stationarity=gate,
            optimizer={'success':bool(opt.success),'message':str(opt.message),'iterations':int(opt.nit),'function_evaluations':int(opt.nfev),'gradient_evaluations':int(opt.njev),'x':opt.x.tolist()})
        if gate['eligible_for_curvature']:
            matrices=[];curvature_points=[]
            try:
                for h in SETTINGS['hessian_steps_radian']:
                    columns=[]
                    for i in range(len(zero)):
                        dq=np.zeros_like(zero);dq[i]=h
                        plus=ev(opt.x+dq,'curvature');minus=ev(opt.x-dq,'curvature')
                        columns.append((np.asarray(plus['gradient_kcal_mol_rad'])-np.asarray(minus['gradient_kcal_mol_rad']))/(2*h))
                        curvature_points.extend([plus['cache_key'],minus['cache_key']])
                    matrices.append(np.asarray(columns).T)
                result['curvature']=curvature_gate(*matrices);result['curvature']['evaluation_keys']=curvature_points
                result['stationary_minimum_qualified']=result['curvature']['status']=='pass'
            except Exception as exc:result['curvature']={'status':'unavailable','reason':str(exc),'evaluation_keys':curvature_points}
    except Exception as exc:result.update(reason=str(exc),traceback=traceback.format_exc())
    result.update(accepted_iterations=accepted,wall_seconds=time.monotonic()-start,successful_new_evaluations=ev.new_calls,cache_hits=ev.cache_hits,
                  numerical_relaxation_free_energy=None,entropy=None,job_id=os.environ['SLURM_JOB_ID'])
    write_new(path,result);return record(path)


def contrast_components(ca,la):
    mace=(ca['MACE_eV']-la['MACE_eV'])*EV_TO_KCAL
    solv=((ca['GFN2_ALPB_hartree']-ca['GFN2_vacuum_hartree'])-(la['GFN2_ALPB_hartree']-la['GFN2_vacuum_hartree']))*HA_TO_KCAL
    return {'native_MACE_R_model_kcal_mol':mace,'solvent_R_kcal_mol':solv,'composite_R_model_kcal_mol':mace+solv}


def collect(manifest,output):
    from accommodation_torsion_analysis import donor_distances
    m=read_json(manifest);out=Path(manifest).parent;endpoints=[];index={}
    for t in m['tasks']:
        rp=out/'endpoints'/t['task_id']/'result.json'
        if rp.exists():
            row=read_json(rp)
            if row['manifest']!=record(manifest):raise InvalidArtifact('endpoint manifest changed')
            row={**row,'receipt':record(rp)}
        else:row={'task_id':t['task_id'],'case_id':t['case_id'],'metal':t['metal'],'status':'unavailable','reason':'endpoint not completed','origin':None,'candidate':None}
        row['reported_donor_contacts_A']={}
        for state in ('origin','candidate'):
            point=row.get(state)
            row['reported_donor_contacts_A'][state]=(donor_distances({'mapping':t['mapping'],'q':point['full_q']}) if point else None)
        endpoints.append(row);index[t['case_id'],t['metal']]=row
    cases=[];primary=read_json(verify(m['primary_result']))
    for case in CASES:
        pair={z:index[case,z] for z in ('Ca','La')};row={'case_id':case,'R0':None,'R_candidate':None,'endpoint_work':{z:pair[z].get('accommodation_work') for z in pair},'both_minima_qualified':all(p.get('stationary_minimum_qualified',False) for p in pair.values())}
        for name,key in (('R0','origin'),('R_candidate','candidate')):
            if all(p[key] is not None for p in pair.values()):row[name]=contrast_components(pair['Ca'][key]['components'],pair['La'][key]['components'])
        row['delta_R_model_kcal_mol']=row['R_candidate']['composite_R_model_kcal_mol']-row['R0']['composite_R_model_kcal_mol'] if row['R_candidate'] and row['R0'] else None
        old={z:next(p for p in primary['points'] if p['case_id']==case and p['metal']==z and p['point']=='origin') for z in ('Ca','La')}
        row['primary_to_Tight_origin_change']={}
        for z in ('Ca','La'):
            if pair[z]['origin'] is not None and old[z]['status']=='complete':
                fields={k:old[z][k] for k in ('MACE_eV','GFN2_vacuum_hartree','GFN2_ALPB_hartree')}
                row['primary_to_Tight_origin_change'][z]=relative_components(pair[z]['origin']['components'],fields)
            else:row['primary_to_Tight_origin_change'][z]=None
        cases.append(row)
    evaluations=[read_json(p) for p in out.glob('endpoints/*/evaluations/*/result.json')]
    native=[read_json(verify(e['artifacts']['MACE'])) for e in evaluations if e['artifacts'].get('MACE')]
    low=[r for e in evaluations if e['artifacts'].get('low_collection')
         for r in read_json(verify(e['artifacts']['low_collection']))['rows'].values()]
    result={'protocol_id':PROTOCOL,'manifest':record(manifest),'analysis_implementation':record(__file__),'endpoints':endpoints,'cases':cases,
        'endpoint_denominator':8,'candidate_endpoints':sum(e.get('candidate') is not None for e in endpoints),
        'qualified_minima':sum(e.get('stationary_minimum_qualified',False) for e in endpoints),
        'actual_evaluations':len(evaluations),'failed_evaluations':sum(e['status']!='complete' for e in evaluations),
        'MACE_evaluation_requests':sum(e['new_MACE_calls'] for e in evaluations),'GFN2_endpoint_requests':sum(e['new_GFN2_calls'] for e in evaluations),
        'observed_MACE_calls_started':sum(r['model_call_started'] for r in native),'MACE_calls_complete':sum(r['status']=='complete' for r in native),
        'observed_ORCA_execution_receipts':sum(r['ORCA_execution_receipt_present'] for r in low),
        'observed_GFN2_SCF_started':sum(r['native_SCF_started'] for r in low),'GFN2_analytic_endpoints_complete':sum(r['status']=='complete' for r in low),
        'new_DFT_calls':0,'entropy':None,'calibrated_decision':None,'baseline_changed':False,'biological_accuracy':None}
    write_new(output,result);return {k:v for k,v in result.items() if k not in ('endpoints','cases')}


def execute(manifest):
    validate(manifest)
    if not os.environ.get('SLURM_JOB_ID'):raise InvalidArtifact('allocation required')
    if int(os.environ['SLURM_CPUS_ON_NODE'])!=32 or int(os.environ['SLURM_NTASKS'])!=32 or int(os.environ['SLURM_MEM_PER_NODE'])!=200000:
        raise InvalidArtifact('allocation differs from declared32CPU/200000MiB layout')
    root=Path(manifest).parent;m=read_json(manifest);lock=(root/'execute.lock').open('a+')
    fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB);start=time.monotonic();gpu=None;receipts=[];error=None
    try:
        gpu=WarmGPU(manifest)
        for task in m['tasks']:
            receipt=optimize_endpoint(manifest,task,gpu);receipts.append(receipt);row=read_json(verify(receipt))
            print(json.dumps({'task_id':task['task_id'],'status':row['status'],'work':row.get('accommodation_work'),'qualified':row['stationary_minimum_qualified']}),flush=True)
    except Exception as exc:error=str(exc)
    finally:
        try:
            if gpu is not None:gpu.close()
        except Exception as exc:error=(error+'; ' if error else '')+str(exc)
        elapsed=time.monotonic()-start;write_new(root/('execution_'+os.environ['SLURM_JOB_ID']+'.json'),{'manifest':record(manifest),'endpoints':receipts,'job_id':os.environ['SLURM_JOB_ID'],'wall_seconds':elapsed,'allocated_core_seconds':elapsed*32,'allocated_GPU_seconds':elapsed,'error':error,'GPU_command':gpu.command if gpu is not None else None})
    result=collect(manifest,root/('result_'+os.environ['SLURM_JOB_ID']+'.json'))
    if error:raise InvalidArtifact(error)
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__);s=p.add_subparsers(dest='operation',required=True)
    q=s.add_parser('prepare')
    for n in ('source','primary-result','agreement','output','cpu-python','gpu-python'):q.add_argument('--'+n,required=True)
    for name in ('validate','execute','gpu-server','collect'):
        q=s.add_parser(name);q.add_argument('--manifest',required=True)
        if name=='collect':q.add_argument('--output',required=True)
    args=vars(p.parse_args());op=args.pop('operation').replace('-','_');result=globals()[op](**args)
    if result is not None:print(json.dumps(result,default=str))

if __name__=='__main__':main()
