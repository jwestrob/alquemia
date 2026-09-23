"""Isolated six-angle warm starts and strict-native five-geometry comparison."""
from __future__ import annotations
import argparse, copy, fcntl, importlib.util, json, os, re, shutil, time, traceback
from pathlib import Path
import numpy as np
from scipy.optimize import minimize
from affordable_common import InvalidArtifact, HA_TO_KCAL, cache_key, paired, read_json, record, verify, write_new, xyz
from accommodation_proposals import MACEProposal
from accommodation_nonlinear import WarmGPU
from adaptive_force_diagnostic import preview, project
from mace_site_kinematics import Kinematics
from mace_hybrid import EV_TO_KCAL, write_xyz
from union_adaptive import snapshot, origin as archived_origin
from slsqp_precision import SETTINGS as FOUR_SETTINGS
from strict_native_pool import recipe, collect_fresh
from compact_solvation import diagnostics
from structure_informed_starts import scf_details
from run_orca_task_manifest import load_manifest_tasks
from nikasha_pool import choose_rows, relative_components
from accommodation_folds_compare import decision
from accommodation_fold_proposals import outcome

PROTOCOL = 'common_six_angular_warm_native_OMOL_strict_GFN2_v1'
SETTINGS = {**FOUR_SETTINGS, 'selected_angular_modes': 6, 'start': 'own_existing_four_angle_proposal'}
CANDIDATES = ('origin', 'adaptive_Ca', 'adaptive_La', 'six_Ca', 'six_La')
FEASIBILITY_SHA = 'd580f4bf6ede0b5fd312602475832b4e29c6b481c7e5475dcbac26702e6b18af'


def data(pin): return read_json(verify(pin))


def private(name):
    spec = importlib.util.spec_from_file_location('_six_private_'+name, Path(__file__).with_name(name+'.py'))
    module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module); return module


def full_q(task, active):
    active = np.asarray(active, dtype=float)
    if active.shape != (6,) or not np.isfinite(active).all(): raise InvalidArtifact('six finite angles required')
    q = np.zeros(task['mode_count']); q[task['active_indices']] = active; return q


angular = private('adaptive_angular_proposals')
angular.full_q = full_q
angular.SETTINGS = copy.deepcopy(SETTINGS)


def cached_point(task, name):
    point = copy.deepcopy(task[name]); mapping = data(task['mapping'])['context']
    native = data(point['MACE']); model = data(task['source_manifest'])['model']
    if 'request' in native:
        request = data(native['request'])
        if native['status'] != 'complete' or native['model'] != model or request['xyz'] != point['coordinate'] or request['charge'] != task['charge'] or request['multiplicity'] != task['multiplicity']:
            raise InvalidArtifact('cached native state receipt differs')
    else:
        archived_origin({'native_MACE_receipt':point['MACE'],'native_MACE_energy_eV':point['MACE_eV'],
                         'charge':task['charge'],'multiplicity':task['multiplicity'],'xyz':point['coordinate']},model)
    if native['forces'] != point['forces'] or native['energy_eV'] != point['MACE_eV']:
        raise InvalidArtifact('cached force/energy differs')
    coords, _, raw, _, _ = project(mapping, point['full_q'], np.load(verify(point['forces']), allow_pickle=False))
    if not np.allclose(coords, [a[1:] for a in xyz(verify(point['coordinate']))], atol=1e-12, rtol=0):
        raise InvalidArtifact('cached geometry differs')
    point.update(active_q_radian=np.asarray(point['full_q'])[task['active_indices']].tolist(),
                 gradient_kcal_mol_rad=raw[task['active_indices']].tolist())
    return point


class SixProposal(MACEProposal):
    def evaluate(self, q, purpose):
        q = np.asarray(q, dtype=float)
        for name in ('origin_cached', 'warm_cached'):
            point = self.task[name]
            if np.array_equal(q, np.asarray(point['full_q'])[self.task['active_indices']]):
                value = cached_point(self.task, name)
                key = cache_key({'manifest': self.manifest, 'task': self.task['task_id'], 'q': q.tolist()})
                value.update(manifest=self.manifest, cache_key=key, reused_scientific_point=name)
                p = self.directory/'evaluations'/key/'result.json'
                if not p.exists(): p.parent.mkdir(parents=True, exist_ok=False); write_new(p, value)
                if read_json(p) != value: raise InvalidArtifact('cached point changed')
                self.requests.append({'purpose': purpose, 'result': record(p), 'reused': True, 'scientific_point': name})
                return value
        return super().evaluate(q, purpose)


def prepare(feasibility, reference, agreement, output):
    if record(feasibility)['sha256'] != FEASIBILITY_SHA: raise InvalidArtifact('frozen feasibility differs')
    f = read_json(feasibility); a = data(f['c5ax_audit']); s = data(f['strict225'])
    strict = {r['case_id']:r for r in s['actual_pools']}; audit = {r['case_id']:r for r in a['rows']}
    labels = {r['case_id']:r['expected_class'] for r in s['rows']}
    ref = read_json(reference)
    if record(reference)['sha256'] != '171bd31d466ff97ef6073ce23286af8519ef23690f7d6ad23bb444cab30ecd5c': raise InvalidArtifact('frozen strict reference differs')
    out = Path(output).resolve(); out.mkdir(parents=True, exist_ok=False); tasks = []; cases = []
    p = data(f['rows'][0]['proposal_manifest'])
    for row in f['rows']:
        cid = row['case_id']; parent = data(row['proposal_manifest'])
        for key in ('model','software','orca','cpu_python','gpu_python'):
            if parent[key] != p[key]: raise InvalidArtifact('native runtime differs across sources')
        old = strict.get(cid) or audit[cid]['strict32']['branches']['fresh']
        cases.append({'case_id':cid, 'expected_class':labels.get(cid,'La'), 'old_matrix':old['matrix'],
                      'old_pool':old['pool'], 'selection':row['selection'], 'prior_source':row['proposal_manifest']})
        for z in ('Ca','La'):
            ep = row['endpoints'][z]; t = copy.deepcopy(ep['task']); r = data(ep['proposal_receipt'])
            t.update(active_indices=row['active_indices'], active_mode_ids=row['selected_ids'],
                     active_roles=['adaptive_physical_angular']*6, selector=row['selection'],
                     original_active_mode_ids=t['active_mode_ids'], origin_cached=r['origin'], warm_cached=r['proposal'],
                     prior_receipt=ep['proposal_receipt'], source_manifest=row['proposal_manifest'])
            tasks.append(t)
    m = {k:p[k] for k in ('model','software','orca','cpu_python','gpu_python','cpu_executable','gpu_executable','optimizer_software')}
    m.update(protocol_id=PROTOCOL, settings=SETTINGS, feasibility=record(feasibility), reference=record(reference),
             agreement=record(agreement), cases=cases, tasks=tasks, declared_case_ids=f['declared_sources'],
             maximum_optimizer_starts=18, maximum_cross_MACE=18, maximum_GFN2=72, reused_strict_cells=108,
             resources={'GPU_cpus':32,'GPUs':1,'GPU_host_mem_MiB':200000,'GFN2_cpus':32,'mpi_ranks':1,'concurrent_tasks':32},
             implementation=snapshot(out/'implementation'), new_q0_calls=0, new_DFT_calls=0)
    mp=out/'manifest.json'; write_new(mp,m); v=validate(mp); write_new(out/'PREFLIGHT.json',v); return v


def validate(manifest):
    m=read_json(manifest); f=data(m['feasibility'])
    if m['protocol_id']!=PROTOCOL or m['settings']!=SETTINGS or m['feasibility']['sha256']!=FEASIBILITY_SHA or m['declared_case_ids']!=f['declared_sources'] or len(m['tasks'])!=18: raise InvalidArtifact('fixed scope differs')
    for key in ('agreement','reference','software','orca','cpu_executable','gpu_executable'): verify(m[key])
    for pin in m['implementation'].values(): verify(pin)
    if [t['task_id'] for t in m['tasks']] != [cid+'__'+z for cid in m['declared_case_ids'] for z in ('Ca','La')]: raise InvalidArtifact('endpoint denominator differs')
    checks=[]
    for c, row in zip(m['cases'],f['rows']):
        ts=[next(t for t in m['tasks'] if (t['case_id'],t['metal'])==(c['case_id'],z)) for z in ('Ca','La')]
        paired(verify(ts[1]['xyz']),verify(ts[0]['xyz']),ts[1]['charge'],ts[0]['charge'])
        if data(ts[0]['mapping'])!=data(ts[1]['mapping']):raise InvalidArtifact('paired map differs')
        origin=[]
        for t in ts:
            ep=row['endpoints'][t['metal']]; original=ep['task']
            if any(t[k]!=original[k] for k in ('xyz','charge','multiplicity','mapping','source_preparation','q0','mode_count')) or t['prior_receipt']!=ep['proposal_receipt'] or t['active_mode_ids']!=row['selected_ids'] or t['active_indices']!=row['active_indices']:raise InvalidArtifact('state/mode mapping differs')
            kin=Kinematics(data(t['mapping'])['context']); names=[v['id'] for v in kin.modes]
            if kin.metal!=0 or xyz(verify(t['xyz']))[0][0]!=t['metal']:raise InvalidArtifact('metal index differs')
            if t['active_mode_ids'] != [names[i] for i in t['active_indices']] or t['active_mode_ids'][:4]!=original['active_mode_ids']:raise InvalidArtifact('six-mode prefix differs')
            old=data(t['prior_receipt'])
            if t['warm_cached']!=old['proposal'] or t['origin_cached']!=old['origin']:raise InvalidArtifact('warm/source point changed')
            for name in ('origin_cached','warm_cached'):
                point=cached_point(t,name); q=np.asarray(point['full_q']); active=q[t['active_indices']]
                if not np.array_equal(q,full_q(t,active)):raise InvalidArtifact('nonselected coordinate nonzero')
                if name=='warm_cached' and np.any(active[4:]!=0):raise InvalidArtifact('new warm coordinates nonzero')
                physical=angular.final_geometry(kin,t,active,[a[0] for a in xyz(verify(t['xyz']))])
                checks.append({'task_id':t['task_id'],'point':name,'geometry':physical})
            origin.append(project(kin.data,t['origin_cached']['full_q'],np.load(verify(t['origin_cached']['forces']),allow_pickle=False)))
        choice=preview(kin.modes,origin[0][1],origin[0][3],origin[1][3],maximum=6)
        if choice!=row['selection'] or any(t['selector']!=choice for t in ts):raise InvalidArtifact('force selector differs')
        if choose_rows(c['old_matrix'],CANDIDATES[:3])!=c['old_pool']:raise InvalidArtifact('old pool replay differs')
    return {'status':'validated','manifest':record(manifest),'sources':9,'searches':18,'maximum_cross_MACE':18,'maximum_GFN2':72,'old_strict_cells':108,'checks':checks,'new_molecular_calls':0}


def optimize(manifest,task,gpu):
    ev=SixProposal(manifest,task,gpu); dest=ev.directory/'result.json'
    if dest.exists(): raise InvalidArtifact('prior attempt remains; no automatic retry')
    began=time.monotonic(); trace=[]; accepted=[]
    result={'task_id':task['task_id'],'case_id':task['case_id'],'metal':task['metal'],'manifest':record(manifest),
            'status':'unavailable','origin':None,'proposal':None,'optimizer':None,'unconstrained_minimum_claimed':False,'composite_gradient':None}
    def evaluate(q,purpose):
        trial={'active_q_radian':np.asarray(q).tolist(),'purpose':purpose,**angular.physical_status(ev.kin,task,q),'MACE_complete':False};trace.append(trial)
        point=ev.evaluate(q,purpose);trial['MACE_complete']=True;return point
    try:
        origin=evaluate(np.zeros(6),'origin');result['origin']=origin
        warmq=np.asarray(task['warm_cached']['full_q'])[task['active_indices']]
        warm=evaluate(warmq,'warm_start');result['warm_start']=warm
        if abs((origin['MACE_eV']-task['q0']['components']['MACE_eV'])*EV_TO_KCAL)>SETTINGS['origin_replay_tolerance_kcal_mol']:raise InvalidArtifact('origin replay failed')
        def objective(q):
            p=evaluate(q,'SLSQP_trial')
            return (p['MACE_eV']-origin['MACE_eV'])*EV_TO_KCAL/HA_TO_KCAL,np.asarray(p['gradient_kcal_mol_rad'])/HA_TO_KCAL
        opt=minimize(objective,warmq,jac=True,method='SLSQP',bounds=[(-.8,.8)]*6,
            constraints=[{'type':'ineq','fun':lambda q:angular.constraints(ev.kin,task,q)[0],'jac':lambda q:angular.constraints(ev.kin,task,q)[1]}],
            callback=lambda q:accepted.append(np.asarray(q).tolist()),options={'maxiter':200,'ftol':1e-8})
        result['optimizer']={'success':bool(opt.success),'message':str(opt.message),'iterations':int(opt.nit),'function_evaluations':int(opt.nfev),'gradient_evaluations':int(opt.njev),'ftol_hartree_equivalent':1e-8}
        check=angular.final_geometry(ev.kin,task,opt.x,[a[0] for a in ev.atoms]);result['final_geometry']=check
        candidate=evaluate(opt.x,'final_candidate');result['proposal']=candidate
        if not opt.success:raise InvalidArtifact('SLSQP did not converge: '+str(opt.message))
        work=candidate['MACE_eV']-origin['MACE_eV']
        if work>SETTINGS['maximum_final_energy_increase_eV']:raise InvalidArtifact('native final energy increased')
        _,_,raw,normed,_=project(ev.kin.data,candidate['full_q'],np.load(verify(candidate['forces']),allow_pickle=False))
        active=set(task['active_indices']);angles=[i for i,v in enumerate(ev.kin.modes) if v['unit']=='radian']
        result.update(status='proposal_available',MACE_proposal_work_kcal_mol=float(work*EV_TO_KCAL),
            MACE_added_work_kcal_mol=float((candidate['MACE_eV']-warm['MACE_eV'])*EV_TO_KCAL),
            boundary_flag=bool(check['boundary_heavy_atoms'] or check['angular_boundary_flag']),
            residual={'force_hamiltonian':'native_vacuum_OMOL','all_mode_ids':[v['id'] for v in ev.kin.modes],
                      'raw_gradient_kcal_per_unit':raw.tolist(),'normalized_gradient_kcal_mol_A':normed.tolist(),
                      'selected_max_load':float(max(abs(normed[i]) for i in active)),
                      'omitted_angular_max_load':float(max(abs(normed[i]) for i in angles if i not in active))})
    except Exception as exc:result.update(reason=str(exc),traceback=traceback.format_exc())
    result.update(trials=trace,accepted_iterations=accepted,requests=ev.requests,
        infeasible_trial_requests=sum(not t['physical_feasible'] for t in trace),
        maximum_trial_heavy_extent_A=max((t['maximum_heavy_displacement_A'] for t in trace),default=None),
        wall_seconds=time.monotonic()-began,job_id=os.environ.get('SLURM_JOB_ID'))
    write_new(dest,result);return record(dest)


def execute(manifest):
    validate(manifest);m=read_json(manifest);root=Path(manifest).resolve().parent
    if not os.environ.get('SLURM_JOB_ID') or int(os.environ.get('SLURM_CPUS_ON_NODE','0'))!=32 or int(os.environ.get('SLURM_MEM_PER_NODE','0'))!=200000:raise InvalidArtifact('existing H200 allocation required')
    with (root/'proposal.lock').open('a+') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB);began=time.monotonic();gpu=None;results=[];cells=[];error=None
        try:
            gpu=WarmGPU(manifest)
            for t in m['tasks']:
                pin=optimize(manifest,t,gpu);r=data(pin);results.append(pin)
                print(json.dumps({'task_id':t['task_id'],'status':r['status']}),flush=True)
                if r['status']!='proposal_available':continue
                for z in ('Ca','La'):
                    own=z==t['metal'];candidate='six_'+t['metal'];other=next(x for x in m['tasks'] if (x['case_id'],x['metal'])==(t['case_id'],z))
                    row={'case_id':t['case_id'],'candidate':candidate,'metal':z,'charge':other['charge'],'multiplicity':1,'status':'unavailable','proposal_receipt':pin}
                    try:
                        if own:xp=r['proposal']['coordinate'];native=data(r['proposal']['MACE']);npin=r['proposal']['MACE']
                        else:
                            d=root/'cross'/(t['case_id']+'__'+candidate+'__'+z);d.mkdir(parents=True,exist_ok=False)
                            atoms=xyz(verify(r['proposal']['coordinate']));atoms[0]=(z,*atoms[0][1:]);write_xyz(d/'core.xyz',atoms);xp=record(d/'core.xyz')
                            req={'manifest':record(manifest),'task_id':other['task_id'],'xyz':xp,'charge':other['charge'],'multiplicity':1}
                            write_new(d/'request.json',req);native,npin=gpu.evaluate(record(d/'request.json'))
                        if native['status']!='complete':raise InvalidArtifact('native cross evaluation failed')
                        row.update(status='complete',xyz=xp,MACE=npin,MACE_eV=native['energy_eV'],native_reused=own)
                    except Exception as exc:row.update(reason=str(exc))
                    cells.append(row)
        except Exception as exc:error=repr(exc)
        finally:
            if gpu is not None:
                try:gpu.close()
                except Exception as exc:error=(error+'; ' if error else '')+repr(exc)
            elapsed=time.monotonic()-began
            receipt={'manifest':record(manifest),'results':results,'cells':cells,'error':error,'wall_seconds':elapsed,'allocated_core_seconds':elapsed*32,'allocated_GPU_seconds':elapsed,'job_id':os.environ['SLURM_JOB_ID'],'GPU_command':gpu.command if gpu else None}
            write_new(root/'GPU_COLLECTION.json',receipt)
    if error:raise InvalidArtifact(error)
    return {'status':'finished','searches':len(results),'native_cells':len(cells)}


def prepare_low(manifest,gpu_collection,output):
    validate(manifest);m=read_json(manifest);g=read_json(gpu_collection)
    if g['manifest']!=record(manifest):raise InvalidArtifact('GPU source differs')
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False);tasks=[];missing=[]
    cells={(r['case_id'],r['candidate'],r['metal']):r for r in g['cells']}
    for c in m['cases']:
        for q in CANDIDATES[3:]:
            for z in ('Ca','La'):
                cell=cells.get((c['case_id'],q,z));task=next(t for t in m['tasks'] if (t['case_id'],t['metal'])==(c['case_id'],z))
                for medium in ('vacuum','alpb'):
                    identity={'case_id':c['case_id'],'candidate':q,'metal':z,'medium':medium};tid='__'.join(identity.values())
                    if not cell or cell['status']!='complete':missing.append({**identity,'task_id':tid,'reason':'new native candidate unavailable'});continue
                    old=c['old_matrix'][z]['origin'];low=old.get('low',old.get('actual_low'))[medium]
                    oldtask=next(t for t in data(low['manifest'])['tasks'] if t['task_id']==low['task_id'])
                    audit=diagnostics(low,oldtask);details=scf_details(verify(low['output']).read_text())
                    d=out/'tasks'/tid;d.mkdir(parents=True);shutil.copyfile(verify(cell['xyz']),d/'core.xyz');(d/'endpoint.inp').write_text(recipe(task['charge'],medium,'fresh'))
                    tasks.append({**identity,'task_id':tid,'charge':task['charge'],'multiplicity':1,'seed_source':None,'gradient_requested':False,
                        'source':{'electron_count':details['electrons'],'parameter_export':audit['parameter_export'],'state_reference':low},
                        'native_cell':cell,'xyz':record(d/'core.xyz'),'input':record(d/'endpoint.inp'),'output_path':str(d/'endpoint.out'),
                        'active_seed_path':str(d/'endpoint.runtime.xtbw'),'active_gbw_path':str(d/'endpoint.runtime.gbw')})
    impl=snapshot(out/'implementation');lm={'protocol_id':PROTOCOL+'_scalar','branch':'fresh','source_manifest':record(manifest),
        'GPU_collection':record(gpu_collection),'tasks':tasks,'missing':missing,'orca':m['orca'],'implementation':impl,'cell_denominator':72,
        'execution_resources':{'mpi_ranks':1,'concurrent_tasks':32},'execution_policy':{'task_runner':impl['run_orca_task_manifest.py'],'runtime_renderer':impl['render_orca_runtime_input.py']}}
    mp=out/'manifest.json';write_new(mp,lm);v=validate_low(mp,True);write_new(out/'PREFLIGHT.json',v);return v


def validate_low(manifest,fresh=False):
    m=read_json(manifest);parent=data(m['source_manifest']);g=data(m['GPU_collection'])
    if m['protocol_id']!=PROTOCOL+'_scalar' or m['branch']!='fresh' or m['execution_resources']!={'mpi_ranks':1,'concurrent_tasks':32}:raise InvalidArtifact('scalar policy differs')
    expected={(c['case_id'],q,z,s) for c in parent['cases'] for q in CANDIDATES[3:] for z in ('Ca','La') for s in ('vacuum','alpb')}
    rows=m['tasks']+m['missing']
    if len(rows)!=72 or {tuple(t[k] for k in ('case_id','candidate','metal','medium')) for t in rows}!=expected:raise InvalidArtifact('72-cell denominator differs')
    for pin in m['implementation'].values():verify(pin)
    load_manifest_tasks(Path(manifest))
    for t in m['tasks']:
        source=next(r for r in g['cells'] if (r['case_id'],r['candidate'],r['metal'])==(t['case_id'],t['candidate'],t['metal']))
        if t['native_cell']!=source or t['charge']!=source['charge'] or t['multiplicity']!=1 or t['seed_source'] is not None or t['gradient_requested'] or verify(t['xyz']).read_bytes()!=verify(source['xyz']).read_bytes() or verify(t['input']).read_text()!=recipe(t['charge'],t['medium'],'fresh'):raise InvalidArtifact('scalar geometry/state/recipe differs')
        if fresh and {p.name for p in Path(t['output_path']).parent.iterdir()}!={'core.xyz','endpoint.inp'}:raise InvalidArtifact('nonfresh scalar directory')
    return {'status':'validated','manifest':record(manifest),'new_calls':len(m['tasks']),'missing':len(m['missing']),'cell_denominator':72}


def execute_low(manifest):
    engine=private('strict_native_pool');engine.validate=validate_low;return engine.execute(manifest)


def collect(manifest,output):
    validate_low(manifest);m=read_json(manifest);mp=Path(manifest).resolve();p=data(m['source_manifest']);g=data(m['GPU_collection'])
    def idx(name):
        path=mp.parent/name;return {r['task_id']:r for r in read_json(path)['rows']} if path.exists() else {}
    before,after=idx('SEEDS_BEFORE.json'),idx('SEEDS_AFTER.json');rows=[]
    for t in m['tasks']:
        r=collect_fresh(mp,t,before.get(t['task_id']),after.get(t['task_id']))
        if r['status']=='complete':
            text=verify(r['actual']['output']).read_text();tol=re.findall(r'Energy Change\s+TolE\s+\.{4}\s+([-+0-9.eE]+)',text)
            if len(tol)!=1 or float(tol[0])!=1e-10 or data(r['actual']['receipt'])['parallelism']['nprocs']!=1:r.update(status='audit_failed',reason='effective tolerance/rank differs',energy_hartree=None)
        rows.append(r)
    rows.extend({**r,'status':'unavailable','energy_hartree':None} for r in m['missing']);lookup={(r['case_id'],r['candidate'],r['metal'],r['medium']):r for r in rows}
    native={(r['case_id'],r['candidate'],r['metal']):r for r in g['cells']};cases=[];ref=data(p['reference'])['branches']['fresh']
    for c in p['cases']:
        matrix=copy.deepcopy(c['old_matrix']);cid=c['case_id']
        for z in ('Ca','La'):
            for q in CANDIDATES[3:]:
                rr={s:lookup[cid,q,z,s] for s in ('vacuum','alpb')};n=native.get((cid,q,z));ok=bool(n and n['status']=='complete' and all(v['status']=='complete' for v in rr.values()))
                matrix[z][q]={'status':'complete' if ok else 'unavailable','xyz':n.get('xyz') if n else None,'MACE':n.get('MACE') if n else None,
                    'components':{'MACE_eV':n['MACE_eV'],'GFN2_vacuum_hartree':rr['vacuum']['energy_hartree'],'GFN2_ALPB_hartree':rr['alpb']['energy_hartree']} if ok else None,
                    'low':{s:r.get('actual') for s,r in rr.items()}}
        pool=choose_rows(matrix,CANDIDATES);calls={};added={}
        for v in ('mathematical','operational'):
            val=pool[v]['composite_R_model_kcal_mol'] if pool['status']=='available' else None;oldval=c['old_pool'][v]['composite_R_model_kcal_mol'];call=decision(val,ref['variants'][v]['bands'])
            calls[v]={'R':val,'delta_R':val-oldval if val is not None else None,'decision':call,'outcome':outcome(call,c['expected_class']),'old_R':oldval,'old_decision':decision(oldval,ref['variants'][v]['bands'])}
        for z in ('Ca','La'):
            n=matrix[z]['six_'+z];old=matrix[z]['adaptive_'+z]
            added[z]=relative_components(n['components'],old['components']) if n['status']=='complete' else None
        cases.append({'case_id':cid,'expected_class':c['expected_class'],'status':pool['status'],'matrix':matrix,'pool':pool,'old_pool':c['old_pool'],'old_band_transfer':calls,'own_added_work':added})
    result={'protocol_id':PROTOCOL,'manifest':record(manifest),'source_manifest':m['source_manifest'],'reference':p['reference'],
        'GPU_collection':m['GPU_collection'],'rows':rows,'cases':cases,'cell_denominator':72,'complete_cells':sum(r['status']=='complete' for r in rows),
        'case_denominator':9,'complete_cases':sum(c['status']=='available' for c in cases),'new_calibration':False,'production_changed':False}
    write_new(output,result);return {k:v for k,v in result.items() if k not in ('rows','cases')}


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);sub=ap.add_subparsers(dest='op',required=True)
    for op,fields in {'prepare':('feasibility','reference','agreement','output'),'validate':('manifest',),'execute':('manifest',),
        'prepare_low':('manifest','gpu_collection','output'),'validate_low':('manifest',),'execute_low':('manifest',),'collect':('manifest','output')}.items():
        cmd=sub.add_parser(op)
        for field in fields:cmd.add_argument('--'+field.replace('_','-'),required=True,type=Path)
    args=vars(ap.parse_args());print(json.dumps(globals()[args.pop('op')](**args),indent=2))
