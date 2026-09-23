"""Source-dihedral starts in the unchanged four-angle physical target space."""
from __future__ import annotations
import argparse
import copy
import json
import fcntl
import os
import re
import shutil
import time
import traceback
from pathlib import Path
import numpy as np
from affordable_common import HA_TO_KCAL, InvalidArtifact, paired, read_json, record, verify, write_new, xyz
from mace_site_kinematics import Kinematics
import adaptive_angular_proposals as angular

PROTOCOL='reference_dihedral_four_angular_alternative_starts_v1'
TEMPLATES={'Ca':'1H4I','La':'4MAE'}
FALLBACKS={'Ca':'q9z4j7-pqq-la_model','La':'a0a3f2yly8-pqq-la_model'}
SETTINGS={'maximum_heavy_displacement_A':.8,'maximum_angle_radian':.8,
          'meaningful_coordinate_change_A':1e-6,'dihedral_replay_tolerance_radian':1e-8,
          'template_oxygen_symmetry':'only_documented_deprotonated_Asp_Glu_terminal',
          'unmapped_chemical_role':'retain_target_q0','clipping_or_rescaling':False,
          'homology_independent_transfer_claimed':False}


def _role(value):
    if isinstance(value,dict):
        return {'chain':value['chain'],'resnum':int(value['resnum']),
                'icode':value.get('icode',value.get('insertion_code','')),'resname':value['resname']}
    match=re.fullmatch(r'([^:]+):([A-Z]{3})(-?\d+)([A-Za-z]?)',value)
    if not match:raise InvalidArtifact('unsupported explicit role selector: '+str(value))
    return dict(chain=match[1],resname=match[2],resnum=int(match[3]),icode=match[4])


def source_record(manifest, cid):
    m=read_json(manifest);c=next(x for x in m['cases'] if x['case_id']==cid)
    c=c.get('source_case',c)
    tasks=[next(t for t in m['tasks'] if t['case_id']==cid and t['metal']==z) for z in ('Ca','La')]
    ca,la=tasks
    paired(verify(la['xyz']),verify(ca['xyz']),la['charge'],ca['charge'])
    if ca['mapping']!=la['mapping'] and read_json(verify(ca['mapping']))!=read_json(verify(la['mapping'])):
        raise InvalidArtifact('paired physical maps differ')
    if ca['active_indices']!=la['active_indices']:raise InvalidArtifact('paired active coordinates differ')
    core=read_json(verify(c['source']['parent']))
    roles={k:_role(v) for k,v in core['fixed_core']['requested_roles'].items()}
    return {'case_id':cid,'biological_group':c['biological_group'],'source_manifest':record(manifest),
            'source_case':c,'tasks':tasks,'roles':roles,
            'fragments':{f['role']:f for f in core['qm_fragments'] if 'role'in f}}


def sources(inputs):
    data=read_json(inputs);out=[]
    for row in data['cases']:
        c=read_json(verify(row['pool_collection']))
        p=next(x for x in c['cases'] if x['case_id']==row['case_id'])
        if p['pool']['status']!='available':raise InvalidArtifact('required old common pool unavailable')
        col=read_json(verify(row['proposal_collection']))
        e=next(x for x in col['endpoints'] if x['case_id']==row['case_id'] and x['metal']=='Ca')
        r=read_json(verify(e['proposal_receipt']))
        out.append(source_record(verify(r['manifest']),row['case_id']) | {'base':row})
    original=out[0]['source_manifest']['path']
    templates={cid:source_record(original,cid) for cid in (*TEMPLATES.values(),*FALLBACKS.values())}
    return data,out,templates


def role_for_mode(src, kin, mode):
    if mode['kind']!='sidechain_torsion':raise InvalidArtifact('template has no peptide-crankshaft equivalent')
    a,b=mode['axis_indices'];meta=kin.data['source_atom_metadata'][a]
    matches=[r for r,v in src['roles'].items() if (v['chain'],v['resnum'],v['icode'])==
             (meta['chain'],meta['resnum'],meta['insertion_code'])]
    if len(matches)!=1:raise InvalidArtifact('selected torsion lacks unique canonical role')
    role=matches[0];res=src['roles'][role]['resname'];chi=mode['id'].rsplit('/',1)[1]
    names={('GLU','chi1'):('N','CA','CB','CG'),('GLU','chi2'):('CA','CB','CG','CD'),
           ('GLU','chi3'):('CB','CG','CD','OE1'),('GLN','chi1'):('N','CA','CB','CG'),
           ('GLN','chi2'):('CA','CB','CG','CD'),('GLN','chi3'):('CB','CG','CD','OE1'),
           ('ASP','chi1'):('N','CA','CB','CG'),('ASP','chi2'):('CA','CB','CG','OD1'),
           ('ASN','chi1'):('N','CA','CB','CG'),('ASN','chi2'):('CA','CB','CG','OD1')}
    if (res,chi) not in names:raise InvalidArtifact('unsupported role chemistry/chi quartet')
    return role,res,chi,names[res,chi]


def indices(src,kin,role,names):
    v=src['roles'][role];meta=kin.data['source_atom_metadata'];result=[]
    for name in names:
        ids=[i for i,a in enumerate(meta) if a and (a['chain'],a['resnum'],a['insertion_code'],a['atom'])==
             (v['chain'],v['resnum'],v['icode'],name)]
        if len(ids)!=1:raise InvalidArtifact('source quartet atom missing/nonunique: '+role+'/'+name)
        result.append(ids[0])
    bonds={frozenset(b) for b in kin.data['bonds']}
    if any(frozenset(b) not in bonds for b in zip(result,result[1:])):
        raise InvalidArtifact('quartet connectivity absent')
    return result


def template_quartet(src,kin,role,names):
    """A template atom need not be inside its old chopped physical map."""
    from openmm import app, unit
    prep=read_json(verify(src['tasks'][0]['source_preparation']))
    source=prep['source'];pdb=app.PDBFile(str(verify(source)))
    atoms=list(pdb.topology.atoms());positions=np.asarray(pdb.positions.value_in_unit(unit.angstrom))
    r=src['roles'][role];found=[]
    for name in names:
        ids=[a.index for a in atoms if (a.residue.chain.id,int(a.residue.id),
              a.residue.insertionCode.strip(),a.name)==(r['chain'],r['resnum'],r['icode'],name)]
        if len(ids)!=1:raise InvalidArtifact('actual template source quartet atom missing/nonunique')
        found.append(ids[0])
    bonds={frozenset((a.index,b.index)) for a,b in pdb.topology.bonds()}
    if any(frozenset(pair) not in bonds for pair in zip(found,found[1:])):
        raise InvalidArtifact('template source topology lacks quartet connectivity')
    for name,index in zip(names,found):
        existing=[i for i,a in enumerate(kin.data['source_atom_metadata']) if a and
                  (a['chain'],a['resnum'],a['insertion_code'],a['atom'])==
                  (r['chain'],r['resnum'],r['icode'],name)]
        if len(existing)>1 or (existing and not np.allclose(positions[index],kin.positions[existing[0]],atol=1e-12,rtol=0)):
            raise InvalidArtifact('actual template source and archived physical coordinates differ')
    return positions[found], {'source':source,'source_indices':found,'quartet_names':names}


def dihedral(p):
    a,b,c,d=np.asarray(p,dtype=float);axis=c-b;axis/=np.linalg.norm(axis)
    v=a-b;v-=axis*np.dot(v,axis);w=d-c;w-=axis*np.dot(w,axis)
    if min(np.linalg.norm(v),np.linalg.norm(w))<1e-10:raise InvalidArtifact('undefined source dihedral')
    return float(np.arctan2(np.dot(np.cross(axis,v),w),np.dot(v,w)))


def periodic_delta(a,b,period=2*np.pi):return float((a-b+period/2)%period-period/2)


def terminal_period(src,template,kin,tkin,role,res,chi):
    if (res,chi) not in {('ASP','chi2'),('GLU','chi3')}:return 2*np.pi
    names=('CG','OD1','OD2') if res=='ASP' else ('CD','OE1','OE2')
    for s,k in ((src,kin),(template,tkin)):
        f=s['fragments'][role]
        if f.get('formal_charge')!=-1 or f.get('carboxyl_hydrogens')!=0:return 2*np.pi
        # Presence only: these names are not a four-atom torsion chain.
        r=s['roles'][role]
        found=[a['atom'] for a in k.data['source_atom_metadata'] if a and
               (a['chain'],a['resnum'],a['insertion_code'])==(r['chain'],r['resnum'],r['icode'])]
        if any(found.count(name)!=1 for name in names):raise InvalidArtifact('carboxylate symmetry atoms missing')
    return np.pi


def transplant(src,template,template_class):
    task=src['tasks'][0];kin=Kinematics(read_json(verify(task['mapping']))['context'])
    tkin=Kinematics(read_json(verify(template['tasks'][0]['mapping']))['context'])
    row={'template_class':template_class,'template_case_id':template['case_id'],
         'template_manifest':template['source_manifest'],'template_mapping':template['tasks'][0]['mapping'],
         'target_case_id':src['case_id'],'mapped':[],'unmapped':[],
         'status':'unsupported','reason':None,'active_q_radian':None,'full_q':None}
    q=np.zeros(4);replay=[]
    try:
        if src['biological_group']==template['biological_group']:raise InvalidArtifact('same biological template group')
        for column,i in enumerate(task['active_indices']):
            mode=kin.modes[i];role,res,chi,names=role_for_mode(src,kin,mode)
            if role not in template['roles'] or template['roles'][role]['resname']!=res:
                row['unmapped'].append({'role':role,'chi':chi,'target_chemistry':res,
                    'template_chemistry':template['roles'].get(role,{}).get('resname'),'delta_radian':0.,
                    'reason':'chemically_nonidentical_or_absent_role_retained_at_target_q0'})
                continue
            ids=indices(src,kin,role,names);tcoords,tmap=template_quartet(template,tkin,role,names)
            if ids[1:3]!=mode['axis_indices']:raise InvalidArtifact('quartet axis differs from physical mode')
            old=dihedral(kin.positions[ids]);desired=dihedral(tcoords)
            period=terminal_period(src,template,kin,tkin,role,res,chi)
            small=np.zeros(len(kin.modes));small[i]=1e-5
            sign=periodic_delta(dihedral(kin.positions_only(small)[ids]),old)/1e-5
            if abs(abs(sign)-1)>1e-7:raise InvalidArtifact('unexpected dihedral/rotation orientation')
            q[column]=periodic_delta(desired,old,period)/float(np.sign(sign))
            row['mapped'].append({'role':role,'chemistry':res,'chi':chi,'mode_id':mode['id'],
                'quartet_names':names,'target_indices':ids,'template_source_mapping':tmap,
                'target_dihedral_radian':old,'template_dihedral_radian':desired,
                'period_radian':float(period),'rotation_sign':float(np.sign(sign)),
                'delta_radian':float(q[column])})
            replay.append((ids,desired,period))
        full=angular.full_q(task,q);p=kin.positions_only(full)
        row.update(active_q_radian=q.tolist(),full_q=full.tolist(),
                   maximum_angle_radian=float(np.max(abs(q))),
                   **angular.physical_status(kin,task,q))
        errors=[abs(periodic_delta(dihedral(p[idx]),desired,period)) for idx,desired,period in replay]
        row['dihedral_replay_max_error_radian']=max(errors,default=None)
        if not replay:raise InvalidArtifact('no chemically matching selected coordinates')
        if max(errors)>SETTINGS['dihedral_replay_tolerance_radian']:raise InvalidArtifact('transplant fails target dihedral replay')
        if np.max(abs(q))>.8+1e-12:raise InvalidArtifact('exact template start exceeds frozen angular domain')
        if not row['physical_feasible']:raise InvalidArtifact('exact template start exceeds frozen heavy displacement domain')
        if row['maximum_heavy_displacement_A']<SETTINGS['meaningful_coordinate_change_A']:
            raise InvalidArtifact('template repeats original source within coordinate tolerance')
        row['geometry_checks']=angular.final_geometry(kin,task,q,[a[0] for a in xyz(verify(task['xyz']))])
        row['status']='start_available'
    except InvalidArtifact as exc:row['reason']=str(exc)
    return row


def prepare(inputs,agreement,output):
    data,srcs,templates=sources(inputs);out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False)
    rows=[]
    for src in srcs:
        starts=[]
        for cls,cid in TEMPLATES.items():
            if templates[cid]['biological_group']==src['biological_group']:cid=FALLBACKS[cls]
            starts.append(transplant(src,templates[cid],cls))
        # Meaningful distinct starts are required geometrically, not by score.
        if all(s['status']=='start_available' for s in starts):
            task=src['tasks'][0];kin=Kinematics(read_json(verify(task['mapping']))['context'])
            diff=float(np.max(np.linalg.norm(kin.positions_only(starts[0]['full_q'])-
                                           kin.positions_only(starts[1]['full_q']),axis=1)))
            if diff<SETTINGS['meaningful_coordinate_change_A']:
                starts[1].update(status='unsupported',reason='La-template duplicates Ca-template start')
        rows.append({'case_id':src['case_id'],'source_manifest':src['source_manifest'],
                     'source_case':src['source_case'],'base':src['base'],'tasks':src['tasks'],'starts':starts})
    count=sum(s['status']=='start_available' for r in rows for s in r['starts'])
    snapshot=out/'preparation_implementation.py';shutil.copyfile(__file__,snapshot)
    result={'protocol_id':PROTOCOL,'settings':SETTINGS,'inputs':record(inputs),'agreement':record(agreement),
            'template_rule':TEMPLATES,'same_biological_group_replacements':FALLBACKS,
            'template_source_manifests':{k:v['source_manifest'] for k,v in templates.items()},
            'case_denominator':len(rows),'start_denominator':2*len(rows),'available_starts':count,
            'maximum_optimizer_starts':2*count,'maximum_new_GFN2_singlepoints':8*count,
            'rows':rows,'new_molecular_calls':0,'production_changed':False,'implementation':record(snapshot)}
    write_new(out/'PREPARATION.json',result)
    return {k:result[k] for k in ('case_denominator','start_denominator','available_starts','maximum_optimizer_starts','maximum_new_GFN2_singlepoints','new_molecular_calls')}|{'result':record(out/'PREPARATION.json')}


def prepare_searches(prepared,output):
    from adaptive_completion import SETTINGS as opt_settings
    p=read_json(prepared);verify(p['implementation']);verify(p['agreement']);verify(p['inputs'])
    if p['protocol_id']!=PROTOCOL or p['settings']!=SETTINGS:raise InvalidArtifact('start preparation policy differs')
    tasks=[];parent=None
    for row in p['rows']:
        old=read_json(verify(row['source_manifest']))
        if parent is None:parent=old
        if any(old[k]!=parent[k] for k in ('model','software','orca','cpu_python','gpu_python')):
            raise InvalidArtifact('source scientific model/runtime mismatch')
        for start in row['starts']:
            if start['status']!='start_available':continue
            for task in row['tasks']:
                t=copy.deepcopy(task)
                t.update(task_id=task['task_id']+'__template_'+start['template_class'],
                         original_task_id=task['task_id'],original_manifest=row['source_manifest'],start=start)
                tasks.append(t)
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False);impl=out/'implementation';impl.mkdir();pins={}
    for path in Path(__file__).parent.glob('*.py'):
        target=impl/path.name;shutil.copyfile(path,target);pins[path.name]=record(target)
    m={'protocol_id':PROTOCOL,'settings':SETTINGS,'optimizer_settings':opt_settings,
       'prepared':record(prepared),'agreement':p['agreement'],'inputs':p['inputs'],
       'implementation':pins,'tasks':tasks,'case_denominator':p['case_denominator'],
       'start_denominator':p['start_denominator'],'declared_optimizer_starts':len(tasks),
       'maximum_new_GFN2_singlepoints':p['maximum_new_GFN2_singlepoints'],
       'resources':{'GPU_cpus':32,'GPUs':1,'GPU_host_mem_MiB':200000},
       **{k:parent[k] for k in ('model','software','orca','cpu_python','gpu_python','cpu_executable','gpu_executable')},
       'new_DFT_calls':0,'new_GFN2_calls_in_this_adapter':0,'production_changed':False}
    path=out/'manifest.json';write_new(path,m);result=validate(path);write_new(out/'PREFLIGHT.json',result);return result


def validate(manifest):
    from adaptive_completion import SETTINGS as opt_settings
    m=read_json(manifest);p=read_json(verify(m['prepared']))
    if m['protocol_id']!=PROTOCOL or m['settings']!=SETTINGS or m['optimizer_settings']!=opt_settings:
        raise InvalidArtifact('frozen start/optimizer policy differs')
    for pin in m['implementation'].values():verify(pin)
    for key in ('agreement','inputs','software','orca','cpu_executable','gpu_executable'):verify(m[key])
    expect=[(r['case_id'],s['template_class'],z) for r in p['rows'] for s in r['starts']
            if s['status']=='start_available' for z in ('Ca','La')]
    if [(t['case_id'],t['start']['template_class'],t['metal']) for t in m['tasks']]!=expect:
        raise InvalidArtifact('admitted paired start membership differs')
    if len(expect)!=m['declared_optimizer_starts'] or len(expect)>32:raise InvalidArtifact('finite search count differs')
    for t in m['tasks']:
        old=read_json(verify(t['original_manifest']))
        origin=next(x for x in old['tasks'] if x['task_id']==t['original_task_id'])
        clean={k:v for k,v in t.items() if k not in ('original_manifest','original_task_id','start')}
        clean['task_id']=t['original_task_id']
        if clean!=origin:raise InvalidArtifact('target source/state/four-angle space changed')
        start=next(s for r in p['rows'] if r['case_id']==t['case_id'] for s in r['starts']
                   if s['template_class']==t['start']['template_class'])
        if t['start']!=start:raise InvalidArtifact('source-informed start changed')
        kin=Kinematics(read_json(verify(t['mapping']))['context'])
        check=angular.final_geometry(kin,t,np.asarray(start['active_q_radian']),[a[0] for a in xyz(verify(t['xyz']))])
        if check!=start['geometry_checks']:raise InvalidArtifact('prepared start geometry does not replay')
        r=read_json(verify(t['origin_reuse']['proposal_receipt']))
        if r['origin']!=t['origin_reuse']['point']:raise InvalidArtifact('origin reuse does not match real receipt')
        native=read_json(verify(t['origin_reuse']['point']['MACE']));req=read_json(verify(native['request']))
        if (native['status']!='complete' or native['model']!=m['model'] or
            req['charge']!=t['charge'] or req['multiplicity']!=t['multiplicity'] or
            xyz(verify(req['xyz']))!=xyz(verify(t['xyz']))):raise InvalidArtifact('origin model/state/coordinates differ')
    return {'status':'prepared_dry_run_pass','manifest':record(manifest),'cases':m['case_denominator'],
            'starts':len(expect),'new_molecular_calls':0,'new_DFT_calls':0}


def optimize(manifest,task,gpu):
    from adaptive_completion import CompletionProposal, SCALE, SETTINGS as opt_settings
    from scipy.optimize import minimize
    from mace_hybrid import EV_TO_KCAL
    ev=CompletionProposal(manifest,task,gpu);dest=ev.directory/'result.json'
    if dest.exists():
        if read_json(dest)['manifest']!=record(manifest):raise InvalidArtifact('existing search source differs')
        return record(dest)
    began=time.monotonic();trace=[];accepted=[]
    result={'task_id':task['task_id'],'case_id':task['case_id'],'metal':task['metal'],
            'template_class':task['start']['template_class'],'manifest':record(manifest),
            'status':'unavailable','origin':None,'start':None,'proposal':None,
            'optimizer':None,'unconstrained_minimum_claimed':False,'composite_gradient':None}
    def evaluate(q,purpose):
        trial={'active_q_radian':np.asarray(q).tolist(),'purpose':purpose,
               **angular.physical_status(ev.kin,task,q),'MACE_complete':False};trace.append(trial)
        point=ev.evaluate(q,purpose);trial['MACE_complete']=True;return point
    try:
        origin=evaluate(np.zeros(4),'origin_reuse');result['origin']=origin
        start=evaluate(np.asarray(task['start']['active_q_radian']),'actual_template_start');result['start']=start
        def objective(q):
            value=evaluate(q,'SLSQP_trial')
            return ((value['MACE_eV']-start['MACE_eV'])*SCALE,
                    np.asarray(value['gradient_kcal_mol_rad'])/HA_TO_KCAL)
        opt=minimize(objective,np.asarray(task['start']['active_q_radian']),jac=True,method='SLSQP',
                     bounds=[(-.8,.8)]*4,
                     constraints=[{'type':'ineq','fun':lambda q:angular.constraints(ev.kin,task,q)[0],
                                   'jac':lambda q:angular.constraints(ev.kin,task,q)[1]}],
                     callback=lambda q:accepted.append(np.asarray(q).tolist()),
                     options={'maxiter':opt_settings['maxiter'],'ftol':opt_settings['optimizer_ftol']})
        result['optimizer']={'success':bool(opt.success),'message':str(opt.message),'iterations':int(opt.nit),
                             'function_evaluations':int(opt.nfev),'gradient_evaluations':int(opt.njev)}
        result['final_geometry']=angular.final_geometry(ev.kin,task,opt.x,[a[0] for a in ev.atoms])
        final=evaluate(opt.x,'final_candidate');result['proposal']=final
        if not opt.success:raise InvalidArtifact('SLSQP did not converge: '+str(opt.message))
        if final['MACE_eV']-start['MACE_eV']>opt_settings['maximum_final_energy_increase_eV']:
            raise InvalidArtifact('candidate native energy increased from its declared start')
        result.update(status='proposal_available',native_work_from_start_kcal_mol=(final['MACE_eV']-start['MACE_eV'])*EV_TO_KCAL,
                      native_work_from_origin_kcal_mol=(final['MACE_eV']-origin['MACE_eV'])*EV_TO_KCAL)
    except Exception as exc:result.update(reason=str(exc),traceback=traceback.format_exc())
    result.update(trials=trace,accepted_iterations=accepted,requests=ev.requests,
                  wall_seconds=time.monotonic()-began,job_id=os.environ.get('SLURM_JOB_ID'))
    write_new(dest,result);return record(dest)


def execute(manifest):
    from accommodation_nonlinear import WarmGPU
    validate(manifest)
    if (not os.environ.get('SLURM_JOB_ID') or int(os.environ.get('SLURM_CPUS_ON_NODE','0'))!=32 or
        int(os.environ.get('SLURM_MEM_PER_NODE','0'))!=200000):raise InvalidArtifact('existing32CPU/200000MiB/GPU allocation required')
    m=read_json(manifest);root=Path(manifest).parent;lock=(root/'proposal.lock').open('a+')
    fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB);began=time.monotonic();gpu=None;error=None;results=[]
    try:
        gpu=WarmGPU(manifest)
        for t in m['tasks']:
            pin=optimize(manifest,t,gpu);results.append(pin)
            print(json.dumps({'task_id':t['task_id'],'status':read_json(verify(pin))['status']}),flush=True)
    except Exception as exc:error=str(exc)
    finally:
        try:
            if gpu is not None:gpu.close()
        except Exception as exc:error=(error+'; ' if error else '')+str(exc)
        elapsed=time.monotonic()-began
        write_new(root/('execution_'+os.environ['SLURM_JOB_ID']+'.json'),
                  {'manifest':record(manifest),'results':results,'error':error,'wall_seconds':elapsed,
                   'allocated_core_seconds':elapsed*32,'allocated_GPU_seconds':elapsed,
                   'job_id':os.environ['SLURM_JOB_ID'],'GPU_command':gpu.command if gpu else None,
                   'new_DFT_calls':0,'new_GFN2_calls':0})
    if error:raise InvalidArtifact(error)
    return {'status':'searches_terminal','completed_search_receipts':len(results)}


def collect(manifest,output):
    m=read_json(manifest);p=read_json(verify(m['prepared']));rows=[];spec=[]
    for case in p['rows']:
        endpoints=[];candidates=[];tasks=[t for t in m['tasks'] if t['case_id']==case['case_id']]
        for t in tasks:
            path=Path(manifest).parent/'proposals'/t['task_id']/'result.json'
            r=read_json(path) if path.exists() else {'status':'unavailable','reason':'search_not_completed'}
            if path.exists() and r['manifest']!=record(manifest):raise InvalidArtifact('search result source mismatch')
            endpoints.append({'task_id':t['task_id'],'status':r['status'],'result':record(path) if path.exists() else None,
                              'reason':r.get('reason')})
            if r['status']=='proposal_available':
                c=r['proposal'];candidates.append({'id':'template_'+t['start']['template_class']+'_search_'+t['metal'],
                    'full_q':c['full_q'],'coordinate':c['coordinate'],'native_reuse':{t['metal']:c['MACE']}})
        status='prepared' if tasks and all(e['status']=='proposal_available' for e in endpoints) else 'unavailable'
        reason=None if status=='prepared' else ('no_admitted_distinct_start' if not tasks else 'required_search_unavailable')
        rows.append({'case_id':case['case_id'],'starts':case['starts'],'endpoints':endpoints,'status':status,'reason':reason})
        spec.append({'case_id':case['case_id'],'status':status,'reason':reason,'candidates':candidates if status=='prepared' else []})
    result={'protocol_id':PROTOCOL,'manifest':record(manifest),'case_denominator':len(rows),
            'available_sources':sum(r['status']=='prepared' for r in rows),'rows':rows,
            'finite_candidate_specification':{'branch':'structure_informed_starts','inputs':m['inputs'],
                'agreement':m['agreement'],'coordinate_limits':{'maximum_angle_radian':.8,'maximum_heavy_displacement_A':.8},
                'maximum_candidates_per_case':4,'cases':spec},'new_molecular_calls_in_collection':0}
    write_new(output,result)
    return {'result':record(output),'available_sources':result['available_sources'],'case_denominator':len(rows)}


def compare(collection,searches,output):
    from collections import Counter
    from accommodation_folds_compare import decision
    from accommodation_nonlinear import relative_components
    from nikasha_pool import choose_rows
    c=read_json(collection);s=read_json(searches);m=read_json(verify(c['manifest']))
    spec=read_json(verify(m['specification']))
    if spec!=s['finite_candidate_specification'] or spec['branch']!='structure_informed_starts':
        raise InvalidArtifact('scored candidate specification differs from actual searches')
    inputs=read_json(verify(spec['inputs']));reference=read_json(verify(inputs['adaptive_reference']))
    if [r['case_id'] for r in c['cases']]!=[r['case_id'] for r in inputs['cases']]:
        raise InvalidArtifact('eight-source score denominator/order differs')
    rows=[]
    def outcome(label,expected):
        return 'correct' if label==expected+'-supported' else ('wrong' if label.endswith('-supported') else label)
    for source,new in zip(inputs['cases'],c['cases']):
        basecol=read_json(verify(source['pool_collection']));base=next(r for r in basecol['cases'] if r['case_id']==source['case_id'])
        prep=next(r for r in s['rows'] if r['case_id']==source['case_id'])
        if base['old_result']['expected_class']!=source['known_class']:raise InvalidArtifact('known class provenance differs')
        expected=source['known_class'];static=base['old_result']['R0']['composite_R_model_kcal_mol']
        row={'case_id':source['case_id'],'known_class':expected,'label_use':'consumed_development',
             'start_statuses':prep['starts'],'branch_status':new['pool']['status'],'reason':new.get('reason'),
             'static':{'R':static,'released_band_call':decision(static,reference['old_frozen_bands'])},
             'variants':{}}
        if new['pool']['status']=='available':
            if choose_rows(new['matrix'],[a['id'] for a in new['candidates']])!=new['pool']:
                raise InvalidArtifact('saved pool algebra does not replay')
        for variant in ('mathematical','operational'):
            old=base['pool'][variant];newscore=new['pool'][variant]
            oldR=old['composite_R_model_kcal_mol'];newR=newscore['composite_R_model_kcal_mol'] if newscore else None
            bands=reference['variants'][variant]['bands'];oldcall=decision(oldR,bands);newcall=decision(newR,bands)
            v={'adaptive_R':oldR,'alternative_start_R':newR,'delta_R_model_kcal_mol':None if newR is None else newR-oldR,
               'adaptive_own_band_call':oldcall,'alternative_start_adaptive_band_transfer':newcall,
               'alternative_start_released_band_transfer':decision(newR,reference['old_frozen_bands']),
               'adaptive_outcome':outcome(oldcall,expected),'alternative_start_transfer_outcome':outcome(newcall,expected),
               'new_protocol_calibration':None,'endpoint_selected_work_from_adaptive':None,
               'selected_candidates':None,'delta_native_R':None,'delta_solvent_R':None}
            if newscore:
                selected={z:new['pool']['rows'][z][variant+'_candidate'] for z in ('Ca','La')}
                v['selected_candidates']=selected
                v['endpoint_selected_work_from_adaptive']={z:relative_components(new['matrix'][z][selected[z]]['components'],
                     base['matrix'][z][base['pool']['rows'][z][variant+'_candidate']]['components']) for z in ('Ca','La')}
                v['delta_native_R']=newscore['native_R_model_kcal_mol']-old['native_R_model_kcal_mol']
                v['delta_solvent_R']=newscore['solvation_delta_R_kcal_mol']-old['solvation_delta_R_kcal_mol']
                if abs(v['delta_native_R']+v['delta_solvent_R']-v['delta_R_model_kcal_mol'])>1e-7:
                    raise InvalidArtifact('native/solvent score-change algebra differs')
            row['variants'][variant]=v
        rows.append(row)
    result={'protocol_id':PROTOCOL,'collection':record(collection),'search_collection':record(searches),
            'inputs':spec['inputs'],'adaptive_reference':inputs['adaptive_reference'],'rows':rows,
            'denominator':len(rows),'available':sum(r['branch_status']=='available' for r in rows),
            'new_calibration':None,'production_changed':False,'new_molecular_calls':0,
            'interpretation':'developmental old-band transfer; no hard-case proposal when source start unsupported',
            'summaries':{v:dict(Counter(r['variants'][v]['alternative_start_transfer_outcome'] for r in rows))
                         for v in ('mathematical','operational')},'implementation':record(__file__)}
    write_new(output,result);return result


def scf_details(text):
    result={}
    for key,label in [('energy','Energy'),('max_density','MAX-Density'),('rms_density','RMS-Density')]:
        matches=re.findall(r'Last '+label+r' change\s+\.\.\.\s+([-+0-9.eE]+)\s+Tolerance\s*:\s*([-+0-9.eE]+)',text)
        if len(matches)!=1:raise InvalidArtifact('missing/nonunique printed convergence diagnostic: '+label)
        actual,tol=map(float,matches[0]);result[key]={'residual':actual,'tolerance':tol,'within_printed_tolerance':abs(actual)<=tol}
    cycles=re.findall(r'SCF CONVERGED AFTER\s+(\d+) CYCLES',text)
    nel=re.findall(r'Number of Electrons\s+NEL\s+\.{4}\s+(\d+)',text)
    mult=re.findall(r'Multiplicity\s+Mult\s+\.{4}\s+(\d+)',text)
    charge=re.findall(r'Total Charge\s+Charge\s+\.{4}\s+(-?\d+)',text)
    if any(len(x)!=1 for x in (cycles,nel,mult,charge)):raise InvalidArtifact('SCF state/convergence audit unavailable')
    block=text.split('ORBITAL ENERGIES')[-1]
    occupations=[float(x) for x in re.findall(r'^\s*\d+\s+([0-9]+\.[0-9]+)\s+[-+0-9.]+\s+[-+0-9.]+\s*$',block,re.M)]
    result.update(cycles=int(cycles[0]),electrons=int(nel[0]),multiplicity=int(mult[0]),charge=int(charge[0]),
                  native_mixer_observed='INFO: Using special xTB SCF mixer' in text,
                  ordinary_restart_warning='MOInp will be ignored' in text,
                  printed_orbital_occupations=occupations,
                  printed_fractional_orbital_count=sum(1e-6<o<2-1e-6 for o in occupations))
    return result


def prepare_numeric(collection,agreement,output):
    from affordable_workflow import dry_run
    c=read_json(collection);m=read_json(verify(c['manifest']))
    row=next(r for r in c['cases'] if r['case_id']=='q88jh5-pqq-la_model')
    old=read_json(verify(row['base_collection']));oldrow=next(r for r in old['cases'] if r['case_id']==row['case_id'])
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False);tasks=[]
    for name,case,candidate in [('old_adaptive',oldrow,'adaptive_La'),('new_template',row,'template_La_search_La')]:
        for medium in ('vacuum','alpb'):
            endpoint=case['matrix']['La'][candidate]['low'][medium]
            om=read_json(verify(endpoint['manifest']));source=next(t for t in om['tasks'] if t['task_id']==endpoint['task_id'])
            text=verify(source['input']).read_text()
            if ('UseXTBMixer true' not in text or 'NoAutostart' not in text or
                'MaxIter 500' not in text or source['charge']!=-2 or source['multiplicity']!=1):
                raise InvalidArtifact('unexpected Q88 primary state/input')
            for precision in ('primary_repeat_1','primary_repeat_2'):
                tid='q88jh5__La__'+name+'__'+medium+'__'+precision
                d=out/'tasks'/tid;d.mkdir(parents=True);xp=d/'core.xyz';shutil.copyfile(verify(source['xyz']),xp)
                ip=d/'endpoint.inp';ip.write_text(text)
                tasks.append({'task_id':tid,'case_id':row['case_id'],'case':row['case_id'],'metal':'La',
                              'geometry':name,'candidate':candidate,'medium':medium,'precision':precision,
                              'charge':-2,'multiplicity':1,'xyz':record(xp),'input':record(ip),'output_path':str(d/'endpoint.out'),
                              'original_task':source,'original_endpoint':endpoint,'restart':None})
    impl=out/'implementation';impl.mkdir();pins={}
    for path in Path(__file__).parent.glob('*.py'):
        dst=impl/path.name;shutil.copyfile(path,dst);pins[path.name]=record(dst)
    result={'protocol_id':PROTOCOL+'_Q88_exact_primary_repeats_v2','collection':record(collection),'agreement':record(agreement),
            'tasks':tasks,'orca':m['orca'],'implementation':pins,
            'execution_resources':{'mpi_ranks':8,'concurrent_tasks':8},
            'execution_policy':{'task_runner':pins['run_orca_task_manifest.py'],'runtime_renderer':pins['render_orca_runtime_input.py']},
            'new_DFT_MACE_calls':0,'restart':None,'endpoint_transfer_tolerance_kcal_mol':.1,
            'contrast_diagnostic_tolerance_kcal_mol':.2,'production_changed':False,
            'native_tightening_status':'unavailable_not_demonstrated_under_native_mixer_override'}
    path=out/'manifest.json';write_new(path,result)
    numerical_validate(path);check=dry_run(path);write_new(out/'PREFLIGHT.json',check);return check


def numerical_validate(manifest):
    m=read_json(manifest)
    if m['protocol_id']!=PROTOCOL+'_Q88_exact_primary_repeats_v2':raise InvalidArtifact('numerical protocol differs')
    if len(m['tasks'])!=8 or {(t['geometry'],t['medium'],t['precision']) for t in m['tasks']}!={
        (g,s,p) for g in ('old_adaptive','new_template') for s in ('vacuum','alpb') for p in ('primary_repeat_1','primary_repeat_2')}:
        raise InvalidArtifact('eight numerical states differ')
    verify(m['agreement']);verify(m['orca']);verify(m['collection'])
    for pin in m['implementation'].values():verify(pin)
    for t in m['tasks']:
        source=t['original_task'];text=verify(source['input']).read_text()
        expected=text
        if (verify(t['input']).read_text()!=expected or verify(t['xyz']).read_bytes()!=verify(source['xyz']).read_bytes() or
            t['charge']!=-2 or t['multiplicity']!=1 or t['restart'] is not None):
            raise InvalidArtifact('exact numerical source/state/recipe changed')
    return {'status':'validated','tasks':8,'molecular_calls':0}


def numerical_collect(manifest,output):
    from compact_solvation import completed, diagnostics
    numerical_validate(manifest);m=read_json(manifest);rows=[]
    for t in m['tasks']:
        r={'task_id':t['task_id'],'geometry':t['geometry'],'medium':t['medium'],'precision':t['precision'],
           'status':'unavailable','energy_hartree':None,'source':t['original_endpoint'],'actual':None}
        pin=completed(manifest,t['task_id'])
        if pin:
            try:
                details=scf_details(verify(pin['output']).read_text());audit=diagnostics(pin,t)
                r.update(actual=pin,details=details,audit=audit)
                if (not details['native_mixer_observed'] or details['electrons']!=438 or details['charge']!=-2 or
                    details['multiplicity']!=1 or audit['charge_sanity_status']!='pass'):
                    raise InvalidArtifact('native precision/state not confirmed')
                r.update(status='complete',energy_hartree=pin['energy_hartree'],
                         difference_from_original_kcal_mol=(pin['energy_hartree']-t['original_endpoint']['energy_hartree'])*HA_TO_KCAL)
            except Exception as exc:r['reason']=str(exc)
        rows.append(r)
    transfers={}
    for geometry in ('old_adaptive','new_template'):
        transfers[geometry]={}
        for precision in ('primary_repeat_1','primary_repeat_2'):
            pair={r['medium']:r for r in rows if r['geometry']==geometry and r['precision']==precision}
            transfers[geometry][precision]=((pair['alpb']['energy_hartree']-pair['vacuum']['energy_hartree'])*HA_TO_KCAL
                if all(r['status']=='complete' for r in pair.values()) else None)
    result={'protocol_id':m['protocol_id'],'manifest':record(manifest),'rows':rows,'denominator':8,
            'complete':sum(r['status']=='complete' for r in rows),'La_transfer_kcal_mol':transfers,
            'full_Ca_La_score':None,'old_results_changed':False,'new_molecular_calls_in_collection':0,
            'native_tightening_status':m['native_tightening_status']}
    write_new(output,result);return result


def main():
    p=argparse.ArgumentParser();s=p.add_subparsers(dest='operation',required=True)
    q=s.add_parser('prepare');q.add_argument('--inputs',required=True);q.add_argument('--agreement',required=True);q.add_argument('--output',required=True)
    q=s.add_parser('prepare-searches');q.add_argument('--prepared',required=True);q.add_argument('--output',required=True)
    for op in ('validate','execute'):
        q=s.add_parser(op);q.add_argument('--manifest',required=True)
    q=s.add_parser('collect');q.add_argument('--manifest',required=True);q.add_argument('--output',required=True)
    q=s.add_parser('compare');q.add_argument('--collection',required=True);q.add_argument('--searches',required=True);q.add_argument('--output',required=True)
    q=s.add_parser('prepare-numeric');q.add_argument('--collection',required=True);q.add_argument('--agreement',required=True);q.add_argument('--output',required=True)
    q=s.add_parser('collect-numeric');q.add_argument('--manifest',required=True);q.add_argument('--output',required=True)
    a=p.parse_args()
    if a.operation=='prepare':r=prepare(a.inputs,a.agreement,a.output)
    elif a.operation=='prepare-searches':r=prepare_searches(a.prepared,a.output)
    elif a.operation=='validate':r=validate(a.manifest)
    elif a.operation=='execute':r=execute(a.manifest)
    elif a.operation=='collect':r=collect(a.manifest,a.output)
    elif a.operation=='compare':r=compare(a.collection,a.searches,a.output)
    elif a.operation=='prepare-numeric':r=prepare_numeric(a.collection,a.agreement,a.output)
    else:r=numerical_collect(a.manifest,a.output)
    print(json.dumps(r,indent=2))

if __name__=='__main__':main()
