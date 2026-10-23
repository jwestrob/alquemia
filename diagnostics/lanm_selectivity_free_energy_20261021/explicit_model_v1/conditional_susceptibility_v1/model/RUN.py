#!/usr/bin/env python3
"""Original-baseline conditional susceptibility: three controls, sixteen fields."""
import argparse
import copy
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import resource
import shutil
import time

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[4]
BOHR=.529177210544
HARTREE=627.509474063
ELECTRIC=332.0637133
ELEFIELD=1439.96455
DEBYE=4.803206802
FIELD_FACTOR=HARTREE*ELEFIELD/(ELECTRIC*BOHR)
FIELDS={'minus_h':-.001,'minus_half':-.0005,'plus_half':.0005,'plus_h':.001}
STATES=tuple(f'{metal}_{kind}' for metal in ('La','Dy') for kind in ('dimer','trimer'))
CONTROL_IDS=('zero_Dy_trimer','ion_La_plus_h','ion_La_minus_h')
THREADS=('OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS')
MAX_CALLS=19

def require(ok,msg):
    if not ok:raise ValueError(msg)
def pin(path):
    p=Path(path).resolve();return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def verify(p):require(pin(p['path'])==p,'Changed pin: '+p['path']);return Path(p['path'])
def read(p):return json.loads(verify(p).read_text())
def save(path,value):
    p=Path(path);tmp=p.with_name(p.name+'.tmp');tmp.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n');os.replace(tmp,p)
def module(p,name):
    s=importlib.util.spec_from_file_location(name,verify(p));m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
def helpers(s):
    return module(s['moment_adapter'],'field_moments'),module(s['induction_helper'],'field_induction'),module(s['native_helper'],'field_native')
def same_input(a,b):return all(verify(a[k]).read_bytes()==verify(b[k]).read_bytes() for k in ('xyz','key','parameters'))
def inventory():
    return [('zero_Dy_trimer','Dy_trimer',0.),('ion_La_plus_h','ion_La',.001),('ion_La_minus_h','ion_La',-.001)]+[(label+'__'+name,label,value) for label in STATES for name,value in FIELDS.items()]
def dot(a,b):return math.fsum(x*y for x,y in zip(a,b))
def unit(v):
    n=math.sqrt(dot(v,v));require(n>1e-12,'Degenerate permanent frame');return [x/n for x in v]

def permanent_moment(row):
    """Native kmpole/rotpole dipoles in the fixed lab frame, restricted to these frames."""
    out=[0.,0.,0.]
    for a in row['frame_signatures']:
        p=a['position_A'];c=a['permanent_record']['components'];q=c[0];mu=c[1:4]
        if a['frame_kind']=='None':require(mu==[0.,0.,0.],'Nonzero unframed dipole');v=[0.,0.,0.]
        else:
            require(a['frame_kind'] in ('Bisector','Z-then-X'),'Unsupported local frame')
            axes=[[x-y for x,y in zip(b,p)] for b in a['axis_positions_A']]
            require(len(axes)==2,'Frame denominator differs');z=unit(axes[0]);x=axes[1]
            if a['frame_kind']=='Bisector':
                x=unit(x);z=unit([u+v for u,v in zip(z,x)])
            proj=dot(x,z);x=unit([u-proj*v for u,v in zip(x,z)])
            y=[z[1]*x[2]-z[2]*x[1],z[2]*x[0]-z[0]*x[2],z[0]*x[1]-z[1]*x[0]]
            v=[BOHR*(mu[0]*x[j]+mu[1]*y[j]+mu[2]*z[j]) for j in range(3)]
        for j in range(3):out[j]+=DEBYE*(q*p[j]+v[j])
    return out

def complete_moment(parsed,row):
    p=row['permanent_lab_Debye'];total=[p[j]+parsed['induced_sum_Debye'][j] for j in range(3)]
    bound=parsed['induced_print_bound_Debye_per_component']+2e-9
    require(max(abs(total[j]-parsed['total_lab_Debye'][j]) for j in range(3))<=.0005+bound,'Permanent/induced versus native M mismatch')
    return dict(parsed,permanent_lab_Debye=p,total_precise_lab_Debye=total,total_precise_lab_au=[x/(DEBYE*BOHR) for x in total],
                total_precise_print_bound_au_per_component=bound/(DEBYE*BOHR))

def render(key,field):
    require(math.isfinite(field) and field in (*FIELDS.values(),0.),'Undeclared field')
    forbidden={'external-field','exfld-freq','polar-pair-response','exchpol','exchange-polar','a-axis','b-axis','c-axis','ewald','solvate','solvation','active','inactive','charge-flux','chgflux'}
    tokens=[line.split('#',1)[0].split() for line in key.splitlines()]
    require(not any(a and a[0].lower() in forbidden for a in tokens),'Nonbaseline or preexisting field key')
    require(sum(bool(a) and a[0].lower()=='save-uinduce' for a in tokens)==1 and sum(bool(a) and a[0].lower()=='polar-print' for a in tokens)==1,'Saved moment reporting differs')
    return key+f'\n# Fixed conditional susceptibility; original baseline\nEXTERNAL-FIELD 0 0 {field*FIELD_FACTOR:.17g}\n'

def source_rows(s,h):
    moment,ind,_=h;mm=read(s['moment_manifest']);mr=read(s['moment_result']);tm=read(s['trimer_manifest']);tr=read(s['trimer_result']);shared=read(s['shared_inventory'])
    require(mr['status']=='zero_field_model_redistribution_complete' and mr['error'] is None and tr['status']=='complete' and tr['error'] is None,'Baseline results unavailable')
    mv=read(s['moment_review']);tv=read(s['trimer_review'])
    require(mv['input_result']==s['moment_result'] and mv['status']=='passed_native_actual_review' and tv['model_result']==s['trimer_result'] and tv['status']=='passed_saved_model_application_admission','Actual baseline review joins differ')
    require(mm['candidate_result']==tm['candidate_result'] and mm['candidate_result']['sha256']=='5001cb5798571245d5dcff7c6afbcc8c482e1525f0e8e09e2c4c1b124eb03832','Candidate differs')
    require(mm['selected_parameters']==tm['selected_parameters'],'Selected parameters differ')
    require(mr['config']['analyze']==tm['model_analyze']==s['analyze'],'Original binary differs')
    rows={}
    for label in (*STATES,'ion_La'):
        metal=label.split('_')[0]
        if label=='ion_La':cid='ion_La';row=copy.deepcopy(mm['rows'][cid]);record=mr['rows'][cid];old=record['moments'];producer=s['moment_result']
        elif label.endswith('dimer'):cid=f'aquo9_central_{metal}__dimer_01';row=copy.deepcopy(mm['rows'][cid]);record=mr['rows'][cid];old=record['moments'];producer=s['moment_result']
        else:cid=f'aquo9_central_{metal}__far__metal_two_water';row=copy.deepcopy(tm['model_rows'][cid]);record=tr['rows'][cid];old=record['parsed_moments'];producer=s['trimer_result']
        require(record['status']=='complete' and same_input(record['runtime_input'],row['input']),'Baseline effective input differs')
        require(len(record['receipts'])==1 and record['receipts'][0]['argv']==[s['analyze']['path'],'system.xyz','EM'] and record['receipts'][0]['returncode']==0,'Baseline native authority differs')
        parent=[line.split() for line in verify(row['source_input']['xyz']).read_text().splitlines()[1:] if line.strip()]
        multipoles=moment.multipoles(verify(row['source_input']['parameters']).read_text(),verify(row['source_input']['key']).read_text())
        require(row['frame_signatures']==[moment.frame_signature(parent,i,multipoles) for i in row['parent_atom_indices']],'Native atom/type/frame signatures differ')
        require([int(parent[i][5]) for i in row['parent_atom_indices']]==row['types'],'Fragment atom types differ')
        parsed=moment.parse_moments(verify(old['analyze_output']),verify(old['induced_output']),row,ind.induction)
        require(parsed==old,'Baseline raw moment replay differs')
        row['permanent_lab_Debye']=permanent_moment(row);row['baseline']=complete_moment(old,row)
        row['source_cell']=cid;row['baseline_authority']={'result':producer,'cell':cid,'receipt':record['receipts'][0]}
        row['energy_replay_kcal_mol']=None
        if label!='ion_La':
            qid=f'pbe0__aquo9__{metal}__'+('water_01__metal_water' if label.endswith('dimer') else 'far__metal_two_water')
            qr=shared['zero_field_records'][qid];qc=qr['cell'];indices=qc['real_parent_atom_indices']
            xyz=[line.split() for line in verify(row['input']['xyz']).read_text().splitlines()[1:] if line.strip()]
            require(row['parent_atom_indices']==indices and [[float(x) for x in a[2:5]] for a in xyz]==[qc['positions_A'][i] for i in indices],'Exact quantum real-atom coordinate join differs')
            row['quantum_zero_field_id']=qid;row['quantum_zero_field_source']=qr['result']
        rows[label]=row
    return rows,mm['selected_parameters']

def prepare(out):
    t=time.monotonic();cpu=time.process_time();s=read(pin(HERE/'SOURCE_INPUTS.json'))
    for p in s.values():verify(p)
    rows,point=source_rows(s,helpers(s));jobs={}
    for jid,label,f in inventory():
        row=rows[label];folder=out/'inputs'/jid;folder.mkdir(parents=True)
        for k,name in (('xyz','system.xyz'),('parameters','amoeba09.prm')):shutil.copyfile(verify(row['input'][k]),folder/name)
        (folder/'system.key').write_text(render(verify(row['input']['key']).read_text(),f))
        jobs[jid]={'id':jid,'state':label,'field_au':[0.,0.,f],'field_MV_cm':[0.,0.,f*FIELD_FACTOR],
                   'input':{k:pin(folder/name) for k,name in (('xyz','system.xyz'),('key','system.key'),('parameters','amoeba09.prm'))}}
    require(len(jobs)==19 and tuple(jobs)[:3]==CONTROL_IDS,'Finite stages differ')
    m={'status':'prepared_unreleased','implementation':pin(__file__),'model_specification':pin(HERE/'MODEL_SPEC.md'),'sources':s,'rows':rows,'jobs':jobs,
       'selected_parameters':point,'controls':list(CONTROL_IDS),'matched_states':list(STATES),'maximum_native_programs':19,'workers':1,'new_native_calls':0,
       'conversion':{'bohr_A':BOHR,'hartree_kcal_mol':HARTREE,'electric':ELECTRIC,'elefield':ELEFIELD,'debye_per_eA':DEBYE,'MV_cm_per_field_au':FIELD_FACTOR},
       'cost':{'wall_seconds':time.monotonic()-t,'CPU_seconds':time.process_time()-cpu,'local_threads':1}}
    save(out/'MANIFEST.json',m);save(out/'CONFIG_RUN_TEMPLATE.json',{'implementation':pin(__file__),'manifest':pin(out/'MANIFEST.json'),'model_specification':pin(HERE/'MODEL_SPEC.md'),'maximum_native_programs':19,'workers':1,'root_release':None})
    return m

def context(cfg):
    require(cfg['implementation']==pin(__file__) and cfg['model_specification']==pin(HERE/'MODEL_SPEC.md') and cfg['maximum_native_programs']==19 and cfg['workers']==1,'Config source/scope differs')
    m=read(cfg['manifest']);require(m['implementation']==cfg['implementation'] and m['model_specification']==cfg['model_specification'] and m['new_native_calls']==0,'Prepared bindings differ')
    for p in m['sources'].values():verify(p)
    require(m['maximum_native_programs']==19 and m['workers']==1 and len(m['jobs'])==19 and list(m['jobs'])[:3]==list(CONTROL_IDS),'Finite inventory differs')
    require([(job['id'],job['state'],job['field_au'][2]) for job in m['jobs'].values()]==inventory(),'Exact nineteen field slots differ')
    for job in m['jobs'].values():
        row=m['rows'][job['state']]
        require(verify(job['input']['key']).read_text()==render(verify(row['input']['key']).read_text(),job['field_au'][2]),'Field key differs')
        for k in ('xyz','parameters'):require(verify(job['input'][k]).read_bytes()==verify(row['input'][k]).read_bytes(),'Physical inputs differ')
        require(job['field_au'][:2]==[0.,0.] and job['field_MV_cm']==[0.,0.,job['field_au'][2]*FIELD_FACTOR],'Field vector differs')
        require(row['permanent_lab_Debye']==permanent_moment(row),'Permanent moment differs')
    return m,helpers(m['sources'])

def release_gate(cfg,m):
    r=read(cfg['root_release'])
    require(r['status']=='approved_original_baseline_conditional_susceptibility' and r['implementation']==cfg['implementation'] and r['manifest']==cfg['manifest'] and r['shared_specification']==m['sources']['shared_specification'] and r['maximum_native_programs']==19 and r['workers']==1 and r['field_forces'] is False and r['fit_authorized'] is False,'Root release differs')
    return r

def parse(output,row,field,h):
    moment,ind,_=h;rr=copy.deepcopy(row);rr['energy_replay_kcal_mol']=None
    if rr['kind']=='ion':rr['kind']='explicit_field_ion' # Replace the zero-field assertion with the analytic control below.
    return complete_moment(moment.parse_moments(output,moment.induced_file(output),rr,ind.induction),row)

def control(job,row,value):
    if job['id']=='zero_Dy_trimer':
        old=row['baseline'];de=value['energy_kcal_mol']-old['energy_kcal_mol'];dm=max(abs(x-y) for a,b in zip(value['induced_Debye'],old['induced_Debye']) for x,y in zip(a,b))
        require(abs(de)<=1e-8 and dm<=1e-6,'Zero-field identity failed')
        return {'passed':True,'energy_difference_kcal_mol':de,'max_induced_difference_Debye':dm,'printed_exact':de==0 and dm==0}
    require(row['types']==[500] and row['COM_A']==[0.,0.,0.] and row['permanent_lab_Debye']==[0.,0.,0.],'Bare ion control definition differs')
    f=job['field_MV_cm'][2]/ELEFIELD;alpha=1.6597;expected_E=-.5*ELECTRIC*alpha*f*f;expected_mu=DEBYE*alpha*f
    key=verify(row['input']['key']).read_text();require(any(a[:3]==['polarize','500','1.659700000000'] or (len(a)>=3 and a[0].lower()=='polarize' and a[1]=='500' and float(a[2])==alpha) for a in (line.split() for line in key.splitlines())),'Frozen ion alpha differs')
    de=value['energy_kcal_mol']-expected_E;du=[value['induced_sum_Debye'][0],value['induced_sum_Debye'][1],value['induced_sum_Debye'][2]-expected_mu]
    require(abs(de)<=1e-7 and max(map(abs,du))<=1e-6,'Bare-ion field sign/scale/energy-factor failed')
    return {'passed':True,'expected_energy_kcal_mol':expected_E,'expected_induced_Debye':[0.,0.,expected_mu],'energy_error_kcal_mol':de,'induced_error_Debye':du}

def state_estimates(label,m,records):
    row=m['rows'][label];e0=row['baseline']['energy_kcal_mol']/HARTREE;mu0=row['baseline']['total_precise_lab_au'];values={}
    for tag,h in [('h',.001),('half',.0005)]:
        p=records[label+'__plus_'+tag]['moments'];n=records[label+'__minus_'+tag]['moments'];ep=p['energy_kcal_mol']/HARTREE;en=n['energy_kcal_mol']/HARTREE
        chiE=-(ep+en-2*e0)/(h*h);chiMu=[(a-b)/(2*h) for a,b in zip(p['total_precise_lab_au'],n['total_precise_lab_au'])];muE=-(ep-en)/(2*h)
        values[tag]={'field_au':h,'chi_E_a0_cubed':chiE,'chi_mu_column_a0_cubed':chiMu,'mu_E_au':muE,'mu0_au':mu0,'mu_E_minus_mu0_au':muE-mu0[2],
                     'even_mu_shift_au':[(a+b)/2-c for a,b,c in zip(p['total_precise_lab_au'],n['total_precise_lab_au'],mu0)],
                     'energy_curvature_print_bound_a0_cubed':2e-8/HARTREE/h**2,'moment_slope_print_bound_a0_cubed':(p['total_precise_print_bound_au_per_component']+n['total_precise_print_bound_au_per_component'])/(2*h)}
    checks={'energy_moment_h':abs(values['h']['chi_E_a0_cubed']-values['h']['chi_mu_column_a0_cubed'][2])<=.5,
            'energy_moment_half':abs(values['half']['chi_E_a0_cubed']-values['half']['chi_mu_column_a0_cubed'][2])<=.5,
            'energy_step':abs(values['h']['chi_E_a0_cubed']-values['half']['chi_E_a0_cubed'])<=.5,
            'moment_step':abs(values['h']['chi_mu_column_a0_cubed'][2]-values['half']['chi_mu_column_a0_cubed'][2])<=.5,
            'odd_h':abs(values['h']['mu_E_minus_mu0_au'])<=1e-4,'odd_half':abs(values['half']['mu_E_minus_mu0_au'])<=1e-4}
    return {'steps':values,'checks':checks,'passed':all(checks.values())}

def summarize(m,records):
    require(set(records)==set(m['jobs']) and all(r['status']=='complete' for r in records.values()),'Incomplete19 denominator')
    states={label:state_estimates(label,m,records) for label in STATES};increments={}
    terms={'La':{'La_trimer':1,'La_dimer':-1},'Dy':{'Dy_trimer':1,'Dy_dimer':-1},'Dy_minus_La':{'Dy_trimer':1,'Dy_dimer':-1,'La_trimer':-1,'La_dimer':1}}
    for name,coeff in terms.items():
        steps={}
        for tag in ('h','half'):
            steps[tag]={key:math.fsum(c*states[label]['steps'][tag][key] for label,c in coeff.items()) for key in ('chi_E_a0_cubed','mu_E_au','mu_E_minus_mu0_au')}
            steps[tag]['chi_mu_column_a0_cubed']=[math.fsum(c*states[label]['steps'][tag]['chi_mu_column_a0_cubed'][j] for label,c in coeff.items()) for j in range(3)]
            for key in ('energy_curvature_print_bound_a0_cubed','moment_slope_print_bound_a0_cubed'):steps[tag][key]=math.fsum(abs(c)*states[label]['steps'][tag][key] for label,c in coeff.items())
        checks={'energy_moment_h':abs(steps['h']['chi_E_a0_cubed']-steps['h']['chi_mu_column_a0_cubed'][2])<=1.,'energy_moment_half':abs(steps['half']['chi_E_a0_cubed']-steps['half']['chi_mu_column_a0_cubed'][2])<=1.,
                'energy_step':abs(steps['h']['chi_E_a0_cubed']-steps['half']['chi_E_a0_cubed'])<=1.,'moment_step':abs(steps['h']['chi_mu_column_a0_cubed'][2]-steps['half']['chi_mu_column_a0_cubed'][2])<=1.}
        increments[name]={'terms':coeff,'steps':steps,'checks':checks,'passed':all(checks.values())}
    return {'states':states,'increments':increments,'numerical_checks_passed':all(v['passed'] for v in (*states.values(),*increments.values())),
            'interpretation':'Conditional increment, not alpha3 or unique mechanistic attribution; consistency gates are not physical error bars.'}

def run(cfg,out):
    tick=time.monotonic();cpu=time.process_time();ch0=resource.getrusage(resource.RUSAGE_CHILDREN);m,h=context(cfg);release=release_gate(cfg,m)
    require(os.environ.get('SLURM_JOB_ID') and os.environ.get('SLURM_CPUS_PER_TASK')=='1' and os.environ.get('SLURM_NTASKS','1')=='1' and all(os.environ.get(k)=='1' for k in THREADS),'One Slurm CPU/task/thread required')
    require(str(out.resolve())==release['execution_root'],'Execution location differs');_,_,native=h
    result={'status':'running','implementation':cfg['implementation'],'config':cfg,'manifest':cfg['manifest'],'release':cfg['root_release'],'job_id':os.environ['SLURM_JOB_ID'],
            'rows':{k:{'status':'not_attempted'} for k in m['jobs']},'native_program_invocations':0,'maximum_native_programs':19,'controls_passed':False,'global_stop':False,'error':None}
    save(out/'RESULT.json',result)
    for jid,job in m['jobs'].items():
        if result['global_stop']:break
        rec={'status':'started','job':job,'receipts':[]};folder=out/'cells'/jid
        try:
            if jid not in CONTROL_IDS:require(result['controls_passed'],'Three initial controls required')
            require(result['native_program_invocations']<19,'Native19 limit reached');folder.mkdir(parents=True,exist_ok=False)
            for k,name in (('xyz','system.xyz'),('key','system.key'),('parameters','amoeba09.prm')):shutil.copyfile(verify(job['input'][k]),folder/name)
            rec['runtime_input']={k:pin(folder/name) for k,name in (('xyz','system.xyz'),('key','system.key'),('parameters','amoeba09.prm'))}
            result['native_program_invocations']+=1;save(out/'CALL_COUNT.json',{'native_program_invocations':result['native_program_invocations'],'maximum':19})
            save(folder/'START.json',dict(rec,implementation=cfg['implementation'],manifest=cfg['manifest'],release=cfg['root_release'],job_id=result['job_id']));rec['start']=pin(folder/'START.json')
            native.invoke([m['sources']['analyze']['path'],'system.xyz','EM'],folder,folder/'analyze.out',rec['receipts'])
            rec['moments']=parse(folder/'analyze.out',m['rows'][job['state']],job['field_au'][2],h)
            if jid in CONTROL_IDS:rec['control']=control(job,m['rows'][job['state']],rec['moments'])
            rec['status']='complete';result['rows'][jid]=rec
            if jid==CONTROL_IDS[-1]:result['controls_passed']=all(result['rows'][k]['control']['passed'] for k in CONTROL_IDS)
        except Exception as exc:result['global_stop']=True;result['error']=repr(exc);rec.update(status='failed',error=repr(exc));result['rows'][jid]=rec
        save(folder/'RESULT.json',rec);save(out/'RESULT.json',result)
    if not result['global_stop']:
        try:
            result['summary']=summarize(m,result['rows']);require(result['summary']['numerical_checks_passed'],'Finite-step consistency failed; cause unresolved')
            result['status']='conditional_model_susceptibility_complete'
        except Exception as exc:result['global_stop']=True;result['error']=repr(exc)
    if result['global_stop']:result['status']='failed_or_numerically_unqualified'
    ch=resource.getrusage(resource.RUSAGE_CHILDREN);result['cost']={'wall_seconds':time.monotonic()-tick,'parent_CPU_seconds':time.process_time()-cpu,'child_user_CPU_seconds':ch.ru_utime-ch0.ru_utime,'child_system_CPU_seconds':ch.ru_stime-ch0.ru_stime,'workers':1}
    save(out/'RESULT.json',result);return result

def admit_result(cfg,m,result):
    bound,h=context(cfg);release_gate(cfg,m);require(bound==m and result['config']==cfg and result['manifest']==cfg['manifest'] and result['release']==cfg['root_release'] and result['implementation']==cfg['implementation'],'Saved binding differs')
    require(result['status']=='conditional_model_susceptibility_complete' and result['error'] is None and not result['global_stop'] and result['native_program_invocations']==19 and result['controls_passed'],'Unavailable model result')
    require(set(result['rows'])==set(m['jobs']),'Full19 slots required')
    for jid,rec in result['rows'].items():
        job=m['jobs'][jid];require(rec['status']=='complete' and rec['job']==job and same_input(rec['runtime_input'],job['input']),'Saved slot/input differs')
        folder=verify(rec['runtime_input']['xyz']).parent;start=read(rec['start'])
        require(start['job']==job and start['runtime_input']==rec['runtime_input'] and start['implementation']==cfg['implementation'] and start['manifest']==cfg['manifest'] and start['release']==cfg['root_release'],'Saved START differs')
        require(len(rec['receipts'])==1,'One analyze per slot required');rr=rec['receipts'][0]
        require(rr['argv']==[m['sources']['analyze']['path'],'system.xyz','EM'] and rr['cwd']==str(folder) and rr['returncode']==0 and all(rr['threads'][k]=='1' for k in THREADS) and rr['output']==pin(folder/'analyze.out'),'Saved native receipt differs')
        value=parse(folder/'analyze.out',m['rows'][job['state']],job['field_au'][2],h);require(value==rec['moments'],'Raw field moment differs')
        if jid in CONTROL_IDS:require(control(job,m['rows'][job['state']],value)==rec['control'],'Saved control differs')
    require(summarize(m,result['rows'])==result['summary'],'Saved curvature arithmetic differs')
    return {'status':'saved_conditional_model_admitted','rows':result['rows'],'summary':result['summary'],'new_native_calls':0}

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=('prepare','run'));p.add_argument('--out',required=True,type=Path);p.add_argument('--config',type=Path);a=p.parse_args();a.out.mkdir(parents=True,exist_ok=False)
    value=prepare(a.out) if a.action=='prepare' else run(json.loads(a.config.read_text()),a.out)
    print(json.dumps({'status':value['status'],'out':str(a.out)}));raise SystemExit(0 if value['status'] in ('prepared_unreleased','conditional_model_susceptibility_complete') else 1)
