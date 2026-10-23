#!/usr/bin/env python3
"""Four fixed consumed D/M sensitivities; no fit, new reference or FD."""
import argparse
import concurrent.futures
import copy
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import resource
import shutil
import threading
import time

HERE=Path(__file__).resolve().parent
PROFILES={'D_plus':(.05,0.),'D_minus':(-.05,0.),'M_plus':(0.,.05),'M_minus':(0.,-.05)}
CENTERS=tuple(f'aquo{n}_central_{m}' for n in (8,9) for m in ('La','Dy'))
TRAIN=tuple(q+'_'+m for q in ('water_minus','water_central','water_plus','aquo8_central','aquo9_central') for m in ('La','Dy'))
MAX_CALLS=236
THREADS=('OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS')

def require(ok,msg):
    if not ok: raise ValueError(msg)
def pin(path):
    p=Path(path).resolve();return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def verify(p):
    require(pin(p['path'])==p,'Changed pin: '+p['path']);return Path(p['path'])
def read(p):return json.loads(verify(p).read_text())
def save(path,x):
    p=Path(path);tmp=p.with_name(p.name+'.tmp');tmp.write_text(json.dumps(x,indent=2,allow_nan=False)+'\n');os.replace(tmp,p)
def module(p,name):
    s=importlib.util.spec_from_file_location(name,verify(p));m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
def same_input(a,b):return all(verify(a[k]).read_bytes()==verify(b[k]).read_bytes() for k in ('xyz','key','parameters'))
def profile_key(key,profile):
    require('polar-pair-response' not in key.lower(),'Unexpected preexisting response')
    dd,dm=PROFILES[profile]
    return key+''.join(f'\nPOLAR-PAIR-RESPONSE {metal} 36 {dd:.17g} {dm:.17g}\n' for metal in (500,501))
def norm(x):return math.sqrt(math.fsum(v*v for v in x))

def helpers(s):
    q=module(s['qualification_adapter'],'sensitivity_qualified');q.PROFILES=dict(PROFILES)
    return q,module(s['moment_adapter'],'sensitivity_moments'),module(s['native'],'sensitivity_native'),module(s['gradient'],'sensitivity_gradient'),module(s['induction'],'sensitivity_induction')

def qualification(s):
    r=read(s['qualification_result']);v=read(s['qualification_review']);c=read(s['qualification_config'])
    require(r['status']=='metal_pair_numerical_qualification_complete' and r['numerical_checks_passed'] is True and r['error'] is None and r['native_program_invocations']==56 and r['internal_FD_energies_completed']==384,'Qualification unavailable')
    require(v['status']=='actual_numerical_qualification_review_passed' and v['result']==s['qualification_result'] and v['config']==s['qualification_config'] and v['new_native_calls']==0 and v['error'] is None,'Independent numerical review unavailable')
    require(r['config']==s['qualification_config'] and c['implementation']==s['qualification_adapter'],'Qualification source join')
    for k in ('prototype_analyze','prototype_testgrad','prototype_source_manifest','source_audit'):
        require(c[k]==s[k],'Qualified binary/source changed');verify(s[k])
    return r,c

def reuse_gradient(s,record,inputs,row,profile,h):
    q,m,n,g,ind=h
    require(record['status']=='complete' and same_input(record['input'],inputs) and len(record['receipts'])==1,'Qualified gradient input differs')
    rec=record['receipts'][0];out=verify(rec['output'])
    require(rec['argv']==[s['prototype_testgrad']['path'],'system.xyz','Y','Y','0.00001'] and rec['returncode']==0 and all(v=='1' for v in rec['threads'].values()),'Qualified gradient invocation differs')
    value=g.gradient(out,row['atoms'],True,n)
    require(value==record['value'] and value['fd_passed'],'Qualified raw gradient replay differs')
    q.response_readback(out.read_text(),profile,row['types'],verify(inputs['key']).read_text());ind.induction(out)
    return {'producer_result':s['qualification_result'],'producer_record':record['job']['id'],'input':record['input'],'receipt':rec,'value':value,'new_native_calls':0,'historical_internal_FD_energies':6*row['atoms']}

def prepare(out):
    tick=time.monotonic();cpu=time.process_time();s=read(pin(HERE/'SOURCE_INPUTS.json'))
    for p in s.values():verify(p)
    qual,qcfg=qualification(s);h=helpers(s);q,m,n,g,ind=h
    data=read(s['training_data']);tm=read(s['trimer_manifest']);mm=read(s['moment_manifest'])
    mr=read(s['moment_result']);tr=read(s['trimer_model_result']);refs=read(s['trimer_comparison']);cr=read(s['center_comparison'])
    require(data['provenance']['candidate_result']['sha256']=='5001cb5798571245d5dcff7c6afbcc8c482e1525f0e8e09e2c4c1b124eb03832' and data['baseline_parameters']==tm['selected_parameters']==mm['selected_parameters'],'Frozen selected point differs')
    require(mr['status']=='zero_field_model_redistribution_complete' and mr['error'] is None and tr['status']=='complete' and tr['error'] is None,'Baseline moment results unavailable')
    require(refs['status']=='complete_saved_trimer_comparison' and len(refs['cycles'])==8 and len(tm['cycles'])==8,'Trimer reference denominator')
    trimers=tuple(k for k,v in tm['model_rows'].items() if v['kind']=='metal_two_water')
    dimers=tuple(sorted({k for v in tm['cycles'].values() for k in v['model_terms'] if '__dimer_' in k}))
    cells=TRAIN+trimers+dimers
    require(len(trimers)==8 and len(dimers)==12 and len(set(cells))==30,'Finite30-cell denominator')
    basemom={k:copy.deepcopy(v['moments']) for k,v in mr['rows'].items()}
    basemom.update({k:copy.deepcopy(v['parsed_moments']) for k,v in tr['rows'].items()})
    static=set(k for v in tm['cycles'].values() for k in v['model_terms'] if k not in cells)
    static.update(k for center in CENTERS for k in m.terms(center,'C') if k not in cells)
    require(len(static)==26,'Unchanged fragment denominator')
    keepmom=set(cells)&set(basemom)|static
    basemom={k:basemom[k] for k in sorted(keepmom)}
    for value in basemom.values():
        verify(value['analyze_output']);verify(value['induced_output'])
    rows={};baseline_force={};resolver=module(s['resolver'],'sensitivity_resolver')
    for cid in cells:
        if cid in mm['rows']:row=copy.deepcopy(mm['rows'][cid])
        elif cid in tm['model_rows']:row=copy.deepcopy(tm['model_rows'][cid])
        else:
            value=data['baseline_model'][cid];src=value['provenance']['effective_input'];header,atoms,coords=q.xyz(verify(src['xyz']))
            prm=verify(src['parameters']).read_text();key=verify(src['key']).read_text();resolved=resolver.parse_parameters(prm,key)
            masses=[resolved['atom'][int(a[5])]['mass'] for a in atoms]
            row={'id':cid,'kind':'published_dimer','metal':cid.rsplit('_',1)[1],'atoms':len(atoms),'types':[int(a[5]) for a in atoms],
                 'names':[a[1] for a in atoms],'native_names':[resolved['atom'][int(a[5])]['symbol'][:3] for a in atoms],
                 'masses_amu':masses,'COM_A':[math.fsum(mass*p[j] for mass,p in zip(masses,coords))/math.fsum(masses) for j in range(3)],
                 'charge':3,'source_input':src,'input':src,'frame_signatures':[m.frame_signature(atoms,i,m.multipoles(prm,key)) for i in range(len(atoms))]}
        folder=out/'base_inputs'/cid;folder.mkdir(parents=True)
        for k,filename in (('xyz','system.xyz'),('key','system.key'),('parameters','amoeba09.prm')):shutil.copyfile(verify(row['input'][k]),folder/filename)
        if cid.startswith('water_'):(folder/'system.key').write_text(m.render((folder/'system.key').read_text(),(folder/'amoeba09.prm').read_text()))
        row['input']={k:pin(folder/f) for k,f in (('xyz','system.xyz'),('key','system.key'),('parameters','amoeba09.prm'))}
        row['energy_replay_kcal_mol']=None
        require('polar-print' in (folder/'system.key').read_text().lower() and 'save-uinduce' in (folder/'system.key').read_text().lower(),'Moment reporting absent')
        _,atoms,coords=q.xyz(folder/'system.xyz');require([int(a[5]) for a in atoms]==row['types'] and len(atoms)==row['atoms'],'Copied atom/type inventory')
        if cid in TRAIN:
            require(coords==data['rows'][cid]['geometry']['positions_A'],'Training geometry changed')
            b=data['baseline_model'][cid];verify(b['native']['gradient_output'])
            baseline_force[cid]={'forces':b['forces'],'energy':b['energy'],'output':b['native']['gradient_output'],'authority':b['provenance']['native_authority']}
        elif 'original__'+cid in qual['rows']:
            b=qual['rows']['original__'+cid];verify(b['receipts'][0]['output'])
            baseline_force[cid]={'forces':b['value']['forces'],'energy':b['value']['energy'],'output':b['receipts'][0]['output'],'authority':s['qualification_result']}
        else:baseline_force[cid]=None
        if cid in basemom and cid in TRAIN:require(basemom[cid]['energy_kcal_mol']==data['baseline_model'][cid]['energy'],'Baseline energy join differs')
        rows[cid]=row
    require(sum(v is not None for v in baseline_force.values())==16,'Baseline force availability differs')
    jobs={};reused={}
    for profile in PROFILES:
        for cid,row in rows.items():
            label=profile+'__'+cid;folder=out/'inputs'/label;folder.mkdir(parents=True)
            for k,fn in (('xyz','system.xyz'),('parameters','amoeba09.prm')):shutil.copyfile(verify(row['input'][k]),folder/fn)
            (folder/'system.key').write_text(profile_key(verify(row['input']['key']).read_text(),profile))
            inp={k:pin(folder/f) for k,f in (('xyz','system.xyz'),('key','system.key'),('parameters','amoeba09.prm'))}
            jobs[label]={'id':label,'profile':profile,'cell':cid,'input':inp}
            old_profile={'D_plus':'direct','M_minus':'mutual'}.get(profile)
            key=(old_profile+'__'+cid) if old_profile else ''
            if key in qual['rows']:reused[label]=reuse_gradient(s,qual['rows'][key],inp,row,profile,h)
    require(len(jobs)==120 and len(reused)==4 and 2*len(jobs)-len(reused)==MAX_CALLS,'Finite program ceiling/reuse differs')
    objective=module(s['objective'],'sensitivity_objective');body=module(s['body'],'sensitivity_body')
    baseline_score=objective.evaluate({'rows':data['rows']},data['baseline_model'],body)
    require(baseline_score==data['baseline_score'],'Exact original25-control score replay differs')
    manifest={'status':'prepared_unreleased','implementation':pin(__file__),'plan':pin(HERE/'PLAN.md'),'sources':s,'profiles':PROFILES,'rows':rows,'jobs':jobs,'reused_gradients':reused,
              'baseline_parameters':data['baseline_parameters'],'training_rows':data['rows'],'baseline_model':data['baseline_model'],'baseline_score':baseline_score,'baseline_forces':baseline_force,
              'baseline_moments':basemom,'unchanged_moment_ids':sorted(static),'trimer_cycles':tm['cycles'],'trimer_reference':refs['cycles'],
              'center_reference_C':{k:cr['parents'][k]['C'] for k in CENTERS},'center_cancellation':{k:mm['permanent_cancellation'][k]['C'] for k in CENTERS},
              'maximum_native_programs':MAX_CALLS,'new_analyze_programs':120,'new_gradient_programs':116,'internal_FD_energies':0,'workers':4,'new_native_calls':0,
              'cost':{'wall_seconds':time.monotonic()-tick,'CPU_seconds':time.process_time()-cpu,'local_threads':1}}
    save(out/'MANIFEST.json',manifest)
    cfg={'implementation':pin(__file__),'manifest':pin(out/'MANIFEST.json'),'plan':pin(HERE/'PLAN.md'),'root_release':None,'maximum_native_programs':MAX_CALLS,'workers':4,'internal_FD_energies':0}
    save(out/'CONFIG_RUN_TEMPLATE.json',cfg)
    return manifest

def summarize(m,records,h):
    objective=module(m['sources']['objective'],'sensitivity_final_objective');body=module(m['sources']['body'],'sensitivity_final_body')
    require(set(records)==set(m['jobs']) and all(r['status']=='complete' for r in records.values()),'Complete120-cell denominator required')
    result={}
    for profile in PROFILES:
        values={cid:records[profile+'__'+cid] for cid in m['rows']}
        score=objective.evaluate({'rows':m['training_rows']},{cid:{'energy':values[cid]['moments']['energy_kcal_mol'],'forces':values[cid]['gradient']['forces']} for cid in TRAIN},body)
        moments=dict(m['baseline_moments']);moments.update({cid:r['moments'] for cid,r in values.items()})
        cycles={}
        for cid,cycle in m['trimer_cycles'].items():
            terms=cycle['model_terms'];ref=m['trimer_reference'][cid]
            energy=math.fsum(c*moments[k]['energy_kcal_mol'] for k,c in terms.items())
            mu=[math.fsum(c*moments[k]['induced_sum_Debye'][j] for k,c in terms.items()) for j in range(3)]
            low=[math.fsum(c*moments[k]['total_lab_Debye'][j] for k,c in terms.items()) for j in range(3)]
            allowance=cycle['M_print_bound_Debye_per_component']+cycle['induced_print_bound_Debye_per_component']+math.fsum(abs(c)*2e-9 for c in terms.values())
            require(max(abs(a-b) for a,b in zip(mu,low))<=allowance,'Trimer permanent/frame/origin consistency failed')
            cycles[cid]={'model_E3_kcal_mol':energy,'baseline_E3_kcal_mol':ref['model_E3_kcal_mol'],'reference_electronic_E3_kcal_mol':ref['reference_electronic_E3_kcal_mol'],
                         'model_minus_reference_E3_kcal_mol':energy-ref['reference_electronic_E3_kcal_mol'],'model_mu3_Debye':mu,'baseline_mu3_Debye':ref['model_mu3_Debye'],
                         'reference_mu3_Debye':ref['reference_mu3_Debye'],'model_minus_reference_mu3_Debye':[a-b for a,b in zip(mu,ref['reference_mu3_Debye'])],
                         'model_mu3_norm_Debye':norm(mu),'corrected_M_mu3_Debye':low,'M_consistency_tolerance_Debye':allowance,
                         'combined_print_bound_Debye_per_component':ref['combined_print_bound_Debye_per_component']}
        centers={}
        for cid in CENTERS:
            proof=m['center_cancellation'][cid];terms=proof['terms'];reference=m['center_reference_C'][cid]
            mu=[math.fsum(c*moments[k]['induced_sum_Debye'][j] for k,c in terms.items()) for j in range(3)]
            centers[cid]={'energy_kcal_mol':moments[cid]['energy_kcal_mol'],'forces':values[cid]['gradient']['forces'],
                          'C_energy_kcal_mol':math.fsum(c*moments[k]['energy_kcal_mol'] for k,c in terms.items()),
                          'C_mu_Debye':mu,'baseline_C_mu_Debye':reference['model_Debye'],'reference_C_mu_Debye':reference['quantum_Debye'],
                          'C_mu_error_Debye':[a-b for a,b in zip(mu,reference['quantum_Debye'])],'combined_print_bound_Debye_per_component':reference['combined_print_bound_Debye_per_component']}
        result[profile]={'increments':PROFILES[profile],'original25_controls':score,'trimer_cycles':cycles,'centers':centers,'full_shell_B_computed':False,'selected_candidate':False}
    return result

class Counter:
    def __init__(self):self.lock=threading.RLock();self.stop=threading.Event();self.calls=0
    def reserve(self):
        with self.lock:
            require(not self.stop.is_set() and self.calls<MAX_CALLS,'Global stop or finite236-call ceiling');self.calls+=1

def context(cfg):
    require(cfg['implementation']==pin(__file__) and cfg['plan']==pin(HERE/'PLAN.md') and cfg['maximum_native_programs']==MAX_CALLS and cfg['workers']==4 and cfg['internal_FD_energies']==0,'Runtime source/scope differs')
    m=read(cfg['manifest']);require(m['implementation']==cfg['implementation'] and m['plan']==cfg['plan'] and m['status']=='prepared_unreleased' and m['new_native_calls']==0 and m['profiles']=={k:list(v) for k,v in PROFILES.items()},'Preparation differs')
    for p in m['sources'].values():verify(p)
    qualification(m['sources']);require(len(m['rows'])==30 and len(m['jobs'])==120 and len(m['reused_gradients'])==4,'Runtime denominator differs')
    for job in m['jobs'].values():
        row=m['rows'][job['cell']]
        require(verify(job['input']['key']).read_text()==profile_key(verify(row['input']['key']).read_text(),job['profile']),'Profile key changed')
        for k in ('xyz','parameters'):require(verify(job['input'][k]).read_bytes()==verify(row['input'][k]).read_bytes(),'Geometry/parameters changed')
    return m,helpers(m['sources'])

def run(cfg,out):
    m,h=context(cfg);q,moment,native,grad,ind=h;s=m['sources'];release=read(cfg['root_release'])
    require(release['status']=='approved_consumed_pair_response_sensitivity' and release['implementation']==cfg['implementation'] and release['manifest']==cfg['manifest'] and release['plan']==cfg['plan'] and release['maximum_native_programs']==MAX_CALLS and release['workers']==4 and release['internal_FD_energies']==0 and release['fit_authorized'] is False,'Root release differs')
    require(read(release['configuration'])==cfg and release['qualification_review']==s['qualification_review'],'Released configuration/review differs')
    require(os.environ.get('SLURM_JOB_ID','').isdigit() and int(os.environ.get('SLURM_NTASKS','1'))==1 and int(os.environ.get('SLURM_CPUS_PER_TASK','0'))==4 and all(os.environ.get(k)=='1' for k in THREADS),'FourCPU single-thread-worker Slurm runtime required')
    tick=time.monotonic();cpu=time.process_time();ch0=resource.getrusage(resource.RUSAGE_CHILDREN);count=Counter()
    result={'status':'running','implementation':cfg['implementation'],'config':cfg,'manifest':cfg['manifest'],'release':cfg['root_release'],'job_id':os.environ['SLURM_JOB_ID'],'rows':{k:{'status':'not_attempted'} for k in m['jobs']},'native_program_invocations':0,'maximum_native_programs':MAX_CALLS,'reused_gradient_programs':4,'internal_FD_energies':0,'global_stop':False,'error':None,'profiles':None}
    def persist():
        with count.lock:result.update(native_program_invocations=count.calls,global_stop=count.stop.is_set());save(out/'RESULT.json',result)
    def worker(label):
        if count.stop.is_set():return
        job=m['jobs'][label];row=m['rows'][job['cell']];folder=out/'cells'/label;r={'id':label,'status':'started','job':job,'receipts':[]}
        try:
            folder.mkdir(parents=True,exist_ok=False)
            for k,fn in (('xyz','system.xyz'),('key','system.key'),('parameters','amoeba09.prm')):shutil.copyfile(verify(job['input'][k]),folder/fn)
            r['runtime_input']={k:pin(folder/f) for k,f in (('xyz','system.xyz'),('key','system.key'),('parameters','amoeba09.prm'))}
            with count.lock:
                count.reserve();save(out/'CALL_COUNT.json',{'native_program_invocations':count.calls,'maximum':MAX_CALLS});save(folder/'START.json',dict(r,implementation=cfg['implementation'],manifest=cfg['manifest'],release=cfg['root_release'],job_id=result['job_id']))
            native.invoke([s['prototype_analyze']['path'],'system.xyz','EM'],folder,folder/'analyze.out',r['receipts'])
            r['moments']=moment.parse_moments(folder/'analyze.out',moment.induced_file(folder/'analyze.out'),row,ind.induction)
            r['analyze_response']=q.response_readback((folder/'analyze.out').read_text(),job['profile'],row['types'],(folder/'system.key').read_text())
            if label in m['reused_gradients']:
                old=m['reused_gradients'][label];rec={'status':'complete','input':old['input'],'receipts':[old['receipt']],'value':old['value'],'job':{'id':old['producer_record']}}
                r['gradient_reuse']=reuse_gradient(s,rec,job['input'],row,job['profile'],h);r['gradient']=old['value']
            else:
                with count.lock:
                    count.reserve();save(out/'CALL_COUNT.json',{'native_program_invocations':count.calls,'maximum':MAX_CALLS})
                native.invoke([s['prototype_testgrad']['path'],'system.xyz','Y','N'],folder,folder/'testgrad.out',r['receipts'])
                r['gradient']=grad.gradient(folder/'testgrad.out',row['atoms'],False,native)
                r['gradient_induction']=ind.induction(folder/'testgrad.out')
                r['gradient_response']=q.response_readback((folder/'testgrad.out').read_text(),job['profile'],row['types'],(folder/'system.key').read_text())
            require(abs(r['gradient']['energy']-r['moments']['energy_kcal_mol'])<=1e-4,'Analyze/gradient energy closure failed')
            r.update(status='complete',start=pin(folder/'START.json'))
        except Exception as exc:
            count.stop.set();r.update(status='failed',error=repr(exc))
        with count.lock:
            result['rows'][label]=r;save(folder/'RESULT.json',r)
            if r['status']=='failed':result['error']=r['error']
            persist()
    persist()
    try:
        with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:list(pool.map(worker,m['jobs']))
        require(not count.stop.is_set() and count.calls==MAX_CALLS,'Incomplete/stopped sensitivity inventory')
        result['profiles']=summarize(m,result['rows'],h);result['status']='consumed_pair_response_sensitivity_complete'
    except Exception as exc:count.stop.set();result.update(status='failed',error=repr(exc))
    ch=resource.getrusage(resource.RUSAGE_CHILDREN);result['cost']={'wall_seconds':time.monotonic()-tick,'parent_CPU_seconds':time.process_time()-cpu,'child_user_CPU_seconds':ch.ru_utime-ch0.ru_utime,'child_system_CPU_seconds':ch.ru_stime-ch0.ru_stime,'workers':4,'native_threads':1};persist();return result

def admit_result(cfg,m,result):
    bound,h=context(cfg);require(bound==m and result['manifest']==cfg['manifest'] and result['config']==cfg,'Result binding differs')
    require(result['status']=='consumed_pair_response_sensitivity_complete' and result['error'] is None and not result['global_stop'] and result['native_program_invocations']==MAX_CALLS and result['internal_FD_energies']==0 and set(result['rows'])==set(m['jobs']),'Result incomplete')
    q,moment,native,grad,ind=h;s=m['sources']
    for label,r in result['rows'].items():
        job=m['jobs'][label];row=m['rows'][job['cell']];require(r['status']=='complete' and r['job']==job and same_input(r['runtime_input'],job['input']),'Raw input/slot differs')
        folder=verify(r['runtime_input']['xyz']).parent;start=read(r['start']);require(start['job']==job and start['runtime_input']==r['runtime_input'] and start['manifest']==cfg['manifest'],'Native start differs')
        require(len(r['receipts'])==(1 if label in m['reused_gradients'] else 2),'Native receipt denominator')
        for rec in r['receipts']:
            require(rec['returncode']==0 and rec['cwd']==str(folder) and all(v=='1' for v in rec['threads'].values()),'Native receipt unavailable');verify(rec['output'])
        require(r['receipts'][0]['argv']==[s['prototype_analyze']['path'],'system.xyz','EM'],'Analyze invocation differs')
        require(moment.parse_moments(folder/'analyze.out',moment.induced_file(folder/'analyze.out'),row,ind.induction)==r['moments'],'Moment raw replay differs')
        q.response_readback((folder/'analyze.out').read_text(),job['profile'],row['types'],(folder/'system.key').read_text())
        if label in m['reused_gradients']:
            old=m['reused_gradients'][label];rec={'status':'complete','input':old['input'],'receipts':[old['receipt']],'value':old['value'],'job':{'id':old['producer_record']}}
            require(reuse_gradient(s,rec,job['input'],row,job['profile'],h)==r['gradient_reuse'] and r['gradient']==old['value'],'Inherited gradient differs')
        else:
            require(r['receipts'][1]['argv']==[s['prototype_testgrad']['path'],'system.xyz','Y','N'],'Gradient invocation differs')
            require(grad.gradient(folder/'testgrad.out',row['atoms'],False,native)==r['gradient'],'Gradient raw replay differs');ind.induction(folder/'testgrad.out')
            q.response_readback((folder/'testgrad.out').read_text(),job['profile'],row['types'],(folder/'system.key').read_text())
        require(abs(r['gradient']['energy']-r['moments']['energy_kcal_mol'])<=1e-4,'Raw energy closure failed')
    require(summarize(m,result['rows'],h)==result['profiles'],'Saved profile arithmetic differs')
    return {'status':'saved_sensitivity_admitted','rows':result['rows'],'profiles':result['profiles'],'new_native_calls':0}

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('action',choices=['prepare','run']);ap.add_argument('--out',required=True,type=Path);ap.add_argument('--config',type=Path);a=ap.parse_args()
    a.out.mkdir(parents=True,exist_ok=False);v=prepare(a.out) if a.action=='prepare' else run(json.loads(a.config.read_text()),a.out)
    print(json.dumps({'status':v['status'],'out':str(a.out)}));raise SystemExit(1 if v['status']=='failed' else 0)
