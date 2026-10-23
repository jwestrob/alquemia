#!/usr/bin/env python3
"""Single49-program numerical qualification. No fitting/cache bridge admission."""
import argparse,ast,concurrent.futures,copy,datetime,hashlib,json,math,os,resource,shutil,threading,time
from pathlib import Path
import numpy as np
import PREPARE as P
import CHECKS as C

def fingerprint(cfg):return hashlib.sha256(json.dumps({k:v for k,v in cfg.items() if k!='root_release'},sort_keys=True,separators=(',',':'),allow_nan=False).encode()).hexdigest()
def context(cfg,runtime=False):
    P.require(cfg['implementation']==P.pin(__file__) and cfg['preparer']==P.pin(P.__file__) and cfg['checks']==P.pin(C.__file__),'Frozen adapter source')
    m=P.read(cfg['manifest']);P.require(m['implementation']==cfg['preparer'] and m['sources']['runtime']==cfg['implementation'] and m['sources']['checks']==cfg['checks'],'Exact preparation/source join')
    P.require(cfg['limits']=={'native_programs':49,'internal_FD_energies':132,'workers':12,'program_timeout_seconds':120},'Fixed49/132 scope')
    [P.verify(q) for q in m['sources'].values()];[P.verify(q) for q in m['radial_inputs'].values()]
    P.require([{k:v for k,v in j.items() if k not in ('input','tangent','runtime_positions_A')} for j in m['programs']]==P.inventory(),'Exact declared jobs')
    for j in m['programs']:
        [P.verify(q) for q in j['input'].values()];P.require(P.verify(j['input']['key']).read_text()==P.render(m['rows'][j['cell']],j,m['units']),'Literal prepared key')
    P.require(sum(6*m['rows'][j['cell']]['atoms'] for j in m['programs'] if j['fd'])==132,'FD denominator')
    old=P.module(m['radial_inputs']['parent_qualification_runner'],'radial_parent_qualification');h=old.helpers(m['helper_sources']);ev=P.module(m['radial_inputs']['inherited_evaluator'],'radial_parent_evaluator')
    if runtime:
        P.require(all(os.environ.get(k)=='1' for k in P.THREADS),'Numerical threads1')
        P.require(os.environ.get('SLURM_JOB_ID','').isdigit() and int(os.environ.get('SLURM_NTASKS','1'))==1 and int(os.environ.get('SLURM_CPUS_PER_TASK','0'))==12,'Actual12CPU single Slurm task')
        ready=P.read(cfg['radial_ready']);sr=P.read(cfg['native_source_review']);P.require(all(sr[k]==v for k,v in ready['source_review_contract'].items()),'Independent native source review')
        b=P.read(cfg['build']);P.require(b['status']=='complete' and b['manifest']==cfg['source_manifest']==ready['source_manifest'] and b['executables']==cfg['executables'],'Actual built native lineage')
        [P.verify(q) for q in cfg['executables'].values()];rr=P.read(cfg['runner_source_review']);P.require(rr['status']=='passed_radial_qualification_source_review' and all(rr[k]==cfg[k] for k in ('implementation','preparer','checks','manifest')),'Exact runner review')
        release=P.read(cfg['root_release']);P.require(release['status']=='released_radial_drive_numerical_qualification' and release['configuration_fingerprint']==fingerprint(cfg) and release['limits']==cfg['limits'] and release['source_manifest']==cfg['source_manifest'] and release['output']==cfg['output'],'Exact root release')
    return m,h,ev

def invocation(h):
    # Reuse the qualified invocation verbatim except its two literal300s limits.
    tree=ast.parse(Path(h['native'].__file__).read_text());fn=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='invoke');fn=copy.deepcopy(fn);fn.name='radial_invoke';changed=0
    for n in ast.walk(fn):
        if isinstance(n,ast.Constant) and n.value==300:n.value=120;changed+=1
    P.require(changed==2,'Exact inherited timeout transformation');ns=dict(vars(h['native']));exec(compile(ast.fix_missing_locations(ast.Module(body=[fn],type_ignores=[])),str(Path(__file__)), 'exec'),ns);return ns['radial_invoke']

def run(cfg,out):
    t=time.monotonic();cpu=time.process_time();ch0=resource.getrusage(resource.RUSAGE_CHILDREN);m,h,ev=context(cfg,True);invoke=invocation(h);lock=threading.RLock();stop=threading.Event()
    result={'status':'running','implementation':P.pin(__file__),'configuration':cfg,'manifest':cfg['manifest'],'source_manifest':cfg['source_manifest'],'build':cfg['build'],'executables':cfg['executables'],'root_release':cfg['root_release'],'job_id':os.environ['SLURM_JOB_ID'],'host':os.uname().nodename,'pid':os.getpid(),'rows':{j['id']:{'status':'not_attempted'} for j in m['programs']},'native_program_invocation_attempts':0,'native_programs_completed':0,'native_programs_returned_zero':0,'internal_FD_energies_observed':0,'internal_FD_energies_requested':0,'internal_FD_energies_completed':0,'operator_checks':{},'directional_checks':{},'field_checks':{},'global_stop':False,'error':None,'numerical_checks_passed':False,'zero_parity_passed':False,'cache_bridge_admitted':False,'physical_accuracy_qualified':False}
    def persist():
        with lock:
            result['global_stop']=stop.is_set();tmp=out/'RESULT.tmp';P.save(tmp,result);os.replace(tmp,out/'RESULT.json')
    def worker(job):
        if stop.is_set():return
        folder=out/'cells'/job['id'];rec={'status':'not_attempted','job':job,'receipts':[]};claimed=False;row=m['rows'][job['cell']];fd=6*row['atoms'] if job['fd'] else 0
        try:
            folder.mkdir(parents=True,exist_ok=False)
            for k,name in P.FILES:shutil.copyfile(P.verify(job['input'][k]),folder/name)
            rec['runtime_input']={k:P.pin(folder/name) for k,name in P.FILES}
            with lock:
                if stop.is_set():return
                P.require(result['native_program_invocation_attempts']<49 and result['internal_FD_energies_requested']+fd<=132,'Hard invocation/FD ceiling')
                claimed=True;rec['status']='started';result['native_program_invocation_attempts']+=1;result['internal_FD_energies_requested']+=fd;result['rows'][job['id']]=copy.deepcopy(rec);persist()
            name=job['program'];args=['system.xyz','Y','Y','0.00001'] if job['fd'] else ['system.xyz','Y','N'] if name=='testgrad' else ['system.xyz','E' if job['stage']=='directions' else 'EM'];output=folder/(name+'.out')
            invoke([cfg['executables'][name]['path']]+args,folder,output,rec['receipts']);rec['actual_output']=P.pin(output)
            with lock:result['native_programs_returned_zero']+=1
            C.parse(rec,job,row,m,h,ev);rec['status']='complete'
            with lock:result['native_programs_completed']+=1;result['internal_FD_energies_completed']+=fd
        except Exception as e:
            stop.set();rec.update(status='failed',error=repr(e))
            with lock:
                if result['error'] is None:result['error']=repr(e)
        finally:
            with lock:
                if folder.exists():
                    for name in ('analyze.out','testgrad.out'):
                        if (folder/name).exists():rec['actual_output']=P.pin(folder/name)
                    P.save(folder/'RESULT.json',rec)
                if claimed and job['fd'] and 'response_readback' in rec:result['internal_FD_energies_observed']+=rec['response_readback']['native_solve_blocks']-1
                if claimed or rec['status']=='failed':result['rows'][job['id']]=copy.deepcopy(rec)
                persist()
    def row(stage,g,program='analyze'):return result['rows']['__'.join((stage,g,program))]
    try:
        persist()
        for stage in ('absent','zero','active','directions','fields','bypass'):
            jobs=[j for j in m['programs'] if j['stage']==stage]
            with concurrent.futures.ThreadPoolExecutor(max_workers=12) as pool:list(pool.map(worker,jobs))
            P.require(not stop.is_set() and all(result['rows'][j['id']]['status']=='complete' for j in jobs),'Failed stage '+stage)
            if stage=='zero':result['zero_parity_passed']=True
            if stage=='active':
                for g in P.METAL:
                    result['operator_checks'][g]=C.operator_checks(row('active',g),row('zero',g),m['rows'][g],m,h);checks={};C.gate(checks,'E_testgrad_analyze_closure',abs(row('active',g)['moments']['energy_kcal_mol']-row('active',g,'testgrad')['value']['energy']),1e-4);result['operator_checks'][g].update(checks)
            if stage=='directions':
                for kind,cells in (('metal_x',P.CENTERS),('water_rotation_x',P.DIMERS)):
                    for g in cells:
                        _,tangent=P.path_geometry(m['rows'][g],kind,0.);gradient=row('active',g,'testgrad')['value']['gradient'];analytic=math.fsum(a*b for aa,bb in zip(gradient,tangent) for a,b in zip(aa,bb));steps=[]
                        for step in (1e-4,5e-5):
                            vals={j['sign']:result['rows'][j['id']]['value']['energy'] for j in jobs if j['cell']==g and j['kind']==kind and j['h']==step};slope=(vals[1]-vals[-1])/(2*step);bound=1e-8/(2*step)+5e-9*sum(abs(x) for v in tangent for x in v);checks={};C.gate(checks,'analytic_vs_FD',abs(slope-analytic),.001+bound);steps.append({'half_step_A':step,'FD':slope,'analytic':analytic,'printing_bound':bound,'checks':checks})
                        checks={};C.gate(checks,'two_step',abs(steps[0]['FD']-steps[1]['FD']),.001+sum(s['printing_bound'] for s in steps));result['directional_checks'][g+'__'+kind]={'steps':steps,'checks':checks}
            if stage=='fields':
                for g in P.DIMERS:
                    f={j['sign']:result['rows'][j['id']] for j in jobs if j['cell']==g};result['field_checks'][g]=C.field_checks(row('active',g),f[-1],f[1],m['rows'][g],m,h)
            persist()
        P.require(result['native_program_invocation_attempts']==result['native_programs_returned_zero']==result['native_programs_completed']==49 and result['internal_FD_energies_requested']==result['internal_FD_energies_observed']==result['internal_FD_energies_completed']==132 and all(v['status']=='complete' for v in result['rows'].values()),'Exact49/132 terminal accounting')
        result.update(status='radial_drive_numerical_qualification_complete',numerical_checks_passed=True)
    except Exception as e:stop.set();result.update(status='radial_drive_numerical_qualification_failed',error=result['error'] or repr(e))
    finally:
        ch=resource.getrusage(resource.RUSAGE_CHILDREN);result['cost']={'wall_seconds':time.monotonic()-t,'controller_CPU_seconds':time.process_time()-cpu,'child_user_CPU_seconds':ch.ru_utime-ch0.ru_utime,'child_system_CPU_seconds':ch.ru_stime-ch0.ru_stime,'peak_child_RSS_KiB':ch.ru_maxrss,'allocated_CPUs':12,'per_receipt_child_deltas_not_summed':True};persist()
    return result
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('mode',choices=('preflight','run'));ap.add_argument('--config',required=True);ap.add_argument('--out',required=True);args=ap.parse_args();cfg=json.loads(Path(args.config).read_text());out=Path(args.out).resolve();out.mkdir(parents=True,exist_ok=False)
    if args.mode=='preflight':
        t=time.monotonic();cpu=time.process_time()
        try:m,h,ev=context(cfg,False);r={'status':'passed_zero_call_preflight','implementation':P.pin(__file__),'configuration':P.pin(args.config),'native_calls':0,'solver_calls':0,'programs':len(m['programs']),'internal_FD_energies':132,'cost':{'wall_seconds':time.monotonic()-t,'CPU_seconds':time.process_time()-cpu}};P.save(out/'RESULT.json',r)
        except Exception as e:P.save(out/'FAILURE.json',{'error':repr(e),'native_calls':0});raise
    else:P.require(str(out)==cfg['output'],'Exact released output');r=run(cfg,out)
    print(json.dumps({'status':r['status']}));raise SystemExit(0 if 'failed' not in r['status'] else 1)
