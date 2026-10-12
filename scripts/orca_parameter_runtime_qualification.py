#!/usr/bin/env python3
"""Finite native GFN2 runtime qualification; isolated artifacts, <=4 CPU children.

Numerical admission is fixed by the supplied October12 task. No installation,
Slurm, protein scoring, retries, or changes to the production scorer.
"""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import argparse, copy, ctypes, datetime, hashlib, importlib.util, json, math
import os, re, shutil, signal, subprocess, sys, threading, time, tomllib

ROOT=Path('/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs')
W=ROOT/'workspaces/orca_parameter_audit_20261012/runtime_qualification_v1'
REPORT=ROOT/'diagnostics/orca_parameter_audit_20261012/runtime'
V=Path('/home/jwestrob/jwestrob/obsidian-vault/projects/lanthanide-binding/orca_parameter_audit_2026-10-12')
SOURCE=Path('/groups/banfield/users/jwestrob/EastRiver/EastRiver_PLM/revision_analysis/2026-10-09_Nikasha_native_parity')
ORCA=Path('/groups/banfield/users/jwestrob/bin/ORCA/orca_6_1_1_linux_x86-64_shared_openmpi418_nodmrg/orca')
DRIVER_SHA='38b5f057452fef275c0a1b98d270ad03d0411dab70453820d1c678bd732f6c83'
NUM=r'[-+]?(?:\d+\.?\d*|\.\d+)(?:[Ee][-+]?\d+)?'
THREAD_KEYS=['OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS','BLIS_NUM_THREADS','GOTO_NUM_THREADS']
ENV={**{k:'1' for k in THREAD_KEYS},'OMP_DYNAMIC':'FALSE','MKL_DYNAMIC':'FALSE','LC_ALL':'C','LANG':'C'}
LOCK=threading.Lock()

def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()

def pin(p):
    p=Path(p).resolve();return {'path':str(p),'sha256':sha(p),'bytes':p.stat().st_size}

def read(p):return json.loads(Path(p).read_text())
def now():return datetime.datetime.now(datetime.timezone.utc).isoformat()

def write(p,d,replace=False):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
    if not replace and p.exists():raise FileExistsError(p)
    tmp=p.with_name(p.name+'.tmp');tmp.write_text(json.dumps(d,indent=2,allow_nan=False)+'\n');tmp.replace(p)

def checked_source(rel):
    m=read(V/'MANIFEST.json');rec=next(x for x in m['files'] if x['path']==rel)
    p=V/rel
    if sha(p)!=rec['sha256']:raise ValueError('snapshot changed: '+rel)
    return p

def upstream():
    spec=importlib.util.spec_from_file_location('saved_native_controls',W/'source/run_native_controls.py')
    m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m

def runtime_binding():
    files={};ldd={}
    for name in ['orca','orca_startup','orca_leanscf']:
        p=ORCA.parent/name;files[str(p.resolve())]=pin(p)
        r=subprocess.run(['/usr/bin/ldd',str(p)],capture_output=True,text=True,check=True)
        ldd[name]=r.stdout
        if 'not found' in r.stdout:raise ValueError('missing runtime library')
        for line in r.stdout.splitlines():
            found=re.search(r'=>\s+(/\S+)',line) or re.search(r'^\s*(/\S+)',line)
            if found:
                q=Path(found[1]).resolve()
                if str(q) not in files:files[str(q)]=pin(q)
    if files[str(ORCA.resolve())]['sha256']!=DRIVER_SHA:raise ValueError('driver differs from historical binary')
    return {'schema':'native_runtime_binding.v1','created_utc':now(),'driver':pin(ORCA),'files':list(files.values()),'ldd':ldd,'environment':{k:os.environ.get(k) for k in ['PATH','LD_LIBRARY_PATH','LD_PRELOAD']},'thread_environment':ENV}

def verify_binding():
    b=read(W/'RUNTIME_BINDING.json')
    for r in b['files']:
        if sha(r['path'])!=r['sha256']:raise ValueError('runtime changed '+r['path'])
    return b

def block(text,name):
    m=re.search(r'^\$'+re.escape(name)+r'\s*\n(.*?)^\$End',text,re.M|re.S)
    return m[1] if m else ''

def scalar(text,name):
    m=re.search(r'&'+re.escape(name)+r'\s*\[.*?\]\s*('+NUM+r')',text)
    return float(m[1]) if m else None

def properties(directory):
    p=Path(directory)/'endpoint.property.txt'
    if not p.exists():return {'charges_e':None,'normal_termination':False,'converged':False}
    text=p.read_text();info=block(text,'Calculation_Info');sp=block(text,'Single_Point_Data')
    chargeblock=block(text,'SCF_Mulliken_Population_Analysis').split('&AtomicCharges')
    rows=[]
    if len(chargeblock)>1:
        tail=re.split(r'^\s*&',chargeblock[1],flags=re.M)[0]
        rows=[(int(i),float(x)) for i,x in re.findall(r'^\s*(\d+)\s+('+NUM+r')\s*$',tail,re.M)]
    return {'energy_hartree':scalar(sp,'FinalEnergy'),'charge':scalar(info,'Charge'),'multiplicity':scalar(info,'Mult'),'valence_electrons':scalar(info,'NumOfElectrons'),'atoms':scalar(info,'NumOfAtoms'),
            'charges_e':[x for _,x in rows],'charge_indices':[i for i,_ in rows],
            'normal_termination':'NORMAL TERMINATION' in block(text,'Calculation_Status'),
            'converged':bool(re.search(r'&Converged\s*\[.*?\]\s*true',sp))}

def compare(a,b,etol=1e-6,qtol=1e-5):
    e=a.get('energy_hartree');ref=b.get('energy_hartree');q=a.get('charges_e');rq=b.get('charges_e')
    de=None if e is None or ref is None else e-ref
    dq=None if not q or not rq or len(q)!=len(rq) else max(abs(x-y) for x,y in zip(q,rq))
    numerical=de is not None and dq is not None and abs(de)<=etol and dq<=qtol
    return {'energy_delta_hartree':de,'energy_delta_kcal_mol':None if de is None else de*627.509474,'max_atom_charge_delta_e':dq,'energy_tolerance_hartree':etol,'charge_tolerance_e':qtol,'numerical_pass':numerical,
            'both_accepted':bool(a.get('accepted') and b.get('accepted')),'accepted_parity':bool(numerical and a.get('accepted') and b.get('accepted'))}

def prepare():
    W.mkdir(parents=True,exist_ok=False);REPORT.mkdir(parents=True,exist_ok=True);(W/'source').mkdir()
    src=SOURCE/'run_native_controls.py';shutil.copyfile(src,W/'source/run_native_controls.py')
    manifest=read(SOURCE/'bundle_manifest.json')
    for job in manifest['jobs']:
        dest=W/job['input_dir'];dest.mkdir(parents=True)
        for name,rec in job['files'].items():
            source=SOURCE/job['input_dir']/name
            if sha(source)!=rec['sha256']:raise ValueError('fixture changed')
            shutil.copyfile(source,dest/name)
        archived=checked_source('la_followup/native_reciprocal_controls/results/'+job['id']+'/receipt.json')
        prop=checked_source('la_followup/native_reciprocal_controls/results/'+job['id']+'/endpoint.property.txt')
        job['archived_receipt']=pin(archived);job['archived_property']=pin(prop)
        job['stage']='La8'
    manifest['remote_root']=str(W);manifest['created_utc']=now();manifest['limits']['workers']=4
    manifest['source_bundle']=pin(SOURCE/'bundle_manifest.json');manifest['source_runner']=pin(src)
    manifest['task_document']=pin(ROOT/'docs/orca_parameter_audit_20261012/runtime_qualification.md')
    manifest['common_rules']=pin(ROOT/'docs/orca_parameter_audit_20261012/COMMON_RULES.md')
    manifest['runtime_helper']=pin(__file__)
    manifest['admission']={'roundtrip_energy_Eh':1e-8,'roundtrip_charge_e':1e-7,'parity_energy_Eh':1e-6,'parity_charge_e':1e-5,'all_printed_density_criteria':True,'no_retries':True}
    write(W/'manifest_La8.json',manifest);write(W/'RUNTIME_BINDING.json',runtime_binding())
    plan={'question':'Does the current exact Linux runtime reproduce the8 archived Ca/La controls, then the finite11-state aquo matrix?','scientific_scope':'implementation qualification only','authorization':'User explicitly permits <=16 total login numerical threads; root allocates <=4 to this task while Slurm is down.','maximum_native_calls':58,'new_other_solver_calls':0,'stages':['serial Ca vacuum scout','remaining7 iff first scout passes','Cu1 aquo6 scout iff La8 passes','remaining44 aquo points'],'settings':'exact archived La inputs; aquo default/roundtrip TolE1e-10, corrected_tight1e-12; no grid/state/geometry/guess edits','prose_clarification':'Root confirmed actual supplied roundtrip1e-10 and corrected_tight1e-12 outrank mistaken tighter-roundtrip summary, before stage2.','failure_policy':'no automatic retries or post-result tolerances; raw agreement and accepted agreement separate','limits':manifest['limits'],'parameter_changes':'La3 arrays; subsequent target metal kcn/shpoly/lgam/refocc only, shell order from pinned schema','source_snapshot':pin(V/'MANIFEST.json')}
    write(W/'PLAN.json',plan);status('prepared',[],next_step='Ca_vacuum_default scout')
    (W/'CHECKPOINT.md').write_text('# Runtime qualification checkpoint\n\nQuestion and frozen scope: PLAN.json. Prepared exact8 fixtures; no molecular calls yet. Current authorization <=4 numerical threads, <=16 across root campaign. First scout precedes other7; no repeats or production changes.\n')
    print(json.dumps({'prepared':str(W),'jobs':8,'runtime_files':len(read(W/'RUNTIME_BINDING.json')['files'])}))

def status(state,owned,**extra):
    with LOCK:write(W/'STATUS.json',{'status':state,'updated_utc':now(),'owned_processes':owned,**extra},replace=True)

def run_job(job,manifest_path):
    module=upstream();jobid=job['id'];d=W/'results'/jobid;d.mkdir(parents=True,exist_ok=False)
    for name,rec in job['files'].items():
        source=W/job['input_dir']/name
        if sha(source)!=rec['sha256']:raise ValueError('prepared input changed')
        shutil.copyfile(source,d/name)
    env=os.environ.copy();env.update(ENV);env['ORCA_TMPDIR']=str(d)
    for key in ['OMP_PLACES','OMP_PROC_BIND']:env.pop(key,None)
    observed={};libraries=set();max_threads=0;max_runnable=0;max_member_threads=0
    def monitor(pgid):
        nonlocal max_threads,max_runnable,max_member_threads
        rss,members=module.group_snapshot(pgid);threads=0;runnable=0
        for member in members:
            p=Path('/proc')/str(member['pid'])
            try:
                st=(p/'status').read_text();n=int(re.search(r'^Threads:\s+(\d+)',st,re.M)[1]);threads+=n
                max_member_threads=max(max_member_threads,n);runnable+=member['state']=='R'
                exe=str((p/'exe').resolve());cmd=(p/'cmdline').read_bytes().replace(b'\0',b' ').decode(errors='replace')
                observed[member['pid']]={'pid':member['pid'],'exe':exe,'command':cmd,'max_threads':max(n,observed.get(member['pid'],{}).get('max_threads',0))}
                for line in (p/'maps').read_text().splitlines():
                    parts=line.split()
                    if len(parts)>5 and parts[-1].startswith('/') and ('.so' in parts[-1] or parts[-1].endswith('/orca')):libraries.add(parts[-1])
            except (FileNotFoundError,ProcessLookupError):pass
        max_threads=max(max_threads,threads);max_runnable=max(max_runnable,runnable)
        if max_member_threads>1:raise RuntimeError('unexpected multithreaded owned process')
        return rss,members
    receipt={'schema':'native_runtime_qualification_receipt.v1','job':job,'manifest_sha256':sha(manifest_path),'binary':pin(ORCA),'runtime_binding':pin(W/'RUNTIME_BINDING.json'),'started_utc':now(),'thread_environment':ENV,'cwd':str(d),'initial_files':sorted(x.name for x in d.iterdir())}
    write(d/'STARTED.json',{'command':[str(ORCA),'endpoint.inp','--bind-to none'],'job':jobid,'start':now()})
    process=module.bounded_process([str(ORCA),'endpoint.inp','--bind-to none'],d,env,240 if job['stage']=='La8' else 60,6442450944,monitor)
    expected=job['expected'];out=module.parse_output(d/'endpoint.out',expected)
    prop=properties(d);q=prop.get('charges_e');checks={'property_normal':prop['normal_termination'],'property_converged':prop['converged'],
        'property_state':all(prop.get(k)==expected[k] for k in ['charge','multiplicity','valence_electrons','atoms']),
        'all_atom_population_count':q is not None and len(q)==expected['atoms'] and prop['charge_indices']==list(range(expected['atoms'])),
        'charge_sum':bool(q) and abs(sum(q)-expected['charge'])<=1e-5,
        'property_stdout_energy':prop.get('energy_hartree') is not None and out.get('energy_hartree') is not None and abs(prop['energy_hartree']-out['energy_hartree'])<=1e-9,
        'version':bool(re.search(r'Program Version\s+6\.1\.1\s', (d/'endpoint.out').read_text()))}
    out.update(charges_e=q,property=prop,property_checks=checks)
    inputs_ok=all(sha(d/name)==rec['sha256'] for name,rec in job['files'].items())
    out['accepted']=bool(process['exit_code']==0 and process['stop_reason'] is None and out['passed'] and all(checks.values()) and inputs_ok)
    out['energy_hartree']=prop.get('energy_hartree')
    if job['stage']=='La8':
        old=read(job['archived_receipt']['path']);a=properties(Path(job['archived_property']['path']).parent);a['accepted']=old['passed']
        receipt['comparison_to_archived']=compare(out,a,1e-8,1e-7)
    else:
        data=read(W/'source/aquo_data.json');rows=[r for r in data['runs'] if r['state']==job['state'] and r['scope']=='aquo' and r['medium']==job['medium']]
        receipt['comparison_to_archived_engines']={engine:compare(out,next(r for r in rows if r['engine']==engine)) for engine in ['xtb','tblite']}
        old=next(r for r in rows if r['engine']=='native' and r['variant']==job['variant']);receipt['comparison_to_archived_Mac_native']=compare(out,old)
    receipt.update(process=process,output=out,inputs_after_match=inputs_ok,passed=out['accepted'],finished_utc=now(),observed_processes=list(observed.values()),observed_library_paths=sorted(libraries),max_process_threads=max_member_threads,max_group_threads_including_idle_driver=max_threads,max_runnable_group_processes=max_runnable)
    receipt['artifacts']={p.name:pin(p) for p in sorted(d.iterdir()) if p.is_file()}
    write(d/'receipt.json',receipt)
    print(json.dumps({'job':jobid,'passed':receipt['passed'],'elapsed_seconds':process['elapsed_seconds'],'energy_hartree':out['energy_hartree'],'max_process_threads':max_member_threads,'old_parity':receipt.get('comparison_to_archived',{}).get('accepted_parity')}),flush=True)
    return receipt

def launch(stage):
    verify_binding();libc=ctypes.CDLL(None,use_errno=True)
    if libc.prctl(36,1,0,0,0)!=0:raise OSError('subreaper unavailable')
    mp=W/'manifest_La8.json';m=read(mp)
    if stage=='scout':jobs=[next(x for x in m['jobs'] if x['id']=='Ca_vacuum_default')]
    elif stage=='La7':
        r=read(W/'results/Ca_vacuum_default/receipt.json')
        if not r['passed'] or not r['comparison_to_archived']['accepted_parity']:raise ValueError('first scout failed')
        jobs=[x for x in m['jobs'] if x['id']!='Ca_vacuum_default']
    else:
        mp=W/'manifest_aquo50.json';m=read(mp)
        if stage=='Cu1':jobs=[x for x in m['jobs'] if x['state']=='Cu1']
        elif stage=='remaining44':
            scout=read(W/'SCOUT_Cu1.json')
            if not scout['corrected_and_roundtrip_pass']:raise ValueError('Cu1 gate failed')
            jobs=[x for x in m['jobs'] if x['state']!='Cu1']
        else:raise ValueError(stage)
    status('running',[],stage=stage,requested_jobs=[j['id'] for j in jobs],runner_pid=os.getpid())
    # Separate dependent default and import stages; at most four serial native workers.
    results=[]
    for group in [[j for j in jobs if j['variant']=='default'],[j for j in jobs if j['variant']!='default']]:
        if not group:continue
        with ThreadPoolExecutor(max_workers=min(4,len(group))) as ex:results.extend(ex.map(lambda job:run_job(job,mp),group))
    write(W/f'EXECUTION_{stage}.json',{'stage':stage,'receipts':[pin(W/'results'/r['job']['id']/'receipt.json') for r in results],'completed':len(results),'accepted':sum(r['passed'] for r in results),'elapsed_sum_seconds':sum(r['process']['elapsed_seconds'] for r in results)})
    status('stage_complete',[],stage=stage,completed=len(results),accepted=sum(r['passed'] for r in results))

def main():
    a=argparse.ArgumentParser();a.add_argument('op',choices=['prepare','scout','La7']);args=a.parse_args()
    if args.op=='prepare':prepare()
    else:launch(args.op)

if __name__=='__main__':main()
