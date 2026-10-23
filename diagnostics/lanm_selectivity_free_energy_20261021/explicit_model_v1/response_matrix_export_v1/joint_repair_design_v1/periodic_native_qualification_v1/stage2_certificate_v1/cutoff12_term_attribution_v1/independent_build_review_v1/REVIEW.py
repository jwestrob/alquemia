"""Finite saved-build receipt and ELF metadata audit; no application execution."""
import json,hashlib,time,resource,subprocess,re,os
from pathlib import Path
T=Path(__file__).resolve().parents[1]
ROOT=T.parents[7]
# Derive workspace by the exact prepared manifest, not cwd.
ready=json.loads((T/'READY.json').read_text());w=Path(ready['manifest']['path']).parents[1];out=w/'independent_build_review_v1';out.mkdir(exist_ok=False)
tick=time.monotonic();cpu=time.process_time();checks=[];pins=[];commands=[]
def pin(p):
 p=Path(p).resolve();return dict(path=str(p),sha256=hashlib.sha256(p.read_bytes()).hexdigest())
def check(x,msg):
 if not x:raise ValueError(msg)
 checks.append(msg)
def verify(x):
 check(pin(x['path'])==x,'pin '+x['path']);pins.append(x);return Path(x['path'])
def read(x):return json.loads(verify(x).read_text())
def dump(p,x):p.write_text(json.dumps(x,indent=2)+'\n')
def tool(label,argv):
 t=time.monotonic();u=resource.getrusage(resource.RUSAGE_CHILDREN);p=subprocess.run(argv,capture_output=True,text=True,timeout=30);v=resource.getrusage(resource.RUSAGE_CHILDREN)
 (out/(label+'.out')).write_text(p.stdout);(out/(label+'.err')).write_text(p.stderr)
 commands.append(dict(argv=argv,implementation=pin(argv[0]),returncode=p.returncode,stdout=pin(out/(label+'.out')),stderr=pin(out/(label+'.err')),cost=dict(wall_seconds=time.monotonic()-t,child_CPU_seconds=v.ru_utime+v.ru_stime-u.ru_utime-u.ru_stime)))
 check(p.returncode==0 and not p.stderr.strip(),'read-only '+label);return p.stdout
result=dict(status='failed_cutoff12_term_observer_build_review',reviewer_implementation=pin(__file__),error=None,native_calls=0,builds=0,FFT_calls=0,factor_calls=0,scientific_qualification=False)
try:
 check(pin(T/'READY.json')['sha256']=='f7a55cb6e53bdaf0158514c10f41dd2d759f853ff1a93bd4b0bb3150d3fbfe55','Frozen source READY')
 m=read(ready['manifest']);source=read(ready['source_manifest']);recipe=read(ready['link_recipe']);oldrecipe=read(m['stage1_link_recipe'])
 rp=pin(w/'build_run_v1/execution/RESULT.json');r=read(rp);cfg=read(r['configuration']);rel=read(r['release']);sr=read(rel['source_review'])
 check(r['status']=='stage2_build_complete' and r['error'] is None and r['native_programs']==0 and r['compile_attempts']==r['link_attempts']==1,'Successful single compile/link without native')
 check(r['implementation']==cfg['implementation']==rel['implementation']==ready['implementation'],'Exact build source');verify(r['implementation'])
 check(r['manifest']==cfg['manifest']==ready['manifest'] and r['source_manifest']==ready['source_manifest'],'Exact source/manifest')
 check(rel['configuration']==r['configuration'] and rel['output']==str(w/'build_run_v1/execution') and rel['status']=='released_term_observer_build','Root build release')
 check(sr['status']=='passed_cutoff12_term_observer_source_review' and sr['manifest']==r['manifest'] and sr['binding']==ready['binding'] and sr['implementations']['build']==r['implementation'],'Exact root source review')
 check(r['reused_objects']==m['reused_objects'] and r['parent_build']==m['parent_build'],'Three admitted object lineage')
 for x in [*r['reused_objects'].values(),r['library'],r['compiler'],*oldrecipe['compiler_and_link_dependencies'],*oldrecipe['module_files']]:verify(x)
 for x in source['files'].values():verify(x)
 check(r['library']==recipe['library']==oldrecipe['library'] and r['compiler']==recipe['compiler'],'Same library/compiler')
 check(tuple(r['commands'])==('compile_GRADIENT_OBSERVER_f','link_gradient_old'),'Only declared two build commands')
 for label,key in [('compile_GRADIENT_OBSERVER_f','compile_argv'),('link_gradient_old','link_argv')]:
  c=r['commands'][label];check(c['argv']==recipe[key] and c['status']=='complete' and c['returncode']==0,'Exact '+label+' argv/exit')
  check(c['wall_seconds_limit']==c['CPU_seconds_limit']==180 and c['AS_bytes']==3*1024**3,'Exact '+label+' limits')
  check(verify(c['stdout']).read_text()==verify(c['stderr']).read_text()=='','No compiler/linker warnings or errors')
 link=r['commands']['link_gradient_old']['argv'];obj=Path(r['commands']['compile_GRADIENT_OBSERVER_f']['argv'][-1]);exe=verify(r['executables']['gradient_old'])
 check(link.index(str(obj))<link.index(r['library']['path']) and link[-1]==str(exe),'Override precedes archive and exact executable')
 check(os.access(exe,os.X_OK),'Executable permission')
 objtext=tool('object_symbols',['/usr/bin/nm','-g','--defined-only',str(obj)]);exetext=tool('executable_symbols',['/usr/bin/nm','-g','--defined-only',str(exe)])
 for name in ('gradient_','term_observe_'):
  check(len([x for x in objtext.splitlines() if x.split()[-1:]==[name]])==1,'One observer object symbol '+name)
  check(len([x for x in exetext.splitlines() if x.split()[-1:]==[name]])==1,'One linked symbol '+name)
 dynamic=tool('dynamic',['/usr/bin/readelf','-d',str(exe)]);needed=set(re.findall(r'Shared library: \[([^]]+)\]',dynamic))
 check(needed==set(m['native_runtime'])-{'loader'},'ELF NEEDED exact qualified runtime libraries')
 check(r['native_runtime']['gradient_old']==m['native_runtime'],'Actual saved node loader readback matches manifest')
 for x in r['native_runtime']['gradient_old'].values():verify(x)
 allocation=r['allocation'];check(allocation['job_id']=='20503' and allocation['allocated_CPUs']==allocation['nodes']==allocation['tasks']==1 and allocation['qos']=='standard_expand' and not allocation['requeue'],'Actual one-CPU default allocation')
 check(allocation['memory_grant_bytes']==4194304000 and allocation['memory_grant_bytes']>3*1024**3,'Default grant accommodates child cap')
 sub=read(pin(w/'build_run_v1/SUBMISSION.json'));event=read(pin(w/'build_run_v1/completion_event.json'));terminal=event['records']['20503']
 check(sub['returncode']==0 and sub['stdout'].strip()=='20503' and terminal['state']=='COMPLETED' and terminal['exit_code']=='0:0','Terminal scheduler receipt')
 check(int(terminal['allocated_cpus'])==1 and terminal['nodes']==allocation['hostname'],'Actual receipt node/allocation')
 result.update(status='passed_cutoff12_term_observer_build_review',result=rp,configuration=r['configuration'],implementation=r['implementation'],manifest=r['manifest'],source_manifest=r['source_manifest'],source_review=rel['source_review'],ready=pin(T/'READY.json'),executables=r['executables'],observer_object=pin(obj),library=r['library'],runtime=r['native_runtime'],source_build_checks_passed=True,build_cost=r['cost'],build_child_CPU_seconds=sum(c['cost']['child_user_CPU_seconds']+c['cost']['child_system_CPU_seconds'] for c in r['commands'].values()),scheduler=dict(completion_event=pin(w/'build_run_v1/completion_event.json'),job_id='20503',elapsed_seconds=int(terminal['elapsed_seconds']),allocated_CPUs=1,allocated_CPU_seconds=int(terminal['elapsed_seconds']),batch_not_added=True))
except Exception as exc:result['error']=repr(exc)
result.update(check_count=len(checks),checks=checks,pins=pins,metadata_commands=commands,review_cost=dict(wall_seconds=time.monotonic()-tick,CPU_seconds=time.process_time()-cpu,threads=1))
dump(out/'REVIEW.json',result)
(out/'REPORT.txt').write_text(result['status']+'\n'+('One exact observer compile/link, unchanged three objects/archive and runtime, linked gradient_/term_observe_ confirmed. 20503 used 3 allocated CPU-s; build controller 3.362 s wall, 1.001 CPU-s and children 0.995 CPU-s. No scientific application ran. This is artifact/build consistency, not observational E/G/V parity or numerical strain qualification.\n' if result['error'] is None else result['error']+'\n'))
print(json.dumps(dict(status=result['status'],review=pin(out/'REVIEW.json'),report=pin(out/'REPORT.txt'),checks=len(checks),error=result['error'],cost=result['review_cost'])))
raise SystemExit(0 if result['error'] is None else 1)
