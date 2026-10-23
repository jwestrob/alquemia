"""Finish only dependency closure/runtime/cost suffix of preserved build review."""
import json,hashlib,time,resource,re,subprocess,os
from pathlib import Path
T=Path(__file__).resolve().parents[1];ready=json.loads((T/'READY.json').read_text());w=Path(ready['manifest']['path']).parents[1]
out=w/'independent_build_review_v1/continuation_v2';out.mkdir(exist_ok=False);tick=time.monotonic();cpu=time.process_time();checks=[]
def pin(p):p=Path(p).resolve();return dict(path=str(p),sha256=hashlib.sha256(p.read_bytes()).hexdigest())
def check(x,s):
 if not x:raise ValueError(s)
 checks.append(s)
def verify(x):check(pin(x['path'])==x,'Pin '+x['path']);return Path(x['path'])
def read(x):return json.loads(verify(x).read_text())
def dump(p,x):p.write_text(json.dumps(x,indent=2)+'\n')
priorpin=pin(w/'independent_build_review_v1/REVIEW.json');prior=read(priorpin)
check(priorpin['sha256']=='ffbec5725f7649d89d5e38ce3d750d24f51590e78cd112a18259bdc43bdce8ce' and prior['check_count']==229 and prior['error']=="ValueError('ELF NEEDED exact qualified runtime libraries')",'Exact completed review prefix and narrow failed assumption')
resultpin=next(x for x in prior['pins'] if x['path']==str(w/'build_run_v1/execution/RESULT.json'));r=read(resultpin);m=read(r['manifest']);recipe=read(ready['link_recipe'])
dynamic=verify(prior['metadata_commands'][-1]['stdout']).read_text();needed=set(re.findall(r'Shared library: \[([^]]+)\]',dynamic))
check(needed<set(m['native_runtime']) and set(m['native_runtime'])-needed=={'loader','libgcc_s.so.1'},'Executable direct dependencies differ only by loader/transitive libgcc_s')
# Read the actual linked libgfortran ELF metadata; no loader/application invocation.
saved_dynamic=pin(w/'independent_build_review_v1/continuation_v1/libgfortran_dynamic.out')
saved_failure=pin(w/'independent_build_review_v1/continuation_v1/FAILURE.json')
transitive=set(re.findall(r'Shared library: \[([^]]+)\]',verify(saved_dynamic).read_text()))
loader_soname=Path(m['native_runtime']['loader']['path']).name
check('libgcc_s.so.1' in transitive and transitive<=((set(m['native_runtime'])-{'loader'})|{loader_soname}),'Qualified transitive dependencies with explicit loader SONAME alias')
check(r['native_runtime']['gradient_old']==m['native_runtime'],'Actual saved node loader closure equals manifest')
for item in r['native_runtime']['gradient_old'].values():verify(item)
a=r['allocation'];check(a['job_id']=='20503' and a['allocated_CPUs']==a['nodes']==a['tasks']==1 and a['qos']=='standard_expand' and not a['requeue'],'Actual one-CPU allocation')
check(a['memory_grant_bytes']==4194304000 and a['memory_grant_bytes']>3*1024**3,'Default grant accommodates child cap')
subpin=pin(w/'build_run_v1/SUBMISSION.json');eventpin=pin(w/'build_run_v1/completion_event.json');sub=read(subpin);event=read(eventpin);terminal=event['records']['20503']
check(sub['returncode']==0 and sub['stdout'].strip()=='20503' and terminal['state']=='COMPLETED' and terminal['exit_code']=='0:0','Terminal scheduler receipt')
check(int(terminal['allocated_cpus'])==1 and terminal['nodes']==a['hostname'],'Receipt node/allocation')
review=dict(status='passed_cutoff12_term_observer_build_review',reviewer_implementation=pin(__file__),result=resultpin,configuration=r['configuration'],implementation=r['implementation'],manifest=r['manifest'],source_manifest=r['source_manifest'],ready=pin(T/'READY.json'),preserved_failed_review=priorpin,completed_prefix_checks=229,remaining_checks=len(checks),checks=checks,exact_original_build_result_unchanged=True,executables=r['executables'],observer_object=pin(r['commands']['compile_GRADIENT_OBSERVER_f']['argv'][-1]),library=r['library'],runtime=r['native_runtime'],source_build_checks_passed=True,original_failure_not_relabelled=True,metadata_correction='Direct ELF NEEDED differs from recursive loader closure. Explicit libgfortran NEEDED proves libgcc_s transitive; actual saved loader/runtime pins exact.',direct_needed=sorted(needed),libgfortran_needed=sorted(transitive),preserved_dependency_alias_failure=saved_failure,reused_libgfortran_dynamic=saved_dynamic,build_cost=r['cost'],build_child_CPU_seconds=sum(c['cost']['child_user_CPU_seconds']+c['cost']['child_system_CPU_seconds'] for c in r['commands'].values()),scheduler=dict(completion_event=eventpin,job_id='20503',elapsed_seconds=int(terminal['elapsed_seconds']),allocated_CPUs=1,allocated_CPU_seconds=int(terminal['elapsed_seconds']),batch_not_added=True),native_calls=0,builds=0,scientific_qualification=False,review_cost=dict(wall_seconds=time.monotonic()-tick,CPU_seconds=time.process_time()-cpu,threads=1),prior_review_cost_separate=prior['review_cost'])
dump(out/'REVIEW.json',review)
(out/'REPORT.txt').write_text('PASSED narrow build review via preserved229-check prefix plus dependency/cost continuation.\nOne exact observer compile/link; unchanged three harness objects/archive; one gradient_ and term_observe_ in binary. All runtime pins match. First review failure was incorrect equality of direct ELF dependencies to recursive loader closure; libgcc_s is explicitly transitive via libgfortran. Original failed receipts remain unchanged, including the narrow loader-alias assertion. The completed suffix uses saved ELF metadata only.\nJob20503:3s x1CPU=3 allocated CPU-s. Build controller3.362s wall/1.001CPU-s; children0.995CPU-s. This review made no scientific calls. Exact five-state E/G/V parity remains a future numerical observation, not established by successful compilation.\n')
print(json.dumps(dict(review=pin(out/'REVIEW.json'),report=pin(out/'REPORT.txt'),status=review['status'],cost=review['review_cost'])))
