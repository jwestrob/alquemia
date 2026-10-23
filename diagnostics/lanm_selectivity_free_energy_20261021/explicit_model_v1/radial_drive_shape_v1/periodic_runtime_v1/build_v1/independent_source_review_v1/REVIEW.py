"""Read source metadata only. No author imports, controls, compiler or scheduler."""
import datetime,hashlib,json,time
from pathlib import Path
HERE=Path(__file__).resolve().parent;P=HERE.parent;PINS={};CHECKS=[]
def pin(p):
 p=Path(p).resolve()
 with p.open('rb') as f:h=hashlib.file_digest(f,'sha256').hexdigest()
 return {'path':str(p),'sha256':h}
def verify(p):
 if p['path'] not in PINS:assert pin(p['path'])==p,p;PINS[p['path']]=p
 else:assert PINS[p['path']]==p
 return Path(p['path'])
def read(p):return json.loads(verify(p).read_text())
def check(n,x):
 assert x,n
 CHECKS.append({'label':n,'passed':True})
def main():
 tick=time.monotonic();cpu=time.process_time();rp=pin(P/'READY.json');r=read(rp)
 check('Exact requested frozen READY',rp['sha256']=='dc2066f4616959fec66ec3f4d6212565472f47070d7f9e2aa6743c94411a8c26')
 for v in r.values():
  if isinstance(v,dict) and set(v)=={'path','sha256'}:verify(v)
 x=read(r['inputs']);h=read(r['build_handoff']);c=read(r['controls']);cfg=read(r['configuration_template']);rel=read(r['release_template']);d=read(r['declaration']);t=read(r['required_review'])
 check('Exact implementation/input/control receipt contract',t['implementation']==r['implementation']==x['implementation']==cfg['implementation'] and t['inputs']==r['inputs']==cfg['inputs'] and t['controls']==r['controls'])
 check('Runtime independent/root admissions exact',x['runtime_source_review']==t['runtime_source_review'] and x['runtime_root_admission']==t['runtime_root_admission']==r['runtime_source_admission'])
 ir=read(t['runtime_source_review']);rr=read(t['runtime_root_admission'])
 check('Admitted runtime READY and source correspondence',ir['ready']==rr['ready']==r['runtime_implementation_ready']==x['runtime_ready'] and rr['independent_review']==t['runtime_source_review'] and ir['native_sources']==h['runtime_native_manifest'])
 check('Declared build handoff and workers',x['handoff']==r['build_handoff']==d['parents']['BUILD_HANDOFF.json'] and cfg['archive_workers']==d['prospective_build']['archive_workers']==12)
 for p in x['sources'].values():verify(p)
 check('Build invocation helper is original self-contained EX COMMON',x['sources']['COMMON.py']['sha256']=='b13a4d20fbd9cc0f17db35b72244c2f9cd4d0d263b65466a7bce07e47decdb6e' and 'import POLICY' not in verify(x['sources']['COMMON.py']).read_text())
 check('30 author controls passed without calls',c['status']=='passed_periodic_v2_build_source_controls' and c['check_count']==len(c['checks'])==30 and c['error'] is None and c['compiler_calls']==c['native_calls']==c['scheduler_calls']==0 and c['inputs']==r['inputs'] and c['implementation']==r['control_source'])
 check('All required failure/finite-scope controls present',{'archive failure stops before wrappers','no retry after archive failure','missing object is terminal','attempt ceiling cannot extend','memory override rejected','no memory/time override'}.issubset(c['checks']))
 radial=read(h['radial_source_manifest']);native=read(h['runtime_native_manifest'])
 check('Complete615 inherited radial source inventory',h['source_inventory']['radial_archive_files']==radial['source_outputs'] and len(radial['source_outputs'])==615)
 check('Exact8 isolated wrapper sources',h['source_inventory']['isolated_objects']==native['files'] and len(native['files'])==len(h['compile_order'])==8)
 for p in native['files'].values():verify(p)
 check('Only initialize/gradient override entire induce object',set(h['link_objects_before_archive'])=={'initialize','gradient_clone','parity'} and all('induce_qualification.f.o' in h['link_objects_before_archive'][role] for role in ('initialize','gradient_clone')) and 'induce_qualification.f.o' not in h['link_objects_before_archive']['parity'])
 b=h['old_baseline_reuse'];ob=read(b['build_result']);oa=read(b['root_admission'])
 check('Exact old cutoff12 executable reused',b['executable']==ob['executables']['gradient_old'] and b['library']==ob['library'] and oa['result']==b['build_result'] and b['ticket_schema']=='STAGE2_BOUND_INPUT_V1')
 check('No candidate or numerical release',cfg['candidate'] is cfg['composition'] is None and cfg['status']=='unbound_periodic_v2_build_configuration' and rel['status']=='unreleased_one_periodic_v2_build' and rel['configuration'] is rel['output'] is rel['build_source_review'] is None)
 check('Exact finite build-only release maxima',rel['maximum_archive_builds']==1 and rel['maximum_wrapper_compiles']==8 and rel['maximum_wrapper_links']==3 and rel['maximum_native_programs']==0)
 check('Declared four loader checks and no old recompilation',r['planned_commands']=={'configure':1,'archive':1,'wrapper_compiles':8,'wrapper_links':3,'old_baseline_compiles_or_links':0,'loader_checks':4})
 wrapper=verify(r['wrapper']).read_text()
 check('12CPU/default-resource/no-requeue wrapper',all(z in wrapper for z in ('--partition=standard\n','--qos=standard_expand\n','--account=banfield\n','--nodes=1\n','--ntasks=1\n','--cpus-per-task=12\n','--no-requeue\n')) and '#SBATCH --mem' not in wrapper and '#SBATCH --time' not in wrapper)
 check('Loader/build overrides cleared',all(z in wrapper for z in ('unset LD_PRELOAD LD_LIBRARY_PATH FFLAGS LDFLAGS CMAKE_GENERATOR CMAKE_TOOLCHAIN_FILE','OPENBLAS_NUM_THREADS=1','OMP_NUM_THREADS=1','PYTHONDONTWRITEBYTECODE=1')))
 report='Independent periodic V2 build adapter source review: PASS\n\nThe immutable handoff fixes615 radial archive files and8 wrapper sources. BUILD adds only declared archive parallelism12; it builds target tinker with OpenMM disabled, then compiles wrappers sequentially against the new radial module directory. Initializer/gradient links place the complete qualification induce object before the new archive; parity uses the normal archive and isolated setup. The old cutoff12 V1 gradient executable is reused with its separate library/module/gate/parser lineage. No molecular executable is invoked.\n\nThe13 build commands have durable attempts and stop on the first failure; four later ldd checks establish loader records only. Slurm12CPU/default memory-time/no-requeue checks and source/tool pins precede compiler commands. Per-process memory/CPU caps and observed child CPU/RSS are distinguished from scheduler-wide accounting. Candidate and configuration/release remain unbound.\n\nRead30 passed author controls without rerun. The initially suspected POLICY import issue is resolved by the frozen INPUTS: the adapter loads original execution_v1 COMMONb13a4d20, which has no POLICY dependency. No correction is required. This is a source receipt, not compile/link success or numerical admission.\n'
 (HERE/'REPORT.txt').write_text(report)
 result=dict(t);result.update(ready=rp,declaration=r['declaration'],configuration_template=r['configuration_template'],report=pin(HERE/'REPORT.txt'),reviewer_implementation=pin(__file__),created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),checks=CHECKS,verified_pins=list(PINS.values()),findings=[],source_only=True,controls_rerun=False,compiler_execution=False,native_admission=False,counts={'saved_control_checks_read':30,'compiler_calls':0,'native_calls':0,'scheduler_calls':0,'new_control_runs':0},cost={'wall_seconds':time.monotonic()-tick,'CPU_seconds':time.process_time()-cpu,'threads':1})
 with (HERE/'REVIEW.json').open('x') as f:json.dump(result,f,indent=2,allow_nan=False);f.write('\n')
 print(json.dumps({'review':pin(HERE/'REVIEW.json'),'checks':len(CHECKS),'pins':len(PINS),'cost':result['cost']}))
if __name__=='__main__':main()
