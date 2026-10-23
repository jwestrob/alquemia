"""Source and existing metadata review only; no binder or runtime execution."""
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
def check(n,v):
 assert v,n
 CHECKS.append({'label':n,'passed':True})
def main():
 wall=time.monotonic();cpu=time.process_time();rp=pin(P/'READY.json');r=read(rp)
 check('Exact requested frozen READY',rp['sha256']=='d5d08408c460dd18527aa303cc0252f286de388030e845a5cbb04aefdbc49266')
 for v in r.values():
  if isinstance(v,dict) and set(v)=={'path','sha256'}:verify(v)
 t=read(r['required_review_template']);x=read(r['inputs']);d=read(r['declaration']);p=read(r['probe_addendum']);c=read(r['controls']);req=read(r['request_template'])
 for q in x['sources'].values():verify(q)
 check('Exact required source/input contract',t['implementation']==r['implementation']==x['implementation']==req['implementation'] and t['inputs']==r['inputs']==req['inputs'] and t['controls']==r['controls'])
 check('20 saved controls passed without scientific calls',c['status']=='passed_central_initializer_binding_source_controls' and c['error'] is None and c['check_count']==len(c['checks'])==20 and c['inputs']==r['inputs'] and c['implementation']==r['control_source'] and all(v==0 for v in c['counts'].values()))
 check('Actual failed/unfinished build and mismatched bindings rejected by saved controls',{'started build rejected','failed build rejected','different job rejected','missing role rejected','wrong root result join rejected','old initializer reuse rejected','release remains unissued'}.issubset(c['checks']))
 m=read(x['central_manifest']);a=read(x['central_input_admission']);ready=read(x['central_ready']);rr=read(x['runtime_root_admission']);ir=read(x['runtime_source_review']);runtime=read(x['runtime_ready'])
 check('Exact central admission manifest and declared terminal candidate',a['status']=='admitted_periodic_v2_central_inputs' and a['manifest']==x['central_manifest'] and a['ready']==x['central_ready'] and a['candidate']==ready['candidate']==m['model']['candidate']==d['candidate'])
 check('Same admitted runtime/native source bundle',rr['ready']==ir['ready']==x['runtime_ready'] and rr['independent_review']==x['runtime_source_review'] and ir['execution_sources']==runtime['execution_sources']==m['execution_sources'] and ir['native_sources']==runtime['native_sources']==m['native_sources'])
 check('Exact NATIVE implementation and settings',runtime['implementations']['native']==x['native_implementation']==r['native_implementation'] and m['model']['settings']==d['settings'])
 check('Approved immutable fixed-probe addendum',p['parent_declaration']==x['declaration']==r['declaration'] and p['status']=='declared_fixed_probe_metadata_binding' and x['probe_addendum']==r['probe_addendum'])
 previous=read(p['source_manifest']);history=read(p['historical_admission']);old=read(p['historical_result'])
 check('Historical evidence authorizes source/data only',history['result']==p['historical_result'] and old['preparation']==p['source_manifest'] and previous['probe_dimension']==p['dimension']==35874 and previous['labels']==p['labels'] and len(p['labels'])==10 and previous['probe_source']==p['probe_source'] and previous['raw_inputs']['probes.dat']==p['probe_data'])
 for name in ('probe_source','probe_data'):verify(p[name])
 cfg=read(x['native_config_template']);release=read(x['native_release_template'])
 check('Compatible producer configuration and no phase prerequisites',cfg['implementation']==x['native_implementation'] and cfg['status']=='unbound_periodic_v2_native_configuration' and all(cfg[k] is None for k in ('certificate','certificate_admission','central_parity_admission','central_scout_review','probe_source','probe_data')))
 check('One native release template remains unissued',release['status']=='unreleased_periodic_v2_native' and release['implementation']==x['native_implementation'] and release['maximum_native_programs']==1)
 check('Actual build/result chain remains unbound',all(r[k] is None for k in ('actual_build','actual_build_admission','bound_runtime_manifest','native_release')) and all(req[k] is None for k in ('build','build_admission','binder_source_review','initializer','certificate','central_parity_admission','central_scout_review','native_output')))
 sub=read(x['build_submission']);check('Declared20739 submission only; no live result read',sub['job_id']=='20739' and x['build_submission']==d['parents']['build_submission'])
 env=read(m['environment']);verify(x['interpreter'])
 check('Exact admitted interpreter path',x['interpreter_path']==str(Path(env['prefix'])/'bin/python') and pin(x['interpreter_path'])==x['interpreter'])
 wrapper=verify(r['wrapper']).read_text()
 check('OneCPU default resources no-requeue wrapper',all(z in wrapper for z in ('--partition=standard\n','--qos=standard_expand\n','--account=banfield\n','--nodes=1\n','--ntasks=1\n','--cpus-per-task=1\n','--no-requeue\n',x['interpreter_path'])) and '#SBATCH --mem' not in wrapper and '#SBATCH --time' not in wrapper)
 check('No automatic downstream operation in binder source',"state_id='central',operation='initialize'" in verify(r['implementation']).read_text() and not any(z in verify(r['implementation']).read_text() for z in ('import subprocess','import numpy','import scipy','sbatch(')))
 report='Independent central initializer binder source review: PASS\n\nThe binder requires the actual successful20739 four-role build and exact root build admission; it cannot substitute pending/failed metadata. It preserves admitted terminal5d992 central inputs, model and cutoff12 settings while creating a separate immutable runtime manifest. The approved ten-probe source/data pins enter that manifest before initialization; no old state, sparse matrix, certificate or parity result is reused. Later per-phase results stay in configurations.\n\nThe produced NATIVE configuration fixes operation initialize/state central and leaves certificate/parity/scout/probe execution fields null. The release template remains unreleased and points to the already admitted runtime source review. NATIVE retains its exact raw build/runtime/state checks and one-call/direct-only denominator. The oneCPU/default-resource wrapper uses the admitted interpreter and single-thread environment. No downstream operation or automatic release exists.\n\nRead20 passed saved controls without rerun. No actual build result, binder execution, native call, probe generation, FFT, assembly, factor or solver was performed. Initialization completion covers header/PME readback; full body/h/H validation remains CERTIFY before later active execution.\n'
 (HERE/'REPORT.txt').write_text(report)
 result=dict(t);result.update(ready=rp,report=pin(HERE/'REPORT.txt'),reviewer_implementation=pin(__file__),created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),checks=CHECKS,verified_pins=list(PINS.values()),findings=[],source_only=True,controls_rerun=False,actual_build_read=False,binder_executed=False,native_admission=False,counts={'saved_control_checks_read':20,'new_control_runs':0,'native_calls':0,'builds':0,'probe_generation':0,'scientific_numerical_calls':0},cost={'wall_seconds':time.monotonic()-wall,'CPU_seconds':time.process_time()-cpu,'threads':1})
 with (HERE/'REVIEW.json').open('x') as f:json.dump(result,f,indent=2,allow_nan=False);f.write('\n')
 print(json.dumps({'review':pin(HERE/'REVIEW.json'),'checks':len(CHECKS),'pins':len(PINS),'cost':result['cost']}))
if __name__=='__main__':main()
