"""Independent source/metadata receipt only; no author imports or numerical replay."""
import datetime,hashlib,json,time
from pathlib import Path
HERE=Path(__file__).resolve().parent;P=HERE.parent;PINS={};CHECKS=[]
def pin(path):
 path=Path(path).resolve()
 with path.open('rb') as f:h=hashlib.file_digest(f,'sha256').hexdigest()
 return {'path':str(path),'sha256':h}
def verify(p):
 if p['path'] not in PINS:assert pin(p['path'])==p,p;PINS[p['path']]=p
 else:assert PINS[p['path']]==p
 return Path(p['path'])
def read(p):return json.loads(verify(p).read_text())
def digest(x):return hashlib.sha256(json.dumps(x,sort_keys=True,separators=(',',':'),allow_nan=False).encode()).hexdigest()
def fp(c):return digest({k:v for k,v in c.items() if k not in ('root_release','source_review','preparation','fingerprint')})
def check(n,x):
 assert x,n
 CHECKS.append({'label':n,'passed':True})
def main():
 wall=time.monotonic();cpu=time.process_time();rp=pin(P/'READY.json');r=read(rp)
 check('Exact requested frozen READY',rp['sha256']=='7b3d2bbb858348922f957baeff4a78a421132e94ef943e61621e802ecf5be2d6')
 m=read(r['manifest']);cfg=read(r['configuration_template']);pcfg=read(r['configuration_prepared']);prep=read(r['preparation']);w=read(r['witness']);decl=read(r['root_declaration']);rel=read(r['root_release_template'])
 t=json.loads((P/'SOURCE_REVIEW_TEMPLATE.json').read_text())
 check('Exact required receipt schema',t==r['source_review_contract'])
 for k in ('implementation','adapter','root_declaration','scheduler_wrapper'):check('Ready/config '+k,r[k]==cfg[k]==pcfg[k]);verify(r[k])
 check('Configuration fingerprint identical across preparation/runtime',fp(cfg)==fp(pcfg)==r['configuration_fingerprint']==prep['configuration_fingerprint']==t['configuration_fingerprint'])
 check('Prepared exact point is retained',prep['witness']==r['witness']==t['witness'] and prep['witness_key']==r['witness_key']==w['key'] and prep['witness_parameters']==r['parameters']==w['parameters'] and digest(w['parameters'])==w['key'])
 check('Prepared fraction and raw solver status preserved',r['fraction']==w['fraction']==decl['policy']['fraction']==.0625 and r['solver_status_retained']==prep['solver_status_retained']=='AlmostSolved' and prep['solver_success'] is False and w['solver_success'] is False)
 check('Author scalar preparation passed without execution',prep['status']=='prepared_boundary_fraction_witness_zero_calls' and prep['eligible'] is True and w['eligible'] is True and prep['error'] is None and prep['native_programs']==prep['solver_calls']==prep['new_derivative_probes']==0 and prep['old28_replayed'] is False and prep['new_cache_constructed'] is False)
 check('One finite declared native scope',r['limits']==cfg['limits']==pcfg['limits']==prep['limits'] and r['limits']['new_parameter_vectors']==1 and r['limits']['new_native_programs']==56 and r['limits']['new_cells']==46 and r['limits']['new_solver_calls']==r['limits']['new_derivative_probes']==0 and r['limits']['workers']==12 and r['limits']['allocated_CPUs']==13)
 check('Original stop and no promotion preserved',prep['original_fit_stop']=='no_actual_safeguarded_improvement' and prep['selected_candidate'] is None and prep['automatic_promotion'] is False and prep['automatic_extension'] is False)
 for p in m['sources'].values():verify(p)
 for k,p in m['unchanged_helpers'].items():check('Config exact inherited helper '+k,p==cfg[k]);verify(p)
 c1,c2=[read(p) for p in r['controls']]
 check('Public-path19 controls read only',c1['status']=='passed_boundary_fraction_manufactured_controls' and c1['checks_passed']==len(c1['checks'])==19 and c1['native_programs']==c1['solver_calls']==c1['scientific_vectors_read']==0 and c1['adapter']==r['adapter'])
 check('Control as-executed source archived by hash',m['sources']['source_history_v1/RUN_controlled_initial.py']['sha256']==c1['run']['sha256'])
 check('Recovery16 controls bind final runner',c2['status']=='passed_fixed_fraction_recovery_controls' and c2['checks_passed']==len(c2['checks'])==16 and c2['run']==r['implementation'] and c2['native_programs']==c2['solver_calls']==c2['scientific_vectors_read']==0)
 check('Final pin and failed-materialization controls present',{'final_ledger_pin_after_materialization','final_saved_result_has_authoritative_pin','unreconciled_terminal_cannot_be_admitted'}.issubset(c2['checks']))
 failure=read(r['original_failure'])
 check('Original preparation failure preserved',failure['status']=='failed_preparation' and failure['native_programs']==failure['solver_calls']==0 and m['failed_preparation']==r['original_failure'])
 check('Runtime preparation/source release still unbound',r['native_release'] is None and cfg['root_release'] is None and cfg['source_review'] is None and cfg['preparation']==r['preparation'])
 check('Root release template binds exact prepared witness',rel['witness']==r['witness'] and rel['configuration_fingerprint']==r['configuration_fingerprint'] and rel['preparation']==r['preparation'] and rel['source_review'] is None)
 wrapper=verify(r['scheduler_wrapper']).read_text()
 check('One13CPU standard/standard_expand no-requeue wrapper',all(s in wrapper for s in ('--partition=standard\n','--qos=standard_expand\n','--account=banfield\n','--cpus-per-task=13\n','--ntasks=1\n','--nodes=1\n','--no-requeue\n','RUN_v2.py run')) and '--mem' not in wrapper and '--time' not in wrapper)
 result=dict(t)
 report='Independent boundary-fraction witness source review: PASS\n\nExactly one external fraction 1/16 is formed from the pinned saved step8 return and restored affine problem. The private raw certificate keeps identical code with only its fraction allowlist replaced; the second original fraction1 call checks the auxiliary-q certificate of that same already-fractioned physical point. Original25 physical gates, fixed feca three electrical norms, signed-W9 ceiling, domain/trust and strict raw-q comparison remain unchanged. AlmostSolved remains explicit.\n\nThe public evaluator accepts only the exact prepared vector, starts a fresh cache with no imports, and enforces1vector/46cells/56programs at the inherited locked invocation boundary. Journal limits are present before initialization. Failure stops subsequent scheduling; final ledger pins follow authoritative materialization, and unavailable final reconciliation removes admission. Slurm13CPUs covers12single-thread workers plus controller; costs retain parent-scheduler accounting separately.\n\nRead19 public-path and16 recovery controls without rerun. Read the saved scalar preparation and exact joins without recalculating its certificate or old28 derivatives. The original fraction-allowlist preparation failure remains preserved. This source receipt does not admit any native result or promote/select a candidate.\n'
 (HERE/'REPORT.txt').write_text(report)
 result.update(ready=rp,manifest=r['manifest'],report=pin(HERE/'REPORT.txt'),reviewer_implementation=pin(__file__),created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),checks=CHECKS,verified_pins=list(PINS.values()),controls=r['controls'],original_failure=r['original_failure'],findings=[],source_only=True,controls_rerun=False,scalar_certificate_recomputed=False,native_candidate_admission=False,selected_candidate=None,counts={'saved_control_checks_read':35,'new_control_runs':0,'native_calls':0,'solver_calls':0,'new_derivative_probes':0,'old_derivative_replays':0},cost={'wall_seconds':time.monotonic()-wall,'CPU_seconds':time.process_time()-cpu,'threads':1})
 with (HERE/'REVIEW.json').open('x') as f:json.dump(result,f,indent=2,allow_nan=False);f.write('\n')
 print(json.dumps({'review':pin(HERE/'REVIEW.json'),'checks':len(CHECKS),'pins':len(PINS),'cost':result['cost']}))
if __name__=='__main__':main()
