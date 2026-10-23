"""Bounded saved joins/arithmetic, reusing admitted full-gradient token parsing."""
from pathlib import Path
from decimal import Decimal,localcontext
import json,hashlib,time,resource
D=Path(__file__).resolve().parents[1];root=D
while root.name!='diagnostics':root=root.parent
W=root.parent/'workspaces'/D.relative_to(root);out=W/'independent_actual_review_v1'
tick=time.monotonic();cpu=time.process_time();checks=[];pins=[]
def pin(p):p=Path(p).resolve();return dict(path=str(p),sha256=hashlib.sha256(p.read_bytes()).hexdigest())
def check(x,s):
 if not x:raise ValueError(s)
 checks.append(s)
def verify(x):check(pin(x['path'])==x,'Pin '+x['path']);pins.append(x);return Path(x['path'])
def read(x):return json.loads(verify(x).read_text())
def dump(p,x):p.write_text(json.dumps(x,indent=2)+'\n')
result=dict(status='failed_term_attribution_actual_review',reviewer_implementation=pin(__file__),error=None,new_native_calls=0,new_FFTs=0,new_factorizations=0,full_gradient_reparses=0,physical_accuracy_qualified=False,original_numerical_gates_overridden=False)
try:
 rp=pin(W/'collector_v1/run_v1/execution/RESULT.json');r=read(rp);cfg=read(r['configuration']);rel=read(r['release']);sr=read(r['source_review']);inp=read(r['inputs']);m=read(r['manifest']);decl=read(m['root_declaration'])
 check(r['status']=='term_attribution_complete_original_strain_failure_preserved' and r['complete'] and r['error'] is None and r['failures']==[],'Complete finite collector')
 check(r['counts']==dict(expected_gradients=5,admitted_gradients=5,terms_per_state=28,total_checks_completed=3),'Exact5/28/3 denominator')
 check(r['implementation']==cfg['implementation']==rel['implementation']==sr['implementation'] and sr['status']=='passed_term_attribution_collector_source_review' and sr['inputs']==r['inputs']==cfg['inputs'],'Frozen collector source/config/review');verify(r['implementation'])
 check(rel['configuration']==r['configuration'] and rel['source_review']==r['source_review'] and rel['status']=='released_term_attribution_collection','Root collection release')
 check(cfg['manifest']==inp['manifest']==r['manifest'] and list(r['rows'])==decl['states']==list(m['rows']),'Exact five declared states')
 build=read(r['build']);br=read(r['build_root_review']);check(br['result']==r['build'] and br['status']=='admitted_cutoff12_term_observer_build' and build['manifest']==r['manifest'],'Admitted actual observer build')
 old=read(decl['cutoff12_collection']);check(r['attribution']['total_original_checks']==old['directions']['strain_1']['cutoff12_original_checks'],'Original failed gates preserved bytewise in saved JSON')
 failed=r['attribution']['total_original_checks'];check(not failed['passed'] and all(not x['check']['passed'] for x in failed['steps']) and not failed['two_step']['passed'],'All original three strain gates remain failed')
 term_energies={};central_vir=None
 for sid,row in r['rows'].items():
  rec=read(row['result']);base=read(row['baseline_result']);prepared=m['rows'][sid]
  check(row['result']==cfg['results'][sid] and rec['status']=='stage2_gradient_complete' and rec['error'] is None and rec['state_id']==sid,'Terminal row '+sid)
  check(rec['manifest']==r['manifest'] and rec['build']==r['build'] and rec['source_manifest']==m['native_sources'] and rec['inputs']==base['inputs']==prepared['inputs'],'Exact state/source/input '+sid)
  check(rec['gradient_output']==row['gradient_output'] and rec['term_output']==row['term_output'] and row['baseline_gradient']==prepared['baseline_gradient']==base['gradient_output'],'Exact raw receipt joins '+sid)
  for x in (row['gradient_output'],row['baseline_gradient'],row['term_output'],rec['configuration'],rec['release'],*rec['inputs'].values()):verify(x)
  check(rec['native_attempts']==1 and rec['counts']==dict(native_programs=1,gradient_calls=1,energy_only_calls=0,extra_operator_applications=0),'One native gradient '+sid)
  c=rec['commands']['native'];check(c['returncode']==0 and c['argv']==[build['executables']['gradient_old']['path'],'system.xyz'],'Actual executable '+sid)
  for x in (c['stdout'],c['stderr']):verify(x)
  ck=row['checks'];check(ck==rec['observer_checks'] and ck['exact_total_token_parity'] and ck['virial_closure'] and ck['post_term_virial_zero'] and ck['energy_sum_error']<=ck['energy_sum_bound'],'Admitted exact parity/closure '+sid)
  # Only the 35-line term trace is independently read here, not full gradients.
  lines=Path(row['term_output']['path']).read_text().splitlines();check(len(lines)==35,'Finite term trace '+sid)
  rows=[l.split() for l in lines[3:33]];check(len(rows)==30 and rows[0][1]=='INITIAL' and rows[-1][1]=='FINAL','Ordered endpoint trace '+sid)
  term_energies[sid]={v[1]:v[3] for v in rows[1:-1]}
  check(list(term_energies[sid])==list(r['attribution']['terms']),'All28 term order '+sid)
  if sid=='central':central_vir={rows[i][1]:(rows[i-1][4:],rows[i][4:]) for i in range(1,29)}
 with localcontext() as ctx:
  ctx.prec=70
  terms=r['attribution']['terms'];summary={}
  for name,t in terms.items():
   before,after=central_vir[name];a=Decimal(after[4])-Decimal(before[4]);check(a==Decimal(t['analytic_decimal']),'Exact yy term increment '+name)
   check(t['central_energy_token']==term_energies['central'][name],'Central term energy '+name)
   fds=[]
   for j,(h,label) in enumerate(((Decimal('.0001'),'1e-04'),(Decimal('.00005'),'5e-05'))):
    lo=term_energies['strain_1__'+label+'__minus'][name];hi=term_energies['strain_1__'+label+'__plus'][name];fd=(Decimal(hi)-Decimal(lo))/(2*h);s=t['steps'][j]
    check(s['minus_energy_token']==lo and s['plus_energy_token']==hi and Decimal(s['FD_decimal'])==fd and s['signed_FD_minus_analytic']==float(fd-a),'Exact saved slope/residual '+name+' '+label);fds.append(fd)
   check(t['two_step_half_minus_full']==float(fds[1]-fds[0]),'Exact step change '+name)
   if t['enabled']:summary[name]=dict(residuals=[float(x-a) for x in fds],step_change=float(fds[1]-fds[0]))
  capture=read(pin(out/'CENTRAL_SUM.json'));check(capture['central_gradient']==r['rows']['central']['gradient_output'] and capture['rows']==11958 and capture['result']==rp,'Existing one-pass central sum authority')
  directions=read(m['directions']);direction=next(x for x in directions['directions'] if x['id']=='strain_1')
  check(direction['anchor_A'][1]==1.8185 and direction['L0_A']==1. and direction['strain_matrix']==[[0.,0.,0.],[0.,1.,0.],[0.,0.,0.]],'Frozen anchored yy direction')
  correction=-Decimal(str(direction['anchor_A'][1]))*Decimal(capture['gradient_sum_tokens'][1])/Decimal(str(direction['L0_A']))
  check(correction==Decimal(capture['correction_token']),'Exact anchored translation correction')
  analytic=sum(Decimal(t['analytic_decimal']) for t in terms.values());diagnostics=[]
  for s in failed['steps']:
   h=Decimal(str(s['half_step_A']));fd=(Decimal(s['energy_plus_token'])-Decimal(s['energy_minus_token']))/(2*h);residual=fd-analytic
   diagnostics.append(dict(half_step_A=float(h),original_residual_decimal=str(residual),anchor_correction_decimal=str(correction),after_anchor_subtraction_decimal=str(residual-correction),after_anchor_subtraction=float(residual-correction),diagnostic_only=True,original_gate_unchanged=True))
  anchor=dict(capture=pin(out/'CENTRAL_SUM.json'),directions=m['directions'],sum_gradient_tokens=capture['gradient_sum_tokens'],correction_decimal=str(correction),steps=diagnostics,truncation_ratio_half_over_full=float(Decimal(diagnostics[1]['after_anchor_subtraction_decimal'])/Decimal(diagnostics[0]['after_anchor_subtraction_decimal'])),scope='Saved chain-rule diagnosis; original failed gates retained; no new qualified test or changed tolerance')
 # Saved scheduler receipt, parent allocations only.
 text=verify(r['scheduler_accounting']).read_text();parents={}
 for l in text.splitlines():
  f=l.split('|')
  if f and '.' not in f[0]:parents[f[0]]=f
 check(set(parents)=={'20513','20514','20515','20516','20517'} and all(x[2:4]==['COMPLETED','0:0'] and x[5]=='1' for x in parents.values()),'Five completed one-CPU parent allocations')
 allocated=sum(int(x[4])*int(x[5]) for x in parents.values());check(allocated==r['new_gradient_scheduler_cost']['allocated_CPU_seconds']==190,'190 allocated CPU-s without batch doublecount')
 source_dir=Path(read(read(pin(D/'INPUTS.json'))['original_source_manifest'])['source_root'])/'source'
 source_pins={name:pin(source_dir/name) for name in ('empole1.f','epolar1.f','pmestuf.f','switch.f','ehal1.f')}
 source_pins['geometry_preparation']=pin(D.parents[1]/'displaced_preparation_v1/PREPARE.py') if False else pin(D.parent.parent/'displaced_preparation_v1/PREPARE.py')
 result.update(status='passed_term_attribution_actual_review',result=rp,implementation=r['implementation'],configuration=r['configuration'],manifest=r['manifest'],source_review=r['source_review'],build=r['build'],all_five_exact_total_token_parities=True,term_count=28,original_three_gates_failed=True,independent_term_arithmetic_passed=True,term_summary=summary,anchor_diagnostic=anchor,source_pins=source_pins,scheduler_allocated_CPU_seconds=allocated,build_cost_separate=r['build_cost_separate'],collector_cost_separate=r['collector_cost'])
except Exception as exc:result['error']=repr(exc)
result.update(check_count=len(checks),checks=checks,pins=pins,review_cost=dict(wall_seconds=time.monotonic()-tick,CPU_seconds=time.process_time()-cpu,threads=1))
dump(out/'REVIEW.json',result);print(json.dumps(dict(status=result['status'],review=pin(out/'REVIEW.json'),checks=len(checks),error=result['error'],anchor=result.get('anchor_diagnostic'),cost=result['review_cost'])))
raise SystemExit(0 if result['error'] is None else 1)
