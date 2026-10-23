"""Saved five-state term attribution; original total gates and exact token parity."""
from pathlib import Path
import sys,json,time,os
from decimal import Decimal,localcontext
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from BINDING import pin,verify,read,save,load,require,THREADS,IDS,validate_manifest
HERE=Path(__file__).resolve().parent

def decomposition(helper,direction,parsed,traces):
 with localcontext() as ctx:
  ctx.prec=70
  total_a,total_print=helper.projection(parsed['central'],direction,None)
  total=helper.slope_checks('strain_1',total_a,total_print,{h:tuple(parsed['strain_1__'+label+'__'+s]['energy_token'] for s in ('minus','plus')) for h,label in zip(helper.STEPS,('1e-04','5e-05'))})
  terms={};central=traces['central']['rows']
  def proj(row):return helper.projection(dict(virial_tokens=row['virial_tokens']),direction,None)
  for i in range(1,29):
   before,bprint=proj(central[i-1]);after,aprint=proj(central[i]);a=after-before;steps=[]
   for h,label in zip(helper.STEPS,('1e-04','5e-05')):
    minus,plus=[traces['strain_1__'+label+'__'+s]['rows'][i]['energy_token'] for s in ('minus','plus')]
    fd=(Decimal(plus)-Decimal(minus))/(2*h);printing=(helper.rounding(plus)+helper.rounding(minus))/(2*h)+bprint+aprint
    steps.append(dict(half_step_A=float(h),FD_kcal_mol_A=float(fd),FD_decimal=str(fd),signed_FD_minus_analytic=float(fd-a),literal_print_allowance=float(printing),minus_energy_token=minus,plus_energy_token=plus))
   terms[central[i]['label']]=dict(enabled=central[i]['enabled'],central_energy_token=central[i]['energy_token'],analytic_kcal_mol_A=float(a),analytic_decimal=str(a),steps=steps,two_step_half_minus_full=float(Decimal(steps[1]['FD_decimal'])-Decimal(steps[0]['FD_decimal'])),scope='Descriptive decomposition; no new per-term physical gate')
  closure=[]
  require(sum(Decimal(x['analytic_decimal']) for x in terms.values())==total_a,'Exact central virial projection closure')
  for k,(h,label) in enumerate(zip(helper.STEPS,('1e-04','5e-05'))):
   termsum=sum(Decimal(x['steps'][k]['FD_decimal']) for x in terms.values())
   minus,plus=[parsed['strain_1__'+label+'__'+s]['energy_token'] for s in ('minus','plus')]
   direct=(Decimal(plus)-Decimal(minus))/(2*h);bound=Decimal(0)
   for sign in ('minus','plus'):
    rows=traces['strain_1__'+label+'__'+sign]['rows'];values=[Decimal(r['energy_token']) for r in rows[1:-1]]
    bound+=(helper.arithmetic(28,sum(map(abs,values)),Decimal(rows[-1]['energy_token']))+sum(helper.rounding(r['energy_token']) for r in rows[1:]))/(2*h)
   require(abs(termsum-direct)<=bound,'Component FD arithmetic closure')
   closure.append(dict(half_step_A=float(h),error=float(abs(termsum-direct)),arithmetic_print_bound=float(bound)))
  return dict(total_original_checks=total,terms=terms,FD_component_closure=closure,central_virial_projection_closure=True,units='q in angstrom, epsilon=q/L0, L0=1 angstrom; dE/dq in kcal/mol/angstrom',dimensionless_strain_units='dE/depsilon in kcal/mol; same numerical value since L0=1 angstrom')

def run(config,release,out):
 out=Path(out);out.mkdir(parents=True,exist_ok=False);tick,cpu=time.monotonic(),time.process_time()
 result=dict(status='term_attribution_collection_incomplete',implementation=pin(__file__),configuration=pin(config),release=pin(release),complete=False,rows={},failures=[],error=None,new_native_calls=0,new_FFTs=0,new_factorizations=0,new_scheduler_queries=0,physical_accuracy_qualified=False,active_response_qualified=False)
 try:
  cfg=json.loads(Path(config).read_text());rel=json.loads(Path(release).read_text());inp=read(cfg['inputs'])
  require(cfg['status']=='bound_term_attribution_collection' and cfg['implementation']==pin(__file__) and cfg['inputs']==pin(HERE/'INPUTS.json'),'Bound collector source')
  require(rel['status']=='released_term_attribution_collection' and rel['implementation']==pin(__file__) and rel['configuration']==pin(config) and Path(rel['output']).resolve()==out.resolve(),'Exact root collection release')
  review=read(rel['source_review']);require(review['status']=='passed_term_attribution_collector_source_review' and review['implementation']==pin(__file__) and review['inputs']==cfg['inputs'],'Collector source review')
  require(all(os.environ.get(k)=='1' for k in THREADS),'One numerical thread')
  require(cfg['manifest']==inp['manifest'],'Prepared observer states');m=read(cfg['manifest']);validate_manifest(m)
  helper=load(inp['numerical_helper'],'term_original_numerics');prior=load(inp['receipt_helper'],'term_original_receipts');observer=load(m['observer'],'term_observer_reader')
  require(inp['numerical_helper']==m['collector'] and inp['numerical_helper']['sha256']=='66c6259d428511b6b399d30aae4dfb9c65ef470135a1324812dc5c0dbc5218a8','Unchanged numerical methods')
  require(float(helper.GATE)==.001 and list(map(float,helper.STEPS))==[1e-4,5e-5],'Original derivative gates')
  build=read(cfg['build']);br=read(cfg['build_root_review'])
  require(build['status']=='stage2_build_complete' and build['error'] is None and build['manifest']==cfg['manifest'] and build['source_manifest']==m['native_sources'] and build['compile_attempts']==build['link_attempts']==1 and build['native_programs']==0,'One isolated observer compile/link')
  require(br['status']=='admitted_cutoff12_term_observer_build' and br['result']==cfg['build'] and br['manifest']==cfg['manifest'],'Root build admission')
  verify(build['executables']['gradient_old']);verify(build['library'])
  subs=read(cfg['submissions']);require(len(subs)==5 and tuple(x['state_id'] for x in subs)==IDS and set(cfg['results'])==set(IDS),'Exact five submissions/results')
  parsed={};traces={}
  for sid,sub in zip(IDS,subs):
   try:
    rec=read(cfg['results'][sid]);row=m['rows'][sid];prior.receipt_identity(rec,sid,row,cfg,m,sub,build)
    require(Path(cfg['results'][sid]['path']).resolve()==Path(row['output']).resolve()/'RESULT.json','Immutable state output')
    nc=read(rec['configuration']);nr=read(rec['release'])
    require(nc['status']=='bound_term_observer_native_configuration' and nc['state_id']==sid and nc['operation']=='gradient' and nc['manifest']==cfg['manifest'] and nc['build']==cfg['build'] and nc['build_root_review']==cfg['build_root_review'],'Native configuration')
    require(nr['status']=='released_term_observer_native' and nr['configuration']==rec['configuration'] and nr['state_id']==sid and nr['maximum_native_programs']==1 and Path(nr['output']).resolve()==Path(row['output']).resolve(),'Native release')
    ns=read(nr['source_review']);require(ns['status']=='passed_cutoff12_term_observer_source_review' and ns['manifest']==cfg['manifest'] and ns['binding']==m['execution_sources']['BINDING.py'] and ns['implementations']['native']==rec['implementation'],'Native source admission')
    for k in ('stdout','stderr'):verify(rec['commands']['native'][k])
    raw=helper.parse_gradient(verify(rec['gradient_output']));old=helper.parse_gradient(verify(row['baseline_gradient']));trace=observer.parse(verify(rec['term_output']))
    require(raw['state_id']==sid and raw['role']=='BYPASS_OLD_ABSENT' and raw['input_binding']==rec['configuration']['sha256'],'Raw output state/config')
    checks=observer.checks(trace,raw,old,helper);require(checks==rec['observer_checks'],'Actual saved per-state observer checks')
    parsed[sid]=raw;traces[sid]=trace
    result['rows'][sid]=dict(status='admitted_exact_term_observation',result=cfg['results'][sid],gradient_output=rec['gradient_output'],term_output=rec['term_output'],baseline_result=row['baseline_result'],baseline_gradient=row['baseline_gradient'],checks=checks,commands=rec['commands'],cost=rec['cost'])
   except Exception as exc:result['failures'].append(dict(state_id=sid,error=repr(exc)))
  result['counts']=dict(expected_gradients=5,admitted_gradients=len(parsed),terms_per_state=28,total_checks_completed=0)
  if len(parsed)==5:
   dirs=read(m['directions']);d=next(x for x in dirs['directions'] if x['id']=='strain_1')
   require(d['kind']=='strain' and d['L0_A']==1 and d['strain_matrix']==[[0.,0.,0.],[0.,1.,0.],[0.,0.,0.]],'Exact yy/L0 projection')
   result['attribution']=decomposition(helper,d,parsed,traces);result['counts']['total_checks_completed']=3
   old=read(read(m['root_declaration'])['cutoff12_collection'])['directions']['strain_1']['cutoff12_original_checks']
   require(result['attribution']['total_original_checks']==old,'All original failed total strain checks reproduced exactly')
  result['scheduler_accounting']=cfg['scheduler_accounting'];result['new_gradient_scheduler_cost']=scheduler(verify(cfg['scheduler_accounting']).read_text(),subs)
  verify(cfg['build_scheduler_accounting']);result['build_cost_separate']=dict(build=cfg['build'],cost=build['cost'],scheduler_accounting=cfg['build_scheduler_accounting'])
  result.update(manifest=cfg['manifest'],build=cfg['build'],build_root_review=cfg['build_root_review'],source_review=rel['source_review'],inputs=cfg['inputs'],historical_failures_preserved=True)
  result['complete']=len(parsed)==5 and not result['failures'] and result['counts']['total_checks_completed']==3
  result['status']='term_attribution_complete_original_strain_failure_preserved' if result['complete'] else 'term_attribution_collection_incomplete'
 except Exception as exc:result['error']=repr(exc)
 result['collector_cost']=dict(wall_seconds=time.monotonic()-tick,CPU_seconds=time.process_time()-cpu,threads=1)
 save(out/'RESULT.json',result);return 0 if result['complete'] else 1

def scheduler(text,submissions):
    expected={str(x['job_id']) for x in submissions};require(len(expected)==5 and len(submissions)==5,'Five distinct submitted allocations including center')
    records=[];parents={}
    for line in text.splitlines():
        if not line.strip():continue
        f=line.split('|');require(len(f)==7,'Seven-column sacct schema')
        job,name,state,code,seconds,cpus,totalcpu=f;root=job.split('.')[0]
        require(root in expected,'No unrelated scheduler records')
        rec=dict(job_id=job,job_name=name,state=state,exit_code=code,elapsed_seconds=int(seconds),allocated_CPUs=int(cpus),allocated_CPU_seconds=int(seconds)*int(cpus),measured_TotalCPU=totalcpu)
        require(rec['elapsed_seconds']>=0 and rec['allocated_CPUs']>=1,'Finite allocation');records.append(rec)
        if job==root:require(job not in parents,'Duplicate parent');parents[job]=rec
    require(set(parents)==expected,'All five parent records')
    require(all(x['state']=='COMPLETED' and x['exit_code']=='0:0' and x['allocated_CPUs']==1 for x in parents.values()),'Five successful one-CPU allocations')
    return {'records':records,'parent_jobs':sorted(parents),'allocation_count':5,'allocated_CPU_seconds':sum(x['allocated_CPU_seconds'] for x in parents.values()),'batch_or_extern_not_added':True,'scope':'Five new term-observer gradient programs including central; build and all historical baseline costs separate.'}

if __name__=='__main__':raise SystemExit(run(*sys.argv[1:]))
