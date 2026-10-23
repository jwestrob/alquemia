"""Saved12 endpoints+one admitted center; raw/corrected original9 gates side by side."""
from BINDING import *
from NUMERICS import corrected_projection,checks,STEPS
import argparse

def run(config,release,out):
 out=Path(out);out.mkdir(parents=True,exist_ok=False);tick,cpu=time.monotonic(),time.process_time()
 result=dict(status='anchored_strain_collection_incomplete',implementation=pin(__file__),configuration=pin(config),release=pin(release),complete=False,rows={},directions={},failures=[],error=None,new_native_calls=0,new_FFTs=0,new_factorizations=0,new_scheduler_queries=0,physical_accuracy_qualified=False,active_response_qualified=False)
 try:
  cfg=json.loads(Path(config).read_text());rel=json.loads(Path(release).read_text());inp=read(cfg['inputs'])
  require(cfg['status']=='bound_anchored_strain_collection' and cfg['implementation']==pin(__file__) and cfg['inputs']==pin(HERE/'COLLECT_INPUTS.json'),'Bound collector source')
  require(rel['status']=='released_anchored_strain_collection' and rel['implementation']==pin(__file__) and rel['configuration']==pin(config) and Path(rel['output']).resolve()==out.resolve(),'Exact root collection release')
  review=read(rel['source_review']);require(review['status']=='passed_anchored_strain_refinement_source_review' and review['manifest']==cfg['manifest'] and review['collector']==pin(__file__) and review['numerics']==pin(HERE/'NUMERICS.py'),'Exact collector review')
  require(all(os.environ.get(k)=='1' for k in THREADS),'One numerical thread')
  require(inp['manifest']==cfg['manifest'],'Exact prepared panel');m=read(cfg['manifest']);validate_manifest(m)
  require(inp['numerics']==pin(HERE/'NUMERICS.py') and inp['numerical_helper']==m['collector'],'Pinned arithmetic')
  helper=load(inp['numerical_helper'],'refinement_original_numerics');prior=load(inp['receipt_helper'],'refinement_original_receipts')
  require(helper.GATE==__import__('decimal').Decimal('.001') and tuple(map(float,STEPS))==tuple(read(m['root_declaration'])['half_steps_A']),'Original gate and declared steps')
  require(cfg['build']==m['build'] and cfg['build_root_review']==m['build_root_review'],'Reused build lineage');build=native_build(cfg,m)
  center_rec=read(m['central_result']);central=helper.parse_gradient(verify(m['central_gradient']));require(central['state_id']=='central' and central['role']=='BYPASS_OLD_ABSENT' and central['input_binding']==center_rec['configuration']['sha256'],'Exact old center readback')
  subs=read(cfg['submissions']);require(len(subs)==12 and tuple(x['state_id'] for x in subs)==IDS and tuple(cfg['results'])==IDS,'All12 fixed slots')
  parsed={}
  for sid,sub in zip(IDS,subs):
   try:
    row=m['rows'][sid];rec=read(cfg['results'][sid]);prior.receipt_identity(rec,sid,row,cfg,m,sub,build)
    require(Path(cfg['results'][sid]['path']).resolve()==Path(row['output']).resolve()/'RESULT.json','Immutable endpoint result')
    nc=read(rec['configuration']);nr=read(rec['release']);ns=read(nr['source_review'])
    require(nc['status']=='bound_anchored_strain_refinement_configuration' and nc['state_id']==sid and nc['operation']=='gradient' and nc['manifest']==cfg['manifest'] and nc['build']==m['build'] and nc['build_root_review']==m['build_root_review'],'Exact endpoint configuration')
    require(nr['status']=='released_anchored_strain_refinement' and nr['configuration']==rec['configuration'] and nr['state_id']==sid and nr['maximum_native_programs']==1 and Path(nr['output']).resolve()==Path(row['output']).resolve(),'Exact endpoint root release')
    require(ns['status']=='passed_anchored_strain_refinement_source_review' and ns['manifest']==cfg['manifest'] and ns['binding']==m['execution_sources']['BINDING.py'] and ns['implementation']==rec['implementation'],'Native source review')
    for x in rec['inputs'].values():verify(x)
    for k in ('stdout','stderr'):verify(rec['commands']['native'][k])
    data=helper.parse_gradient(verify(rec['gradient_output']));require(data['state_id']==sid and data['role']=='BYPASS_OLD_ABSENT' and data['input_binding']==rec['configuration']['sha256'],'Native output state/config')
    parsed[sid]=data;result['rows'][sid]=dict(result=cfg['results'][sid],gradient_output=rec['gradient_output'],energy_token=data['energy_token'],polarization_token=data['polarization_token'],commands=rec['commands'],cost=rec['cost'],job_id=sub['job_id'])
   except Exception as exc:result['failures'].append(dict(state_id=sid,error=repr(exc)))
  completed=0
  for d in read(m['directions'])['directions']:
   name=d['id']
   if name not in DIRECTIONS:continue
   needed=[name+'__'+h+'__'+s for h in STEP_LABELS for s in ('minus','plus')]
   if any(s not in parsed for s in needed):continue
   projection=corrected_projection(helper,central,d)
   energies={step:tuple(parsed[name+'__'+label+'__'+s]['energy_token'] for s in ('minus','plus')) for step,label in zip(STEPS,STEP_LABELS)}
   raw=checks(helper,name,projection['raw'],projection['raw_print'],energies)
   corrected=checks(helper,name,projection['corrected'],projection['corrected_print'],energies)
   result['directions'][name]=dict(raw_virial_only_checks=raw,corrected_anchored_checks=corrected,projection={k:([str(x) for x in v] if isinstance(v,list) else str(v)) for k,v in projection.items()},units='q in angstrom; dE/dq in kcal/mol/angstrom; L0=1angstrom',original_historical_failures_preserved=True);completed+=3
  result['counts']=dict(expected_new_gradients=12,admitted_new_gradients=len(parsed),reused_central_gradients=1,corrected_checks_completed=completed,raw_checks_completed=completed)
  result['scheduler_accounting']=cfg['scheduler_accounting'];result['new_gradient_scheduler_cost']=scheduler(verify(cfg['scheduler_accounting']).read_text(),subs)
  result.update(manifest=cfg['manifest'],build=m['build'],build_root_review=m['build_root_review'],central_result=m['central_result'],central_gradient=m['central_gradient'],source_review=rel['source_review'],root_declaration=m['root_declaration'],build_cost='Historical reused build; no new build cost')
  result['complete']=len(parsed)==12 and not result['failures'] and completed==9
  result['numerical_checks_passed']=result['complete'] and all(x['corrected_anchored_checks']['passed'] for x in result['directions'].values())
  result['status']='anchored_strain_checks_passed' if result['numerical_checks_passed'] else ('anchored_strain_numerical_failed' if result['complete'] else 'anchored_strain_collection_incomplete')
 except Exception as exc:result['error']=repr(exc)
 result['collector_cost']=dict(wall_seconds=time.monotonic()-tick,CPU_seconds=time.process_time()-cpu,threads=1);save(out/'RESULT.json',result);return 0 if result['complete'] else 1

def scheduler(text,submissions):
    expected={str(x['job_id']) for x in submissions};require(len(expected)==12 and len(submissions)==12,'Twelve distinct submitted allocations with central reused')
    records=[];parents={}
    for line in text.splitlines():
        if not line.strip():continue
        f=line.split('|');require(len(f)==7,'Seven-column sacct schema')
        job,name,state,code,seconds,cpus,totalcpu=f;root=job.split('.')[0]
        require(root in expected,'No unrelated scheduler records')
        rec=dict(job_id=job,job_name=name,state=state,exit_code=code,elapsed_seconds=int(seconds),allocated_CPUs=int(cpus),allocated_CPU_seconds=int(seconds)*int(cpus),measured_TotalCPU=totalcpu)
        require(rec['elapsed_seconds']>=0 and rec['allocated_CPUs']>=1,'Finite allocation');records.append(rec)
        if job==root:require(job not in parents,'Duplicate parent');parents[job]=rec
    require(set(parents)==expected,'All twelve parent records')
    require(all(x['state']=='COMPLETED' and x['exit_code']=='0:0' and x['allocated_CPUs']==1 for x in parents.values()),'Twelve successful one-CPU allocations')
    return {'records':records,'parent_jobs':sorted(parents),'allocation_count':12,'allocated_CPU_seconds':sum(x['allocated_CPU_seconds'] for x in parents.values()),'batch_or_extern_not_added':True,'scope':'Twelve new declared refinement gradients; central and historical builds are reused, with no new central/build cost.'}

if __name__=='__main__':raise SystemExit(run(*sys.argv[1:]))
