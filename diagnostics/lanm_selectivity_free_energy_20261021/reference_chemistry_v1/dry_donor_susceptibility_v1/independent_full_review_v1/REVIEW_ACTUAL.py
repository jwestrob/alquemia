#!/usr/bin/env python3
"""Full10 saved-only raw response/cost review, preserving the admitted scout."""
import hashlib,importlib.util,json,math,re,time
from decimal import Decimal as D,getcontext
from pathlib import Path
getcontext().prec=60
ROOT=Path.cwd();REL=Path('reference_chemistry_v1/dry_donor_susceptibility_v1');DB=ROOT/'diagnostics/lanm_selectivity_free_energy_20261021'/REL;WB=ROOT/'workspaces/lanm_selectivity_free_energy_20261021'/REL
HELPER=DB/'independent_scout_review_v1/REVIEW_ACTUAL.py'
assert hashlib.sha256(HELPER.read_bytes()).hexdigest()=='c188bee7587c761fe8ec695f5ca5caa21b9a2ba05248cd031e5dfe60d62ebda0'
sp=importlib.util.spec_from_file_location('preserved_scout_reader',HELPER);h=importlib.util.module_from_spec(sp);sp.loader.exec_module(h)
pin,read,save,serial,ck,verify,dd,hb,one,section,array,close=h.pin,h.read,h.save,h.serial,h.ck,h.verify,h.dd,h.hb,h.one,h.section,h.array,h.close
N=h.N;COLLECTION=WB/'full_collection_v1/result/RESULT.json';OUT=WB/'independent_full_review_v1'

def main():
 if not COLLECTION.exists():
  print(json.dumps({'status':'not_terminal_review_not_run','expected_collection':str(COLLECTION),'new_native_calls':0}));return
 OUT.mkdir(exist_ok=False);tick=time.monotonic();cpu=time.process_time();r=read(COLLECTION);ready=read(DB/'READY.json');base=dict(implementation=ready['implementation'],panel=ready['panel'],collection=pin(COLLECTION),review_implementation=pin(__file__),helper=pin(HELPER),new_native_calls=0)
 try:
  ck(r['implementation']==ready['implementation'] and r['panel']==ready['panel'] and r['stage']=='full','frozen full10 source/panel')
  sp=importlib.util.spec_from_file_location('frozen_quantum10_full',verify(r['implementation']));mod=importlib.util.module_from_spec(sp);sp.loader.exec_module(mod);p=mod.load(verify(r['panel']))
  ck(r['cell_ids']==[c['cell_id'] for c in p['cells']] and r['declared']==r['complete']==10,'complete exact10')
  scoutpath=WB/'independent_scout_review_v1/REVIEW.json';scout=read(scoutpath);ck(pin(scoutpath)['sha256']=='2b396c26c6fbfa0cc83e963ac9b42fb23210b8307e655c3646e023d4eafff567' and scout['status']=='passed_dry_donor_susceptibility_scout_review' and scout['implementation']==r['implementation'] and scout['panel']==r['panel'],'preserved admitted La5 review')
  oldrows=read(verify(scout['collection']))['rows'];raw={m:{} for m in ['La','Dy']};native={};snapshots=[]
  for c in p['cells']:
   cid=c['cell_id'];metal=c['metal'];rec=r['rows'][cid];admitted=mod.admit_result(verify(rec['result']),c,p,r['panel']['path']);ck(admitted==rec,cid+'/raw native admission')
   if metal=='La':ck(rec==oldrows[cid],cid+'/exact inherited La row')
   text=verify(rec['output']).read_text();pt=verify(rec['properties']['property_file']).read_text();g=section(pt,'Geometry');ck('&Units "Bohr"' in g,cid+'/Bohr positions');xyz=re.findall(r'^\s*([A-Za-z]+)\s+('+N+r')\s+('+N+r')\s+('+N+r')\s*$',g,re.M);ck(len(xyz)==8 and [v[0] for v in xyz]==c['symbols'],cid+'/eight real sites')
   block=section(pt,'SCF_Dipole_Moment');vs={k:array(block,name,3) for k,name in [('mu','dipoleTotal'),('nuc','dipoleNucContrib'),('ele','dipoleElecContrib')]};z=[11,6,6,8,8,1,1,1]
   for k in range(3):
    predicted=sum((D(z[i])*dd(xyz[i][k+1]) for i in range(8)),D(0));bound=sum((D(abs(z[i]))*hb(xyz[i][k+1]) for i in range(8)),D(0))+hb(vs['nuc'][k]);close(predicted,dd(vs['nuc'][k]),bound,cid+'/nuclear lab origin'+str(k));close(dd(vs['mu'][k]),dd(vs['nuc'][k])+dd(vs['ele'][k]),sum((hb(vs[t][k]) for t in vs),D(0)),cid+'/dipole components'+str(k))
   origin=one(r'^Position of electric origin\s*\.\.\.\s*('+N+r')\s+('+N+r')\s+('+N+r')\s*$',text,cid+'/native origin');ck(list(map(dd,origin))==[0,0,0],cid+'/lab zero')
   et=array(section(pt,'SCF_Energy'),'totalEnergy',1)[0];ft=one(r'&FinalEnergy\s+\[&Type "Double"\]\s+('+N+r')',section(pt,'Single_Point_Data'),cid+'/property final');dt=one(r'^Dispersion correction\s+('+N+r')\s*$',text,cid+'/D4');nt=one(r'^FINAL SINGLE POINT ENERGY\s+('+N+r')\s*$',text,cid+'/native final');ns=one(r'^Total Energy\s*:\s*('+N+r')\s+Eh',text,cid+'/native electronic')
   close(dd(et),dd(ns),hb(et)+hb(ns),cid+'/electronic stage');close(dd(ft),dd(nt),hb(ft)+hb(nt),cid+'/final stage');close(dd(ft),dd(et)+dd(dt),hb(ft)+hb(et)+hb(dt),cid+'/electronic plus D4')
   mu=sum((D(str(t))*dd(v) for t,v in zip(p['field_direction_lab'],vs['mu'])),D(0));mb=sum((abs(D(str(t)))*hb(v) for t,v in zip(p['field_direction_lab'],vs['mu'])),D(0))
   raw[metal][c['field_tag']]=dict(e=dd(et),final=dd(ft),d4=dd(dt),mu=mu,mu3=[dd(t) for t in vs['mu']],eb=hb(et),fb=hb(ft),db=hb(dt),mb=mb,energy_tokens=dict(electronic=et,final=ft,D4=dt),dipole_tokens=vs,field_au=c['field_au'])
   native[cid]={k:v for k,v in rec['parsed'].items() if any(x in k.lower() for x in ['scf','grid','cosx','criteria','energy','d4'])};snapshots.append(dict(cell_id=cid,result=rec['result'],output=rec['output'],property=rec['properties']['property_file'],start=rec['start'],elapsed_seconds=rec['elapsed_seconds'],allocated_elapsed_CPU_seconds_proxy=rec['allocated_elapsed_CPU_seconds_proxy']))
  values={};gates={};d4={}
  for metal in ['La','Dy']:
   v=raw[metal];ck(set(v)=={'zero','minus_h','minus_half','plus_half','plus_h'},metal+'/complete five fields');zero=v['zero'];values[metal]={};d4[metal]=all(abs(q['d4']-zero['d4'])<=q['db']+zero['db'] for q in v.values());ck(d4[metal]==r['response']['D4_field_independence'][metal],metal+'/D4 equality')
   for tag,step in [('h',D('.001')),('half',D('.0005'))]:
    a,b=v['minus_'+tag],v['plus_'+tag];est={'chi_energy':-(b['e']+a['e']-2*zero['e'])/step**2,'chi_moment':(b['mu']-a['mu'])/(2*step),'odd_energy_moment':-(b['e']-a['e'])/(2*step),'zero_moment':zero['mu'],'even_moment_numerator':b['mu']+a['mu']-2*zero['mu']};actual=r['response']['per_metal'][metal]['estimates'][tag]
    for k in est:close(est[k],actual[k],actual['roundoff_bounds'][k],metal+'/'+tag+'/'+k)
    est['total_chi_energy']=-(b['final']+a['final']-2*zero['final'])/step**2;total=r['response']['per_metal'][metal]['final_total_estimates'][tag];close(est['total_chi_energy'],total['chi_energy'],total['roundoff_bounds']['chi_energy'],metal+'/'+tag+'/total curvature')
    est['decimal_print_bounds']={'chi_energy':(a['eb']+b['eb']+2*zero['eb'])/step**2,'chi_moment':(a['mb']+b['mb'])/(2*step),'odd_energy_moment':(a['eb']+b['eb'])/(2*step),'zero_moment':zero['mb']};values[metal][tag]=est;gates[metal+'__energy_vs_moment_'+tag]=abs(est['chi_energy']-est['chi_moment'])<=D('.5');gates[metal+'__odd_slope_'+tag]=abs(est['odd_energy_moment']-est['zero_moment'])<=D('.0001')
   for key in ['chi_energy','chi_moment']:gates[metal+'__'+key+'_step']=abs(values[metal]['h'][key]-values[metal]['half'][key])<=D('.5')
  paired={}
  for tag in ['h','half']:
   paired[tag]={k:values['Dy'][tag][k]-values['La'][tag][k] for k in ['chi_energy','chi_moment','odd_energy_moment','zero_moment','even_moment_numerator']};actual=r['response']['Dy_minus_La']['estimates'][tag]
   for k,v in paired[tag].items():close(v,actual[k],actual['roundoff_bounds'][k],'paired/'+tag+'/'+k)
   paired[tag]['decimal_print_bounds']={k:values['La'][tag]['decimal_print_bounds'][k]+values['Dy'][tag]['decimal_print_bounds'][k] for k in values['La'][tag]['decimal_print_bounds']};gates['Dy_minus_La__energy_vs_moment_'+tag]=abs(paired[tag]['chi_energy']-paired[tag]['chi_moment'])<=1
  for key in ['chi_energy','chi_moment']:gates['Dy_minus_La__'+key+'_step']=abs(paired['h'][key]-paired['half'][key])<=1
  ck(gates==r['response']['checks'] and len(gates)==16,'independent original16 gates');passed=all(gates.values()) and all(d4.values());ck(passed==r['numerical_checks_passed']==r['response']['passed'],'full producer qualification replay');ck(mod.score(r['rows'],p,'full')==r['response'],'saved full response replay')
  # Independently rederive just the two model center quartets/zeros, using the previously admitted40 authority.
  inv=read(verify(p['source_inventory']));mr=read(verify(p['model_actual_review']));model=read(verify(inv['model40_result']));mm=read(verify(model['manifest']));ck(mr['status']=='passed_dry_acetate_electrical_actual_review' and mr['result']==inv['model40_result'],'exact model40 actual review');K=D(str(mm['conversion']['HARTREE']));mvalues={};comparison={}
  def model_energy(pinrecord):return one(r'^\s*Total Potential Energy\s*:\s*('+N+r')\s+Kcal/mole',verify(pinrecord).read_text(),'model center energy')
  for metal in ['La','Dy']:
   gid=metal+'_central';zero=model_energy(mm['rows'][gid]['zero_energy']['native']['analyze_output']);mvalues[metal]={};comparison[metal]={}
   ck(mm['rows'][gid]['charge_origin']['origin_A']==[0,0,0] and dd(mm['rows'][gid]['charge_origin']['charge_e'])==2,'model charged lab origin '+metal)
   for tag,step in [('h',D('.001')),('half',D('.0005'))]:
    ep=model_energy(model['rows'][gid+'__plus_'+tag]['receipts'][0]['output']);en=model_energy(model['rows'][gid+'__minus_'+tag]['receipts'][0]['output']);chi=-(dd(ep)+dd(en)-2*dd(zero))/(step**2*K);mu=-(dd(ep)-dd(en))/(2*step*K);q=model['summary']['states'][gid]['steps'][tag];close(chi,q['chi_E_a0_cubed'],D('1e-12'),'model raw chi '+metal+'/'+tag);close(mu,q['mu_E_au'],D('1e-14'),'model raw mu '+metal+'/'+tag);mvalues[metal][tag]={'chi_energy':chi,'mu_energy':mu,'chi_print_bound':(hb(ep)+hb(en)+2*hb(zero))/(step**2*K),'mu_print_bound':(hb(ep)+hb(en))/(2*step*K)};comparison[metal][tag]={'model_chi':chi,'quantum_chi_E':values[metal][tag]['chi_energy'],'quantum_chi_mu':values[metal][tag]['chi_moment'],'model_minus_quantum_E':chi-values[metal][tag]['chi_energy'],'model_minus_quantum_mu':chi-values[metal][tag]['chi_moment'],'model_mu':mu,'quantum_zero_mu':raw[metal]['zero']['mu'],'model_minus_quantum_zero_mu':mu-raw[metal]['zero']['mu']}
    for est in ['chi_energy','chi_moment']:close(chi-values[metal][tag][est],r['model_comparison']['model_minus_quantum'][metal][tag][est],r['response']['per_metal'][metal]['estimates'][tag]['roundoff_bounds'][est]+q['chi_print_bound_a0_cubed'],'paired table '+metal+'/'+tag+'/'+est)
  comparison['paired']={}
  for tag in ['h','half']:
   v=mvalues['Dy'][tag]['chi_energy']-mvalues['La'][tag]['chi_energy'];comparison['paired'][tag]={'model_chi':v,'quantum_chi_E':paired[tag]['chi_energy'],'quantum_chi_mu':paired[tag]['chi_moment'],'model_minus_quantum_E':v-paired[tag]['chi_energy'],'model_minus_quantum_mu':v-paired[tag]['chi_moment'],'model_mu':mvalues['Dy'][tag]['mu_energy']-mvalues['La'][tag]['mu_energy'],'quantum_zero_mu':raw['Dy']['zero']['mu']-raw['La']['zero']['mu']}
   for est in ['chi_energy','chi_moment']:close(v-paired[tag][est],r['model_comparison']['paired_model_minus_quantum'][tag][est],r['response']['Dy_minus_La']['estimates'][tag]['roundoff_bounds'][est]+sum(model['summary']['states'][m+'_central']['steps'][tag]['chi_print_bound_a0_cubed'] for m in ['La','Dy']),'paired table/'+tag+'/'+est)
  # Retain stable all-job cost snapshots, including both collectors and inherited La work.
  events={};jobs={}
  for name,path in [('La_native',WB/'run_v1/completion_event.json'),('La_collector',WB/'scout_collection_v1/completion_event.json'),('Dy_native',WB/'run_v1/remaining_submission_v1/completion_event.json'),('full_collector',WB/'full_collection_v1/completion_event.json')]:
   ev=read(path);events[name]={'path':str(path),'records':ev.get('records',{})};jobs.update(ev.get('records',{}))
  expected={str(i) for i in list(range(18832,18838))+list(range(18839,18845))};ck(set(jobs)==expected,'ten native plus two collector scheduler records');ck(all(v['state']=='COMPLETED' and v['exit_code']=='0:0' for v in jobs.values()),'all12 scheduler jobs complete')
  costs={jid:dict(elapsed_seconds=int(v['elapsed_seconds']),allocated_CPUs=int(v['allocated_cpus']),allocated_CPU_seconds=int(v['elapsed_seconds'])*int(v['allocated_cpus'])) for jid,v in jobs.items()};save(OUT/'SCHEDULER_SNAPSHOT.json',events);save(OUT/'RAW_REPLAY.json',serial({'raw':raw,'Decimal_estimates':values,'paired':paired,'model_center_replay':mvalues,'comparison':comparison,'numerical_checks':gates,'D4_field_independence':d4,'native_metadata':native,'native_snapshots':snapshots}));save(OUT/'CHECKS.json',h.checks)
  result=dict(base,status='passed_dry_donor_susceptibility_full_review' if passed else 'failed_dry_donor_susceptibility_full_review',scout_review=pin(scoutpath),model_actual_review=p['model_actual_review'],model_result=inv['model40_result'],check_count=len(h.checks),checks=pin(OUT/'CHECKS.json'),raw_replay=pin(OUT/'RAW_REPLAY.json'),complete_cells=10,saved_native_admissions=10,numerical_checks=16,numerical_checks_passed=passed,D4_field_independence=d4,cost_reviewed=True,scheduler_snapshot=pin(OUT/'SCHEDULER_SNAPSHOT.json'),scheduler_costs=costs,total_scheduler_CPU_seconds=sum(v['allocated_CPU_seconds'] for v in costs.values()),native_elapsed_CPU_seconds_proxy=sum(x['allocated_elapsed_CPU_seconds_proxy'] for x in snapshots),maximum_independent_derived_difference=max(h.diffs),review_cost={'wall_seconds':time.monotonic()-tick,'CPU_seconds':time.process_time()-cpu,'threads':1},limitations=['Numerical resolutions and printing/arithmetic bounds are not physical accuracy thresholds.','Native EnergyCheck and other failed printed SCF criteria remain explicit.','Global uniform-tangent response does not isolate local induction or charge transfer.','No parameter change, molecular call, new target or promotion by this review.'])
  save(OUT/'REVIEW.json',result);report='# Independent full10 dry-donor susceptibility review\n\n'+('All ten raw admissions and16 numerical checks pass. ' if passed else 'Full10 is numerically unqualified. ')+'The La scout receipt remains unchanged. Zero native calls by this review.\n\n| Quantity | h, au | Model χE | Quantum χE | Quantum χμ | Model−quantum χE | Model−quantum χμ |\n|---|---:|---:|---:|---:|---:|---:|\n'
  for kind in ['La','Dy','paired']:
   for tag,step in [('h','.001'),('half','.0005')]:
    v=comparison[kind][tag];report+=f"| {'Dy−La' if kind=='paired' else kind} | {step} | {float(v['model_chi']):.9f} | {float(v['quantum_chi_E']):.9f} | {float(v['quantum_chi_mu']):.9f} | {float(v['model_minus_quantum_E']):.9f} | {float(v['model_minus_quantum_mu']):.9f} |\n"
  report+='\nSusceptibilities are in a0³. Both estimators/steps, raw tokens and propagated bounds are retained.\n\n| Central projection, ea0 | Model μE at half step | Quantum zero μ |\n|---|---:|---:|\n'
  for kind in ['La','Dy','paired']:
   v=comparison[kind]['half'];report+=f"| {'Dy−La' if kind=='paired' else kind} | {float(v['model_mu']):.9f} | {float(v['quantum_zero_mu']):.9f} |\n"
  report+=f"\nAll nuclear/electronic/total dipole, laboratory-origin, state/basis, field and electronic/final/D4 closures pass. Fixed-geometry D4 is field invariant. Full scheduler cost including ten native jobs and two collectors: {result['total_scheduler_CPU_seconds']} CPU-s; exact per-job breakdown retained. Numerical consistency does not establish physical reference accuracy; failed native printed SCF criteria and frozen-f/basis/global-projection limitations remain. No unique mechanism or parameter repair is inferred.\n";(OUT/'REPORT.md').write_text(report)
  print(json.dumps({'review':pin(OUT/'REVIEW.json'),'status':result['status'],'scheduler_CPU_seconds':result['total_scheduler_CPU_seconds']},indent=2))
 except Exception as exc:
  save(OUT/'CHECKS.json',h.checks);save(OUT/'REVIEW.json',dict(base,status='incomplete_or_failed_dry_donor_susceptibility_full_review',error=repr(exc),checks=pin(OUT/'CHECKS.json'),cost_reviewed=False));raise
if __name__=='__main__':main()
