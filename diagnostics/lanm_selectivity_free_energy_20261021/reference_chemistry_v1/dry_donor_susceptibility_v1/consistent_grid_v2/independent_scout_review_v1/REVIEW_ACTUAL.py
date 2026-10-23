#!/usr/bin/env python3
"""Saved matching-grid Dy5 raw/Decimal review adapted from the admitted La5 reviewer. Zero native."""
import argparse,ast,hashlib,importlib.util,json,math,re,time
from decimal import Decimal as D,getcontext
from pathlib import Path
getcontext().prec=60
ROOT=Path.cwd();REL=Path('reference_chemistry_v1/dry_donor_susceptibility_v1/consistent_grid_v2');DB=ROOT/'diagnostics/lanm_selectivity_free_energy_20261021'/REL;WB=ROOT/'workspaces/lanm_selectivity_free_energy_20261021'/REL
COLLECTION=WB/'scout_collection_v1/result/RESULT.json';OUT=WB/'independent_scout_review_v1/actual_v1'
SOURCE_READY={'path': '/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/diagnostics/lanm_selectivity_free_energy_20261021/reference_chemistry_v1/dry_donor_susceptibility_v1/consistent_grid_v2/READY.json', 'sha256': '73c729953e7854af542613361836e56ba159f3037fb1f0bcb93240c40737598d'}
SOURCE_REVIEW={'path': '/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/lanm_selectivity_free_energy_20261021/reference_chemistry_v1/dry_donor_susceptibility_v1/consistent_grid_v2/independent_review_v1/REVIEW.json', 'sha256': 'b1453891f22b268c6f88afa6fe07a8aa7069245128d2e412b8f4816ee8fb7729'}
PRIOR_REVIEWER={'path': '/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/diagnostics/lanm_selectivity_free_energy_20261021/reference_chemistry_v1/dry_donor_susceptibility_v1/independent_scout_review_v1/REVIEW_ACTUAL.py', 'sha256': 'c188bee7587c761fe8ec695f5ca5caa21b9a2ba05248cd031e5dfe60d62ebda0'}
N=r'[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][-+]?\d+)?'
def pin(p):
 p=Path(p).resolve();return dict(path=str(p),sha256=hashlib.sha256(p.read_bytes()).hexdigest())
def read(p):return json.loads(Path(p).read_text())
def save(p,v):Path(p).write_text(json.dumps(v,indent=2,allow_nan=False)+'\n')
def serial(x):
 if isinstance(x,D):return str(x)
 if isinstance(x,dict):return {k:serial(v) for k,v in x.items()}
 if isinstance(x,list):return [serial(v) for v in x]
 return x
checks=[];diffs=[]
def ck(ok,label):
 checks.append(dict(check=label,passed=bool(ok)));assert ok,label
def verify(p):ck(pin(p['path'])==p,'pin:'+p['path']);return Path(p['path'])
def dd(t):return D(t.replace('D','E').replace('d','e'))
def hb(t):return D('.5')*D(10)**dd(t).as_tuple().exponent

def one(pattern,text,name):
 a=re.findall(pattern,text,re.M);ck(len(a)==1,name+'/single token');return a[0]
def section(text,name):return one(r'^\$'+re.escape(name)+r'\n([\s\S]*?)^\$End\s*$',text,'property '+name)
def array(text,name,n):
 a=one(r'&'+name+r'[^\n]*\n([\s\S]*?)(?=\n\s*&|\Z)',text,'array '+name)
 vals=re.findall(r'^\s*(\d+)\s+('+N+r')\s*$',a,re.M);ck([int(i) for i,_ in vals]==list(range(n)),'array ordered '+name);return [v for _,v in vals]
def close(x,y,bound,name):
 delta=abs(float(x)-float(y));diffs.append(delta);ck(delta<=float(bound)+8*math.ulp(max(abs(float(x)),abs(float(y)),1.)),name)

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--config',required=True);args=ap.parse_args();cfg=read(args.config)
 ck(cfg['review_implementation']==pin(__file__) and cfg['source_ready']==SOURCE_READY and cfg['source_review']==SOURCE_REVIEW,'frozen review/source configuration')
 if any(cfg[k] is None for k in ('collection','native_event','collector_event')):
  print(json.dumps({'status':'unbound_review_not_run','expected_collection':str(COLLECTION),'new_native_calls':0}));return
 ck(Path(cfg['collection']['path']).resolve()==COLLECTION.resolve(),'exact prospective collection path')
 r=read(verify(cfg['collection']))
 if r.get('status') not in ('dry_donor_susceptibility_scout_complete','dry_donor_susceptibility_scout_unqualified'):
  print(json.dumps({'status':'not_terminal_review_not_run','expected_collection':str(COLLECTION),'new_native_calls':0}));return
 ready=read(verify(SOURCE_READY));source_review=read(verify(SOURCE_REVIEW));verify(PRIOR_REVIEWER)
 OUT.mkdir(parents=True,exist_ok=False);tick=time.monotonic();cpu=time.process_time();base=dict(implementation=ready['implementation'],panel=ready['panel'],collection=pin(COLLECTION),review_implementation=pin(__file__),source_review=SOURCE_REVIEW,prior_reviewer=PRIOR_REVIEWER,authorship='Singer authored the matching-grid adapter and adapted Arendt original Decimal reviewer; source was separately reviewed by Arendt. This receipt is saved raw/independent-arithmetic evidence, not a claim of independent source authorship.',new_native_calls=0)
 try:
  ck(r['implementation']==ready['implementation'] and r['panel']==ready['panel'] and r['stage']=='scout','exact final scout source/panel')
  src=verify(r['implementation']);sp=importlib.util.spec_from_file_location('dry_scout_frozen',src);mod=importlib.util.module_from_spec(sp);sp.loader.exec_module(mod)
  p=mod.load(verify(r['panel']));ck(r['cell_ids']==p['scout_cell_ids'] and r['declared']==5 and r['complete']==5,'complete exact matching-grid Dy5')
  ck(source_review['status']=='passed_dry_donor_susceptibility_source_review' and source_review['implementation']==r['implementation'] and source_review['panel']==r['panel'],'actual independent source review bound')
  raw={};native={};snapshots=[]
  for cid in p['scout_cell_ids']:
   rec=r['rows'][cid];c=next(c for c in p['cells'] if c['cell_id']==cid);q=mod.admit_result(verify(rec['result']),c,p,r['panel']['path']);ck(q==rec,cid+'/single raw native admission')
   text=verify(rec['output']).read_text()
   ck(c['metal']=='Dy' and c['expected_grid_fingerprints']==dict(c['original_grid_fingerprints'],grid2=c['original_grid_fingerprints']['grid3']),'only expected Dy grid2 changes')
   ck(mod.grid.grid_fingerprints(text)==c['expected_grid_fingerprints'] and 'Recomputing exchange energy using gridx3' in text,cid+'/actual matching full physical grids')
   proptext=verify(rec['properties']['property_file']).read_text();g=section(proptext,'Geometry');ck('&Units "Bohr"' in g,cid+'/Bohr property units');xyz=re.findall(r'^\s*([A-Za-z]+)\s+('+N+r')\s+('+N+r')\s+('+N+r')\s*$',g,re.M);ck(len(xyz)==8 and [v[0] for v in xyz]==c['symbols'],cid+'/real property sites')
   block=section(proptext,'SCF_Dipole_Moment');vs={k:array(block,name,3) for k,name in [('mu','dipoleTotal'),('nuc','dipoleNucContrib'),('ele','dipoleElecContrib')]};z=[11,6,6,8,8,1,1,1]
   for k in range(3):
    predicted=sum((D(z[i])*dd(xyz[i][k+1]) for i in range(8)),D(0));bound=sum((D(abs(z[i]))*hb(xyz[i][k+1]) for i in range(8)),D(0))+hb(vs['nuc'][k]);close(predicted,dd(vs['nuc'][k]),bound,cid+'/effective nuclear origin'+str(k))
    close(dd(vs['mu'][k]),dd(vs['nuc'][k])+dd(vs['ele'][k]),sum((hb(vs[t][k]) for t in vs),D(0)),cid+'/dipole components'+str(k))
   origin=one(r'^Position of electric origin\s*\.\.\.\s*('+N+r')\s+('+N+r')\s+('+N+r')\s*$',text,cid+'/native origin');ck(list(map(dd,origin))==[0,0,0],cid+'/laboratory zero')
   et=array(section(proptext,'SCF_Energy'),'totalEnergy',1)[0];ft=one(r'&FinalEnergy\s+\[&Type "Double"\]\s+('+N+r')',section(proptext,'Single_Point_Data'),cid+'/final property');dt=one(r'^Dispersion correction\s+('+N+r')\s*$',text,cid+'/D4')
   nt=one(r'^FINAL SINGLE POINT ENERGY\s+('+N+r')\s*$',text,cid+'/native final');ns=one(r'^Total Energy\s*:\s*('+N+r')\s+Eh',text,cid+'/native electronic')
   close(dd(et),dd(ns),hb(et)+hb(ns),cid+'/electronic property/output');close(dd(ft),dd(nt),hb(ft)+hb(nt),cid+'/final property/output');close(dd(ft),dd(et)+dd(dt),hb(ft)+hb(et)+hb(dt),cid+'/D4 energy stage closure')
   mu=sum((D(str(t))*dd(v) for t,v in zip(p['field_direction_lab'],vs['mu'])),D(0));mub=sum((abs(D(str(t)))*hb(v) for t,v in zip(p['field_direction_lab'],vs['mu'])),D(0))
   raw[c['field_tag']]=dict(e=dd(et),final=dd(ft),d4=dd(dt),mu=mu,mu3=[dd(t) for t in vs['mu']],eb=hb(et),fb=hb(ft),db=hb(dt),mb=mub,energy_tokens=dict(electronic=et,final=ft,D4=dt,native_final=nt,native_electronic=ns),dipole_tokens=vs,field_au=c['field_au'],field_vector_au=c['field_vector_au'])
   native[cid]={k:rec['parsed'].get(k) for k in ['native_scf_mode','native_scf_stop_signals','printed_scf_criteria','native_diagnostic_lines','dft_xc_grid_metadata','cosx_metadata','physical_grid_fingerprints','grid_parse_error','energy_hartree','electronic_energy_hartree','D4_energy_hartree']};snapshots.append(dict(cell_id=cid,result=rec['result'],output=rec['output'],property=rec['properties']['property_file'],start=rec['start'],elapsed_seconds=rec['elapsed_seconds'],allocated_elapsed_CPU_seconds_proxy=rec['allocated_elapsed_CPU_seconds_proxy']))
  ck(set(raw)=={'zero','minus_h','minus_half','plus_half','plus_h'},'five signed raw slots')
  zero=raw['zero'];values={};gates={};d4_ok=all(abs(v['d4']-zero['d4'])<=v['db']+zero['db'] for v in raw.values());ck(d4_ok==r['response']['D4_field_independence']['Dy'],'D4 field equality replay')
  for label,h in [('h',D('.001')),('half',D('.0005'))]:
   a,b=raw['minus_'+label],raw['plus_'+label];est={'chi_energy':-(b['e']+a['e']-2*zero['e'])/h**2,'chi_moment':(b['mu']-a['mu'])/(2*h),'odd_energy_moment':-(b['e']-a['e'])/(2*h),'zero_moment':zero['mu'],'even_moment_numerator':b['mu']+a['mu']-2*zero['mu']};est['total_chi_energy']=-(b['final']+a['final']-2*zero['final'])/h**2
   actual=r['response']['per_metal']['Dy']['estimates'][label]
   for k in ['chi_energy','chi_moment','odd_energy_moment','zero_moment','even_moment_numerator']:close(est[k],actual[k],actual['roundoff_bounds'][k],'Decimal derived '+label+'/'+k)
   close(est['total_chi_energy'],r['response']['per_metal']['Dy']['final_total_estimates'][label]['chi_energy'],r['response']['per_metal']['Dy']['final_total_estimates'][label]['roundoff_bounds']['chi_energy'],'Decimal total curvature '+label)
   est['decimal_printing_bounds']={'chi_energy':(a['eb']+b['eb']+2*zero['eb'])/h**2,'chi_moment':(a['mb']+b['mb'])/(2*h),'odd_energy_moment':(a['eb']+b['eb'])/(2*h),'zero_moment':zero['mb']};values[label]=est
   gates['Dy__energy_vs_moment_'+label]=abs(est['chi_energy']-est['chi_moment'])<=D('.5');gates['Dy__odd_slope_'+label]=abs(est['odd_energy_moment']-est['zero_moment'])<=D('.0001')
  for key in ['chi_energy','chi_moment']:gates['Dy__'+key+'_step']=abs(values['h'][key]-values['half'][key])<=D('.5')
  ck(gates==r['response']['checks'] and len(gates)==6,'independent original six gates');passed=all(gates.values()) and d4_ok
  ck(passed==r['numerical_checks_passed'] and passed==r['response']['passed'],'producer numerical qualification replay')
  ck(mod.score(r['rows'],p,'scout')==r['response'],'single saved response replay')
  # Snapshot all owned scheduler records, including collector, rather than infer allocation from admitted rows.
  events={};jobs={}
  for name in ('native','collector'):
   path=verify(cfg[name+'_event']);ev=read(path);events[name]={'source':cfg[name+'_event'],'records':ev.get('records',{})}
   for jid,v in ev.get('records',{}).items():
    ck(jid not in jobs,'unique scheduler job record');jobs[jid]=v
  native_ids={str(r['rows'][cid]['job_id']) for cid in p['scout_cell_ids']};collector_ids=set(events['collector']['records'])
  ck(len(native_ids)==5 and set(events['native']['records'])==native_ids and len(collector_ids)==1 and len(jobs)==6,'exact five producer jobs and one collector')
  ck(all(jobs[j]['state']=='COMPLETED' and jobs[j]['exit_code']=='0:0' and int(jobs[j]['allocated_cpus'])==4 for j in native_ids),'five native jobs completed with four CPUs')
  cj=next(iter(collector_ids));ck((jobs[cj]['state'],jobs[cj]['exit_code'])==(('COMPLETED','0:0') if passed else ('FAILED','1:0')) and int(jobs[cj]['allocated_cpus'])==1,'collector exit retains numerical failure if present')
  costs={jid:dict(elapsed_seconds=int(v['elapsed_seconds']),allocated_CPUs=int(v['allocated_cpus']),allocated_CPU_seconds=int(v['elapsed_seconds'])*int(v['allocated_cpus'])) for jid,v in jobs.items()};save(OUT/'SCHEDULER_SNAPSHOT.json',events)
  save(OUT/'RAW_REPLAY.json',serial({'raw':raw,'Decimal_estimates':values,'numerical_checks':gates,'D4_field_independence':d4_ok,'native_SCf_metadata':native,'native_snapshots':snapshots}));save(OUT/'CHECKS.json',checks)
  result=dict(base,status='passed_dry_donor_susceptibility_scout_review' if passed else 'failed_dry_donor_susceptibility_scout_review',review_passed=True,failed_numerical_checks=[k for k,v in gates.items() if not v],check_count=len(checks),checks=pin(OUT/'CHECKS.json'),raw_replay=pin(OUT/'RAW_REPLAY.json'),complete_cells=5,declared_cells=5,saved_native_admissions=5,numerical_checks=6,numerical_checks_passed=passed,D4_field_independence=d4_ok,cost_reviewed=True,scheduler_snapshot=pin(OUT/'SCHEDULER_SNAPSHOT.json'),scheduler_costs=costs,total_scheduler_CPU_seconds=sum(v['allocated_CPU_seconds'] for v in costs.values()),native_elapsed_CPU_seconds_proxy=sum(x['allocated_elapsed_CPU_seconds_proxy'] for x in snapshots),maximum_independent_derived_difference=max(diffs),review_cost={'wall_seconds':time.monotonic()-tick,'CPU_seconds':time.process_time()-cpu,'threads':1},limitations=['Printing/arithmetic bounds and numerical gates do not qualify physical accuracy.','Full native SCF stopping/printed criteria and final integration stages retained.','Uniform tangent susceptibility does not isolate local charge-transfer or induction mechanisms.','Old EnGrad center remains contextual; matched SP zero used in every estimator.'])
  save(OUT/'REVIEW.json',result)
  title='passes' if passed else 'fails numerical qualification';report=f'# Matching-grid Dy5 susceptibility scout raw/Decimal review\n\nThe complete five-cell scout {title}; all raw outputs/properties are admitted and all six original gates independently replayed. Zero native calls by this review.\n\n'
  report+='| Field size, au | chi_E, a0³ | chi_mu, a0³ | odd energy moment, ea0 | zero moment, ea0 |\n|---|---:|---:|---:|---:|\n'
  for tag,h in [('h','.001'),('half','.0005')]:
   v=values[tag];report+=f"| {h} | {float(v['chi_energy']):.9f} | {float(v['chi_moment']):.9f} | {float(v['odd_energy_moment']):.9f} | {float(v['zero_moment']):.9f} |\n"
  report+=f"\nD4 field-independence: {d4_ok}. Electronic/final/D4 and nuclear/electronic/total dipole closures, explicit lab origin and all state/basis/field joins are retained. Scheduler cost across five native jobs plus collector: {result['total_scheduler_CPU_seconds']} CPU-s; per-job breakdown in REVIEW.json. Native stopping criteria and printing bounds are not physical accuracy guarantees. No automatic La release is performed. The original ten-cell unqualified panel is preserved; the paired result requires the same new recipe for both metals.\n"
  (OUT/'REPORT.md').write_text(report);print(json.dumps({'review':pin(OUT/'REVIEW.json'),'status':result['status'],'cost':result['total_scheduler_CPU_seconds']},indent=2))
 except Exception as exc:
  save(OUT/'CHECKS.json',checks);save(OUT/'REVIEW.json',dict(base,status='incomplete_or_failed_dry_donor_susceptibility_scout_review',error=repr(exc),checks=pin(OUT/'CHECKS.json'),cost_reviewed=False));raise
if __name__=='__main__':main()
