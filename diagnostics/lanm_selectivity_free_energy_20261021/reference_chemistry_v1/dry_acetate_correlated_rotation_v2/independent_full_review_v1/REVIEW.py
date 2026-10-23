#!/usr/bin/env python3
"""Terminal canonical8 saved readback; original three gates and explicit cost scopes."""
import argparse,csv,hashlib,importlib.util,json,math,os,re,time
from decimal import Decimal,getcontext
from pathlib import Path
HERE=Path(__file__).resolve().parent
sp=importlib.util.spec_from_file_location('canonical_rotation_v2',HERE.parent/'REFERENCE.py');a=importlib.util.module_from_spec(sp);sp.loader.exec_module(a)
SOURCE_SHA='007ad5b22ffb09fd3ed700c64d5ff1adf76b579b71fadaa5cf9bdda4b1620d4a'
D=Decimal;getcontext().prec=40;NUMBER=r'[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][-+]?\d+)?'
THREADS=('OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','BLIS_NUM_THREADS')

def decimals(text):
 native=text.split('****END OF INPUT****')[-1];values={};tokens={}
 for key,pattern in [('CCSD_T',r'^\s*FINAL SINGLE POINT ENERGY\s+'),('CCSD',r'^\s*E\(CCSD\)\s*\.{2,}\s*'),('triples',r'^\s*Triples Correction \(T\)\s*\.{2,}\s*')]:
  hits=re.findall(pattern+'('+NUMBER+')',native,re.M);a.b.check(len(hits)==1,'Missing/duplicate canonical '+key);tokens[key]=hits[0].replace('D','E').replace('d','e');values[key]=D(hits[0].replace('D','E').replace('d','e'));a.b.check(values[key].is_finite(),'Nonfinite canonical energy')
 a.b.check(abs(values['CCSD_T']-values['CCSD']-values['triples'])<D('1e-8'),'Unchanged CCSD/triples closure')
 return values,tokens

def cc_trace(native):
 header=re.findall(r'^Iter\s+E\(tot\)\s+E\(Corr\)\s+Delta-E\s+Residual\s+Time\s*$',native,re.M)
 a.b.check(len(header)==1,'Unique canonical residual table required')
 table=native.split(header[0],1)[1].split('COUPLED CLUSTER ENERGY',1)[0]
 records=[[float(D(v.replace('D','E').replace('d','e'))) for v in x] for x in re.findall(r'^\s*('+NUMBER+r')\s+('+NUMBER+r')\s+('+NUMBER+r')\s+('+NUMBER+r')\s+('+NUMBER+r')\s*$',table,re.M)]
 tol=[D(x.replace('D','E').replace('d','e')) for x in re.findall(r'Convergence tolerance \(max\. residuum\)\s*\.{2,}\s*('+NUMBER+r')',native)]
 a.b.check(tol==[D('1e-8')] and records and all(math.isfinite(v) for x in records for v in x) and abs(records[-1][3])<=1e-8,'Original canonical residual gate')
 a.b.check(bool(re.search(r'---\s*The\s+(?:Coupled[ -]Cluster|CCSD)\s+iterations\s+have\s+converged\s*---',native,re.I)),'Actual canonical convergence marker')
 return records

def half_token(value):return D('.5')*D(10)**value.as_tuple().exponent

def accounting(path,ids):
 # Root supplies an exact sacct -P table with a header; scheduler step rows are preserved but not double-counted.
 rows=list(csv.DictReader(Path(path).read_text().splitlines(),delimiter='|'));required={'JobIDRaw','State','ElapsedRaw','AllocCPUS','ExitCode'}
 a.b.check(rows and required<=set(rows[0]),'Scheduler header requires JobIDRaw|State|ElapsedRaw|AllocCPUS|ExitCode')
 native={};support={};steps=[]
 for row in rows:
  jid=row['JobIDRaw'].strip()
  if not jid.isdigit():steps.append(row);continue
  a.b.check(jid not in native and jid not in support,'Duplicate scheduler parent job')
  elapsed=int(row['ElapsedRaw']);cpus=int(row['AllocCPUS']);a.b.check(elapsed>=0 and cpus>=0,'Invalid scheduler allocation')
  item=dict(raw=row,elapsed_seconds=elapsed,allocated_CPUs=cpus,allocated_CPU_seconds=elapsed*cpus)
  if jid in ids:
   a.b.check(row['State']=='COMPLETED' and row['ExitCode']=='0:0' and cpus==4,'Canonical native scheduler not complete4CPU');native[jid]=item
  else:support[jid]=item
 a.b.check(set(native)==ids and len(native)==8,'Exact eight native scheduler jobs required')
 return dict(source=a.b.pin(path),native_jobs=native,support_jobs_from_root_table=support,step_rows_not_double_counted=steps,
  native_allocated_CPU_seconds=sum(v['allocated_CPU_seconds'] for v in native.values()),support_allocated_CPU_seconds=sum(v['allocated_CPU_seconds'] for v in support.values()),meaning='ElapsedRaw times AllocCPUS; CPU utilization is not inferred. Root-supplied non-native parent rows are separate support accounting.')


def main(args):
 a.b.check(all(os.environ.get(k)=='1' for k in THREADS),'One saved-analysis thread required');a.b.check(args.collection.is_file(),'Terminal collection absent');a.b.check(not args.out.exists(),'Preserve independent review output')
 tick,cpu=time.monotonic(),time.process_time();rp=a.b.pin(args.collection);r=a.read(args.collection);checks={}
 def check(name,v):checks[name]=bool(v);a.b.check(v,name)
 check('exact_source',a.b.pin(a.__file__)['sha256']==SOURCE_SHA)
 check('full8_terminal',r['status']=='canonical_rotation_complete' and r['admitted_cells']==r['declared_cells']==8)
 p=a.load(a.b.verify(r['panel']));check('collector_source',r['implementation']==a.b.pin(a.__file__))
 inv=a.read(a.b.verify(r['inventory']));check('all8_slots',set(inv)==set(r['slots'])==set(r['native_admissions'])==set(a.ALL))
 rows={};energies={};raw={};jobids=set();nativeproxy=0.
 for c in p['cells']:
  cid=c['cell_id'];q=a.admit(p,inv[cid]);check(cid+'/exact_raw_admission',q==r['native_admissions'][cid] and all(q['checks'].values()))
  check(cid+'/slot_pin',r['slots'][cid]==dict(status='complete',result=a.b.pin(inv[cid])))
  expected=p['_recovery_inputs']['legacy_implementation'] if cid==a.SCOUT else a.b.pin(a.__file__);check(cid+'/producer_lineage',q['adapter']==expected)
  text=a.b.verify(q['output']).read_text();native=text.split('****END OF INPUT****')[-1];e,tokens=decimals(text);energies[cid]=e
  check(cid+'/printed_energy_tokens',tokens==q['energy_decimal_tokens_hartree'])
  trace=cc_trace(native);check(cid+'/CC_residual_replay',trace==q['CC_residual_trace'] and q['CC_printed_tolerance']==[1e-8])
  check(cid+'/RHF_NoFrozenCore',bool(re.search(r'Reference Wavefunction\s*\.{2,}\s*RHF\b',native)) and bool(re.search(r'Frozen core treatment\s*\.{2,}\s*NO frozen core\b',native)))
  for name,pattern,target in [('electrons',r'Number of [Ee]lectrons\s+(?:NEL\s+)?\.{2,}\s*(\d+)',40),('correlated',r'Number of correlated electrons\s*\.{2,}\s*(\d+)',40),('AO',r'Number of basis functions\s*\.{2,}\s*(\d+)',189)]:
   values=[int(x) for x in re.findall(pattern,native,re.I)];check(cid+'/'+name,bool(values) and all(x==target for x in values))
  check(cid+'/ECP',re.findall(r'Atom\s+(\d+)\s*[A-Za-z]+\s+ECP group\s*=>',native)==['0'] and re.findall(r'replacing\s+(\d+)\s+core electrons',native)==[str(46 if c['metal']=='La' else 55)])
  jobids.add(str(q['job_id']));nativeproxy+=q['allocated_cpu_seconds'];rows[cid]=q
  raw[cid]=dict(result=a.b.pin(inv[cid]),input=q['input'],output=q['output'],start=q['start'],producer=q['adapter'],release=q['release'],energy_tokens=tokens,
   SCF_stop_signals=q['native_scf_stop_signals'],SCF_check_mode=q['native_scf_mode'],printed_SCF_criteria=q['printed_scf_criteria'],CC_residual_trace=trace,CC_printed_tolerance=q['CC_printed_tolerance'],final_CC_residual=trace[-1][3],T1_diagnostic=q['T1_diagnostic'],largest_printed_amplitudes=q['largest_printed_amplitudes'],native_checks=q['checks'],elapsed_seconds=q['elapsed_seconds'],allocated_CPU_seconds_proxy=q['allocated_cpu_seconds'],job_id=q['job_id'])
 check('exact8_native_jobs',jobids=={'18774','18861','18862','18863','18864','18865','18866','18867'})
 HA=D(str(a.b.HA));der={};bounds={};gate={};arithmetic={}
 for metal in ('La','Dy'):
  der[metal]={};bounds[metal]={}
  for step,label in [('0.01','010'),('0.005','005')]:
   plus=energies['cc__'+metal+'__rotation__plus_'+label+'__complex'];minus=energies['cc__'+metal+'__rotation__minus_'+label+'__complex']
   der[metal][step]={k:(plus[k]-minus[k])*HA/(2*D(step)) for k in plus}
   bounds[metal][step]={k:(half_token(plus[k])+half_token(minus[k]))*HA/(2*D(step)) for k in plus}
 der['paired']={step:{k:der['Dy'][step][k]-der['La'][step][k] for k in der['La'][step]} for step in ('0.01','0.005')}
 bounds['paired']={step:{k:bounds['Dy'][step][k]+bounds['La'][step][k] for k in bounds['La'][step]} for step in ('0.01','0.005')}
 qm=a.read(a.b.verify(p['comparison_authorities']['reference_derivative_result']));mm=a.read(a.b.verify(p['comparison_authorities']['model_derivative_result']))
 check('original_comparison_authorities',r['comparison_authorities']==p['comparison_authorities'] and qm['admitted_cells']==16 and qm['stage_numerical_checks_passed'] is True and mm['complete_cells']==16 and mm['error'] is None)
 for metal in ('La','Dy','paired'):
  delta=der[metal]['0.005']['CCSD_T']-der[metal]['0.01']['CCSD_T'];gate[metal]=abs(delta)<=D('.1')
  check(metal+'/original_step_gate',gate[metal]==r['scoring']['checks'][metal]['passed'] and r['scoring']['checks'][metal]['resolution_kcal_mol_A']==.1 and abs(float(delta)-r['scoring']['checks'][metal]['difference_kcal_mol_A'])<=1e-8 and abs(float(delta)-r['scoring']['groups'][metal]['smaller_minus_larger_step_kcal_mol_A'])<=1e-8)
  arithmetic[metal]={'smaller_minus_larger_step':float(delta),'steps':{}}
  for step in ('0.01','0.005'):
   for component,v in der[metal][step].items():
    stored=r['scoring']['groups'][metal]['steps'][step][component];check(metal+'/'+step+'/'+component,abs(float(v)-stored['derivative_kcal_mol_A'])<=1e-8 and abs(float(bounds[metal][step][component])-stored['print_bound_kcal_mol_A'])<=1e-15)
   cc=float(der[metal][step]['CCSD_T']);ref=qm['scoring']['directions']['rotation'][metal]['derivatives_kcal_mol_A'][step];model=mm['comparison']['directions']['rotation'][metal]['steps'][step]['derivatives_kcal_mol_A']
   expected=dict(canonical_CCSD_T=cc,PBE0_electronic=ref['electronic'],PBE0_D4=ref['PBE0_D4'],D4_only=ref['D4'],fixed_model=model['Total'],model_minus_CC=model['Total']-cc,PBE0_electronic_minus_CC=ref['electronic']-cc,PBE0_D4_minus_CC=ref['PBE0_D4']-cc)
   check(metal+'/'+step+'/comparison',all(abs(v-r['comparison'][metal][step][k])<=1e-8 for k,v in expected.items()) and r['comparison'][metal][step]['fixed_model_components']==model and r['comparison'][metal][step]['cross_method_accuracy_gate'] is None)
   arithmetic[metal]['steps'][step]=dict(components_kcal_mol_A={k:float(v) for k,v in der[metal][step].items()},print_bounds_kcal_mol_A={k:float(v) for k,v in bounds[metal][step].items()},comparison=dict(expected,model_label='selected 17621 AMOEBA'))
 check('unchanged3_checks',len(r['scoring']['checks'])==r['scoring']['declared_checks']==3 and all(gate.values())==r['scoring']['all_checks_passed'])
 check('native_proxy_sum',abs(nativeproxy-r['native_allocated_CPU_seconds'])<=1e-6);check('no_accuracy_or_CP_force_promotion',not r['physical_state_qualified'] and not r['method_accuracy_qualified'] and not r['CP_force_qualified'] and r['new_native_calls']==0)
 costs=accounting(args.scheduler,jobids)
 # Separate saved recovery costs; these are already-consumed support CPU, never added to native endpoint totals.
 # Explicit sibling workspace from the bound original panel avoids depending on native task depth.
 ws=Path(r['panel']['path']).parents[2]/'dry_acetate_correlated_rotation_v2'
 support={}
 original=Path(r['panel']['path']).parents[1]
 paths={'original_preparation':original/'preparation_v1/PREFLIGHT.json','original_saved_controls':original/'source_controls_v1/RESULT.json','original_root_source_review':original/'root_source_review_v1/REVIEW.json','recovery_source_review':ws/'independent_source_review_v1/REVIEW.json','failed_original_scout_review':Path(p['_recovery_inputs']['failed_original_review']['path']),
  'recovery_preparation':ws/'preparation_v1/RESULT.json','recovery_controls':ws/'saved_controls_v1/RESULT.json','corrected_scout_review':ws/'independent_actual_review_v1/actual_v1/REVIEW.json'}
 for name,path in paths.items():
  q=a.read(path);support[name]=dict(source=a.b.pin(path),status=q['status'],cost=q.get('review_cost',q.get('cost')),new_native_calls=q['new_native_calls'])
  check(name+'/zero_native',q['new_native_calls']==0)
 check('collector_unchanged',a.b.pin(args.collection)==rp);args.out.mkdir(parents=True,exist_ok=False)
 out=dict(status='passed_canonical_rotation_full8_actual_review',implementation=a.b.pin(__file__),native_adapter=a.b.pin(a.__file__),collection=rp,panel=r['panel'],legacy_inputs=a.b.pin(HERE.parent/'INPUTS.json'),checks=checks,check_count=len(checks),admitted_cells=8,old_scout_reused=1,new_native_endpoints=7,original_three_checks=gate,all_three_step_checks_passed=all(gate.values()),independent_arithmetic=arithmetic,raw_native=raw,scheduler_cost=costs,native_allocated_CPU_seconds_proxy=nativeproxy,collector_saved_cost=r['cost'],saved_recovery_support_costs=support,new_native_calls=0,
  limitations=['Original three 0.1 kcal/mol/Angstrom checks qualify only sampled two-step numerical consistency. Cross-method columns have no added accuracy gate.','Model column is exact selected 17621 AMOEBA, not MACE. CCSD/triples components are accounting, not causal attribution. D4 stays with PBE0 and is never added to canonical CC.','Raw SCF stopping criteria, T1 diagnostics and largest amplitudes are retained. Canonical CCSD(T), finite basis, frozen-f ECP and single-reference/physical-state limits remain. No CP-corrected force or physical truth claim.','Scheduler allocation cost, native elapsed-times-ranks proxies and local recovery/reporting CPU are separate scopes; step rows are not double-counted.'],review_cost=dict(wall_seconds=time.monotonic()-tick,CPU_seconds=time.process_time()-cpu,threads=1))
 receipt=a.b.write(args.out/'REVIEW.json',out)
 lines=['# Canonical rotation full-eight independent review','',str(sum(gate.values()))+'/3 original step checks pass; all 8 raw endpoints admitted.','', '| Metal | q | Canonical CCSD(T) | PBE0 electronic | PBE0-D4 | selected 17621 AMOEBA |','|---|---:|---:|---:|---:|---:|']
 for metal,v in arithmetic.items():
  for step,z in v['steps'].items():
   c=z['comparison'];lines.append('| '+metal+' | '+step+' | '+' | '.join('%.9f'%c[k] for k in ('canonical_CCSD_T','PBE0_electronic','PBE0_D4','fixed_model'))+' |')
 lines+=['','Units kcal/mol/Angstrom.','']+['- '+x for x in out['limitations']]
 (args.out/'REPORT.md').write_text('\n'.join(lines)+'\n');print(json.dumps(dict(review=receipt,all_three_checks=all(gate.values()),native_scheduler_CPU_seconds=costs['native_allocated_CPU_seconds'],support_scheduler_CPU_seconds=costs['support_allocated_CPU_seconds']),indent=2))
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--collection',type=Path,required=True);ap.add_argument('--scheduler',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);args=ap.parse_args()
 preexisted=args.out.exists();tick,cpu=time.monotonic(),time.process_time()
 try:main(args)
 except Exception as exc:
  if not preexisted:
   args.out.mkdir(parents=True,exist_ok=True)
   failure=dict(status='failed_canonical_rotation_full8_actual_review',implementation=a.b.pin(__file__),collection=a.b.pin(args.collection) if args.collection.is_file() else None,scheduler=a.b.pin(args.scheduler) if args.scheduler.is_file() else None,error=repr(exc),new_native_calls=0,review_cost=dict(wall_seconds=time.monotonic()-tick,CPU_seconds=time.process_time()-cpu,threads=1))
   a.b.write(args.out/'FAILURE.json',failure)
  raise
