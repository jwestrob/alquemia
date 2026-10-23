#!/usr/bin/env python3
"""One saved replay of all18; independent Decimal CP/component/gate arithmetic."""
import argparse,importlib.util,json,math,os,re,time
from decimal import Decimal,getcontext
from pathlib import Path
HERE=Path(__file__).resolve().parent
s=importlib.util.spec_from_file_location('tighter18_admission',HERE.parent/'REFERENCE.py');a=importlib.util.module_from_spec(s);s.loader.exec_module(a)
SOURCE_SHA='81c106cf256dc2d0e233f6ed8a668d87ecf1bc4bdf267c8872230ea0c8447591'
D=Decimal;getcontext().prec=40
F=r'[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][-+]?\d+)?'

def tokens(pin):
 text=a.b.verify(pin).read_text().split('****END OF INPUT****')[-1];e={};evidence={}
 for name,pattern in {'final':r'^\s*FINAL SINGLE POINT ENERGY\s+('+F+r')\s*$', 'total':r'^\s*E\(CCSD\(T\)\)\s*\.{2,}\s*('+F+r')', 'CCSD':r'^\s*E\(CCSD\)\s*\.{2,}\s*('+F+r')', 'triples':r'^\s*Triples Correction \(T\)\s*\.{2,}\s*('+F+r')', 'RHF':r'^\s*E\(0\)\s*\.{2,}\s*('+F+r')'}.items():
  hits=re.findall(pattern,text,re.M)
  if name=='RHF' and not hits:hits=re.findall(r'^\s*Reference energy\s*\.{2,}\s*('+F+r')',text,re.M)
  a.b.check(hits and len(set(hits))==1,'Missing/ambiguous raw '+name);e[name]=D(hits[0].replace('D','E').replace('d','e'));a.b.check(e[name].is_finite(),'Nonfinite raw energy');evidence[name]=hits
 a.b.check(abs(e['final']-e['total'])<=D('1e-8') and abs(e['final']-e['CCSD']-e['triples'])<=D('1e-8'),'Raw CCSD/triples/final closure')
 return e,evidence


def main():
 ap=argparse.ArgumentParser();ap.add_argument('--result',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);ap.add_argument('--scheduler-accounting',type=Path);args=ap.parse_args()
 a.b.check(all(os.environ.get(k)=='1' for k in a.THREADS),'One local thread required');a.b.check(args.result.is_file(),'Actual collector absent; no partial review');a.b.check(not args.out.exists(),'Preserve review output')
 tick,cpu=time.monotonic(),time.process_time();rp=a.b.pin(args.result);r=a.read(rp);checks={}
 def check(name,value):checks[name]=bool(value);a.b.check(value,name)
 check('exact_frozen_source',a.b.pin(a.__file__)['sha256']==SOURCE_SHA)
 check('complete18_15plus3',r['declared']==r['complete']==18 and r['new_declared']==r['new_complete']==15 and r['reuse_declared']==r['reused_complete']==3)
 check('terminal_scientific_status',r['status'] in ('tighter_local_controls_passed','tighter_local_controls_failed'))
 p=a.load(a.b.verify(r['panel']));check('collector_source',r['implementation']==p['adapter'] and r['parser']==p['parser'])
 ids={c['cell_id'] for c in p['cells']};check('exact18_rows',set(r['rows'])==ids);inv=a.read(r['inventory']);check('exact15_paths',set(inv)==set(p['new_cell_ids']))
 reused=a.reused_rows(a.read(p['configuration']),p);raw={};rhf={};groups={};aux={};cost_new=0.;cost_reuse=0.;jobs={};max_closure=0.
 HA=D(str(a.b.HA))
 for c in p['cells']:
  cid=c['cell_id'];q=reused[cid] if cid in a.REUSED else a.admit_result(inv[cid],c,p,r['panel']['path'])
  check(cid+'/raw_native_replay',q==r['rows'][cid] and q['status']=='complete' and all(q['checks'].values()))
  local,lt=tokens(q['output']);canonical,ct=tokens(p['controls'][c['canonical_id']]['output'])
  check(cid+'/energy_tokens',abs(local['final']-D(str(q['energy_hartree'])))<=D('1e-10') and abs(canonical['final']-D(str(p['controls'][c['canonical_id']]['energy_hartree'])))<=D('1e-10'))
  check(cid+'/raw_CCSD',float(local['CCSD'])==q['CCSD_energies_hartree'][0] and float(local['triples'])==q['triples_hartree'][0])
  rhf[cid]=local['RHF']-canonical['RHF'];key=c['donor']+'__'+c['metal'];sign=D(1 if c['role']=='complex' else -1)
  g=groups.setdefault(key,{k:[D(0),D(0)] for k in ('final','RHF','CCSD','triples')})
  for k in g:g[k][0]+=sign*local[k]*HA;g[k][1]+=sign*canonical[k]*HA
  aux.setdefault(key,[]).append(q['aux_C']);closure=float(local['final']-local['CCSD']-local['triples']);max_closure=max(max_closure,abs(closure))
  start=a.read(q['start']);jobs[cid]=dict(job_id=start['job_id'],ranks=start['ranks'],nodes=start.get('nodes'),elapsed_seconds=q['elapsed_seconds'],allocated_CPU_seconds_proxy=q['allocated_cpu_seconds'],result=q['result'],output=q['output'],start=q['start'])
  check(cid+'/cost_scope',q['elapsed_seconds']>=0 and abs(q['allocated_cpu_seconds']-q['elapsed_seconds']*start['ranks'])<=1e-8)
  if cid in a.REUSED:cost_reuse+=q['allocated_cpu_seconds']
  else:cost_new+=q['allocated_cpu_seconds'];check(cid+'/one_node8',start['ranks']==8 and start['nodes']==1)
  raw[cid]=dict(local_tokens=lt,canonical_tokens=ct,component_closure_hartree=closure,CC_terminal_residuals=q['CC_terminal_residuals'],triples_terminal_deltas=q['triples_terminal_deltas'],actual_triples_settings=q['actual_triples_settings'],native_checks=q['checks'])
 independent={'RHF__'+cid:abs(v)<=D('1e-6') for cid,v in rhf.items()};cp={};components={};paired={};exchange={}
 for key,g in groups.items():
  cp[key]=g['final'];independent['CP__'+key]=abs(cp[key][0]-cp[key][1])<=D('.25');independent['AuxC__'+key]=all(x==aux[key][0] for x in aux[key])
  components[key]={k:float(v[0]-v[1]) for k,v in g.items()};check(key+'/CP_float_closure',max(abs(float(v)-r['per_metal_CP_local_canonical_kcal_mol'][key][i]) for i,v in enumerate(cp[key]))<=1e-7)
  check(key+'/CP_error_closure',abs(float(cp[key][0]-cp[key][1])-r['CP_errors_kcal_mol'][key])<=1e-7)
 for donor in ('water','acetate_bi','acetamide'):
  paired[donor]=[cp[donor+'__Dy'][i]-cp[donor+'__La'][i] for i in (0,1)];independent['paired_CP__'+donor]=abs(paired[donor][0]-paired[donor][1])<=D('.25')
  check(donor+'/paired_closure',max(abs(float(v)-r['paired_CP_local_canonical'][donor][i]) for i,v in enumerate(paired[donor]))<=1e-7)
 for donor in ('acetate_bi','acetamide'):
  exchange[donor]=[paired[donor][i]-paired['water'][i] for i in (0,1)];independent['paired_exchange__'+donor]=abs(exchange[donor][0]-exchange[donor][1])<=D('.25')
  check(donor+'/exchange_closure',max(abs(float(v)-r['paired_exchange_local_canonical'][donor][i]) for i,v in enumerate(exchange[donor]))<=1e-7)
 check('all35_unchanged_flags',len(independent)==35 and independent==r['checks']);check('scientific_verdict',all(independent.values())==r['all_controls_passed'])
 check('scientific_status_agrees',r['status']==('tighter_local_controls_passed' if all(independent.values()) else 'tighter_local_controls_failed'))
 check('no_cluster_or_physical_release',r['cluster_release_authorized'] is False and r['physical_accuracy_qualified'] is False and r['new_native_calls']==0)
 check('new_cost_sum',abs(cost_new-r['new_native_allocated_CPU_seconds'])<=1e-6);check('reuse_cost_sum',abs(cost_reuse-r['reused_native_allocated_CPU_seconds'])<=1e-6)
 failed=a.read(p['failed_tight_controls']);check('failed_Tight18_preserved',failed['complete']==18 and failed['status']=='local_controls_failed' and failed['all_controls_passed'] is False)
 accounting=None if args.scheduler_accounting is None else dict(source=a.b.pin(args.scheduler_accounting),record=json.loads(args.scheduler_accounting.read_text()))
 a.b.check(a.b.pin(args.result)==rp,'Collector changed during review');args.out.mkdir(parents=True,exist_ok=False)
 result=dict(status='passed_tighter_controls_actual_review',implementation=a.b.pin(__file__),native_adapter=p['adapter'],result=rp,panel=r['panel'],checks=checks,check_count=len(checks),admitted_cells=18,new_admitted=15,reused_admitted=3,qualification_checks=independent,all35_passed=all(independent.values()),failed_qualification_ids=[k for k,v in independent.items() if not v],
  CP_errors_kcal_mol={k:float(v[0]-v[1]) for k,v in cp.items()},paired_CP_errors_kcal_mol={k:float(v[0]-v[1]) for k,v in paired.items()},paired_exchange_errors_kcal_mol={k:float(v[0]-v[1]) for k,v in exchange.items()},
  local_minus_canonical_CP_components_kcal_mol=components,raw_energy_and_convergence=raw,maximum_raw_CCSD_triples_closure_hartree=max_closure,failed_Tight_authority=p['failed_tight_controls'],
  actual_cost=dict(new15_native_allocated_CPU_seconds_proxy=cost_new,reused3_native_allocated_CPU_seconds_proxy=cost_reuse,historical_abandoned_scheduler_CPU_seconds=r['historical_abandoned_scheduler_CPU_seconds'],collector_wall_seconds=r['wall_seconds'],collector_CPU_seconds=r['CPU_seconds'],scheduler_accounting=accounting,jobs=jobs),
  limitations=['Qualification is against this fixed18 canonical control panel at unchanged tolerances; no W8/W9 cluster, CP/B or physical-state accuracy release.', 'RHF and all-electron-with-ECP/frozen-f conventions are retained. Local CCSD/triples differences are accounting components, not uniquely identified physical causes.', 'Iterative triples use the native reported tighter/default triples settings and unchanged tolerance; no claim that every triples threshold was separately tightened.', 'Native elapsed times times allocated ranks are proxies; scheduler accounting is separate when supplied. Original abandoned attempt costs and failed Tight18 remain preserved.'],
  new_native_calls=0,cost=dict(wall_seconds=time.monotonic()-tick,CPU_seconds=time.process_time()-cpu,threads=1))
 out=a.b.write(args.out/'REVIEW.json',result)
 (args.out/'REPORT.md').write_text('# Tighter local-correlation controls: independent actual review\n\nAll18 raw roles admitted (15 new +3 reused). Unchanged35 qualification checks: '+str(sum(independent.values()))+'/35.\n\nFailed checks: '+(', '.join(result['failed_qualification_ids']) or 'none')+'.\n\n'+ '\n'.join('- '+x for x in result['limitations'])+'\n')
 print(json.dumps(dict(review=out,all35_passed=result['all35_passed'],CP_errors=result['CP_errors_kcal_mol'],paired=result['paired_CP_errors_kcal_mol'],exchange=result['paired_exchange_errors_kcal_mol'],cost=result['cost']),indent=2))
if __name__=='__main__':main()
