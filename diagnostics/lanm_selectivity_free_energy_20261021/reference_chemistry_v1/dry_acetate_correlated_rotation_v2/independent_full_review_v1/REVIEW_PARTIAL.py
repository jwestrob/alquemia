#!/usr/bin/env python3
"""Saved incomplete-panel account; immutable full-eight failure remains authoritative."""
import csv, importlib.util, json, math, os, re, time
from decimal import Decimal as D
from pathlib import Path

HERE=Path(__file__).resolve().parent
sp=importlib.util.spec_from_file_location('frozen_full8', HERE/'REVIEW.py')
f=importlib.util.module_from_spec(sp);sp.loader.exec_module(f)
a=f.a;b=a.b
ready=a.read(HERE/'READY.json')
assert b.pin(f.__file__)==ready['implementation']
assert all(os.environ.get(k)=='1' for k in f.THREADS)
tick,cpu=time.monotonic(),time.process_time()
out=Path(ready['output_directory']).parent/'partial_account_v1'
assert not out.exists()
collection=b.pin(ready['collection_path']);r=a.read(b.verify(collection))
failure=b.pin(Path(ready['output_directory'])/'FAILURE.json')
assert a.read(b.verify(failure))['error']=="ValueError('full8_terminal')"
assert collection==a.read(b.verify(failure))['collection']
assert r['status']=='canonical_rotation_incomplete' and r['admitted_cells']==4 and r['declared_cells']==8
p=a.load(b.verify(r['panel']));inv=a.read(b.verify(r['inventory']))
assert set(inv)==set(a.ALL)==set(r['slots'])
assert r['implementation']==b.pin(a.__file__)
assert set(r['native_admissions'])=={c for c in a.ALL if '__La__' in c}
assert all(v is None for v in r['scoring']['groups'].values())
assert all(v['passed'] is None and v['resolution_kcal_mol_A']==.1 for v in r['scoring']['checks'].values())
assert r['scoring']['declared_checks']==3 and r['scoring']['all_checks_passed'] is False
records={};energies={};costs={};raw_parse_replay_count=0
for cid,filename in inv.items():
 q=a.read(filename);c=next(c for c in p['cells'] if c['cell_id']==cid)
 assert r['slots'][cid]['result']==b.pin(filename) and q['cell']==c and q['panel']==r['panel']
 for key,name in [('start','START.json'),('input','endpoint.inp'),('output','endpoint.out'),('stderr','endpoint.err')]:
  assert q[key]==b.pin(Path(filename).parent/name)
 s=a.read(b.verify(q['start']))
 assert all(s[k]==q[k] for k in ('adapter','cell','panel','release','input','job_id'))
 assert s['ranks']==4 and s['nodes']==1 and s['memory_MiB']==16000 and s['maxcore_MiB']==3000
 assert s['command']==[p['orca']['path'],'endpoint.inp','--bind-to none']
 assert b.verify(q['input']).read_text()==a.deck(p,c,4,3000)
 release=a.read(b.verify(q['release']))
 assert release['task_output_directories'][cid]==str(Path(filename).resolve().parent)
 expected=p['_recovery_inputs']['legacy_implementation'] if cid==a.SCOUT else b.pin(a.__file__)
 assert q['adapter']==expected and q['new_native_calls']==1
 if cid!=a.SCOUT: assert q['release']==ready['remaining_release']
 native=b.verify(q['output']).read_text().split('****END OF INPUT****')[-1]
 stderr=b.verify(q['stderr']).read_text()
 table=native.split('Iter       E(tot)',1)[1]
 trace=[[D(v) for v in row] for row in re.findall(r'^\s*('+f.NUMBER+r')\s+('+f.NUMBER+r')\s+('+f.NUMBER+r')\s+('+f.NUMBER+r')\s+('+f.NUMBER+r')\s*$',table,re.M)]
 assert trace and all(v.is_finite() for row in trace for v in row)
 tolerance=re.findall(r'Convergence tolerance \(max\. residuum\)\s*\.{2,}\s*('+f.NUMBER+')',native)
 assert list(map(D,tolerance))==[D('1e-8')]
 common=dict(result=b.pin(filename),input=q['input'],output=q['output'],stderr=q['stderr'],start=q['start'],producer=q['adapter'],release=q['release'],job_id=q['job_id'],returncode=q['returncode'],resources={k:s[k] for k in ('ranks','nodes','memory_MiB','maxcore_MiB')},CC_row_count=len(trace),last_CC_row=[str(v) for v in trace[-1]],CC_tolerance=1e-8,normal_termination='ORCA TERMINATED NORMALLY' in native,SCF_criteria_lines=[l.strip() for l in native.splitlines() if l.strip().startswith('Last ') and 'Tolerance' in l],SCF_energy_stop='Energy Check signals convergence' in native,SCF_mode='Total+1el-Energy' if 'Total+1el-Energy' in native else None)
 if '__La__' in cid:
  replay=a.admit(p,filename);raw_parse_replay_count+=1
  assert replay==r['native_admissions'][cid] and all(replay['checks'].values()) and r['slots'][cid]['status']=='complete'
  e,tokens=f.decimals(native);energies[cid]=e
  assert tokens==q['energy_decimal_tokens_hartree'] and f.cc_trace(native)==q['CC_residual_trace']
  closure=e['CCSD']+e['triples']-e['CCSD_T'];bound=sum(f.half_token(v) for v in e.values())
  assert abs(closure)<=bound
  common.update(status='admitted',energy_tokens=tokens,CCSD_plus_T_minus_final_Eh=str(closure),closure_print_bound_Eh=str(bound),T1=q['T1_diagnostic'],printed_SCF_criteria=q['printed_scf_criteria'],native_checks=q['checks'])
 else:
  assert q['status']=='failed_admission' and q['error']=="ValueError('Missing/ambiguous finite CCSD_T energy')" and r['slots'][cid]['status']=='unavailable'
  assert len(trace)==100 and trace[-1][3]>D('1e-8')
  assert 'The Coupled-Cluster iterations have NOT converged' in native and 'ORCA finished by error termination in MDCI' in native
  assert not common['normal_termination'] and not re.search(r'^\s*FINAL SINGLE POINT ENERGY\s+',native,re.M)
  assert not re.search(r'^\s*Triples Correction \(T\)\s*\.',native,re.M)
  try:a.parse(native,c,q['returncode'])
  except ValueError as exc: assert repr(exc)==q['error'];raw_parse_replay_count+=1
  else:raise AssertionError('Failed Dy unexpectedly admitted')
  common.update(status='unavailable',parser_error=q['error'],cause='100 CC iterations exhausted without residual convergence; no triples/final CCSD(T) energy',last20_residual_range=[float(min(x[3] for x in trace[-20:])),float(max(x[3] for x in trace[-20:]))],MPI_segmentation_signal_reported='Segmentation fault' in stderr,canonical_energy=None)
 assert math.isfinite(q['elapsed_seconds']) and q['elapsed_seconds']>=0 and q['allocated_cpu_seconds']==4*q['elapsed_seconds']
 costs[cid]=dict(job_id=q['job_id'],status=q['status'],elapsed_seconds=q['elapsed_seconds'],CPU_seconds_wall_times_ranks=q['allocated_cpu_seconds'])
 records[cid]=common

# The producer score requires all eight and stays null. This is a separate La-only account requested by root.
qm=a.read(b.verify(p['comparison_authorities']['reference_derivative_result']))
mm=a.read(b.verify(p['comparison_authorities']['model_derivative_result']))
assert r['comparison_authorities']==p['comparison_authorities']
assert qm['admitted_cells']==16 and qm['stage_numerical_checks_passed'] and mm['complete_cells']==16 and mm['error'] is None
ha=D(str(b.HA));steps={}
for step,label in [('0.01','010'),('0.005','005')]:
 plus=energies['cc__La__rotation__plus_'+label+'__complex'];minus=energies['cc__La__rotation__minus_'+label+'__complex']
 slope={k:(plus[k]-minus[k])*ha/(2*D(step)) for k in plus}
 bounds={k:(f.half_token(plus[k])+f.half_token(minus[k]))*ha/(2*D(step)) for k in plus}
 ref=qm['scoring']['directions']['rotation']['La']['derivatives_kcal_mol_A'][step]
 model=mm['comparison']['directions']['rotation']['La']['steps'][step]['derivatives_kcal_mol_A']
 original=r['comparison']['La'][step]
 assert original['canonical_CCSD_T'] is None and original['PBE0_electronic']==ref['electronic'] and original['PBE0_D4']==ref['PBE0_D4'] and original['fixed_model_components']==model
 cc=float(slope['CCSD_T'])
 steps[step]=dict(derivatives_kcal_mol_A={k:float(v) for k,v in slope.items()},decimal_derivatives={k:str(v) for k,v in slope.items()},printing_bounds_kcal_mol_A={k:float(v) for k,v in bounds.items()},CCSD_plus_triples_derivative_minus_total=float(slope['CCSD']+slope['triples']-slope['CCSD_T']),comparison=dict(canonical_CCSD_T=cc,PBE0_electronic=ref['electronic'],PBE0_D4=ref['PBE0_D4'],D4_only=ref['D4'],selected_17621_AMOEBA=model['Total'],AMOEBA_minus_CC=model['Total']-cc,PBE0_electronic_minus_CC=ref['electronic']-cc,PBE0_D4_minus_CC=ref['PBE0_D4']-cc,cross_method_accuracy_threshold=None))
delta=D(steps['0.005']['decimal_derivatives']['CCSD_T'])-D(steps['0.01']['decimal_derivatives']['CCSD_T'])
proxy=sum(c['CPU_seconds_wall_times_ranks'] for c in costs.values())
assert abs(proxy-r['native_allocated_CPU_seconds'])<1e-7 and len(r['native_costs'])==8
assert {str(x['job_id']) for x in costs.values()}==set(map(str,ready['native_job_ids']))
sub=a.read(b.verify(ready['remaining_submission']));assert {str(j['job_id']):j['cell_id'] for j in sub['jobs']}=={str(v['job_id']):k for k,v in costs.items() if k!=a.SCOUT}
collector_sub=a.read(b.verify(ready['collector_submission']));assert collector_sub['job_id']==18868
scheduler_path=Path(ready['scheduler_path']);scheduler={'source':None,'limitation':'Expected root-supplied SACCT.psv absent; no scheduler query made. Native wall-times-ranks is a proxy, not measured utilization.'}
if scheduler_path.is_file():
 sr=list(csv.DictReader(scheduler_path.read_text().splitlines(),delimiter='|'));parents=[x for x in sr if x['JobIDRaw'].isdigit()]
 native=[x for x in parents if x['JobIDRaw'] in set(map(str,ready['native_job_ids']))]
 assert len(native)==8 and all(int(x['AllocCPUS'])==4 for x in native)
 scheduler=dict(source=b.pin(scheduler_path),parent_rows=parents,step_rows_not_double_counted=[x for x in sr if not x['JobIDRaw'].isdigit()],native_allocated_CPU_seconds=sum(int(x['ElapsedRaw'])*int(x['AllocCPUS']) for x in native),support_allocated_CPU_seconds=sum(int(x['ElapsedRaw'])*int(x['AllocCPUS']) for x in parents if x not in native))
support={}
ws=Path(ready['collection_path']).parents[2]
paths={'original_failed_scout_review':b.verify(p['_recovery_inputs']['failed_original_review']),'recovery_preparation':ws/'preparation_v1/RESULT.json','recovery_controls':ws/'saved_controls_v1/RESULT.json','corrected_scout_review':ws/'independent_actual_review_v1/actual_v1/REVIEW.json'}
for name,path in paths.items():
 q=a.read(path);assert q['new_native_calls']==0
 support[name]=dict(source=b.pin(path),status=q['status'],cost=q.get('review_cost',q.get('cost')))
assert b.pin(ready['collection_path'])==collection
result=dict(status='reviewed_incomplete_canonical_rotation_panel',implementation=b.pin(__file__),frozen_reviewer=ready['implementation'],frozen_full8_failure=failure,collection=collection,panel=r['panel'],remaining_release=ready['remaining_release'],remaining_submission=ready['remaining_submission'],collector_submission=ready['collector_submission'],declared_cells=8,admitted_cells=4,unavailable_cells=4,full_panel_qualified=False,collector_scoring_unchanged=r['scoring'],raw_parse_replay_count=raw_parse_replay_count,endpoints=records,La_only=dict(steps=steps,smaller_minus_larger_CCSD_T_derivative=float(delta),original_resolution_kcal_mol_A=.1,sampled_two_step_check_passed=abs(delta)<=D('.1'),scope='Separate descriptive La-only arithmetic; no change to incomplete collector or full-eight gate.'),Dy=None,paired=None,comparison_authorities=p['comparison_authorities'],native_costs=costs,native_wall_times_ranks_CPU_seconds=proxy,La_proxy_CPU_seconds=sum(v['CPU_seconds_wall_times_ranks'] for k,v in costs.items() if '__La__' in k),failed_Dy_proxy_CPU_seconds=sum(v['CPU_seconds_wall_times_ranks'] for k,v in costs.items() if '__Dy__' in k),scheduler=scheduler,collector_cost=r['cost'],saved_support=support,new_native_calls=0,limitations=['No Dy or paired canonical derivative exists; unconverged CCSD iterates are never substituted for CCSD(T).','La two-step agreement is sampled numerical consistency only. No cross-method accuracy gate, physical-state or CP-force qualification.','Selected 17621 AMOEBA is the model authority. D4 stays with PBE0; it is never added to canonical CC.','SCF density/DIIS/orbital-rotation failures and T1 diagnostics remain explicit. Frozen-f ECP, finite basis and single-reference limitations persist.','Centroid rotation changes metal-donor distances as well as orientation. Component derivatives do not identify a uniquely faulty physical term.','No automatic retry, iteration extension, new solver settings or source recovery is authorized by this account.'],review_cost=dict(wall_seconds=time.monotonic()-tick,CPU_seconds=time.process_time()-cpu,threads=1))
out.mkdir(parents=True,exist_ok=False);rp=b.write(out/'REVIEW.json',result)
lines=['Incomplete canonical rotation review','', '4/8 admitted (La); all four Dy results unavailable. Frozen full-eight reviewer failed unchanged and remains preserved.','', 'La-only slopes, kcal/mol/Angstrom:','q | CCSD(T) | CCSD | triples | PBE0 electronic | PBE0-D4 | selected17621 AMOEBA']
for step,z in steps.items():
 c=z['comparison'];e=z['derivatives_kcal_mol_A'];lines.append(step+' | '+' | '.join('%.9f'%v for v in [e['CCSD_T'],e['CCSD'],e['triples'],c['PBE0_electronic'],c['PBE0_D4'],c['selected_17621_AMOEBA']]))
lines += ['',f'La smaller-minus-larger step: {float(delta):.9f}; original resolution .1; pass={abs(delta)<=D(".1")}. Full panel remains unqualified.',f'All-eight native wall-times-ranks proxy: {proxy:.6f} CPU-s. Scheduler accounting: '+str(scheduler.get('native_allocated_CPU_seconds','unavailable')),'']+result['limitations']
(out/'REPORT.txt').write_text('\n'.join(lines)+'\n')
print(json.dumps(dict(review=rp,report=b.pin(out/'REPORT.txt'),La=result['La_only'],cost=result['review_cost']),indent=2))
