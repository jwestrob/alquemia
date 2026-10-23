"""Finite source/metadata joins; no imports of campaign or numerical modules."""
import ast
import copy
import hashlib
import json
from pathlib import Path
import time

B=Path(__file__).resolve().parents[8]
D=Path(__file__).resolve().parent.parent
W=Path(str(D).replace('/diagnostics/','/workspaces/'))/'independent_source_review_v1'
checks=[];verified=set();tick=time.perf_counter();cpu=time.process_time()

def pin(p):
    p=Path(p).resolve();return dict(path=str(p),sha256=hashlib.sha256(p.read_bytes()).hexdigest())

def verify(x):
    assert pin(x['path'])==x,x['path'];verified.add(x['path']);return Path(x['path'])

def read(x):return json.loads(verify(x).read_text())

def check(name,ok):
    assert ok,name;checks.append(name)

def digest(x):return hashlib.sha256(json.dumps(x,sort_keys=True,separators=(',',':'),allow_nan=False).encode()).hexdigest()

ready=json.loads((D/'READY.json').read_text())
for x in ready.values():
    if isinstance(x,dict) and set(x)=={'path','sha256'}:verify(x)
cfg=read(ready['configuration_template']);prep=read(ready['preparation']);inv=read(ready['reuse_inventory']);controls=read(ready['controls'])
check('unchanged fingerprint',digest({k:v for k,v in cfg.items() if k not in ('root_release','source_review','preparation')})==ready['configuration_fingerprint']==prep['configuration_fingerprint'])
for key in ('implementation','adapter','core','reuse_implementation','journal_implementation','declaration','execution_context','reuse_inventory'):
    check('configuration binding '+key,cfg[key]==ready[key])
check('zero-call preparation',prep['status']=='prepared_physical_gate_component_fit_zero_native' and prep['native_programs']==prep['real_proposal_calls']==0 and prep['error'] is None)
check('867 exact metadata rows',len(inv['registry'])==867 and all(len(r['parameters'])==12 and digest(r['parameters'])==k for k,r in inv['registry'].items()))
keys=[inv['origin_key'],*inv['derivative_keys'],inv['pilot_key'],inv['full_pilot_key']]
check('27 distinct seeds and exact24 derivative identities',len(set(keys))==27 and len(inv['derivative_keys'])==24 and set(prep['vectors'])==set(keys))
check('690 imported cells recorded by preparation',prep['import_counts']==ready['counts'] and prep['import_counts']['unique_cells']==690)
check('origin remains incumbent',prep['selected_key']==prep['origin_key']==inv['origin_key'] and prep['candidate_frozen'] is False)
check('native full-pilot failure distinct',inv['full_pilot_key']==ready['full_pilot']['key'] and prep['vectors'][inv['full_pilot_key']]['import_role']=='failed_full_pilot')
full=read(inv['authorities']['full_pilot_root']);check('full pilot24/25 preserved',full['status']=='reviewed_complete_native_tradeoff_physical_gate_failure_preserved' and full['physical_gates_passed']==24)
check('all24 focused saved controls pass',controls['status']=='passed_physical_gate_reuse_controls' and controls['checks_passed']==len(controls['checks'])==24 and controls['new_native_calls']==controls['solver_calls']==0)
for filename,key in [('DRIVER.py','implementation'),('ADAPTER.py','adapter'),('REUSE.py','reuse_implementation'),('JOURNAL.py','journal_implementation')]:
    check('controls bound to '+filename,controls['sources'][filename]==ready[key])
core_review=read(ready['core_review']);check('independent CORE review retained',core_review['core']==ready['core'] if 'core' in core_review else ready['core'] in core_review.values())
parent=D.parent/'persistence_performance_v1/future_adapter_v1'
check('journal unchanged bytewise',(D/'JOURNAL.py').read_bytes()==(parent/'JOURNAL.py').read_bytes())
def method(path,name):return next(n for n in ast.walk(ast.parse(path.read_text())) if isinstance(n,ast.FunctionDef) and n.name==name)
for name in ('_invoke','authorize','_stop'):
    check('unchanged '+name,ast.dump(method(D/'ADAPTER.py',name))==ast.dump(method(parent/'ADAPTER.py',name)))
new=copy.deepcopy(method(D/'ADAPTER.py','_cell'))
lock=next(x for x in new.body if isinstance(x,ast.With) and ast.unparse(x.items[0].context_expr)=='lock')
check('only cell delta is lazy exact import',ast.unparse(lock.body[1]).startswith('saved_key = self._saved_metal') and ast.unparse(lock.body[2]).startswith('if token not in self.cache and saved_key is not None:'))
del lock.body[1:3]
check('remaining cell path unchanged',ast.dump(new)==ast.dump(method(parent/'ADAPTER.py','_cell')))
declaration=read(ready['declaration']);check('limits exact',cfg['limits']==ready['limits']==declaration['limits'])
check('physical gates retained and old ceilings descriptive',declaration['actual_acceptance']['all25_original_physical_gates'] is True and declaration['actual_acceptance']['all_four_old_block_scores_at_or_below_fixed_origin24_ceiling'] is False)
check('unbound future release',cfg['root_release'] is None and cfg['source_review'] is None and ready['execution_released'] is False)
report={**ready['source_review_required'], 'reviewer_source':pin(__file__),'ready':pin(D/'READY.json'),
    'reuse_implementation':ready['reuse_implementation'],'journal_implementation':ready['journal_implementation'],
    'configuration_template':ready['configuration_template'],'preparation':ready['preparation'],'controls':ready['controls'],
    'reuse_inventory':ready['reuse_inventory'],'core_review':ready['core_review'],'parent_journal_review':ready['parent_journal_review'],
    'findings':[
      'Immutable terminal aggregate admissions precede exact metadata indexing; other vector bodies load only on an exact full-vector or effective-metal identity.',
      'All27 seeds retain original scores and native provenance. The failed full pilot is available native data, not an initial or automatically feasible incumbent.',
      'First saved proposal requires exact origin, probe identities, affine arrays, retained cones and settings. Driver allows zero-call import only at step1.',
      'Actual selection checks all25 energy/force gates and strict q improvement. Four old-score ceilings are descriptive; electrical losses remain reported.',
      'Only complete eligible trial batches can select; derivative-role requests alone never select. Current and earlier-fraction ties are preserved.',
      'Eight outer steps,200 new vector keys and6496 invocation intents remain enforced. Imported costs are separate; no old mutable ledger is read.',
      'Native invocation, immediate stop, attempt-intent accounting and journal code are unchanged. No retries or partial-batch selection were introduced.'
    ],
    'checks_passed':len(checks),'checks':checks,'exact_pins_verified':len(verified),
    'author_controls_rerun':False,'core_tests_rerun':False,'solver_calls':0,'native_calls':0,
    'limitations':['Source and saved metadata review only; no actual concurrent run, new optimizer proposal or new scientific result was exercised.',
                   'CORE numerical formulation relies on separate independent review037e41a6; prior reviewer authorship of reused pilot arithmetic and conic helpers is disclosed.',
                   'Prior numerical admissions are reused, not repeated; no new electrical physical acceptance threshold or fresh-validation release follows.'],
    'cost':dict(wall_seconds=time.perf_counter()-tick,CPU_seconds=time.process_time()-cpu,numerical_threads=1)}
report.pop('recommend_additional_pins',None)
W.mkdir(parents=True,exist_ok=True)
p=W/'REVIEW.json';assert not p.exists();p.write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(dict(review=pin(p),checks=len(checks),pins=len(verified),cost=report['cost'])))
