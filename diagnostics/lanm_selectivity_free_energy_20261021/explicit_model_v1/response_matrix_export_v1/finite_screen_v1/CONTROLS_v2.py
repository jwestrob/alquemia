#!/usr/bin/env python3
"""Narrow JSON-list replay controls; inherit unchanged kernel's16 controls."""
import ast
import copy
import importlib.util
import json
from pathlib import Path
import time

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[4]
W=ROOT/'workspaces/lanm_selectivity_free_energy_20261021/explicit_model_v1/response_matrix_export_v1/finite_screen_v1'
s=importlib.util.spec_from_file_location('screen_v3_controls',HERE/'SCREEN_v3.py');m=importlib.util.module_from_spec(s);s.loader.exec_module(m)
tick,cpu=time.monotonic(),time.process_time();checks=[]
a=ast.parse((HERE/'SCREEN_v2.py').read_text());b=ast.parse((HERE/'SCREEN_v3.py').read_text())
assert [x.name for x,y in zip(a.body,b.body) if ast.dump(x)!=ast.dump(y)]==['zero_replay'];checks.append('only JSON-array zero-replay function changed')
prep=json.loads((W/'preparation_v3/MANIFEST.json').read_text());target=m.read(prep['targets']);g=target['groups']
observed={'original_energy_checks':[dict(r,candidate_kcal_mol=r['selected17676_model_kcal_mol']) for r in g['original_energy_force']['energy_checks']],
          'pair_dipoles':{r['id']:{'change_from_baseline':[0.,0.,0.]} for r in g['pair_and_collective_dipoles']['pair_cycles']},
          'CPB_dipoles':{parent:{k:{'change_from_baseline':[0.,0.,0.]} for k in ('C','P','B')} for parent in g['pair_and_collective_dipoles']['parents']},
          'trimer_E3_mu3':{cid:{'energy_kcal_mol':{'change_from_baseline':0.},'mu_Debye':{'change_from_baseline':[0.,0.,0.]}} for cid in g['trimer_energy_and_dipoles']['cycles']},'field_rows':[None]*28}
baseline={cid:{'mu_print_bound':1e-5} for cid in (*prep['rows'],*prep['static'])}
m.zero_replay(observed,target,baseline);checks.append('all declared JSON-list zero-cycle formats accepted')
for kind in ('pair','CPB','trimer'):
    bad=copy.deepcopy(observed)
    if kind=='pair':next(iter(bad['pair_dipoles'].values()))['change_from_baseline']=[1.,0.,0.]
    elif kind=='CPB':next(iter(bad['CPB_dipoles'].values()))['B']['change_from_baseline']=[1.,0.,0.]
    else:next(iter(bad['trimer_E3_mu3'].values()))['mu_Debye']['change_from_baseline']=[1.,0.,0.]
    try:m.zero_replay(bad,target,baseline)
    except ValueError:checks.append(kind+' discrepancy rejected')
    else:raise AssertionError(kind+' discrepancy incorrectly accepted')
old=m.read(m.pin(W/'source_controls_v1/RESULT.json'));assert old['status']=='passed' and old['count']==16
out=W/'source_controls_v2';out.mkdir(exist_ok=False)
result={'status':'passed','implementation':m.pin(__file__),'screen':m.pin(HERE/'SCREEN_v3.py'),'preparation':m.pin(W/'preparation_v3/MANIFEST.json'),'inherited_unchanged_kernel_controls':m.pin(W/'source_controls_v1/RESULT.json'),'inherited_count':16,'new_checks':checks,'new_count':len(checks),'actual_candidate_screen_executed':False,'new_native_calls':0,'cost':{'wall_seconds':time.monotonic()-tick,'CPU_seconds':time.process_time()-cpu,'threads':1}}
m.save(out/'RESULT.json',result);print(json.dumps({'result':m.pin(out/'RESULT.json'),'new_count':len(checks),'cost':result['cost']}))
