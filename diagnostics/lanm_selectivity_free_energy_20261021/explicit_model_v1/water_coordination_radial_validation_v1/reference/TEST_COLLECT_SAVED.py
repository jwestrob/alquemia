#!/usr/bin/env python3
"""One saved-only exception-boundary regression; no native/candidate data."""
import importlib.util
import json
import os
from pathlib import Path
import time
from types import SimpleNamespace

HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('radial_ref',HERE/'RUN.py')
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
assert all(os.environ.get(k)=='1' for k in m.THREADS)
wall,cpu=time.monotonic(),time.process_time()
out=m.W/'explicit_model_v1/water_coordination_radial_validation_v1/reference/collector_controls_v2'
out.mkdir(parents=True,exist_ok=False)
cells=[{'cell_id':cid,'geometry_id':gid,'geometry':{'geometry_id':gid}} for cid,gid in zip(m.CELLS,m.GIDS)]
geometry={'rows':{gid:{'input':{},'bodies':{}} for gid in m.GIDS}}
paths={};directories={}
for i,cid in enumerate(m.CELLS):
    folder=out/'synthetic'/cid;folder.mkdir(parents=True)
    path=folder/'RESULT.json'
    path.write_text('{"truncated":' if i==0 else json.dumps({'status':'failed','allocated_cpu_seconds':i}))
    paths[cid]=str(path);directories[cid]=str(folder)
inventory=out/'INVENTORY.json';m.save(inventory,paths)
m.gate=lambda *args: ({},SimpleNamespace(),{},geometry,cells,{'anchor_reference':{}},{'task_output_directories':directories})
receipt=m.collect(None,None,None,inventory,out/'COLLECTED.json')
r=m.read(receipt);checks=[]
def check(name,value):
    checks.append({'name':name,'passed':bool(value)})
    assert value,name
check('exact_four_slots_retained',r['declared']==4 and set(r['rows'])==set(m.GIDS))
check('malformed_cell_unavailable_with_raw_authority',r['rows'][m.GIDS[0]]['reference'] is None and 'JSONDecodeError' in r['rows'][m.GIDS[0]]['reason'] and r['rows'][m.GIDS[0]]['raw_result']==m.pin(paths[m.CELLS[0]]))
check('all_three_subsequent_cells_processed',all('Native reference admission/provenance failed' in r['rows'][gid]['reason'] for gid in m.GIDS[1:]))
check('no_partial_pass_or_fake_values',r['available']==0 and r['status']=='unavailable' and all(v['reference'] is None for v in r['rows'].values()))
check('subsequent_cost_evidence_retained',r['allocated_native_CPU_seconds']==6)
check('zero_native_calls',r['native_calls_this_collection']==0)
result={'status':'collector_exception_boundary_controls_passed','implementation':m.pin(HERE/'RUN.py'),'test_source':m.pin(__file__),'collection':receipt,'checks':checks,'check_count':len(checks),'native_calls':0,'candidate_or_new_validation_outcomes_read':False,'local_threads':1,'wall_seconds':time.monotonic()-wall,'CPU_seconds':time.process_time()-cpu}
print(json.dumps(m.save(out/'RESULT.json',result)))
