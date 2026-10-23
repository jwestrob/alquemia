#!/usr/bin/env python3
"""Only the added runtime binary-pin guard; no native calls or broad rerun."""
import ast
import importlib.util
import json
from pathlib import Path
import tempfile
import time

t=time.monotonic();cpu=time.process_time();here=Path(__file__).resolve().parent
s=importlib.util.spec_from_file_location('orientation_pin_control',here/'RUN.py');m=importlib.util.module_from_spec(s);s.loader.exec_module(m)
new=(here/'RUN.py').read_text();old=(here/'draft_v1/RUN.py').read_text()
added="    for p in cfg['executables'].values():verify(p)\n"
at=old.index("    require(set(bound['rows'])==set(CELLS)",old.index('def bound_inputs'))
assert old[:at]+added+old[at:]==new
tree=ast.parse(new);fn=next(x for x in tree.body if isinstance(x,ast.FunctionDef) and x.name=='bound_inputs')
loop=next(x for x in fn.body if isinstance(x,ast.For) and 'executables' in ast.unparse(x))
code=compile(ast.Module(body=[loop],type_ignores=[]),'<actual runtime pin loop>','exec')
checks={'only_added_runtime_guard':True,'actual_guard_extracted':True}
with tempfile.TemporaryDirectory() as td:
    p=Path(td)/'nonexecutable_fixture';p.write_text('saved binary identity stand-in\n');cfg={'executables':{'analyze':m.pin(p)}}
    exec(code,{'cfg':cfg,'verify':m.verify});checks['unchanged_pin_passes']=True
    p.write_text('changed fixture bytes\n')
    try:exec(code,{'cfg':cfg,'verify':m.verify})
    except ValueError:checks['changed_pin_rejected']=True
    else:raise AssertionError('Changed runtime pin accepted')
w=m.WORK/'paired_water_orientation_validation_v1/model';prep=m.read(m.pin(w/'preparation_v2/PREPARATION.json'));bound=m.read(prep['bound'])
assert prep['status']=='saved_preparation_passed' and prep['native_calls']==0 and bound['implementation']==m.pin(here/'RUN.py')
checks['actual_revised_preparation_passed']=True
prior=m.read(m.pin(w/'source_controls_v1/RESULT.json'));assert prior['status']=='passed' and prior['check_count']==70 and prior['implementation']['sha256']==m.pin(here/'draft_v1/RUN.py')['sha256'];checks['unchanged_prior_controls_bound']=True
result={'status':'passed','implementation':m.pin(here/'RUN.py'),'control_source':m.pin(__file__),'preparation':m.pin(w/'preparation_v2/PREPARATION.json'),
    'prior_controls':m.pin(w/'source_controls_v1/RESULT.json'),'preserved_draft_source':m.pin(here/'draft_v1/RUN.py'),'checks':checks,'check_count':len(checks),
    'native_calls':0,'reserved_outcomes_opened':False,'local_threads':1,'elapsed_seconds':time.monotonic()-t,'CPU_seconds':time.process_time()-cpu}
out=w/'runtime_pin_controls_v2';out.mkdir(parents=True,exist_ok=False);m.save(out/'RESULT.json',result);print(json.dumps({'result':m.pin(out/'RESULT.json')}))
