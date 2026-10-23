#!/usr/bin/env python3
"""Four source-review binding controls; no repeated scoring tests."""
import copy
import importlib.util
import json
from pathlib import Path
import time
HERE=Path(__file__).resolve().parent
s=importlib.util.spec_from_file_location('collector_join',HERE/'RUN.py');m=importlib.util.module_from_spec(s);s.loader.exec_module(m)
start,cpu=time.monotonic(),time.process_time();out=m.WORK/'collection/review_join_controls_v2';out.mkdir(parents=True,exist_ok=False)
cfg=m.load(HERE/'CONFIG_SOURCE.json');rr=m.read(cfg['reference_source_review']);original=m.read
checks=[{'name':'actual_source_configuration_review_join','passed':True}]
for field in ('implementation','configuration','prepared_manifest'):
    bad=copy.deepcopy(rr);bad[field]=None
    m.read=lambda p,bad=bad:bad if p==cfg['reference_source_review'] else original(p)
    try:m.load(HERE/'CONFIG_SOURCE.json')
    except ValueError as exc:assert str(exc)=='Reference review source/preparation differs';checks.append({'name':'reject_changed_review_'+field,'passed':True})
    else:raise AssertionError(field)
m.read=original
result={'status':'review_join_controls_passed','implementation':m.pin(HERE/'RUN.py'),'configuration':m.pin(HERE/'CONFIG_SOURCE.json'),'test_source':m.pin(__file__),'checks':checks,'check_count':len(checks),'new_native_calls':0,'fit_or_fresh_outcomes_read':False,'local_threads':1,'wall_seconds':time.monotonic()-start,'CPU_seconds':time.process_time()-cpu}
print(json.dumps(m.save(out/'RESULT.json',result)))
