"""Focused exponent-boundary controls; no actual fixture or physical calculation."""
import decimal as dc
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
import time

p=Path(__file__).resolve().with_name('TAIL.py')
s=importlib.util.spec_from_file_location('tail_edge',p)
t=importlib.util.module_from_spec(s);sys.modules[s.name]=t;s.loader.exec_module(t)
start=time.process_time(); ar=t.Arithmetic(40,emin=-10,emax=10)
checks=[]
x=ar.point('-23.025850929940456840179914546843642076011')
center=ar.near.exp(x.lo)
assert center == dc.Decimal('1e-10')
checks.append({'name':'correctly_rounded_center_is_normal_minimum','passed':True})
assert ar.near.next_minus(center).adjusted() < ar.near.Emin
checks.append({'name':'outward_lower_padding_is_subnormal','passed':True})
try: ar.exp(x)
except t.BoundError: pass
else: raise AssertionError('subnormal padding was not rejected')
checks.append({'name':'padding_exponent_rejection_after_valid_input_and_exp','passed':True})
# A valid normal result away from the boundary must remain admitted.
iv=ar.exp(ar.point('-20'))
assert iv.lo>0 and iv.lo<iv.hi
checks.append({'name':'nearby_normal_positive_exp_is_enclosed','passed':True})
out={'status':'passed_focused_exponent_boundary_controls','implementation':{'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()},
     'controls_implementation':{'path':str(Path(__file__).resolve()),'sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},
     'checks':checks,'count':len(checks),'cpu_seconds':time.process_time()-start,'new_native_calls':0,'actual_box_evaluations':0,
     'prior_control_note':'CONTROLS.py reject_subnormal_padding rejects its out-of-domain input before sqrt; this focused control separately exercises actual post-exp padding rejection.'}
dest=Path(sys.argv[1]);dest.parent.mkdir(parents=True,exist_ok=True)
with dest.open('x') as f:json.dump(out,f,indent=2,sort_keys=True);f.write('\n')
print(json.dumps({'status':out['status'],'count':len(checks),'cpu_seconds':out['cpu_seconds']}))
