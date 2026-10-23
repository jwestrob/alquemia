"""One finite stdlib source/mock check; never imports or calls a numeric backend."""
import ast
import copy
import importlib.util
import json
from pathlib import Path
import sys
import time
import traceback

HERE=Path(__file__).resolve().parent
tick,cpu=time.monotonic(),time.process_time()
spec=importlib.util.spec_from_file_location('nominal_worker',HERE/'RUN.py')
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
checks=[]

def check(name, fn):
    try:
        fn();checks.append(dict(name=name,passed=True))
    except Exception as exc:
        checks.append(dict(name=name,passed=False,error=repr(exc),traceback=traceback.format_exc()))

def require(condition):
    assert condition

def reject(fn):
    try:fn()
    except (AssertionError,ValueError,KeyError):return
    raise AssertionError('invalid fixture accepted')

base=dict(U_B=3.,U_R=7.,delta_odd=.25,e_app=None,application_error_bound_established=False)
check('nominal formula; odd term outside relative upper',lambda:require(m.nominal_parameters(base,base)['tau']==1e-6*10+.25))
for key,value in [('e_app',0.),('application_error_bound_established',True),('delta_odd',-.1),('U_R',-1.),('U_B',float('nan')),('delta_odd',float('inf')),('U_B',-7.)]:
    bad=dict(base,**{key:value});check('reject '+key+' '+str(value),lambda bad=bad:reject(lambda:m.nominal_parameters(bad,bad)))
check('reject mismatched reviewed odd bound',lambda:reject(lambda:m.nominal_parameters(base,dict(base,delta_odd=.2))))

def p(name):return dict(path=name,sha256=name)
cfg={key:p(key) for key in ('nominal_contract','pme_result','pme_root_review','pme_bound_review','composition_manifest','pme_export','campaign_result','campaign_configuration','candidate','terminal_fit_review')}
cfg['nominal_contract']['sha256']=m.NOMINAL_SHA;cfg['metal']='La'
norms=dict(arithmetic='ordinary_double_not_outward',q_negative_count=0,q_min=0.,q_max=.1,q_abs_max=.1,q_conjugate_difference_max=.01,U_R_nominal=7.,evenness_perturbation_nominal=.25)
row=dict(status='complete',site_checks=dict(passed=True),export=cfg['pme_export'],inputs={'system.xyz':p('geometry')},operator_norms=norms)
candidate={key:cfg[source] for key,source in [('result','campaign_result'),('configuration','campaign_configuration'),('selected','candidate'),('actual_review','terminal_fit_review')]}
objects=dict(nominal_contract=dict(recommended_contract=dict(tau_nominal='1e-6*(U_B+U_R)+delta_odd',e_app=None,application_error_bound_established=False)),
    pme_root_review=dict(status='admitted_full_box_operator_initialization_nominal_data',execution=cfg['pme_result'],independent_actual_review=cfg['pme_bound_review'],composition=cfg['composition_manifest'],e_app=None,application_error_bound_established=False,cells={'La':copy.deepcopy(row)}),
    pme_result=dict(status='full_box_operator_initialization_complete',error=None,composition=cfg['composition_manifest'],candidate=candidate,cells={'La':copy.deepcopy(row)}),
    pme_bound_review=dict(result=cfg['pme_result'],cells={'La':copy.deepcopy(row)}))
inputs={key:cfg[key] for key in ('nominal_contract','pme_result','pme_root_review','pme_bound_review','metal')}
inputs.update(base);inputs.update(geometry=p('geometry'),operator_producer=p('producer'),operator_preparation_result=p('preparation'))
review=copy.deepcopy(inputs)
original_read,original_verify=m.read,m.verify
m.read=lambda pin:objects[pin['path']]
m.verify=lambda pin:require(isinstance(pin,dict) and set(pin)=={'path','sha256'})
check('synthetic matched nominal provenance',lambda:m.nominal_provenance(cfg,inputs,review))
for label,path,value in [
    ('wrong snapshot',('pme_root_review','cells','La','export'),p('other')),
    ('wrong candidate',('pme_result','candidate','selected'),p('other')),
    ('wrong geometry',('pme_result','cells','La','inputs','system.xyz'),p('other')),
    ('negative q',('pme_result','cells','La','operator_norms','q_negative_count'),1),
    ('failed initialization',('pme_result','error'),'failure'),
    ('invented application bound',('pme_root_review','e_app'),0.),
    ('wrong independent result',('pme_bound_review','result'),p('other'))]:
    parent=objects
    for part in path[:-1]:parent=parent[part]
    prior=parent[path[-1]];parent[path[-1]]=value
    check('reject '+label,lambda:reject(lambda:m.nominal_provenance(cfg,inputs,review)))
    parent[path[-1]]=prior
m.read,m.verify=original_read,original_verify

old=ast.parse((HERE.parent/'static_factor_worker_v1/RUN.py').read_text());new=ast.parse((HERE/'RUN.py').read_text())
def functions(tree):return {x.name:x for x in tree.body if isinstance(x,ast.FunctionDef)}
oldf,newf=functions(old),functions(new)
for name in ('pin','verify','read','load','supervise'):
    check('unchanged '+name,lambda name=name:require(ast.dump(oldf[name])==ast.dump(newf[name])))
def calls(node,tail):
    return [ast.dump(x) for x in ast.walk(node) if isinstance(x,ast.Call) and isinstance(x.func,ast.Attribute) and x.func.attr==tail]
for name in ('setrlimit','factor_once','audit_factor','ten_probes','validate_lower'):
    check('unchanged '+name+' calls',lambda name=name:require(calls(oldf['worker'],name)==calls(newf['worker'],name)))
check('compile source',lambda:compile((HERE/'RUN.py').read_text(),str(HERE/'RUN.py'),'exec'))
check('no numerical module imported',lambda:require('numpy' not in sys.modules and 'scipy' not in sys.modules and 'sksparse' not in sys.modules))
result=dict(status='passed_nominal_worker_source_controls' if all(x['passed'] for x in checks) else 'failed_nominal_worker_source_controls',
    implementation=m.pin(HERE/'RUN.py'),test_source=m.pin(__file__),specification=m.pin(HERE/'SPEC.md'),checks=checks,
    cost=dict(wall_seconds=time.monotonic()-tick,CPU_seconds=time.process_time()-cpu,threads=1),
    native_calls=0,matrix_assemblies=0,matrix_array_reads=0,numeric_factorizations=0,solves=0)
out=Path(sys.argv[1]);out.parent.mkdir(parents=True,exist_ok=True)
with out.open('x') as handle:handle.write(json.dumps(result,indent=2,allow_nan=False)+'\n')
print(json.dumps(dict(status=result['status'],checks=len(checks),cost=result['cost'])))
raise SystemExit(0 if all(x['passed'] for x in checks) else 1)
