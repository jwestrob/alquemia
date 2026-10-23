"""Eight-parameter solver and actual donor driver clipping/selection controls; no engines."""
import ast
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
import time
import numpy as np

start=time.monotonic();cpu=time.process_time();here=Path(__file__).resolve().parent
def pin(p):
    p=Path(p);return {'path':str(p.resolve()),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def load(p,name):
    s=importlib.util.spec_from_file_location(name,p);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
directions=here.parent/'paired_water_six_parameter_v1/DIRECTIONS.py';new=load(directions,'donor_directions');oldpath=here.parent/'paired_water_repair_v1/joint_refinement_v1/REFINE.py';old=load(oldpath,'four_directions')
checks={}
def check(name,ok):checks[name]=bool(ok);new.require(ok,name)
def function(p,name):
    t=p.read_text();n=next(x for x in ast.parse(t).body if isinstance(x,ast.FunctionDef)and x.name==name);return ast.get_source_segment(t,n)
check('exact_frozen_generic_directions',pin(directions)['sha256']=='741167a8c5ce3c9e07af6a94baf40894c277b5c47950f31b3c526da125ae11ca')
for n in (8,):
    J=np.vstack((np.diag(np.arange(1,n+1)),np.full((1,n),.3)));r=np.linspace(-.2,.2,n+1)
    for label,z in [('interior',np.full(n,.5)),('clipped',np.linspace(0,1,n))]:
        scaled,bounded,diag=new.directions(J,r,z);lo=np.maximum(-.025,-z);hi=np.minimum(.025,1-z)
        check(f'{n}_{label}_feasible',all(np.all(d>=lo-1e-12)and np.all(d<=hi+1e-12)for d in(scaled,bounded)))
        check(f'{n}_{label}_bounded_no_worse',np.linalg.norm(r+J@bounded)**2<=np.linalg.norm(r+J@scaled)**2+1e-12)
        gradient=J.T@(r+J@bounded);kkt=[]
        for j,v in enumerate(bounded):
            kkt.append(abs(gradient[j])<1e-10 if lo[j]+1e-10<v<hi[j]-1e-10 else gradient[j]>=-1e-10 if abs(v-lo[j])<1e-10 else gradient[j]<=1e-10)
        check(f'{n}_{label}_bounded_KKT',all(kkt))
        if n==4:
            a,b,c=old.directions(J,r,z);check(f'four_dimensional_exact_regression_{label}',np.array_equal(a,scaled)and np.array_equal(b,bounded)and c==diag)
tree=ast.parse((here/'REFINE.py').read_text());main=next(x for x in tree.body if isinstance(x,ast.FunctionDef)and x.name=='main')
loop=next(x for x in ast.walk(main)if isinstance(x,ast.For)and isinstance(x.target,ast.Name)and x.target.id=='j'and any(isinstance(y,ast.Call)and isinstance(y.func,ast.Attribute)and isinstance(y.func.value,ast.Name)and y.func.value.id=='points'and y.func.attr=='append'for y in ast.walk(x)))
initial=np.array([3.6318,.1518,3.3,.2987,3.347237658604727,.139,3.2657502966294456,.139]);lower=np.array([.85*v if j%2==0 else .05 for j,v in enumerate(initial)]);upper=np.array([1.15*v if j%2==0 else .5 for j,v in enumerate(initial)]);width=upper-lower
for label,z in [('interior',np.full(8,.5)),('clipped',np.array([0,.002,.998,1,0,1.,.003,.997]))]:
    center=lower+z*width;env={'np':np,'z':z,'center':center,'width':width,'lower':lower,'upper':upper,'points':[]}
    exec(compile(ast.fix_missing_locations(ast.Module(body=[loop],type_ignores=[])),'actual_driver_derivatives','exec'),env)
    points=env['points'];check(label+'_sixteen_derivative_slots',len(points)==16)
    for j in range(8):
        a=np.asarray(points[2*j]['parameters']);b=np.asarray(points[2*j+1]['parameters']);span=(b[j]-a[j])/width[j]
        check(f'{label}_{j}_actual_clipped_width',span>0 and abs(span-(min(1,z[j]+.005)-max(0,z[j]-.005)))<1e-12)
        check(f'{label}_{j}_only_one_parameter_changed',all(a[k]==center[k]and b[k]==center[k]for k in range(8)if k!=j))
selection=next(x for x in ast.walk(main)if isinstance(x,ast.For)and isinstance(x.target,ast.Name)and x.target.id=='value'and isinstance(x.iter,ast.Name)and x.iter.id=='completed')
current={'score':{'objective':1.}};first={'score':{'objective':.5}};env={'selected':current,'completed':[{'score':{'objective':1.}}]}
code=compile(ast.fix_missing_locations(ast.Module(body=[selection],type_ignores=[])),'actual_selection','exec');exec(code,env);check('current_wins_equal_objective',env['selected']is current)
env={'selected':current,'completed':[first,{'score':{'objective':.5}},{'score':{'objective':.7}}]};exec(code,env);check('first_strict_minimum_wins_trial_ties',env['selected']is first)
oldsource=here.parent/'paired_water_six_parameter_v1/REFINE.py';oldmain=next(x for x in ast.parse(oldsource.read_text()).body if isinstance(x,ast.FunctionDef)and x.name=='main')
def batch_text(source,node):return ast.get_source_segment(source.read_text(),next(x for x in ast.walk(node)if isinstance(x,ast.FunctionDef)and x.name=='batch'))
expected=batch_text(oldsource,oldmain).replace('<=384','<=480').replace('obj.evaluate(ev.data,m,body)','obj.evaluate(ev.data,m,body,ev.water_baseline)')
check('batch_only_changes_vector_cap_and_explicit_water_argument',batch_text(here/'REFINE.py',main)==expected)
out=Path(sys.argv[1]).resolve();out.mkdir(parents=True,exist_ok=False)
result={'status':'passed','source':pin(__file__),'driver':pin(here/'REFINE.py'),'directions':pin(directions),'original_directions':pin(oldpath),'checks':checks,'synthetic_numerical_controls_not_molecular_results':True,'new_native_calls':0,'maximum_local_numerical_threads':1,'elapsed_seconds':time.monotonic()-start,'parent_cpu_seconds':time.process_time()-cpu}
(out/'RESULT.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({'status':'passed','checks':len(checks),'result':pin(out/'RESULT.json')}))
