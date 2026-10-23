#!/usr/bin/env python3
"""Bounded geometry/guard controls; never read a candidate or reserved outcome."""
import ast
import copy
import importlib.util
import json
from pathlib import Path
import time

start=time.monotonic();cpu=time.process_time();here=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('stageA_model_source',here/'RUN.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
cfg=m.default_config();data,geometries=m.inputs(cfg);checks={}
def check(name,ok):
    checks[name]=bool(ok);assert ok,name
def rejects(name,fn):
    try:fn()
    except (ValueError,AssertionError):checks[name]=True
    else:raise AssertionError(name)
for gid,g in geometries.items():
    parent=m.parent_id(gid);row=data['rows'][parent];text=m.typed_xyz(g,row);lines=text.splitlines()[1:]
    check('typed_coordinates_'+gid,[[float(v) for v in x.split()[2:5]] for x in lines]==g['positions_A'])
    check('types_'+gid,[int(x.split()[5]) for x in lines]==row['input']['types'])
    check('bonds_'+gid,[[int(v)-1 for v in x.split()[6:]] for x in lines]==row['input']['bonds_zero_based'])
    peer=gid[:-2]+('Dy' if gid.endswith('La') else 'La')
    check('paired_coordinates_'+gid,g['positions_A']==geometries[peer]['positions_A'])
    check('parent_key_source_'+gid,parent in data['rows'])
check('reserved16',len(m.RESERVED)==len(set(m.RESERVED))==16)
check('consumed12',len(m.CHECKS)==len(set(m.CHECKS))==12)
check('FD1344',6*sum(data['rows'][g]['geometry']['atoms'] for g in m.CHECKS)==1344)
check('two_reserved_programs_per_cell',2*len(m.RESERVED)==cfg['maximum_reserved_programs']==32)
first=next(iter(geometries));g=copy.deepcopy(geometries[first]);row=copy.deepcopy(data['rows'][m.parent_id(first)])
bad=copy.deepcopy(g);bad['positions_A'][1][0]+=.01
rejects('changed_reserved_coordinate_rejected',lambda:m.typed_xyz(bad,row))
badrow=copy.deepcopy(row);badrow['input']['types'][1]+=1
rejects('changed_type_rejected',lambda:m.typed_xyz(g,badrow))
tree=ast.parse((here/'RUN.py').read_text());bind=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='bind')
guard=next(n.value.args[0] for n in bind.body if isinstance(n,ast.Expr) and isinstance(n.value,ast.Call) and getattr(n.value.func,'id',None)=='require')
compiled=compile(ast.Expression(guard),'<actual candidate-entry condition>','eval')
base={'status':'bounded_paired_donor_complete','error':None,'candidate_frozen':True,'development_checks_passed':True}
check('positive_entry_condition_only',eval(compiled,{'candidate':base}))
for key,value in (('status','running'),('error','failure'),('candidate_frozen',False),('development_checks_passed',False)):
    bad=dict(base);bad[key]=value;check('entry_rejects_'+key,not eval(compiled,{'candidate':bad}))
old={'energy':0.,'forces':[[0.,0.,0.]]};good={'fd_passed':True,'energy':1e-4,'forces':[[1e-3,0.,0.]]};m.replay_check(good,old);check('unchanged_numerical_limits_pass',True)
for name,change in (('Cartesian',{'fd_passed':False}),('energy',{'energy':1.0001e-4}),('force',{'forces':[[1.0001e-3,0.,0.]]})):
    bad=dict(good);bad.update(change);rejects('numerical_rejects_'+name,lambda bad=bad:m.replay_check(bad,old))
check('run_source_AST',True)
result={'status':'passed','implementation':m.pin(here/'RUN.py'),'control_source':m.pin(__file__),'specification':m.pin(here/'SPEC.md'),
    'reserved_panel':cfg['reserved_panel'],'training_result':cfg['training_result'],'candidate_config':cfg['candidate_config'],
    'checks':checks,'check_count':len(checks),'candidate_bound':False,'candidate_result_opened':False,'reserved_outcomes_opened':False,
    'native_calls':0,'local_threads':1,'elapsed_seconds':time.monotonic()-start,'CPU_seconds':time.process_time()-cpu}
out=m.WORK/'paired_donor_validation_stageA_v1/model/source_controls_v2';out.mkdir(parents=True,exist_ok=False);m.save(out/'RESULT.json',result)
print(json.dumps({'result':m.pin(out/'RESULT.json'),'checks':len(checks),'elapsed_seconds':result['elapsed_seconds']}))
