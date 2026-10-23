#!/usr/bin/env python3
"""Exercise actual binder loop on saved initial rows; create no bound candidate."""
import argparse
import ast
import copy
import importlib.util
import json
from pathlib import Path
import time

HERE=Path(__file__).resolve().parent
s=importlib.util.spec_from_file_location('continuation_transfer_loop_control',HERE/'RUN.py');r=importlib.util.module_from_spec(s);s.loader.exec_module(r)


def main(out):
    started=time.monotonic();out.mkdir(parents=True,exist_ok=False)
    data=r.loaded(r.verify(r.PINS['pins']['prepared32']));cfg=data['config']
    parent=r.read(r.PINS['pins']['continuation_prior_result']);prep=r.read(r.PINS['pins']['continuation_preparation'])
    assert prep['status']=='prepared' and prep['new_native_program_invocations']==0
    selected=prep['initial'];assert selected['parameters']==parent['selected']['parameters']
    ecfg=r.read(r.PINS['pins']['continuation_evaluator_config']);previous=r.read(r.PINS['pins']['previous_transfer_result'])
    native=r.module(cfg['native_helper'],'loop_saved_native');calls=[]
    def forbidden(*a,**k):calls.append(True);raise AssertionError('No native calls permitted')
    native.invoke=forbidden
    def execute(source,directory):
        tree=ast.parse(source.read_text());bind=next(n for n in tree.body if isinstance(n,ast.FunctionDef)and n.name=='bind')
        loop=next(n for n in bind.body if isinstance(n,ast.For)and isinstance(n.iter,ast.Attribute)and n.iter.attr=='CELLS')
        loop=copy.deepcopy(loop);loop.iter=ast.Name(id='test_cells',ctx=ast.Load());program=ast.fix_missing_locations(ast.Module(body=[loop],type_ignores=[]))
        seen=[]
        def tracked(gid,value,row,inp,econfig,path,prior_transfer,immutable_parent):
            assert immutable_parent is parent,'Stable17589 parent was shadowed by prior reuse'
            branch=r.selected_authority(gid,value,row,inp,econfig,path,prior_transfer,immutable_parent)
            assert branch=='inherited17589_native_cache';seen.append(gid);return branch
        namespace=dict(r.__dict__,data=data,cfg=cfg,values=selected['parameters'],selected=selected,ecfg=ecfg,
            candidate_path=directory/'NO_CANDIDATE_CREATED.json',previous=previous,continuation_parent=parent,prior=parent,
            inputs={},reusable={},authority=[],branches={},out=directory,native=native,test_cells=r.WATER,selected_authority=tracked)
        try:exec(compile(program,str(source),'exec'),namespace)
        except AssertionError as exc:return {'passed':False,'error':str(exc),'rows_admitted_before_failure':seen}
        assert namespace['continuation_parent']is parent and len(namespace['reusable'])==8
        assert all(namespace['reusable'][k]['value']['native']==selected['model'][k]['native']for k in r.WATER)
        return {'passed':True,'rows_admitted':seen,'branches':namespace['branches']}
    old=execute(HERE/'draft_v1/RUN.py',out/'old_loop_reproduction')
    assert not old['passed'] and old['error']=='Stable17589 parent was shadowed by prior reuse' and len(old['rows_admitted_before_failure'])==1
    new=execute(HERE/'RUN.py',out/'corrected_loop')
    assert new['passed'] and set(new['rows_admitted'])==set(r.WATER)and not calls
    for p in (r.PINS['pins']['continuation_prior_result'],parent['evaluator']['ledger'],r.PINS['pins']['continuation_preparation']):r.verify(p)
    assert not list(out.rglob('BOUND.json'))
    result={'status':'passed','implementation':r.pin(__file__),'source':r.pin(HERE/'RUN.py'),'old_source':r.pin(HERE/'draft_v1/RUN.py'),
        'parent':r.PINS['pins']['continuation_prior_result'],'saved_preparation':r.PINS['pins']['continuation_preparation'],
        'old_loop_reproduction':old,'corrected_loop':new,'new_native_calls':0,'candidate_bound':False,
        'continuation_execution_outcome_opened':False,'reserved_validation_opened':False,
        'elapsed_seconds':time.monotonic()-started,'maximum_local_threads':1}
    r.save(out/'RESULT.json',result);print(json.dumps({'status':'passed','result':r.pin(out/'RESULT.json')}))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--out',required=True,type=Path);a=p.parse_args();main(a.out)
