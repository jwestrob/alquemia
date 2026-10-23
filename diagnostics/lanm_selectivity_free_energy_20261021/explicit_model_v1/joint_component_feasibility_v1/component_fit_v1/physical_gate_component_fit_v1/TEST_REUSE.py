"""Focused exact-reuse controls. No optimizer, molecular helper, or native call."""
import ast
import copy
import importlib.util
import json
import os
from pathlib import Path
import tempfile
import threading
import time
from types import SimpleNamespace

HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('reuse_test_adapter',HERE/'ADAPTER.py')
A=importlib.util.module_from_spec(spec);spec.loader.exec_module(A)


def run(out):
    tick,cpu=time.monotonic(),time.process_time();checks=[]
    def check(name,value):
        if not value:raise AssertionError(name)
        checks.append(name)
    def rejects(name,fn):
        try:fn()
        except ValueError:checks.append(name);return
        raise AssertionError('No rejection: '+name)
    check('one thread',all(os.environ.get(k)=='1' for k in (*A.THREADS,'BLIS_NUM_THREADS')))
    parent=HERE.parent/'persistence_performance_v1/future_adapter_v1'
    check('unchanged journal bytes',(HERE/'JOURNAL.py').read_bytes()==(parent/'JOURNAL.py').read_bytes())
    def method(path,name):
        return ast.dump(next(n for n in ast.walk(ast.parse(path.read_text())) if isinstance(n,ast.FunctionDef) and n.name==name),include_attributes=False)
    check('unchanged native invocation/locking/durable attempt path',method(HERE/'ADAPTER.py','_invoke')==method(parent/'ADAPTER.py','_invoke'))
    check('unchanged explicit role authorization',method(HERE/'ADAPTER.py','authorize')==method(parent/'ADAPTER.py','authorize'))
    check('unchanged stop persistence',method(HERE/'ADAPTER.py','_stop')==method(parent/'ADAPTER.py','_stop'))
    check('limits preserved',A.LIMITS=={'outer_steps':8,'new_parameter_vectors':200,'new_native_programs':6496,'reused_origin_vectors':1,'reused_initial_derivative_vectors':24,'workers':12,'numerical_threads_per_process':1,'aggregate_login_threads':16,'native_timeout_seconds':300})
    # Manufactured46-cell identity schema, with five force programs per metal.
    base=SimpleNamespace(Evaluator=object,effective=lambda values,metal,domain:[values[0 if metal=='La' else 1]])
    Cls=A.evaluator_class(base)
    rows={metal+str(i):{'metal':metal,'needs_gradient':i<5} for metal in ('La','Dy') for i in range(23)}
    source_values=[1.,2.];source_key=A.digest(source_values);loads=[]
    def identity(values,cid):
        metal=rows[cid]['metal'];ident={'cell_id':cid,'effective':base.effective(values,metal,None),'science':'fixed'}
        return A.digest(ident),ident,'unused',{'effective':ident['effective']}
    def vector(values):
        cells={}
        for cid in rows:
            key,ident,_,meta=identity(values,cid)
            cells[cid]={'status':'complete','cache_key':key,'cache_identity':ident,'request_provenance':{'cache_key':key,'cache_identity':ident,'render':meta,'vector':values[:]}}
        return {'status':'complete','parameters':values[:],'model':{},'response':{},'score':{'original_EF_checks_passed':False},'cells':cells}
    original=vector(source_values)
    def instance():
        ev=Cls.__new__(Cls);ev.lock=threading.RLock();ev.evaluation_lock=threading.Lock();ev.stop_event=threading.Event();ev.domain={};ev.m={'rows':rows};ev.cache={};ev.cell_locks={}
        ev.registry={source_key:{'parameters':source_values[:]}};ev.complete_vectors={};ev.loaded_saved_vectors={};ev.imported_cell_keys=set();ev.initial_vector=[0.,0.]
        ev.saved_metal_index={m:{A.digest(base.effective(source_values,m,None)):source_key} for m in ('La','Dy')}
        ev.cache_identity=identity;ev.checked_values=lambda v:v[:];ev._gate=lambda:None;ev.allowed={source_key:{'parameters':source_values[:],'role':'trial','outer_step':1}};ev.executor_identity={'mode':'mock'}
        ev._invoke=lambda *args:(_ for _ in ()).throw(AssertionError('Forbidden native helper'))
        return ev
    old_load=A.R.load
    try:
        def load(entry):loads.append(entry);return copy.deepcopy(original)
        A.R.load=load
        ev=instance()
        check('exact saved full vector requires zero programs',ev.required_programs(source_values)==0)
        check('one changed metal requires28 programs',ev.required_programs([1.001,2.])==28)
        check('two changed metals require56 programs',ev.required_programs([1.001,2.001])==56)
        check('nearby vector is not an exact reuse',ev.required_programs([1.+1e-12,2.])==28)
        check('no body loaded for metadata budget',len(loads)==0)
        answer=ev.evaluate(source_values)
        check('exact full vector imported once',len(loads)==1 and len(ev.imported_cell_keys)==46)
        check('failed physical gate retained',answer['score']['original_EF_checks_passed'] is False)
        check('saved source bytes untouched',answer==original)
        ev.evaluate(source_values)
        check('second exact request does not reread authority',len(loads)==1)
        rejects('explicit trial/probe role remains required',lambda:ev.evaluate([1.001,2.]))
        loads.clear();ev=instance();answer=ev._cell([1.,2.001],'La0')
        check('unchanged metal imported lazily',len(loads)==1 and answer['cache_key']==identity([1.,2.001],'La0')[0])
        check('new request provenance preserves current whole vector',answer['request_provenance']['vector']==[1.,2.001])
        check('stored source request unchanged',original['cells']['La0']['request_provenance']['vector']==source_values)
        ev._cell([1.,2.001],'La1');check('second cell reuses loaded native authority',len(loads)==1)
        malformed=copy.deepcopy(original);malformed['cells']['La0']['cache_identity']={'corrupted':True}
        rejects('changed saved scientific identity rejected',lambda:instance()._import_saved(source_key,malformed))
        malformed=copy.deepcopy(original);malformed['cells']['La0']['request_provenance']['render']={'corrupted':True}
        rejects('changed saved rendered identity rejected',lambda:instance()._import_saved(source_key,malformed))
        malformed=copy.deepcopy(original);malformed['parameters']=[1.001,2.]
        rejects('different saved parameter vector rejected',lambda:instance()._import_saved(source_key,malformed))
    finally:A.R.load=old_load
    with tempfile.TemporaryDirectory() as folder:
        p=Path(folder)/'source.json';p.write_text('{}');authority=A.pin(p);p.write_text('{"changed":true}')
        rejects('changed immutable artifact rejected',lambda:A.R.read(authority))
    out=Path(out);out.mkdir(parents=True,exist_ok=False)
    result={'status':'passed_physical_gate_reuse_controls','sources':{k:A.pin(HERE/k) for k in ('ADAPTER.py','DRIVER.py','REUSE.py','JOURNAL.py')},'implementation':A.pin(__file__),'checks':checks,'checks_passed':len(checks),'new_native_calls':0,'solver_calls':0,'scientific_parameter_evaluations':0,'actual_worker_concurrency_tested':False,'cost':{'wall_seconds':time.monotonic()-tick,'CPU_seconds':time.process_time()-cpu,'numerical_threads':1}}
    A.save(out/'RESULT.json',result);print(json.dumps({'result':A.pin(out/'RESULT.json'),'checks':len(checks),'cost':result['cost']}))


if __name__=='__main__':
    import sys
    run(sys.argv[1])
