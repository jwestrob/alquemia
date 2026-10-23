"""One artificial read-only reproduction and repaired symbolic constructor only."""
import ast,hashlib,importlib.util,json,os
from pathlib import Path
import sys,time,traceback,warnings
import numpy as np
from scipy.sparse import csc_array
from sksparse import cholmod

HERE=Path(__file__).resolve().parent
out=Path(sys.argv[1]);out.mkdir(parents=True,exist_ok=False)
start,cpu=time.monotonic(),time.process_time();checks=[]
def pin(p):
    p=Path(p).resolve();return dict(path=str(p),sha256=hashlib.sha256(p.read_bytes()).hexdigest())
def check(name,value):
    checks.append(dict(name=name,passed=bool(value)))
    if not value:raise ValueError(name)
def payload(a):return hashlib.sha256(memoryview(a).cast('B')).hexdigest()
result=dict(status='started',checks=checks,implementation=pin(__file__),runner=pin(HERE/'RUN.py'),declaration=pin(HERE/'DECLARATION.md'),
            constructor_attempts=0,completed_symbolic_constructors=0,numeric_factorizations=0,solve_calls=0,
            full_matrix_reads=0,native_calls=0,error=None)
try:
    check('one thread',all(os.environ.get(k)=='1' for k in ('OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','BLIS_NUM_THREADS')))
    old=(HERE.parent/'RUN.py').read_text();new=(HERE/'RUN.py').read_text()
    before='k=csc_array((np.array(data,copy=True),idx,ptr),shape=(N,N),copy=False);backend.validate_lower(k)'
    after="k=csc_array((np.array(data,copy=True),np.array(idx,copy=True,order='C'),\n                     np.array(ptr,copy=True,order='C')),shape=(N,N),copy=False);backend.validate_lower(k)"
    check('exact sole runner delta',old.count(before)==1 and old.replace(before,after)==new)
    cfg=json.loads((HERE/'CONFIG_UNRELEASED.json').read_text());inputs=json.loads((HERE/'INPUTS.json').read_text());parent=json.loads(Path(inputs['parent_configuration']['path']).read_text())
    check('exact sole config delta',dict(parent,implementation=pin(HERE/'RUN.py'))==cfg)
    bpin=cfg['backend'];check('backend source pin',pin(bpin['path'])==bpin)
    s=importlib.util.spec_from_file_location('qualified_backend',bpin['path']);backend=importlib.util.module_from_spec(s);s.loader.exec_module(backend)
    values=dict(data=np.array([4.,1.,.5,3.,.25,2.],dtype=np.float64),indices=np.array([0,1,2,1,2,2],dtype=np.int32),indptr=np.array([0,3,5,6],dtype=np.int32))
    for name,a in values.items():np.save(out/(name+'.npy'),a,allow_pickle=False)
    filepins={name:pin(out/(name+'.npy')) for name in values}
    arrays={name:np.load(p['path'],mmap_mode='r',allow_pickle=False) for name,p in filepins.items()}
    data,idx,ptr=(arrays[k] for k in ('data','indices','indptr'))
    bad=csc_array((np.array(data,copy=True),idx,ptr),shape=(3,3),copy=False);bad.data[bad.indptr[:-1]]-=.125
    check('serialized read-only index views',not bad.indices.flags.writeable and not bad.indptr.flags.writeable and bad.data.flags.writeable)
    result['constructor_attempts']+=1
    try:
        obj=cholmod.CholeskyFactor(bad,**backend.SETTINGS)
        result['completed_symbolic_constructors']+=1
        result['old_exception']=None
    except Exception as exc:result['old_exception']=dict(type=type(exc).__name__,message=str(exc))
    check('old error reproduced exactly',result['old_exception']==dict(type='ValueError',message='buffer source array is read-only'))
    good=csc_array((np.array(data,copy=True),np.array(idx,copy=True,order='C'),np.array(ptr,copy=True,order='C')),shape=(3,3),copy=False)
    good.data[good.indptr[:-1]]-=.125
    check('all repaired buffers writable contiguous',all(a.flags.writeable and a.flags.c_contiguous for a in (good.data,good.indices,good.indptr)))
    check('all payloads and dtypes unchanged',all(payload(getattr(bad,k))==payload(getattr(good,k)) and getattr(bad,k).dtype==getattr(good,k).dtype for k in values))
    result['constructor_attempts']+=1
    with warnings.catch_warnings(record=True) as seen:
        warnings.simplefilter('always');obj=cholmod.CholeskyFactor(good,**backend.SETTINGS)
    result['completed_symbolic_constructors']+=1;result['state']=backend.metadata(obj)
    result['warnings']=[dict(type=type(x.message).__name__,message=str(x.message)) for x in seen]
    check('repaired symbolic only succeeds',obj.N==3 and obj.is_super and not obj.is_numeric and not seen)
    check('serialized files unchanged',all(pin(p['path'])==p for p in filepins.values()))
    result.update(status='passed_writable_sparse_index_constructor_recovery',files=filepins,payload_hashes={k:payload(getattr(good,k)) for k in values},settings=backend.SETTINGS,artificial_shift=.125)
except Exception as exc:
    result.update(status='failed_writable_sparse_index_constructor_recovery',error=dict(type=type(exc).__name__,message=str(exc),traceback=traceback.format_exc()))
result['cost']=dict(wall_seconds=time.monotonic()-start,CPU_seconds=time.process_time()-cpu,numerical_threads=1)
(out/'RESULT.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
print(json.dumps(dict(status=result['status'],result=pin(out/'RESULT.json'),checks=len(checks),error=result['error'],cost=result['cost'])))
raise SystemExit(0 if result['error'] is None else 1)
