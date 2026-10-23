"""Exactly the declared tiny symbolic graphs and disposable limit test."""
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import resource
import subprocess
import sys
import time
import traceback

import numpy as np
from scipy.sparse import csc_array


def pin(p):
    p=Path(p).resolve();return dict(path=str(p),sha256=hashlib.sha256(p.read_bytes()).hexdigest())


def graph_counts(n, edges, perm):
    position={old:new for new,old in enumerate(perm)};adj=[set() for _ in range(n)]
    for a,b in edges:
        a,b=position[a],position[b];adj[a].add(b);adj[b].add(a)
    counts=[]
    for j in range(n):
        higher=sorted(k for k in adj[j] if k>j);counts.append(1+len(higher))
        for a in higher:
            for b in higher:
                if a!=b:adj[a].add(b)
    return counts


def main():
    config=Path(sys.argv[1]);out=Path(sys.argv[2]);cfg=json.loads(config.read_text())
    assert pin(__file__)==cfg['controls'];assert pin(cfg['implementation']['path'])==cfg['implementation']
    assert pin(cfg['declaration']['path'])==cfg['declaration']
    env=json.loads(Path(cfg['environment']['path']).read_text());assert pin(cfg['environment']['path'])==cfg['environment']
    assert Path(sys.prefix).resolve()==Path(env['prefix'])
    assert all(os.environ.get(k)=='1' for k in env['threads']['environment'])
    for p in env['runtime_files']:assert pin(p['path'])==p
    spec=importlib.util.spec_from_file_location('symbolic_helper',cfg['implementation']['path'])
    helper=importlib.util.module_from_spec(spec);spec.loader.exec_module(helper)
    out.mkdir(parents=True,exist_ok=False);tick,cpu=time.perf_counter(),time.process_time()
    result=dict(status='started',configuration=pin(config),implementation=cfg['implementation'],
                declaration=cfg['declaration'],environment=cfg['environment'],rows=[],checks=[],
                symbolic_analyses=0,numeric_factorizations=0,native_calls=0,error=None)
    def ck(label,passed):result['checks'].append(dict(name=label,passed=bool(passed)))
    def persist():
        result['cost']=dict(wall_seconds=time.perf_counter()-tick,CPU_seconds=time.process_time()-cpu,numerical_threads=1)
        (out/'RESULT.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    persist()
    try:
        n=6;path=[(i,i+1) for i in range(5)];cycle=path+[(0,5)]
        graphs=[('diagonal',[]),('path',path),('star',[(0,i) for i in range(1,6)]),
                ('cycle',cycle),('dense',[(i,j) for i in range(6) for j in range(i)]),('indefinite_cycle',cycle)]
        for name,edges in graphs:
            a=np.eye(n)*6
            for i,j in edges:a[max(i,j),min(i,j)]=-.25
            if name=='indefinite_cycle':a[np.diag_indices(n)]=-1.
            a=csc_array(a);saved=(a.data.copy(),a.indices.copy(),a.indptr.copy())
            f,report=helper.analyze_only(a,max_numeric_estimate_bytes=1<<20,reserve_bytes=4096)
            result['symbolic_analyses']+=1
            expected=graph_counts(n,edges,report['permutation'])
            ck(name+'/exact_graph_column_counts',expected==report['column_counts'])
            ck(name+'/symbolic_only',f is not None and not f.is_numeric and report['rcond_before_numeric']==-1)
            ck(name+'/input_unchanged',all(np.array_equal(x,y) for x,y in zip(saved,(a.data,a.indices,a.indptr))))
            ck(name+'/ample_plan',report['numeric_memory_allowed'])
            denied=helper.memory_plan(n,report['column_counts'],report['index_bytes'],report['input_CSC_bytes'],max_numeric_estimate_bytes=0,reserve_bytes=4096)
            ck(name+'/zero_budget_denies',not denied['numeric_memory_allowed'])
            report['name']=name;report['independent_graph_column_counts']=expected;result['rows'].append(report)
            del f;persist()
        ck('same_pattern_indefinite_not_spd_proof',result['rows'][3]['column_counts']==result['rows'][5]['column_counts'] and not result['rows'][5]['is_numeric'])
        limits_before=resource.getrlimit(resource.RLIMIT_AS)
        script='''import importlib.util,json,sys
spec=importlib.util.spec_from_file_location("s",sys.argv[1]);s=importlib.util.module_from_spec(spec);spec.loader.exec_module(s)
before=s.virtual_bytes();limit=before+8*1024*1024;r=s.set_worker_address_space_limit(limit)
small=bytearray(1024*1024);r["small_allocation_passed"]=len(small)==1024*1024
try:
 large=bytearray(16*1024*1024);r["large_refused"]=False
except MemoryError:
 r["large_refused"]=True
print(json.dumps(r))
'''
        t=time.perf_counter();u0=resource.getrusage(resource.RUSAGE_CHILDREN)
        child=subprocess.run([sys.executable,'-c',script,cfg['implementation']['path']],capture_output=True,text=True,check=False)
        u1=resource.getrusage(resource.RUSAGE_CHILDREN)
        (out/'LIMIT_CHILD.stdout').write_text(child.stdout);(out/'LIMIT_CHILD.stderr').write_text(child.stderr)
        result['limit_child']=dict(exit_code=child.returncode,stdout=pin(out/'LIMIT_CHILD.stdout'),stderr=pin(out/'LIMIT_CHILD.stderr'),
                wall_seconds=time.perf_counter()-t,CPU_seconds=(u1.ru_utime+u1.ru_stime)-(u0.ru_utime+u0.ru_stime),max_RSS_KiB=u1.ru_maxrss)
        payload=json.loads(child.stdout) if child.returncode==0 else {}
        result['limit_child']['readback']=payload
        ck('child_limit_bound',child.returncode==0 and payload.get('small_allocation_passed') and payload.get('large_refused'))
        ck('parent_limit_unchanged',resource.getrlimit(resource.RLIMIT_AS)==limits_before)
        ck('exact_scope',result['symbolic_analyses']==6 and result['numeric_factorizations']==0)
        result['status']='passed_tiny_symbolic_memory_controls' if all(x['passed'] for x in result['checks']) else 'failed_tiny_symbolic_memory_controls'
    except Exception as e:
        result['status']='failed_tiny_symbolic_memory_controls';result['error']=dict(type=type(e).__name__,message=str(e),traceback=traceback.format_exc())
    persist();print(json.dumps(dict(status=result['status'],result=pin(out/'RESULT.json'),cost=result['cost'])))
    return 0 if result['status']=='passed_tiny_symbolic_memory_controls' else 1


if __name__=='__main__':raise SystemExit(main())
