"""Node import/runtime metadata only: no box, operator, factor or molecular call."""
import ctypes
import hashlib
import json
import os
from pathlib import Path
import resource
import sys
import time

def main(authority,out):
    tick,cpu=time.monotonic(),time.process_time();audit=json.loads(Path(authority).read_text())
    if Path(sys.prefix).resolve()!=Path(audit['prefix']).resolve():raise ValueError('Python prefix')
    for item in audit['runtime_files']:
        with Path(item['path']).open('rb') as stream:sha=hashlib.file_digest(stream,'sha256').hexdigest()
        if sha!=item['sha256']:raise ValueError('Runtime file pin')
    import numpy as np
    import scipy
    import scipy.sparse
    versions=dict(python='.'.join(map(str,sys.version_info[:3])),numpy=np.__version__,scipy=scipy.__version__)
    if any(versions[k]!=audit['versions'][k] for k in versions):raise ValueError('Runtime versions')
    prefix=Path(sys.prefix)
    blas=ctypes.CDLL(str(prefix/'lib/libopenblasp-r0.3.34.so'));omp=ctypes.CDLL(str(prefix/'lib/libgomp.so.1'))
    blas.openblas_get_num_threads.restype=ctypes.c_int;omp.omp_get_max_threads.restype=ctypes.c_int
    actual=dict(openblas=blas.openblas_get_num_threads(),openmp=omp.omp_get_max_threads())
    if actual!={'openblas':1,'openmp':1}:raise ValueError('Actual numerical threads')
    if os.environ.get('LD_PRELOAD') or os.environ.get('LD_LIBRARY_PATH'):raise ValueError('Loader override')
    result=dict(status='passed_stage1_node_runtime_preflight',versions=versions,actual_threads=actual,
                prefix=str(prefix),address_space_limit=list(resource.getrlimit(resource.RLIMIT_AS)),
                native_calls=0,operator_applications=0,factorizations=0,
                cost=dict(wall_seconds=time.monotonic()-tick,CPU_seconds=time.process_time()-cpu,
                          peak_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
    with Path(out).open('x') as stream:json.dump(result,stream,indent=2);stream.write('\n')

if __name__=='__main__':main(*sys.argv[1:])
