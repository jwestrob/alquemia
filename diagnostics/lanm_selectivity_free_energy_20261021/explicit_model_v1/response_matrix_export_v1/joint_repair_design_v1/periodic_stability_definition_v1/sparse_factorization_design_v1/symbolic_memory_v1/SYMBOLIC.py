"""Public symbolic-only fill report. Does not call any numeric factor method."""
import os
import resource
from pathlib import Path

import numpy as np
from scipy.sparse import csc_array, issparse
from sksparse import cholmod

SETTINGS=dict(lower=True,order='amd',sym_kind='sym',supernodal_mode='simplicial')


def virtual_bytes():
    return int(Path('/proc/self/statm').read_text().split()[0])*os.sysconf('SC_PAGE_SIZE')


def set_worker_address_space_limit(limit_bytes):
    """Call only inside the dedicated disposable worker, before assembly."""
    if type(limit_bytes) is not int or limit_bytes<=0:
        raise ValueError('positive integer address-space ceiling required')
    previous=resource.getrlimit(resource.RLIMIT_AS)
    if previous[0]!=resource.RLIM_INFINITY and limit_bytes>previous[0]:
        raise ValueError('requested ceiling would relax inherited soft limit')
    if previous[1]!=resource.RLIM_INFINITY and limit_bytes>previous[1]:
        raise ValueError('requested ceiling exceeds inherited hard limit')
    current=virtual_bytes()
    if current>=limit_bytes:
        raise MemoryError('existing virtual mappings already exceed requested ceiling')
    resource.setrlimit(resource.RLIMIT_AS,(limit_bytes,previous[1]))
    actual=resource.getrlimit(resource.RLIMIT_AS)
    if actual!=(limit_bytes,previous[1]):
        raise RuntimeError('address-space limit did not bind')
    return dict(previous=list(previous),current=list(actual),mapped_before_bytes=current,
                ceiling_kind='total virtual address space; not RSS')


def memory_plan(n, counts, index_bytes, input_bytes, *,
                max_numeric_estimate_bytes, reserve_bytes):
    if any(type(v) is not int or v<0 for v in (max_numeric_estimate_bytes,reserve_bytes)):
        raise ValueError('explicit nonnegative integer planning limits required')
    counts=[int(x) for x in counts]
    if len(counts)!=n or index_bytes not in (4,8) or any(c<1 or c>n-j for j,c in enumerate(counts)):
        raise ValueError('invalid symbolic counts/index width')
    nnz=sum(counts)
    # Perm/ColCount/p/nz/next/prev: 6*n+5 integer entries. IPerm optional.
    packed=(8+index_bytes)*nnz+index_bytes*(6*n+5)
    iperm=index_bytes*n
    estimate=2*(packed+iperm)+2*input_bytes+128*n+reserve_bytes
    return dict(symbolic_L_nnz=nnz,work_proxy_sum_squared_column_counts=sum(c*c for c in counts),
                work_proxy_is_not_measured_common_fl=True,packed_factor_array_bytes=packed,
                optional_inverse_permutation_bytes=iperm,input_CSC_bytes=input_bytes,
                forecast=dict(factor_payload_copies=2,input_additional_copies=2,
                              workspace_bytes_per_coordinate=128,caller_reserve_bytes=reserve_bytes,
                              estimated_additional_bytes=estimate,ceiling_bytes=max_numeric_estimate_bytes,
                              guaranteed_peak_bound=False),
                numeric_memory_allowed=estimate<=max_numeric_estimate_bytes)


def analyze_only(lower, *, max_numeric_estimate_bytes, reserve_bytes):
    """Consumes sparse_assembly_v1's lower CSC; returns (factor-or-None,report).

    No numeric arrays are requested from CHOLMOD. An over-budget plan returns
    None so the caller cannot accidentally use this returned object to factor.
    A caller must also install the independent worker address-space ceiling.
    """
    if not issparse(lower) or lower.format!='csc' or lower.dtype!=np.dtype('float64'):
        raise ValueError('real float64 lower CSC required')
    a=csc_array(lower,copy=False);n=a.shape[0]
    if a.ndim!=2 or a.shape[1]!=n or n==0 or not a.has_canonical_format:
        raise ValueError('square canonical lower CSC required')
    for j in range(n):
        lo,hi=a.indptr[j:j+2];rows=a.indices[lo:hi];values=a.data[lo:hi]
        if np.any(rows<j) or np.sum(rows==j)!=1 or not np.all(np.isfinite(values)):
            raise ValueError('lower triangle with exactly one explicit finite diagonal required')
    mapped_before=virtual_bytes()
    f=cholmod.CholeskyFactor(a,**SETTINGS)
    counts=np.asarray(f.colcount).copy();perm=np.asarray(f.perm).copy()
    if f.is_numeric or f.is_super or f.N!=n or sorted(perm.tolist())!=list(range(n)):
        raise RuntimeError('unexpected symbolic state/permutation')
    width=np.dtype(f.itype).itemsize
    p=memory_plan(n,counts,width,int(a.data.nbytes+a.indices.nbytes+a.indptr.nbytes),
                  max_numeric_estimate_bytes=max_numeric_estimate_bytes,reserve_bytes=reserve_bytes)
    report=dict(dimension=n,settings=SETTINGS,column_counts=counts.tolist(),permutation=perm.tolist(),
                index_dtype=str(np.dtype(f.itype)),index_bytes=width,input_index_bytes=a.indices.dtype.itemsize,
                is_numeric=bool(f.is_numeric),is_super=bool(f.is_super),rcond_before_numeric=float(f.rcond),
                mapped_virtual_bytes_before_analysis=mapped_before,
                mapped_virtual_bytes_after_analysis=virtual_bytes(),
                raw_common_lnz=None,raw_common_fl=None,raw_common_memory_usage=None,
                raw_common_memory_inuse=None,raw_minor=None,raw_status=None,
                common_fields_unavailable='not exposed by this public Python API',**p)
    try:
        report['public_nnz']=int(f.nnz)
    except Exception as e:
        report['public_nnz']=None;report['public_nnz_error']=repr(e)
    report['status']='symbolic_complete_numeric_stage_not_executed' if p['numeric_memory_allowed'] else 'symbolic_complete_numeric_memory_rejected'
    return (f if p['numeric_memory_allowed'] else None),report
