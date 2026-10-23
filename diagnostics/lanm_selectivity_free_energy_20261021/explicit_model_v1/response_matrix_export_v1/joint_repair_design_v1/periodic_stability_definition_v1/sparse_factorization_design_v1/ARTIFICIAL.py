#!/usr/bin/env python3
"""Tiny manufactured factor/certificate checks, never a molecular matrix."""
from fractions import Fraction as F
import hashlib
import json
import math
import os
from pathlib import Path
import time
import numpy as np
from scipy.sparse import csc_matrix
from scipy.sparse.linalg import splu


def pin(path):
    path=Path(path).resolve()
    return {'path':str(path),'sha256':hashlib.sha256(path.read_bytes()).hexdigest()}


def rational_certificate(K,Q):
    """Exact rational verification of the *stored binary* tiny triangular factor."""
    n=len(K); a=[[F(float(v)) for v in row] for row in K]; q=[[F(float(v)) for v in row] for row in Q]
    assert all(q[i][i]>0 and all(q[i][j]==0 for j in range(i+1,n)) for i in range(n))
    residual=[[a[i][j]-sum(q[i][k]*q[j][k] for k in range(n)) for j in range(n)] for i in range(n)]
    error=max(sum(abs(x) for x in row) for row in residual)
    # M = diag(Q)-abs(strict_lower(Q)); M^-1 dominates abs(Q^-1).
    y=[F(0)]*n;z=[F(0)]*n
    for i in range(n):y[i]=(1+sum(abs(q[i][j])*y[j] for j in range(i)))/q[i][i]
    for i in range(n-1,-1,-1):z[i]=(1+sum(abs(q[j][i])*z[j] for j in range(i+1,n)))/q[i][i]
    inv2=max(y)*max(z); margin=1/inv2-error
    return {'positive_definite_certified_for_stored_tiny_K':margin>0,'exact_margin':str(margin),
            'residual_infinity_upper_exact':str(error),'inverse_two_norm_squared_upper_exact':str(inv2)}


def main(out):
    assert all(os.environ.get(k)=='1' for k in ('OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','BLIS_NUM_THREADS'))
    tick,cpu=time.monotonic(),time.process_time(); checks=[];results={}
    def check(name,value):
        if not value:raise AssertionError(name)
        checks.append(name)
    K=np.full((3,3),.6);np.fill_diagonal(K,1.)
    check('negative_Gershgorin_does_not_mean_indefinite', min(np.diag(K)-(np.sum(abs(K),axis=1)-abs(np.diag(K))))<0 and np.linalg.eigvalsh(K)[0]>0)
    lu=splu(csc_matrix(K),permc_spec='MMD_AT_PLUS_A',diag_pivot_thresh=0.,options={'Equil':False,'SymmetricMode':True})
    check('artificial_symmetric_permutations',np.array_equal(lu.perm_r,lu.perm_c))
    Pr=np.eye(3)[np.argsort(lu.perm_r),:]
    check('documented_SuperLU_permutation',np.max(abs(Pr@K@Pr.T-lu.L.toarray()@lu.U.toarray()))<1e-14)
    D=lu.U.diagonal();check('positive_candidate_diagonal',np.all(D>0))
    Q=lu.L.toarray()*np.sqrt(D)[None,:]
    results['exact_tiny_factor_verification']=rational_certificate(Pr@K@Pr.T,Q)
    check('independent_rational_factor_certificate',results['exact_tiny_factor_verification']['positive_definite_certified_for_stored_tiny_K'])
    wrong=Q.copy();wrong[2,0]+=2
    check('bad_factor_rejected_independently',not rational_certificate(Pr@K@Pr.T,wrong)['positive_definite_certified_for_stored_tiny_K'])
    B=np.diag([-1.,2.]);Z=np.array([[math.sqrt(2.)],[0.]])
    check('B_indefinite_full_operator_positive',np.linalg.eigvalsh(B)[0]<0 and np.linalg.eigvalsh(B+Z@Z.T)[0]>0)
    S=np.eye(1)+Z.T@np.linalg.solve(B,Z)
    check('signed_small_Schur_inertia_recovery',S[0,0]<0)
    badR=2*np.ones((2,2));badB=np.diag([-1.,1.]);A=badB+badR
    check('negative_subspace_diagonal_only_is_insufficient',A[0,0]>0 and np.linalg.eigvalsh(A)[0]<0)
    base=np.diag([.5,2.,3.,4.]);E=np.eye(4)[:,:3]
    for beta,expected in ((-.25,True),(-.75,False),(.5,True)):
        S=np.eye(3)+beta*E.T@np.linalg.solve(base,E)
        check('rank3_update_'+str(beta),(np.linalg.eigvalsh(S)[0]>0)==expected and (np.linalg.eigvalsh(base+beta*E@E.T)[0]>0)==expected)
    A=np.diag([1e-5,100.]);tau=1e-6*100.
    check('positive_definite_alone_fails_relative_threshold',np.linalg.eigvalsh(A)[0]>0 and np.linalg.eigvalsh(A-tau*np.eye(2))[0]<0)
    result={'status':'passed_tiny_artificial_sparse_design_controls','source':pin(__file__),'checks':checks,'count':len(checks),'results':results,
            'largest_dimension':4,'sparse_factorizations':1,'box_matrices_read_or_built':0,'native_molecular_calls':0,'dependency_installs':0,
            'rigorous_scope':'Exact Fraction arithmetic certifies only the explicitly stored tiny binary matrix/factor above; all other tiny eigenvalue checks are ordinary doubles.',
            'wall_seconds':time.monotonic()-tick,'CPU_seconds':time.process_time()-cpu,'numerical_threads':1}
    out=Path(out);out.mkdir(parents=True,exist_ok=False);(out/'RESULT.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'result':pin(out/'RESULT.json'),'checks':len(checks),'wall_seconds':result['wall_seconds'],'CPU_seconds':result['CPU_seconds']}))


if __name__=='__main__':
    import sys
    main(sys.argv[1])
