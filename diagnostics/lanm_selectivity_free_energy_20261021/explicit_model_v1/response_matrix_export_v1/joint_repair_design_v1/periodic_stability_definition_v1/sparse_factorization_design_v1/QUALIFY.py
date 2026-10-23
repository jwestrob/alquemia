"""Finite manufactured CHOLMOD checks. Never reads a saved molecular matrix."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import sys
import time
import traceback
import warnings

import numpy as np
from scipy.sparse import csc_array
from sksparse import cholmod

EPS = np.finfo(float).eps
TINY = np.finfo(float).tiny
THRESHOLD = 1e-6
SETTINGS = dict(lower=True, order='amd', sym_kind='sym', supernodal_mode='simplicial')


def pin(p):
    p = Path(p).resolve()
    return dict(path=str(p), sha256=hashlib.sha256(p.read_bytes()).hexdigest())


def verify(p):
    assert pin(p['path']) == p
    return json.loads(Path(p['path']).read_text())


def guard(a):
    a = np.asarray(a)
    if a.ndim != 2 or a.shape[0] != a.shape[1]:
        raise ValueError('nonsquare')
    if np.iscomplexobj(a):
        raise ValueError('complex unsupported')
    a = np.asarray(a, dtype=float)
    if not np.all(np.isfinite(a)):
        raise ValueError('nonfinite')
    if np.max(np.abs(a-a.T)) > 1e-10:
        raise ValueError('asymmetry above inherited bound')
    return a


def norm(a):
    return float(np.linalg.norm(a, np.inf))


def dense(a):
    e, v = np.linalg.eigh(a)
    residual = float(np.max(np.abs(a@v-v*e)))
    orth = float(np.max(np.abs(v.T@v-np.eye(len(a)))))
    return dict(eigenvalues=e.tolist(), residual=residual, orthogonality=orth,
                passed=residual <= 1e-10*max(1., float(e[-1])) and orth <= 1e-10)


def metadata(f):
    r = dict(minor=None, status=None, unavailable_reason='not exposed by public scikit-sparse 0.5.0 API')
    for k in ('N', 'is_numeric', 'is_ll', 'is_super', 'is_lower', 'order'):
        try:
            x = getattr(f, k)
            r[k] = x.item() if isinstance(x, np.generic) else x
        except Exception as e:
            r[k] = repr(e)
    return r


def exception(e):
    return dict(type=type(e).__name__, message=str(e), attributes={k: repr(v) for k, v in vars(e).items()})


def factor_checks(f, k):
    n = len(k)
    p = f.get_perm()
    l = f.get_factor(kind='LL', lower=True).toarray()
    b = np.arange(1, n+1, dtype=float)
    x = np.asarray(f.solve(b, system='A'))
    ferr = norm(k[np.ix_(p,p)]-l@l.T)/max(1.,norm(k))
    serr = norm(k@x-b)/max(TINY,norm(k)*norm(x)+norm(b))
    tol = 100*n*EPS
    positive = bool(np.all(np.isfinite(l)) and np.all(np.diag(l)>0))
    return dict(permutation=p.tolist(), L=l.tolist(), rhs=b.tolist(), solution=x.tolist(),
                factor_backward_error=ferr, solve_backward_error=serr, arithmetic_bound=tol,
                finite_positive_factor=positive,
                passed=positive and np.all(np.isfinite(x)).item() and ferr <= tol and serr <= tol)


def manufactured():
    h = np.array([[1,1,1,1],[1,-1,1,-1],[1,1,-1,-1],[1,-1,-1,1]],float)/2
    b = np.array([[-.2,.1],[.1,1.]])
    return [
        ('diagonal6_spd',np.diag([.5,1,2,3,4,6]),6.),
        ('coupled3_negative_gershgorin',np.full((3,3),.9)+np.eye(3)*.1,3.),
        ('hadamard4_spd',h@np.diag([.2,.7,1.1,2.])@h.T,2.1),
        ('threshold_above',np.diag([1.01e-6,.1,.5,1.]),1.),
        ('threshold_below',np.diag([.99e-6,.1,.5,1.]),1.),
        ('threshold_equality',np.diag([1e-6,1.]),1.),
        ('threshold_plus_machine_epsilon',np.diag([1e-6+EPS,1.]),1.),
        ('singular2',np.diag([0.,1.]),1.),
        ('indefinite2',np.diag([-.2,1.]),1.),
        ('B_indefinite2',b,1.2),
        ('B_plus_PSD_R_positive2',b+np.diag([.5,0.]),1.2)]


def run(config, outdir):
    start, cpu = time.perf_counter(), time.process_time()
    cfg = json.loads(Path(config).read_text())
    assert pin(__file__) == cfg['implementation']
    declaration = verify(cfg['declaration'])
    env = verify(cfg['environment'])
    assert env['status'] == 'passed_exact_isolated_environment_audit'
    assert Path(sys.prefix).resolve() == Path(env['prefix'])
    for p in env['runtime_files']:
        assert pin(p['path']) == p
    assert all(os.environ.get(k) == '1' for k in env['threads']['environment'])
    assert declaration['maximum_dimension'] == 6 and declaration['threshold']['epsilon'] == THRESHOLD
    outdir = Path(outdir); outdir.mkdir(parents=True, exist_ok=False)
    result = dict(status='started', implementation=pin(__file__), configuration=pin(config),
                  declaration=cfg['declaration'], environment=cfg['environment'], settings=SETTINGS,
                  manufactured=[], rank3=[], input_errors=[], factorizations=0, rank_mutations=0,
                  native_calls=0, saved_molecular_matrix_reads=0, error=None)
    checks = []
    def check(name, passed):
        checks.append(dict(name=name, passed=bool(passed)))
    def persist():
        result['checks'] = checks
        result['cost'] = dict(wall_seconds=time.perf_counter()-start,
                              CPU_seconds=time.process_time()-cpu, numerical_threads=1)
        p = outdir / 'RESULT.json'
        tmp = outdir / 'RESULT.tmp'
        tmp.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
        tmp.replace(p)
    def factor(k):
        assert len(k) <= declaration['maximum_dimension']
        result['factorizations'] += 1
        assert result['factorizations']+result['rank_mutations'] <= declaration['maximum_numerical_factorization_attempts']
        f = cholmod.CholeskyFactor(csc_array(k), **SETTINGS)
        error = None
        with warnings.catch_warnings(record=True) as ww:
            warnings.simplefilter('always')
            try:
                f.factorize(csc_array(k), ldl=False, beta=0.)
            except Exception as e:
                error = exception(e)
        return f,error,[dict(category=type(w.message).__name__,message=str(w.message)) for w in ww]
    persist()
    try:
        for name,a,u in manufactured():
            a = guard(a); k = a-THRESHOLD*u*np.eye(len(a))
            dr, kr = dense(a),dense(k)
            margin = kr['eigenvalues'][0]
            resolution = 100*len(a)*EPS*max(1.,norm(a))
            expected = 'numerically_unresolved' if abs(margin)<=resolution else ('positive' if margin>0 else 'not_positive')
            f,error,ww = factor(k)
            rec = dict(name=name,A=a.tolist(),U=u,K=k.tolist(),dense=dr,shifted_dense=kr,
                       gershgorin_lower=float(np.min(np.diag(a)-(np.sum(np.abs(a),axis=1)-np.abs(np.diag(a))))),
                       roundoff_resolution=resolution,decision=expected,backend_metadata=metadata(f),exception=error,warnings=ww)
            if error is None:
                rec['factor_checks'] = factor_checks(f,k)
            result['manufactured'].append(rec)
            check(name+'/dense_arithmetic',dr['passed'] and kr['passed'])
            check(name+'/supplied_upper_bound',dr['eigenvalues'][-1] <= u)
            if expected == 'positive':
                check(name+'/positive_factor',error is None and rec['factor_checks']['passed'])
            elif expected == 'not_positive':
                check(name+'/nonpositive_exception',error is not None and error['type']=='CholmodNotPositiveDefiniteError')
            else:
                check(name+'/unresolved_not_admitted',rec['decision']=='numerically_unresolved')
            persist()
        a0=np.diag(declaration['rank3_cases']['base_diagonal']); u=6.
        k0=a0-THRESHOLD*u*np.eye(6); e=np.eye(6)[:,:3]
        for beta in declaration['rank3_cases']['beta_values']:
            f,error,ww=factor(k0)
            if error is not None:
                raise RuntimeError('declared rank3 base factor failed')
            xe=np.asarray(f.solve(e,system='A'))
            small=np.eye(3)+beta*(e.T@xe)
            k=k0+beta*(e@e.T)
            dr,sr=dense(k),dense(small)
            expected=dr['eigenvalues'][0]>0
            rec=dict(beta=beta,K0=k0.tolist(),E=e.tolist(),K_updated=k.tolist(),small=small.tolist(),
                     small_dense=sr,updated_dense=dr,base_factor_checks=factor_checks(f,k0),
                     three_solve_backward_errors=[norm(k0@xe[:,i]-e[:,i])/max(TINY,norm(k0)*norm(xe[:,i])+norm(e[:,i])) for i in range(3)],
                     expected_positive=bool(expected),exception=None,warnings=ww,failed_object_reused=False)
            f.change_factor(kind='LDL')
            with warnings.catch_warnings(record=True) as ww2:
                warnings.simplefilter('always')
                try:
                    if beta:
                        result['rank_mutations']+=1
                        c=csc_array(np.sqrt(abs(beta))*e)
                        (f.update if beta>0 else f.downdate)(c)
                    l,d=f.get_factor(kind='LDL',lower=True)
                    l,d=l.toarray(),d.toarray()
                    p=f.get_perm(); ld_error=norm(k[np.ix_(p,p)]-l@d@l.T)/max(1.,norm(k))
                    positive=bool(np.all(np.isfinite(d)) and np.all(np.diag(d)>0))
                    rec.update(LDL_diagonal=np.diag(d).tolist(),LDL_backward_error=ld_error,
                               backend_positive_diagonal=positive)
                    if positive:
                        rec['updated_factor_checks']=factor_checks(f,k)
                    else:
                        rec['rejection']='nonpositive LDL diagonal; factor discarded without solve or LL conversion'
                except Exception as err:
                    rec['exception']=exception(err)
                    rec['rejection']='backend exception; disposable factor discarded'
                rec['warnings'] += [dict(category=type(w.message).__name__,message=str(w.message)) for w in ww2]
            rec['backend_metadata']=metadata(f)
            result['rank3'].append(rec)
            check(str(beta)+'/schur_equivalence',bool(sr['eigenvalues'][0]>0)==bool(expected))
            check(str(beta)+'/dense_arithmetic',dr['passed'] and sr['passed'])
            check(str(beta)+'/base_and_three_solves',rec['base_factor_checks']['passed'] and max(rec['three_solve_backward_errors'])<=600*EPS)
            if expected:
                check(str(beta)+'/updated_factor',rec['exception'] is None and rec.get('updated_factor_checks',{}).get('passed',False))
            else:
                check(str(beta)+'/failed_update_not_promoted',rec['exception'] is not None or rec.get('backend_positive_diagonal') is False)
            if 'LDL_backward_error' in rec:
                check(str(beta)+'/explicit_updated_matrix',rec['LDL_backward_error']<=600*EPS)
            del f
            persist()
        bad=[('nonsquare',np.ones((2,3))),('asymmetric',np.array([[1,.01],[0,1]])),
             ('NaN',np.array([[np.nan,0],[0,1]])),('infinity',np.diag([np.inf,1])),('complex',np.eye(2,dtype=complex))]
        for name,a in bad:
            try:
                guard(a); err=None
            except ValueError as er:
                err=exception(er)
            result['input_errors'].append(dict(name=name,exception=err,numeric_attempts=0))
            check(name+'/guard_rejection',err is not None)
        check('complete_inventory',len(result['manufactured'])==11 and len(result['rank3'])==4 and len(result['input_errors'])==5)
        check('finite_cap',result['factorizations']+result['rank_mutations']<=64)
        result['status']='passed_manufactured_backend_qualification' if all(c['passed'] for c in checks) else 'failed_manufactured_backend_qualification'
    except Exception as er:
        result['status']='failed_manufactured_backend_qualification'
        result['error']=dict(**exception(er),traceback=traceback.format_exc())
    finally:
        persist()
    print(json.dumps(dict(status=result['status'],result=pin(outdir/'RESULT.json'),cost=result['cost'])))
    return 0 if result['status']=='passed_manufactured_backend_qualification' else 1


if __name__=='__main__':
    p=argparse.ArgumentParser(); p.add_argument('--config',required=True);p.add_argument('--out',required=True)
    a=p.parse_args();raise SystemExit(run(a.config,a.out))
