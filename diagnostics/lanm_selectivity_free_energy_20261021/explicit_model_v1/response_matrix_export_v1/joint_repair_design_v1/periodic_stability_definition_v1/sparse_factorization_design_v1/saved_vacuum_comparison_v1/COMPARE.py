"""One released, saved-only six-matrix numerical comparison; no native calls."""
import argparse
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import sys
import time
import traceback
import warnings

import numpy as np
from scipy.sparse import csc_array
from sksparse import cholmod


def pin(p):
    p=Path(p).resolve()
    return dict(path=str(p),sha256=hashlib.sha256(p.read_bytes()).hexdigest())


def verify(p):
    assert pin(p['path'])==p, 'Changed pin '+p['path']
    return Path(p['path'])


def read(p):
    return json.loads(verify(p).read_text())


def load(p,name):
    spec=importlib.util.spec_from_file_location(name,verify(p))
    m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
    return m


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--config',type=Path,required=True)
    ap.add_argument('--release',type=Path,required=True);ap.add_argument('--out',type=Path,required=True)
    args=ap.parse_args();tick,cpu=time.perf_counter(),time.process_time()
    cfg=json.loads(args.config.read_text());release=json.loads(args.release.read_text())
    assert cfg['implementation']==pin(__file__)
    assert release['status']=='released_six_saved_vacuum_backend_comparison'
    assert release['configuration']==pin(args.config) and release['implementation']==pin(__file__)
    assert release['maximum_factorizations']==6 and release['numerical_threads']==1
    assert Path(release['output']).resolve()==args.out.resolve()
    assert release['backend_review']==cfg['backend_review']
    backend_review=read(cfg['backend_review'])
    ready=read(cfg['backend_ready']);env=read(ready['environment_audit'])
    assert backend_review['status']=='passed_manufactured_backend_source_and_result_review'
    assert backend_review['backend_ready']==cfg['backend_ready']
    assert backend_review['implementation']==ready['implementation'] and backend_review['result']==ready['result']
    assert ready['status']=='passed_isolated_cholmod_manufactured_qualification'
    assert read(ready['result'])['status']=='passed_manufactured_backend_qualification'
    assert Path(sys.prefix).resolve()==Path(env['prefix'])
    for p in env['runtime_files']:verify(p)
    assert all(os.environ.get(k)=='1' for k in env['threads']['environment'])
    sources=read(cfg['parent_inputs']);old=read(sources['qualification_result'])
    review=read(sources['qualification_review']);manifest=read(sources['qualification_manifest'])
    assert review['status']=='passed_local_response_native_qualification'
    assert review['result']==sources['qualification_result'] and old['numerical_checks_passed']
    assert old['manifest']==sources['qualification_manifest'] and old['native_program_invocations']==57
    verify(cfg['declaration']);verify(sources['source_manifest']);verify(sources['build_receipt'])
    helper=load(sources['matrix_parser'],'saved_vacuum_parser')
    backend=load(ready['implementation'],'qualified_tiny_backend')
    assert backend.THRESHOLD==1e-6
    snapshots=sources['snapshots']
    assert len(snapshots)==len({x['id'] for x in snapshots})==len({x['export']['path'] for x in snapshots})==6
    args.out.mkdir(parents=True,exist_ok=False)
    result=dict(status='started',implementation=pin(__file__),configuration=pin(args.config),
                release=pin(args.release),backend_ready=cfg['backend_ready'],backend_review=cfg['backend_review'],
                parent_inputs=cfg['parent_inputs'],declaration=cfg['declaration'],
                rows={x['id']:dict(status='untouched',source=x['export']) for x in snapshots},
                factorization_attempts=0,dense_eigensystems=0,solves=0,new_native_calls=0,
                periodic_qualification=False,error=None)
    def persist():
        result['cost']=dict(wall_seconds=time.perf_counter()-tick,CPU_seconds=time.process_time()-cpu,numerical_threads=1)
        tmp=args.out/'RESULT.tmp'
        tmp.write_text(json.dumps(result,indent=2,allow_nan=False,default=lambda x:x.item() if isinstance(x,np.generic) else (_ for _ in ()).throw(TypeError(type(x).__name__)))+'\n')
        tmp.replace(args.out/'RESULT.json')
    persist()
    try:
        for item in snapshots:
            cid=item['id'];rec=old['rows'][cid];row=manifest['rows'][rec['job']['cell']]
            assert rec['status']=='complete' and item['export']==rec['exports'][0]
            assert item['xyz']==rec['runtime_input']['xyz'] and Path(item['export']['path']).name=='response_matrix_001.dat'
            verify(item['xyz']);parsed=helper.parse_export(item['export'],row,manifest['units'])
            saved=rec['matrix_qualification']['exports'][0]
            assert parsed['sequence']==1 and all(parsed[k]==saved[k] for k in ('A','vectors','sites','checks','restored','peek'))
            a=np.asarray(parsed['A'],float);alpha=np.repeat([s['alpha_A3'] for s in parsed['sites']],3)
            assert a.shape==(len(alpha),len(alpha)) and 0<len(alpha)<=84 and np.all(np.isfinite(alpha)) and np.all(alpha>0)
            root=np.sqrt(alpha);b=backend.guard(root[:,None]*a*root[None,:])
            u=max(math.fsum(max(abs(float(b[i,j])),abs(float(b[j,i]))) for j in range(len(b))) for i in range(len(b)))
            assert np.isfinite(u) and u>0
            k=b-1e-6*u*np.eye(len(b))
            result['dense_eigensystems']+=1
            dr=backend.dense(b);actual=np.asarray(saved['spectrum']['normalized_A'])
            assert len(actual)==len(b)
            ev=np.asarray(dr['eigenvalues']);arith=helper.arithmetic(len(b),b,actual,ev)
            diff=float(np.max(abs(actual-ev)));sym=float(np.max(abs(b-b.T)))
            resolution=100*len(b)*np.finfo(float).eps*max(1.,backend.norm(b))
            margin=float(ev[0]-1e-6*u)
            rr=dict(status='started',source=item['export'],geometry=item['xyz'],profile=rec['job']['profile'],
                    dimension=len(b),sites=len(parsed['sites']),normalization='sqrt(alpha) A sqrt(alpha)',
                    symmetry_max=sym,dense=dr,saved_dense_min=float(actual[0]),saved_dense_max=float(actual[-1]),
                    dense_replay_max_difference=diff,inherited_arithmetic_bound=arith,
                    original_dense_criterion_passed=bool(ev[0]>1e-6*ev[-1]>0),
                    saved_dense_criterion_passed=bool(actual[0]>1e-6*actual[-1]>0),
                    sufficient_U=u,raw_matrix_infinity_norm=backend.norm(b),shift=1e-6*u,
                    shifted_dense_margin=margin,roundoff_resolution=resolution,
                    sufficient_decision='not_evaluated',exception=None,warnings=[],
                    numerical_integrity_passed=bool(dr['passed'] and diff<=arith and sym<=1e-10 and ev[-1]<=u+arith))
            result['rows'][cid]=rr;persist()
            assert rr['numerical_integrity_passed'],'Dense/normalization integrity failed'
            result['factorization_attempts']+=1;assert result['factorization_attempts']<=6
            f=None
            with warnings.catch_warnings(record=True) as ww:
                warnings.simplefilter('always')
                try:
                    f=cholmod.CholeskyFactor(csc_array(k),**backend.SETTINGS)
                    f.factorize(csc_array(k),ldl=False,beta=0.)
                    result['solves']+=1
                    rr['factor_checks']=backend.factor_checks(f,k)
                    rr['numerical_integrity_passed'] &= bool(rr['factor_checks']['passed'])
                    if abs(margin)<=resolution:
                        rr['sufficient_decision']='numerically_unresolved'
                    elif margin>0 and rr['factor_checks']['passed']:
                        rr['sufficient_decision']='passed_numerical_sufficient_test'
                    else:
                        rr['sufficient_decision']='numerical_integrity_failure'
                        rr['numerical_integrity_passed']=False
                except cholmod.CholmodNotPositiveDefiniteError as e:
                    rr['exception']=backend.exception(e)
                    rr['sufficient_decision']='numerically_unresolved' if abs(margin)<=resolution else 'not_certified_by_sufficient_test'
                    if margin>resolution:rr['numerical_integrity_passed']=False
                except Exception as e:
                    rr['exception']=backend.exception(e);rr['numerical_integrity_passed']=False
                    rr['sufficient_decision']='numerical_integrity_failure'
                rr['warnings']=[dict(category=type(w.message).__name__,message=str(w.message)) for w in ww]
            rr['backend_metadata']=backend.metadata(f) if f is not None else None
            rr['status']='complete';del f;persist()
        result['summary']=dict(total=6,completed=sum(x['status']=='complete' for x in result['rows'].values()),
                               sufficient_passes=sum(x.get('sufficient_decision')=='passed_numerical_sufficient_test' for x in result['rows'].values()),
                               dense_passes=sum(x.get('original_dense_criterion_passed',False) for x in result['rows'].values()),
                               numerical_integrity_passed=all(x.get('numerical_integrity_passed',False) for x in result['rows'].values()))
        result['status']='completed_saved_vacuum_backend_comparison' if result['summary']['numerical_integrity_passed'] else 'failed_saved_vacuum_backend_comparison'
    except Exception as e:
        result['status']='failed_saved_vacuum_backend_comparison'
        result['error']=dict(type=type(e).__name__,message=str(e),traceback=traceback.format_exc())
    persist()
    print(json.dumps(dict(status=result['status'],result=pin(args.out/'RESULT.json'),cost=result['cost'])))
    return 0 if result['status']=='completed_saved_vacuum_backend_comparison' else 1


if __name__=='__main__':
    raise SystemExit(main())
