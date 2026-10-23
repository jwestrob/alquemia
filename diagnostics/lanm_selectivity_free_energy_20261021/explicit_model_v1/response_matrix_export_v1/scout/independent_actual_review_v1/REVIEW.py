#!/usr/bin/env python3
"""One saved actual review of scout18564; no native invocation."""
import hashlib
import importlib.util
import json
from decimal import Decimal, localcontext
from pathlib import Path
import re
import time
import numpy as np

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[5]
BASE=ROOT/'workspaces/lanm_selectivity_free_energy_20261021/explicit_model_v1/response_matrix_export_v1/scout'
OUT=BASE/'independent_actual_review_v1'

def pin(p):
    p=Path(p);return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}

def verify(p):
    assert pin(p['path'])==p,p
    return Path(p['path'])

def sites(p):
    t=[s.split() for s in verify(p).read_text().splitlines() if s.strip()]
    assert int(t[0][0])==4 and len(t)==5
    assert [(int(s[0]),s[1]) for s in t[1:]]==[(1,'La3'),(2,'O'),(3,'H'),(4,'H')]
    a=np.array([[float(v) for v in s[2:5]] for s in t[1:]])
    assert np.isfinite(a).all()
    return a

def raw(p):
    t=[s.split() for s in verify(p).read_text().splitlines() if s.strip()]
    assert len(t)==36 and t[0]==['TINKER_RESPONSE_MATRIX_V1'] and t[-1]==['END_TINKER_RESPONSE_MATRIX_V1']
    assert t[1][2:]==['4','4','12'] and t[6]==['RESTORED','T']
    ss=[v for v in t if v[0]=='SITE'];vv=[v for v in t if v[0]=='VECTOR'];aa=[v for v in t if v[0]=='A']
    assert len(ss)==4 and len(vv)==len(aa)==12
    assert [int(v[1]) for v in vv]==[int(v[1]) for v in aa]==list(range(1,13))
    assert [(int(v[1]),int(v[2]),int(v[3]),v[4]) for v in ss]==[(1,1,500,'La3'),(2,2,36,'O'),(3,3,37,'H'),(4,4,37,'H')]
    A=np.array([[float(x) for x in v[2:]] for v in aa]);V=np.array([[float(x) for x in v[2:]] for v in vv])
    assert A.shape==(12,12) and V.shape==(12,4) and np.isfinite(A).all() and np.isfinite(V).all()
    return t,ss,vv,A,V

def main():
    tick,cpu=time.monotonic(),time.process_time()
    result_path=BASE/'run_v1/execution/RESULT.json';r=json.loads(result_path.read_text());cfg=r['configuration']
    assert cfg==json.loads((BASE/'run_v1/CONFIG_RUN.json').read_text())
    m=json.loads(verify(cfg['manifest']).read_text())
    spec=importlib.util.spec_from_file_location('frozen_scout',verify(cfg['implementation']));mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
    admission=mod.admit_result(cfg,m,r)
    assert admission['status']=='admitted_response_matrix_scout' and admission['new_native_calls']==0
    assert len(admission['qualification']['checks'])==37 and all(c['passed'] for c in admission['qualification']['checks'].values())
    for k in ('source_manifest','analyze','build_receipt','root_release','declaration'):verify(cfg[k])
    build=json.loads(verify(cfg['build_receipt']).read_text());assert build['job_id']=='18563' and build['manifest']==cfg['source_manifest'] and build['executables']['analyze']==cfg['analyze']
    for cmd in build['commands']:assert cmd['returncode']==0;verify(cmd['log'])
    verify(build['cache'])
    off=r['rows']['off'];on=r['rows']['on'];old=m['old_zero']
    assert verify(off['moments']['analyze_output']).read_bytes()==verify(old['analyze_output']).read_bytes()
    u0=sites(old['induced_output']);uoff=sites(off['moments']['induced_output']);uon=sites(on['moments']['induced_output'])
    assert np.array_equal(u0,uoff) and np.array_equal(uoff,uon)
    off_text=verify(off['moments']['analyze_output']).read_text();on_text=verify(on['moments']['analyze_output']).read_text()
    stripped=''.join(s for s in on_text.splitlines(True) if 'Response Matrix Written To' not in s)
    assert stripped==off_text
    q=[float(v) for v in re.findall(r'Induced Dipoles\s*:\s*Iterations\s+\d+\s+RMS Residual\s+(\S+)',on_text)]
    assert q==[0.,0.]
    exports=[raw(p) for p in on['exports']]
    assert exports[0][0][1][1]=='1' and exports[1][0][1][1]=='2'
    assert np.array_equal(exports[0][3],exports[1][3]) and np.array_equal(exports[0][4],exports[1][4])
    units=m['units'];K=units['electric'];D=units['debye_per_eA'];B=units['bohr_A'];H=units['hartree_kcal_mol'];n=4;dim=12
    A=exports[0][3];ed,ep,u,v=exports[0][4].T
    alpha=np.repeat([float(s[5]) for s in exports[0][1]],3)
    assert np.array_equal(alpha,np.repeat([1.6597,.837,.496,.496],3))
    assert np.array_equal(ed,ep) and np.array_equal(u,v)
    ar=128*dim*np.finfo(float).eps*max(1.,np.max(abs(A)),np.max(abs(ed)),np.max(abs(u)))
    sym=float(np.max(abs(A-A.T)));assert sym<=ar
    eig=np.linalg.eigvalsh(A);normeig=np.linalg.eigvalsh(np.sqrt(alpha)[:,None]*A*np.sqrt(alpha)[None,:])
    assert eig[0]>0 and normeig[0]>1e-6*normeig[-1]
    L=np.linalg.cholesky(A)
    def solve(b):return np.linalg.solve(L.T,np.linalg.solve(L,b))
    dense=solve(ed);res=float(np.linalg.norm(A@u-ed));amp=float(np.linalg.norm(np.eye(dim)-A@np.diag(alpha),2));an=np.sqrt(dim)*ar
    cgbound=amp*np.sqrt(n)*.5e-10/D+an
    assert res<=cgbound and np.linalg.norm(dense-u)<=(res+np.linalg.norm(A@dense-ed)+2*an)/eig[0]
    token=re.findall(r'^\s*Polarization\s+(\S+)\s+\d+\s*$',on_text,re.M);assert len(token)==1
    with localcontext() as ctx:
        ctx.prec=50
        contraction=-Decimal(str(K))/2*sum(Decimal(vv[4])*Decimal(vv[3]) for vv in exports[0][2])
        energy_error=float(abs(contraction-Decimal(token[0])))
    assert energy_error<=5e-9+128*dim*np.finfo(float).eps*max(1.,abs(float(contraction)))
    site_error=float(np.max(abs(D*u.reshape(4,3)-uon)));assert site_error<=.5e-6+D*ar
    z=np.tile([0.,0.,1.],4);coupling=H/(K*B);chi=solve(z).reshape(4,3).sum(axis=0)*coupling/B
    field_reports={};oldsites={};fieldbounds={}
    for label,target in m['old_fields'].items():
        moment=target['moments'];f=target['field_au'];oldsites[label]=sites(moment['induced_output'])
        txt=verify(moment['analyze_output']).read_text();qs=[float(x) for x in re.findall(r'Induced Dipoles\s*:\s*Iterations\s+\d+\s+RMS Residual\s+(\S+)',txt)]
        assert len(qs)==2 and all(np.isfinite(qs)) and min(qs)>=0 and max(qs)<=1e-8
        pred=solve(ed+coupling*f*z)
        rb=amp*np.sqrt(n)*(max(qs)+.5e-10)/D+an
        errbound=(rb+np.linalg.norm(A@pred-ed-coupling*f*z)+an)/eig[0]
        error=float(np.max(abs(D*pred.reshape(4,3)-oldsites[label])));bound=.5e-6+D*errbound
        assert error<=bound
        fieldbounds[label]=errbound;field_reports[label]={'site_max_error_Debye':error,'bound_Debye':bound}
    slopes={}
    for tag,h in (('h',.001),('half',.0005)):
        obs=(oldsites['plus_'+tag]-oldsites['minus_'+tag]).sum(axis=0)/(2*h*D*B)
        bound=n*.5e-6/(h*D*B)+np.sqrt(n)*(fieldbounds['plus_'+tag]+fieldbounds['minus_'+tag])/(2*h*B)+128*dim*np.finfo(float).eps*max(1.,np.max(abs(chi)))
        error=float(np.max(abs(chi-obs)));assert error<=bound
        slopes[tag]={'raw_column_a0_cubed':obs.tolist(),'max_error_a0_cubed':error,'bound_a0_cubed':bound}
    receipt={'status':'passed_baseline_response_matrix_scout_review','result':pin(result_path),'implementation':cfg['implementation'],'reviewer_implementation':pin(__file__),'configuration':pin(BASE/'run_v1/CONFIG_RUN.json'),'manifest':cfg['manifest'],'source_manifest':cfg['source_manifest'],'build_receipt':cfg['build_receipt'],'executable':cfg['analyze'],'root_release':cfg['root_release'],
      'complete_native_slots':2,'raw_export_files':on['exports'],'replayed_numerical_checks':37,'all_saved_admission_checks_passed':True,
      'independent_arithmetic':{'baseline_analyze_bytes_exact':True,'on_after_removing_two_announcements_bytes_exact':True,'old_off_on_induced_sites_exact':True,'two_export_matrices_RHS_solutions_exact':True,'dual_RHS_solutions_exact':True,'matrix_symmetry_max':sym,'A_eigenvalues':eig.tolist(),'normalized_A_eigenvalues':normeig.tolist(),'condition_A':float(eig[-1]/eig[0]),'condition_normalized_A':float(normeig[-1]/normeig[0]),'native_residual_norm':res,'post_peek_residual_bound':cgbound,'native_vs_Cholesky_solution_norm':float(np.linalg.norm(dense-u)),'polarization_energy_Decimal_kcal_mol':str(contraction),'polarization_print_error_kcal_mol':energy_error,'native_site_rounding_error_Debye':site_error,'predicted_chi_column_a0_cubed':chi.tolist(),'old_field_sites':field_reports,'old_field_columns':slopes},
      'actual_native_cost':r['cost'],'review_cost':{'wall_seconds':time.monotonic()-tick,'CPU_seconds':time.process_time()-cpu,'local_threads':1},'new_native_calls':0,'new_reference_calls':0,'tests_repeated':0,'error':None,
      'limitations':['Qualifies original baseline operator on this scout only; remaining53 require their own raw qualification.','No changed response physics, candidate forces, fit, reference promotion, GPU or periodic qualification.']}
    OUT.mkdir(exist_ok=False);(OUT/'REVIEW.json').write_text(json.dumps(receipt,indent=2,allow_nan=False)+'\n')
    print(json.dumps({'review':pin(OUT/'REVIEW.json'),'independent_arithmetic':receipt['independent_arithmetic'],'review_cost':receipt['review_cost']},indent=2))

if __name__=='__main__':main()
