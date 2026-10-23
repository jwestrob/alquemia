#!/usr/bin/env python3
"""Manufactured controls only: no real saved vectors, Jacobians or solver calls."""
import importlib.util,json,math,time,sys,tempfile,hashlib
from pathlib import Path
import numpy as np
from scipy import sparse
P=Path(__file__).resolve().parent;spec=importlib.util.spec_from_file_location('tested_runtime',P/'RUN.py');R=importlib.util.module_from_spec(spec);spec.loader.exec_module(R)
start=time.monotonic();cpu=time.process_time();checks=[]
def check(name,ok):
    if not ok:raise AssertionError(name)
    checks.append(name)
def rejects(name,f):
    try:f()
    except (ValueError,FileExistsError,RuntimeError):checks.append(name);return
    raise AssertionError(name)
def forbidden(*a,**k):raise AssertionError('Real solver construction forbidden in controls')
R.clarabel.DefaultSolver=forbidden
# Load source helpers/pins only, never vectors, LINEARIZATION or CONIC_DATA.
inv=R.read(R.pin(P.parent/'INVENTORY.json'));basepin=R.pin(Path(inv['sources']['settings']['path']).with_name('COMPARE.py'))
h,B,edits=R.helpers(inv,basepin);check('dimension_only_scalar_helper',set(edits)=={12,13} and all(edits.values()))
groups=[]
for i,name in enumerate(R.NAMES):
 r=[1.] if name.startswith('E3/') else [1.,0.,0.];J=np.zeros((len(r),14));J[0,0]=-2.;J[0,12]=1e-15
 groups.append(dict(name=name,r=r,J=J,scale=1.,reference=[0.]*len(r),model=r))
A=np.zeros((6,15));A[0,14]=-1;A[2,0]=1;A[3,14]=-1;A[4,14]=1;A[5,1]=1
b=np.array([0.,.25,.025,0.,.5,.1]);layout=[dict(name='block/P_energy',kind='SOC',start=0,stop=2),dict(name='trust/upper/0',kind='nonnegative',start=2,stop=3),dict(name='epigraph/lower',kind='nonnegative',start=3,stop=4),dict(name='epigraph/upper',kind='nonnegative',start=4,stop=5),dict(name='signed_W9_z/upper',kind='nonnegative',start=5,stop=6)]
m=R.adapt(sparse.csc_matrix(A),b,layout,groups,.0625)
check('only_two_epigraph_rows_removed',m['removed']==['epigraph/lower','epigraph/upper'])
check('six_cones_and_two_bounds_added',m['A'].shape==(24,15) and len(m['cones'])==8)
x=np.zeros(15);x[0]=.001;x[1]=.002;x[14]=.01
slack=m['b']-m['A']@x;old=b-A@np.r_[x[:14],.25]
check('inherited_rows_exact_fixed_s_substitution',all(np.array_equal(slack[z['start']:z['stop']],old[z['original_rows']]) for z in m['layout'] if z['original_rows'] is not None))
check('retained_t_column_zero',all(m['A'][z['start']:z['stop'],14].nnz==0 for z in m['layout'] if z['original_rows'] is not None))
for g,z in zip(groups,[z for z in m['layout'] if z['name'].startswith('trimer/')]):
 check('new_cone_sign_'+g['name'],np.array_equal(slack[z['start']:z['stop']],np.r_[1-x[-1],np.asarray(g['r'])+g['J']@x[:14]]))
check('finite_radial_column_retained',m['A'][:,12].nnz==6)
check('zero_conic_margins_nonnegative',min(B['cone_margins'](m,m['b']).values())>=0)
bad=[dict(g) for g in groups];bad[0]['scale']=0;rejects('zero_norm_scale_rejected',lambda:R.adapt(sparse.csc_matrix(A),b,layout,bad,.0625))
rejects('missing_group_rejected',lambda:R.adapt(sparse.csc_matrix(A),b,layout,groups[:-1],.0625))
rejects('wrong_row_map_rejected',lambda:R.adapt(sparse.csc_matrix(A),b,layout[:-1],groups,.0625))
p=dict(raw_center_q=.0625,q_current=.06250001,center=[.5]*14,width=[1.]*14,lower=[0.]*14,upper=[1.]*14,bounds_lower=[-.025]*14+[0.],bounds_upper=[.025]*14+[.06250001])
p['blocks']={k:dict(r=[.25],J=np.zeros((1,14)),denominator=1.) for k in ('old_EF','pair_mu','C_mu','chi','P_energy','B_energy')}
p['energy']=dict(names=[str(i) for i in range(10)],r=[0.]*10,J=np.zeros((10,14)))
p['forces']={str(i):dict(r=[0.]*27,J=np.zeros((27,14)),weights=[1/27]*27) for i in range(15)}
p['signed_W9_z']=dict(r_Debye=0.,J_Debye=[0.]*14,absolute_ceiling_Debye=.1)
z=R.certificate(p,groups,np.zeros(15),h['certify_fraction']);check('scalar_zero_has_all25_and_exact_boundaries',z['exact_sign_feasible'] and len(z['original25'])==25 and z['q']==p['raw_center_q'])
fraction=np.zeros(15);fraction[-1]=.1;check('unimproved_residual_rejects_positive_t',not R.certificate(p,groups,fraction,h['certify_fraction'])['exact_sign_feasible'])
qbad=dict(p,raw_center_q=np.nextafter(.0625,0.));check('one_bit_q_violation_not_padded',R.certificate(qbad,groups,np.zeros(15),h['certify_fraction'])['margins']['literal_q_cap']<0)
qhead=dict(p,blocks={k:dict(v) for k,v in p['blocks'].items()});qhead['blocks']['P_energy']['r']=[.1];check('max_q_preserves_nonlimiting_headroom',R.certificate(qhead,groups,np.zeros(15),h['certify_fraction'])['exact_sign_feasible'])
rejects('nonfinite_point_rejected',lambda:R.certificate(p,groups,[float('nan')]*15,h['certify_fraction']))
# A constructor failure consumes one durable attempt; a second call cannot reach it.
with tempfile.TemporaryDirectory() as tmp:
 t=Path(tmp);cfg=dict(attempt_path=str(t/'LOCK.json'),configuration_fingerprint='manufactured',implementation=R.pin(P/'RUN.py'));result=dict(solver_calls=0,solver_constructions=0)
 def fail(*a):raise RuntimeError('manufactured constructor failure before solve')
 rejects('constructor_failure_preserved',lambda:R.once(cfg,m,None,t,result,B,factory=fail))
 check('durable_attempt_precedes_constructor',(t/'LOCK.json').exists() and result['solver_constructions']==1 and result['solver_calls']==0)
 rejects('same_lock_forbids_second_attempt',lambda:R.once(cfg,m,None,t,dict(solver_calls=0,solver_constructions=0),B,factory=forbidden))
 # Release binding is acyclic and preserves mode/output constraints.
 cfg.update(output_directory=str(t/'prepare'),source_review={'path':'manufactured','sha256':'review'},runtime_declaration={'path':'manufactured','sha256':'declaration'})
 release={k:cfg[k] for k in ('implementation','output_directory','source_review','runtime_declaration')};release.update(status='released_trimer_local_span_prepare',mode='prepare',configuration_fingerprint='manufactured',maximum_solver_calls=0)
 R.save(t/'RELEASE.json',release);cfg['root_release']=R.pin(t/'RELEASE.json');R.release_gate(cfg,'prepare','manufactured');checks.append('acyclic_exact_release_binding')
 rejects('release_mode_mismatch_rejected',lambda:R.release_gate(cfg,'solve','manufactured'))
 changed=dict(cfg,output_directory=str(t/'other'));rejects('release_output_mismatch_rejected',lambda:R.release_gate(changed,'prepare','manufactured'))
 # Cycles helper receives manufactured seven-term values, no actual cells.
 terms={'fixture':{'model_terms':{'a':1,'b':-1},'reference_electronic_E3_kcal_mol':0.,'model_E3_kcal_mol':0.,'reference_mu3_Debye':[0.,0.,0.],'model_mu3_Debye':[0.,0.,0.],'reference_print_bound_Debye_per_component':0.,'model_print_bound_Debye_per_component':0.}}
 vals={'a':{'energy':3.,'mu_Debye':[1.,-2.,3.]},'b':{'energy':1.,'mu_Debye':[.5,-1.,2.]}}
 c=h['cycles'](terms,vals)['fixture'];check('exact_saved_cycle_helper_signed_arithmetic',c['candidate_E3_kcal_mol']==2. and c['candidate_mu3_Debye']==[.5,-1.,1.])
out=Path(sys.argv[1]);out.mkdir(parents=True,exist_ok=False);R.save(out/'RESULT.json',dict(status='passed_manufactured_runtime_controls',implementation=R.pin(P/'RUN.py'),test_source=R.pin(__file__),checks=checks,count=len(checks),actual_solver_calls=0,actual_stencil_reads=0,actual_Jacobians_constructed=0,native_calls=0,cost=dict(wall_seconds=time.monotonic()-start,CPU_seconds=time.process_time()-cpu,threads=1)))
print(len(checks),R.pin(out/'RESULT.json'))
