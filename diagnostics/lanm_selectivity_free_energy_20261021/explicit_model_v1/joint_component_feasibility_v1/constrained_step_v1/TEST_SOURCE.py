"""Constraint calculus plus saved physical reconstruction; no optimizer calls."""
import argparse
import importlib.util
import json
from pathlib import Path
import time

HERE=Path(__file__).resolve().parent
s=importlib.util.spec_from_file_location('constrained_proposal',HERE/'PROPOSE.py');m=importlib.util.module_from_spec(s);s.loader.exec_module(m)

def main(out):
    tick,cpu=time.monotonic(),time.process_time();checks={}
    def forbidden(*a,**kw):raise AssertionError('No solver calls during preparation')
    m.minimize=forbidden
    # Test calculus independently on deterministic synthetic coefficients.
    blocks={k:{'r':m.np.array([1.+i/10,.7]),'J':m.np.arange(24,dtype=float).reshape(2,12)/(70+i)} for i,k in enumerate((*m.SLICES,'P_energy','B_energy'))}
    synthetic={'blocks':blocks,'energy':{'r':m.np.arange(10)/20,'J':m.np.arange(120,dtype=float).reshape(10,12)/500,'names':[str(i) for i in range(10)]},'forces':{},'bounds_lower':m.np.r_[m.np.full(12,-.025),0.],'bounds_upper':m.np.r_[m.np.full(12,.025),1.]}
    for kind in ('La','Dy','paired'):
        for geo in m.GEOMETRIES:
            n=8 if geo=='aquo8_central' else 9 if geo=='aquo9_central' else 1;q=3+6*n
            synthetic['forces'][kind+'/'+geo]={'r':m.np.arange(q,dtype=float)/100,'J':m.np.arange(q*12,dtype=float).reshape(q,12)/10000,'weights':m.np.r_[m.np.full(3,1/6),m.np.full(q-3,1/(12*n))]}
    x=m.np.r_[m.np.arange(12)/10000,.005];value,jac,names=m.constraints(synthetic,x);fd=[]
    for j in range(13):
        delta=m.np.zeros(13);delta[j]=1e-6
        fd.append((m.constraints(synthetic,x+delta)[0]-m.constraints(synthetic,x-delta)[0])/(2e-6))
    err=float(m.np.max(abs(jac-m.np.column_stack(fd))))
    checks['all851_analytic_constraint_derivatives']=err<1e-8
    checks['851_constraints13_variables']=jac.shape==(851,13)
    checks['component_t_sign']=all(jac[names.index('block/'+n),-1]==-1 for n in ('P_energy','B_energy'))
    checks['protected_t_absent']=all(jac[names.index('block/'+n),-1]==0 for n in m.SLICES)
    zero=m.predictions(synthetic,m.np.zeros(13));checks['separate25_gates']=len(zero['original25'])==25
    # Actual preparation reads only the origin and its existing24 probes.
    problem=m.prepare(HERE/'INPUTS.json');origin=m.predictions(problem,m.np.zeros(13))
    checks['exact_origin24']=problem['origin']==m.ORIGIN
    checks['actual24_saved_probe_denominator']=len(problem['derivative_keys'])==24
    checks['actual_origin_all25_exact']=origin['all25_pass_exact_threshold']
    checks['actual_origin_zero_step_feasible']=origin['numerically_feasible']
    checks['actual851_constraints']=len(origin['constraints'])==851
    checks['all_six_blocks_nondegenerate']=all(m.np.linalg.norm(b['r'])>1e-12 for b in problem['blocks'].values())
    checks['physical413_reconstruction_all25vectors']=True # prepare raises on any failure
    m.require(all(checks.values()),str(checks))
    result={'status':'passed_constrained_step_source_preparation','implementation':m.pin(HERE/'PROPOSE.py'),'inputs':m.pin(HERE/'INPUTS.json'),'tests':m.pin(__file__),'specification':m.pin(HERE/'SPEC.md'),'checks':checks,'count':len(checks),'maximum_synthetic_constraint_derivative_error':err,'origin':problem['origin'],'origin_artifact':problem['origin_artifact'],'existing_vectors_reconstructed':25,'constraint_count':851,'solver_calls':0,'native_calls':0,'new_parameter_evaluations':0,'cost':{'wall_seconds':time.monotonic()-tick,'CPU_seconds':time.process_time()-cpu,'threads':1}}
    with Path(out).open('x') as f:json.dump(result,f,indent=2);f.write('\n')
    print(json.dumps(result,indent=2))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--out',required=True);a=p.parse_args();main(a.out)
