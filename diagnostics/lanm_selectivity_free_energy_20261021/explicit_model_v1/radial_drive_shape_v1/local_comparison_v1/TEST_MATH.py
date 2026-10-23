#!/usr/bin/env python3
"""Manufactured dimension/cone/certificate correspondence; no solver/native."""
import argparse,copy,json,math,subprocess,time
from pathlib import Path
import numpy as np
import clarabel
import RUN as R
import MATH as M
import SOLVE as S

def main(cfg,out):
    tick,cpu=time.monotonic(),time.process_time();checks=[]
    def check(name,ok):R.require(ok,name);checks.append({'name':name,'passed':True})
    def forbidden(*a,**kw):raise AssertionError('No solver/native call permitted')
    oldsolver,oldsub=clarabel.DefaultSolver,subprocess.run;clarabel.DefaultSolver=forbidden;subprocess.run=forbidden
    try:
        sources=R.read(cfg['inputs'])['sources'];old=R.module(cfg['core'],'radial_math_control_core');oldattr=R.module(sources['attribution_math'],'radial_math_control_attribution');core=M.core14(old);attr=M.attribution14(oldattr)
        def columns(n):return np.arange(n*12,dtype=float).reshape(n,12)*1e-6
        blocks={}
        for k,(a,b) in old.SLICES.items():
            r=np.full(b-a,.01);blocks[k]={'r':r,'J':columns(len(r)),'denominator':float(r@r)*2}
        for k in ('P_energy','B_energy'):blocks[k]={'r':np.array([1.,2.]),'J':columns(2),'denominator':100.}
        forces={k+'/'+geo:{'r':np.full(3+6*n,.1),'J':columns(3+6*n),'weights':np.full(3+6*n,1/(3+6*n))} for k in ('La','Dy','paired') for geo,n in [('aquo8_central',8),('aquo9_central',9),('water_minus',1),('water_central',1),('water_plus',1)]}
        p={'center':np.zeros(12),'lower':-np.ones(12),'upper':np.ones(12),'width':np.full(12,2.),'z':np.full(12,.5),'q_current':.1,'raw_center_q':.1,'blocks':blocks,'energy':{'r':np.full(10,.1),'J':columns(10),'names':['g'+str(i) for i in range(10)]},'forces':forces,'bounds_lower':np.r_[np.full(12,-.025),0.],'bounds_upper':np.r_[np.full(12,.025),.1],'mapping_S':1.}
        p14=copy.deepcopy(p)
        for k in ('center','lower','upper','width','z'):p14[k]=np.r_[p[k],np.full(2,p[k][0])]
        p14['bounds_lower']=np.r_[np.full(14,-.025),0.];p14['bounds_upper']=np.r_[np.full(14,.025),.1]
        for block in [*p14['blocks'].values(),p14['energy'],*p14['forces'].values()]:block['J']=np.column_stack([block['J'],np.zeros((len(block['r']),2))])
        d=np.arange(12,dtype=float)*.0001;d14=np.r_[d,0.,0.]
        c12=old.certify_fraction(p,np.r_[d,.09],.5,False);c14=core.certify_fraction(p14,np.r_[d14,.09],.5,False)
        check('all851_scalar_margins_equal_at_zero_new_columns',c12['certificate']['residual_constraints']==c14['certificate']['residual_constraints'])
        check('all25_physical_metrics_equal',c12['certificate']['original25']==c14['certificate']['original25'])
        check('first12_physical_parameters_equal',c12['parameters']==c14['parameters'][:12] and c14['parameters'][12:]==[0.,0.])
        check('851_and25_denominators_preserved',len(c14['certificate']['residual_constraints'])==851 and len(c14['certificate']['original25'])==25)
        p14['signed_W9_z']={'r_Debye':.2,'J_Debye':np.r_[np.full(12,.01),.03,-.04],'absolute_ceiling_Debye':.2,'reference_Dy_minus_La_z_Debye':-.1}
        channels={'chi_states_a0_cubed':(['La_dimer','La_trimer','Dy_dimer','Dy_trimer'],np.zeros(4)),'C_dipoles_Debye':(['aquo8_central_La','aquo8_central_Dy','aquo9_central_La','aquo9_central_Dy'],np.zeros((4,3))),'pair_dipoles_Debye':(['p'],np.zeros((1,3)))}
        p14['electrical_unweighted']={k:{'names':names,'r':r,'J':np.zeros((*r.shape,14))} for k,(names,r) in channels.items()};p14['electrical_targets']={'chi_a0_cubed':{k:0. for k in channels['chi_states_a0_cubed'][0]},'C_Debye':{k:[0.,0.,0.] for k in channels['C_dipoles_Debye'][0]}}
        models={case:S.assemble(core,attr,p14,case) for case in S.CASES}
        check('same_rhs_and_layout',np.array_equal(models['original12']['b'],models['radial14']['b']) and models['original12']['layout']==models['radial14']['layout'])
        check('original12_is_literal_column_restriction',np.array_equal(models['original12']['A'].toarray(),models['radial14']['A'][:,list(range(12))+[14]].toarray()))
        for case in S.CASES:
            model=models[case];x=np.r_[d,math.sqrt(.09)] if case=='original12' else np.r_[d14,math.sqrt(.09)];full=S.expand(x,case);slack=model['b']-model['A']@x
            check(case+'_20_SOC',sum(a['kind']=='SOC' for a in model['layout'])==20)
            check(case+'_890_affine_sides',sum(a['kind']=='nonnegative' for a in model['layout'])==890)
            cert=core.certify_fraction(p14,np.r_[full[:14],full[-1]**2],1.,False)['certificate']['residual_constraints'];maxerr=0.
            for item in model['layout']:
                name=item['name'];s=slack[item['start']:item['stop']]
                if name in cert:maxerr=max(maxerr,abs(float(s[0]**2-s[1:]@s[1:] if item['kind']=='SOC' else s[0])-cert[name]))
            check(case+'_cone_scalar_correspondence',maxerr<1e-12)
            require_names={a['name'] for a in model['layout']};check(case+'_no_EF_ceiling_same3electrical','block/old_EF' not in require_names and all('block/'+k in require_names for k in ('pair_mu','C_mu','chi')))
            tail=slack[-2:];res=.2+math.fsum(float(a)*float(b) for a,b in zip(p14['signed_W9_z']['J_Debye'],full[:14]));check(case+'_signed_guard_sides',np.max(abs(tail-np.array([.2-res,.2+res])))<1e-14)
        report=S.certificate(core,attr,p14,np.r_[d14,.09],.5,False)
        check('signed_violation_cannot_be_admitted',not report['signed_W9_safeguard_pass'] and not report['retained_exact_sign_feasible'])
        check('auxiliary_q_only_physical_parameters_unchanged',report['physical_fraction_changed'] is False and report['parameters']==c14['parameters'])
        check('raw_solver_status_not_upgraded',report['raw_certificate']['certificate']['solver_success'] is False and not report['model_admitted'] and not report['native_trial_authorized'])
        check('source_only_zero_solver_native',True)
    finally:clarabel.DefaultSolver=oldsolver;subprocess.run=oldsub
    r={'status':'passed_radial_comparison_manufactured_controls','implementation':R.pin(__file__),'sources':{n:R.pin(Path(__file__).parent/(n+'.py')) for n in ('RUN','ADAPTER','MATH','SOLVE')},'checks':checks,'count':len(checks),'actual_jacobians':0,'native_calls':0,'solver_constructions':0,'solver_calls':0,'cost':{'wall_seconds':time.monotonic()-tick,'CPU_seconds':time.process_time()-cpu,'threads':1}};R.save(out/'RESULT.json',r);print(json.dumps({'status':r['status'],'checks':len(checks)}))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--config',required=True);p.add_argument('--out',required=True);a=p.parse_args();out=Path(a.out).resolve();out.mkdir(parents=True,exist_ok=False)
    try:main(json.loads(Path(a.config).read_text()),out)
    except Exception as e:R.save(out/'FAILURE.json',{'error':repr(e),'native_calls':0,'solver_calls':0});raise
