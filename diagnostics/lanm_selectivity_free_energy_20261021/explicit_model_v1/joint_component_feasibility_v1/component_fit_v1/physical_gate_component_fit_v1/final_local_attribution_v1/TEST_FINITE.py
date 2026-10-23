"""Focused artificial wrapper/cone controls; no scientific solve or native call."""
import copy
import json
import math
from pathlib import Path
import socket
import sys
import threading
import time
from types import SimpleNamespace
import numpy as np
import RUN
import SOLVE

def main():
    tick,cpu=time.monotonic(),time.process_time();out=Path(sys.argv[1]);out.mkdir(parents=True,exist_ok=False)
    cfg=json.loads((RUN.HERE/'CONFIG_TEMPLATE.json').read_text())
    core=RUN.module(cfg['core'],'test_finite_core');old=RUN.module(cfg['original_adapter'],'test_finite_parent')
    calls=[]
    def forbidden(*a,**kw):calls.append('forbidden');raise AssertionError('No solver or native calls authorized')
    SOLVE.clarabel.DefaultSolver=forbidden
    checks=[]
    def check(name,ok):
        checks.append({'name':name,'passed':bool(ok)});assert ok,name
    def rejects(name,fn):
        try:fn()
        except ValueError:check(name,True)
        else:check(name,False)
    try:
        blocks={}
        for name,n in [('old_EF',413),('pair_mu',102),('C_mu',12),('chi',4),('P_energy',2),('B_energy',2)]:
            r=np.full(n,.1);J=np.zeros((n,12));J[:,0]=.01
            blocks[name]={'r':r,'J':J,'denominator':1. if name.endswith('energy') else math.fsum(float(x*x) for x in r)}
        forces={}
        for kind in ('La','Dy','paired'):
            for geo in core.GEOMETRIES:
                n=8 if geo=='aquo8_central' else 9 if geo=='aquo9_central' else 1
                forces[kind+'/'+geo]={'r':np.zeros(3+6*n),'J':np.zeros((3+6*n,12)),
                    'weights':np.r_[np.full(3,1/6),np.full(6*n,1/(12*n))]}
        ch=['La_dimer','La_trimer','Dy_dimer','Dy_trimer'];cn=[f'aquo{n}_central_{m}' for n in (8,9) for m in ('La','Dy')]
        p={'blocks':blocks,'forces':forces,'energy':{'r':np.zeros(10),'J':np.zeros((10,12)),'names':[str(i) for i in range(10)]},
            'center':np.zeros(12),'lower':np.full(12,-1.),'upper':np.ones(12),'width':np.full(12,2.),'z':np.full(12,.5),
            'q_current':.02+SOLVE.CUSHION,'raw_center_q':.02,'bounds_lower':np.r_[np.full(12,-.025),0.],
            'bounds_upper':np.r_[np.full(12,.025),.02+SOLVE.CUSHION],'mapping_S':math.log(4),
            'electrical_unweighted':{'chi_states_a0_cubed':{'names':ch,'r':np.array([1.,2.,3.,4.]),'J':np.zeros((4,12))},
                'C_dipoles_Debye':{'names':cn,'r':np.zeros((4,3)),'J':np.zeros((4,3,12))},
                'pair_dipoles_Debye':{'names':['test_pair'],'r':np.zeros((1,3)),'J':np.zeros((1,3,12))}},
            'electrical_targets':{'chi_a0_cubed':dict(zip(ch,[.5,1.5,2.5,3.5])),'C_Debye':{k:[0.,0.,0.] for k in cn}}}
        point=np.r_[np.full(12,-.001),math.sqrt(.02)]
        original,_,names=core.constraints(p,np.r_[point[:12],point[-1]**2]);lookup=dict(zip(names,original))
        for case in SOLVE.CASES:
            m=SOLVE.assemble(core,p,case);layout=m['layout'];slack=m['b']-m['A']@point
            check(case+'_SOC_count',sum(a['kind']=='SOC' for a in layout)=={'A':17,'B':20,'C':17,'D':2}[case])
            check(case+'_old_EF_never_retained',not any(a['name']=='block/old_EF' for a in layout))
            check(case+'_electrical_norms_only_B',sum(a['name'] in ('block/pair_mu','block/C_mu','block/chi') for a in layout)==(3 if case=='B' else 0))
            check(case+'_physical_family',sum(a['name'].startswith('energy/') for a in layout)==(0 if case=='D' else 20))
            check(case+'_box_and_trust',sum(a['name'].startswith('box/') for a in layout)==(0 if case=='C' else 24) and sum(a['name'].startswith('trust/') for a in layout)==24)
            differences=[]
            for a in layout:
                if a['name'] not in lookup:continue
                s=slack[a['start']:a['stop']];v=s[0]**2-float(s[1:]@s[1:]) if a['kind']=='SOC' else s[0]
                differences.append(abs(v-lookup[a['name']]))
            check(case+'_inherited_scalar_cone_identity',max(differences)<1e-12)
        broken=copy.deepcopy(p);broken['energy']['r'][0]=2.
        d=np.r_[np.zeros(12),.019]
        row=SOLVE.certificate(core,broken,d,1.,'D',False)
        check('D_retains_physical_failure_descriptively',not row['all25_predicted'] and row['retained_numerically_feasible_1e_8'])
        row=SOLVE.certificate(core,broken,d,1.,'A',False)
        check('A_rejects_same_physical_failure',not row['retained_exact_sign_feasible'])
        edge=copy.deepcopy(p);edge['center'][0]=.999;edge['z'][0]=.9995
        # Make this control isolate only the box; all norm slopes zero.
        for b in edge['blocks'].values():b['J'][:]=0
        point=np.r_[np.array([.001]+[0.]*11),.019]
        c=SOLVE.certificate(core,edge,point,1.,'C',False);a=SOLVE.certificate(core,edge,point,1.,'A',False)
        check('C_diagnostic_outside_box_allowed_only_in_affine_scope',not c['literal_original_box_pass'] and c['retained_numerically_feasible_1e_8'] and not c['native_trial_authorized'])
        check('A_retains_same_box_face',not a['retained_exact_sign_feasible'])
        check('raw_epigraph_failure_preserved',c['raw_epigraph_margin']<0 and c['constructed_auxiliary_certificate']['certificate']['minimum_epigraph_margin']>0)
        check('auxiliary_does_not_change_physical_vector',c['raw_certificate']['parameters']==c['parameters']==c['constructed_auxiliary_certificate']['parameters'])
        check('no_CORE_eligible_used_as_case_admission',c['raw_certificate']['eligible'] is False and c['retained_numerically_feasible_1e_8'])
        class FakeEvaluator:
            def _parse(self,*a):return forbidden(*a)
            def _assemble(self,*a):return forbidden(*a)
        fake=SimpleNamespace(Evaluator=FakeEvaluator);parent=old.evaluator_class(fake);child=RUN.evaluator_class(old,fake)
        for name in ('_cell','_invoke','_import_saved','required_programs'):
            check('inherited_'+name,getattr(parent,name).__code__ is getattr(child,name).__code__)
        check('parse_and_assemble_inherited',child._parse is FakeEvaluator._parse and child._assemble is FakeEvaluator._assemble)
        obj=child.__new__(child);obj.mode='run';obj.stop_event=threading.Event();obj.ledger={'native_program_invocations':672};obj.batch={'hostname':socket.gethostname()};obj._gated=True
        rejects('cap672_fails_before_native',obj._gate)
        obj.ledger['native_program_invocations']=671;obj._gate();check('cap671_boundary_allowed',True)
        obj.mode='prepare';rejects('prepare_cannot_enter_native_gate',obj._gate)
        obj.allowed={};rejects('unknown_vector_cannot_evaluate',lambda:obj.evaluate([0.]*12))
        rejects('allowlist_cannot_expand',lambda:obj.authorize([0.]*12,'trial',1))
        check('zero_solver_native_calls',calls==[])
        status='passed_focused_final_local_controls';error=None
    except Exception as exc:status='failed_focused_final_local_controls';error=repr(exc)
    result={'status':status,'error':error,'implementation':RUN.pin(__file__),'runner':RUN.pin(RUN.__file__),'solver':RUN.pin(SOLVE.__file__),
        'checks':checks,'solver_calls':0,'native_calls':0,'artificial_only':True,'cost':{'wall_seconds':time.monotonic()-tick,'CPU_seconds':time.process_time()-cpu,'threads':1}}
    (out/'RESULT.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n');print(json.dumps({'result':RUN.pin(out/'RESULT.json'),'status':status,'error':error,'checks':len(checks)}))
    raise SystemExit(bool(error))

if __name__=='__main__':main()
