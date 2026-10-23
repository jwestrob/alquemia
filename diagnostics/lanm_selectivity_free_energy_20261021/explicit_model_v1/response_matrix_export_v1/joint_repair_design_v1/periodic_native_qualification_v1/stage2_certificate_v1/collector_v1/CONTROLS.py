"""Manufactured parser, unequal-dual and finite-step controls; no physical files."""
import copy
from decimal import Decimal
import importlib.util
import json
from pathlib import Path
import sys
import time


def main(out):
    tick,cpu=time.monotonic(),time.process_time();out=Path(out);out.mkdir(parents=True,exist_ok=False)
    spec=importlib.util.spec_from_file_location('collector',Path(__file__).with_name('COLLECT.py'))
    m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);checks=[]
    def check(name,value):
        checks.append(dict(label=name,passed=bool(value)))
        if not value:raise AssertionError(name)
    def rejects(name,fn):
        try:fn()
        except (ValueError,KeyError):check(name,True)
        else:check(name,False)
    def f(value):return format(Decimal(str(value)),'.17E')
    def output(role='ACTIVE'):
        lines=['STAGE2_ONE_ANALYTIC_GRADIENT_V1','UNITS E=kcal/mol G=positive_dE/dr_kcal/mol/angstrom VIRIAL=kcal/mol',
               'TOTAL_ENERGY '+f(10),'POLARIZATION_ENERGY '+f('-3.34'),
               'CONTRACTION_CONSTANTS electric dielec '+f(2)+' '+f(1),
               'TOTAL_GRADIENT 1 '+' '.join(map(f,(2,3,4)))]
        for i,row in enumerate(((5,0,0),(0,6,0),(0,0,7)),1):lines.append('TOTAL_VIRIAL '+str(i)+' '+' '.join(map(f,row)))
        lines+=['ROLE '+role,'STATE_ID '+('central' if role=='ACTIVE' else role),'INPUT_BINDING '+'a'*64]
        if role=='BYPASS_OLD_ABSENT':lines+=['INTERNAL_COUNTERS unavailable_unchanged_old_archive']
        else:
            lines+=['COUNTS induced_solve direct mutual 1 1 7','COUNTS extra_diagnostic_operator_applications=0 internal_FD=0',
                    'PCG_PRE_PEEK 3 '+' '.join(map(f,('1e-10','1e-22','2e-22','1e-9')))]
            if role=='ACTIVE':lines+=['SHARED_H_READBACK_MAX_ERROR '+f(0),
                'RAW_DIRECT 1 '+' '.join(map(f,(1,2,3,4,5,6))),
                'CORRECTED_RHS 1 '+' '.join(map(f,('1.1','2.2','3.3','4.1','5.2','6.3'))),
                'FINAL_POST_PEEK 1 '+' '.join(map(f,(1,'.1','.2','.3','.4','.5','.6'))),
                'CORRECTION 1 '+' '.join(map(f,('.4','.1','.2','.3')))]
        return '\n'.join(lines+['END gradient_calls=1 energy_only_calls=0 internal_FD=0',''])
    def parse(text):
        path=out/'manufactured.dat';path.write_text(text);return m.parse_gradient(path,n=1)
    try:
        text=output();center=parse(text)
        check('Unequal channels retained',center['polarization_contraction']['maximum_u_v_difference']==.3)
        check('Stationary dual contraction passes',center['polarization_contraction']['passed'])
        check('Both shared-h channels pass',center['shared_h']['passed'])
        check('Actual mutual count retained',center['counters']['mutual_applications']==7)
        check('Half literal decimal unit',m.rounding('1.2300E-02')==Decimal('0.0000005'))
        check('Zero token precision retained',m.rounding('0.00000000000000000E+00')==Decimal('5e-18'))
        wrong=text.replace('POLARIZATION_ENERGY '+f('-3.34'),'POLARIZATION_ENERGY '+f('-1.54'))
        check('Incorrect equal-channel energy rejected',parse(wrong)['polarization_contraction']['passed'] is False)
        wrong=text.replace('CORRECTED_RHS 1 '+f('1.1'),'CORRECTED_RHS 1 '+f('1.11'))
        check('Wrong shared direct channel rejected',parse(wrong)['shared_h']['passed'] is False)
        rejects('Duplicate/missing atom order',lambda:parse(text.replace('TOTAL_GRADIENT 1','TOTAL_GRADIENT 2')))
        rejects('Iteration-cap termination rejected',lambda:parse(text.replace('PCG_PRE_PEEK 3','PCG_PRE_PEEK 100')))
        rejects('Unconverged original polar-eps rejected',lambda:parse(text.replace(f('1e-10'),f('2e-8'))))
        rejects('Trailing displaced state rejected',lambda:parse(text+'EXTRA\n'))
        rejects('NaN rejected',lambda:parse(text.replace('TOTAL_ENERGY '+f(10),'TOTAL_ENERGY NaN')))
        rejects('Unknown role rejected',lambda:parse(text.replace('ROLE ACTIVE','ROLE OTHER')))
        parsed_bypasses={role:parse(output(role)) for role in m.BYPASSES}
        check('Old archive counters unknown',parsed_bypasses[m.BYPASSES[0]]['counters'] is None)
        check('Three bypass tokens exact',len({v['E_G_virial_token_sha256'] for v in parsed_bypasses.values()})==1)
        analytic,b=m.projection(center,dict(kind='translation'),[1,0,0])
        check('Positive gradient sign',analytic==2)
        strain=dict(kind='strain',strain_matrix=[[0,0,0],[0,1,0],[0,0,0]],L0_A=1.)
        check('Strain uses virial instead of Cartesian-only slope',m.projection(center,strain,[1,0,0])[0]==6)
        energies={s:(f(10-2*s),f(10+2*s)) for s in m.STEPS}
        slopes=m.slope_checks('manufactured',analytic,b,energies)
        check('Both steps and comparison pass linear energy',slopes['passed'])
        check('Reversed analytic sign fails',not m.slope_checks('sign',-analytic,b,energies)['passed'])
        broken=dict(energies);s=m.STEPS[1];broken[s]=(f(10-Decimal('2.01')*s),f(10+Decimal('2.01')*s))
        check('Two-step inconsistency retained',not m.slope_checks('step',analytic,b,broken)['two_step']['passed'])
        expected={'central':dict(role='ACTIVE')};parsed={'central':center,**parsed_bypasses};directions=[]
        for role in m.BYPASSES:expected[role]=dict(role=role)
        for i,name in enumerate(m.DIRECTION_IDS):
            d=dict(id=name,kind='linear',coefficients_dr_dq=None);slope=Decimal(2)
            if name.startswith('strain'):
                axis=int(name[-1]);matrix=[[0.,0.,0.] for _ in range(3)];matrix[axis][axis]=1.
                d.update(kind='strain',strain_matrix=matrix,L0_A=1.);slope=Decimal(5+axis)
            directions.append(d)
            for step in m.STEPS:
                for sign in (-1,1):
                    key=f'{name}:{step}:{sign}';expected[key]=dict(role='ACTIVE',direction=name,step=step,sign=sign)
                    value=copy.deepcopy(center);value['energy_token']=f(10+sign*step*slope);parsed[key]=value
        m.load_tangent=lambda item:[1,0,0]
        full=m.summarize_panel(parsed,expected,dict(directions=directions),'full')
        check('All36 cells yield exactly24 FD checks',full['complete'] and full['passed'] and full['FD_checks_completed']==24)
        missing=dict(parsed);del missing[next(k for k in missing if k.startswith('metal_0:'))]
        partial=m.summarize_panel(missing,expected,dict(directions=directions),'full')
        check('Missing endpoint cannot pass or omit denominator',not partial['complete'] and not partial['passed'] and partial['FD_checks_completed']==21)
        scout=m.summarize_panel({k:parsed[k] for k in ('central',)+m.BYPASSES},
            {k:expected[k] for k in ('central',)+m.BYPASSES},dict(directions=directions),'central_scout')
        check('Central scout does not claim finite differences',scout['passed'] and scout['FD_checks_completed']==0)
        bad=copy.deepcopy(parsed);bad[m.BYPASSES[2]]['E_G_virial_token_sha256']='changed'
        check('Bypass token mismatch fails full panel',not m.summarize_panel(bad,expected,dict(directions=directions),'full')['passed'])
        status='passed_stage2_collector_manufactured_controls';error=None
    except Exception as exc:status='failed_stage2_collector_manufactured_controls';error=repr(exc)
    result=dict(status=status,error=error,implementation=m.pin(__file__),collector=m.pin(Path(__file__).with_name('COLLECT.py')),
                checks=checks,check_count=len(checks),new_native_calls=0,new_FFT_replays=0,new_factorizations=0,
                actual_physical_rows_read=0,cost=dict(wall_seconds=time.monotonic()-tick,CPU_seconds=time.process_time()-cpu,threads=1))
    (out/'RESULT.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n');print(json.dumps(result))
    return 0 if error is None else 1


if __name__=='__main__':sys.exit(main(sys.argv[1]))
