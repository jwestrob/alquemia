#!/usr/bin/env python3
"""Focused saved/synthetic controls; this harness never invokes a native program."""
import argparse
import copy
import importlib.util
import json
import math
from pathlib import Path
import re
import time

HERE=Path(__file__).resolve().parent
s=importlib.util.spec_from_file_location('field_test_source',HERE/'REFERENCE.py')
r=importlib.util.module_from_spec(s);s.loader.exec_module(r)


def main(panel,out):
    tick,cpu=time.monotonic(),time.process_time();out.mkdir(parents=True,exist_ok=False)
    p=r.load(panel);checks={}
    def check(name,ok):
        checks[name]=bool(ok)
        if not ok:raise AssertionError(name)
    def rejects(name,fn):
        try:fn()
        except (ValueError,AssertionError,KeyError):check(name,True)
        else:check(name,False)
    precision={'energy_hartree':1e-14,'moment_z_au':1e-12}
    zero=dict(energy_hartree=-100.,dipole_au=[0,0,.25],precision=precision)
    rr={tag:dict(energy_hartree=-100.-.25*f-10*f*f,dipole_au=[0,0,.25+20*f],precision=precision)
        for tag,f in r.FIELDS.items()}
    est=r.estimate(rr,zero);check('quadratic_all_six_gates',r.consistency(est,.5,True)['passed'])
    check('curvature_sign_and_factor',max(abs(est[k]['chi_energy']-20) for k in est)<1e-7)
    bad=copy.deepcopy(rr);bad['plus_h']['dipole_au'][2]+=.01
    check('energy_moment_disagreement_fails',not r.consistency(r.estimate(bad,zero),.5,True)['passed'])
    bad=copy.deepcopy(rr)
    for tag,f in r.FIELDS.items():bad[tag]['energy_hartree']+=2e6*f**4
    check('quartic_step_mismatch_fails',not r.consistency(r.estimate(bad,zero),.5,True)['passed'])
    bad=copy.deepcopy(rr)
    for tag,f in r.FIELDS.items():
        bad[tag]['energy_hartree']-=5000*f**3/6;bad[tag]['dipole_au'][2]+=2500*f*f
    check('nonlinear_odd_slope_is_retained_and_fails',not r.consistency(r.estimate(bad,zero),.5,True)['passed'])
    bad=copy.deepcopy(rr);bad['plus_h']['energy_hartree']=math.nan
    rejects('nonfinite_response_rejected',lambda:r.estimate(bad,zero))
    combined=r.combine({'a':est,'b':est},{'a':1,'b':-1})
    check('signed_difference_zero',combined['h']['chi_energy']==0)
    check('bounds_use_absolute_coefficients',combined['h']['roundoff_bounds']['odd_minus_zero']==2*est['h']['roundoff_bounds']['odd_minus_zero'])
    check('decimal_precision',r.token_roundoff('-1.2340e+02')>=.005)
    rows={gid+'__Fz_'+tag:dict(properties=copy.deepcopy(v)) for gid in r.BASES for tag,v in rr.items()}
    fake=copy.deepcopy(p)
    for gid in r.BASES:fake['zero_replays'][gid]=zero
    full=r.score(rows,fake,'full');check('all_four_and_paired_quadratic_pass',full['numerical_checks_passed'])
    check('conditional_and_paired_zero',full['increments']['La']['estimates']['h']['chi_energy']==full['paired']['estimates']['h']['chi_energy']==0)
    rows.pop(next(iter(rows)));missing=r.score(rows,fake,'full')
    check('missing_retains_denominator',len(missing['per_state'])==4 and not missing['numerical_checks_passed'])
    c=p['cells'][0];z=p['zero_replays'][c['zero_id']];record=p['zero_field_records'][c['zero_id']]
    old=r.b.verify(record['output']).read_text();property_text=r.b.verify(record['property_file']).read_text()
    # Synthetic origin-zero property/output fixture, derived only from admitted zero data.
    replacements=[('dipoleTotal','total_dipole_au','Total Dipole Moment'),
                  ('dipoleNucContrib','nuclear_dipole_au','Nuclear contribution'),
                  ('dipoleElecContrib','electronic_dipole_au','Electronic contribution')]
    for field,key,label in replacements:
        vals=z['moment'][key];head,tail=property_text.split('&'+field+' ',1)
        block,sep,rest=tail.partition('\n   &')
        for i,value in enumerate(vals):block=re.sub(r'^\s*'+str(i)+r'\s+('+r.moment.NUMBER+r')\s*$',f'{i} {value:.17e}',block,flags=re.M)
        property_text=head+'&'+field+' '+block+sep+rest
        old=re.sub(r'^\s*'+re.escape(label)+r'\s*:.*$',label+' : '+' '.join(f'{x:.9f}' for x in vals),old,flags=re.M)
    old=re.sub(r'^Position of electric origin.*$','Position of electric origin ... 0.000000 0.000000 0.000000',old,flags=re.M)
    echo=f'| 900>  EField 0.0,0.0,{c["field_au"]:.7f}\n| 901>  EFieldOrigin 0.0,0.0,0.0\n| 902>  Origin 0.0,0.0,0.0\n'
    old=old.replace('****END OF INPUT****','\n'+echo+'****END OF INPUT****',1)
    op=out/'synthetic_origin_zero.out';pp=out/'synthetic_origin_zero.property.txt';op.write_text(old);pp.write_text(property_text)
    parsed=r.parse(old,c,0);check('synthetic_declared_field_echo',parsed['status']=='complete')
    field=r.field_properties(op,pp,c,parsed);check('synthetic_origin_zero_full_precision',r.moment.difference(field['dipole_au'],z['dipole_au'])<1e-12)
    wrong=old.replace(f'EField 0.0,0.0,{c["field_au"]:.7f}','EField 0.0,0.0,0.0010000',1)
    check('wrong_field_sign_rejected',not r.parse(wrong,c,0)['checks']['declared_field_echo'])
    badc=copy.deepcopy(c);badc['expected_aux_J_basis_functions']+=1
    check('wrong_aux_count_rejected',not r.parse(old,badc,0)['checks']['full_parent_aux_J_basis'])
    op.write_text(old.replace('Position of electric origin ... 0.000000','Position of electric origin ... 0.100000'))
    rejects('wrong_property_origin_rejected',lambda:r.field_properties(op,pp,c,parsed));op.write_text(old)
    pp.write_text(property_text.replace('&Charge [&Type "Integer"] 3','&Charge [&Type "Integer"] 0'))
    rejects('wrong_property_charge_rejected',lambda:r.field_properties(op,pp,c,parsed));pp.write_text(property_text)
    pp.write_text(property_text.replace('&FinalEnergy [&Type "Double"]', '&MissingEnergy [&Type "Double"]'))
    rejects('missing_property_energy_rejected',lambda:r.field_properties(op,pp,c,parsed));pp.write_text(property_text)
    root=out/'synthetic_claims';cid=r.SCOUTS[0];ids=[x['cell_id'] for x in p['cells']]
    r.claim(root,root/'tasks'/cid/'execution',cid,ids)
    rejects('duplicate_attempt_rejected',lambda:r.claim(root,root/'tasks'/cid/'execution',cid,ids))
    r.stop(root,cid,'synthetic_failure');nextid=r.SCOUTS[1]
    rejects('failure_stops_new_claims',lambda:r.claim(root,root/'tasks'/nextid/'execution',nextid,ids))
    cap=out/'synthetic_cap';cap.mkdir();r.save(cap/'NATIVE_LEDGER.json',ids)
    rejects('sixteen_native_ceiling',lambda:r.claim(cap,cap/'tasks'/cid/'execution',cid,ids))
    # The original COM zero is intentionally invalid for the new explicit-origin parser.
    rejects('old_COM_not_silently_accepted',lambda:r.field_properties(r.b.verify(record['output']),r.b.verify(record['property_file']),c,z['native']))
    result=dict(status='passed_saved_controls',implementation=r.b.pin(r.__file__),test_implementation=r.b.pin(__file__),
                panel=r.b.pin(panel),checks=checks,passed=len(checks),new_native_calls=0,native_subprocess_invocations=0,
                synthetic_fixtures_not_native_evidence=True,local_threads=1,
                cost=dict(wall_seconds=time.monotonic()-tick,CPU_seconds=time.process_time()-cpu))
    print(json.dumps(r.save(out/'RESULT.json',result)))


if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('--panel',required=True);a.add_argument('--out',type=Path,required=True);x=a.parse_args();main(x.panel,x.out)
