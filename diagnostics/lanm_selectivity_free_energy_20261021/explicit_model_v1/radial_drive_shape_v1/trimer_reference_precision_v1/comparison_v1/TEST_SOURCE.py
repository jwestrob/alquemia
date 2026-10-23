#!/usr/bin/env python3
"""Manufactured comparison controls; no real reference arithmetic or subprocesses."""
import argparse
import copy
import json
import os
from pathlib import Path
import subprocess
import time
import COMPARE as c


def main(out):
    out = Path(out).resolve(); out.mkdir(parents=True, exist_ok=False)
    tick, cpu = time.monotonic(), time.process_time(); checks = []
    def check(ok, name):
        c.require(ok, name); checks.append(name)
    def rejects(fn, name):
        try:
            fn()
        except (ValueError, KeyError, TypeError):
            checks.append(name); return
        raise ValueError('Accepted corrupt control: '+name)
    def no_process(*a, **kw):
        raise AssertionError('Real process forbidden during manufactured controls')
    subprocess.Popen = no_process
    result = {'status': 'failed_reference_precision_comparison_controls', 'error': None,
              'implementation': c.pin(__file__), 'comparator': c.pin(c.__file__),
              'new_native_reference_solver_calls': 0, 'actual_reference_arithmetic': False}
    try:
        inputs = c.read(c.pin(c.HERE/'INPUTS.json'))
        for p in inputs.values():
            if isinstance(p, dict) and set(p) == {'path', 'sha256'}:
                c.verify(p)
        checks.append('all_frozen_input_pins')
        inv = c.read(inputs['inventory']); m = c.read(inputs['runtime_manifest'])
        check([r['cell_id'] for r in inv['rows']] == [r['cell_id'] for r in m['rows']], 'fourteen_metadata_order')
        check(all(a['old'] == b for a,b in zip(m['rows'],inv['rows'])), 'literal_old_metadata_identity')
        template = c.read(c.pin(c.HERE/'CONFIG_TEMPLATE.json'))
        check(all(template[k] is None for k in ('scout_result','all14_result','scout_root_review','all14_root_review')), 'future_result_admissions_unbound')
        helper = c.load(inputs['cycle_helper'])
        ids = [str(i) for i in range(14)]; records = {k:{'status':'complete'} for k in ids}
        c.denominator(ids,records); checks.append('complete_fourteen')
        rejects(lambda:c.denominator(ids,{k:v for k,v in records.items() if k!='0'}),'missing_state')
        bad=copy.deepcopy(records);bad['0']['status']='failed'
        rejects(lambda:c.denominator(ids,bad),'failed_state')
        rejects(lambda:c.denominator(ids[:-1]+[ids[0]],records),'duplicate_state')
        rejects(lambda:c.denominator(ids,{**records,'extra':{'status':'complete'}}),'unexpected_state')
        phase={'status':'reference_precision_all14_complete','error':None,'manifest':'M','implementation':'I',
               'phase':'remaining','native_attempts_this_phase':12,'native_attempts_total':14,'scheduler_job_id':None,'cells':records}
        c.phase(phase,ids,'M','I','remaining',12);checks.append('remaining_phase_counts')
        for key,value in [('native_attempts_this_phase',14),('native_attempts_total',13),('status','started'),('scheduler_job_id','fake'),('error','failed')]:
            bad=copy.deepcopy(phase);bad[key]=value
            rejects(lambda b=bad:c.phase(b,ids,'M','I','remaining',12),'reject_phase_'+key)
        names=['Energy change','MAX-Density change','RMS-Density change','Orbital Gradient','Orbital Rotation']
        limits=['1e-11','1e-8','1e-9','1e-8','1e-8']
        cur={'status':'passed_current_scf_residuals','all_required_passed':True,'mode':'All-Criteria','algorithm':'SOSCF',
             'required_checks':{k:True for k in names},'summary':{k:{'value_token':'0','tolerance_token':v} for k,v in zip(names,limits)}}
        c.current_gate(cur);checks.append('saved_current_residual_pass')
        for key,value in [('Energy change','-1.01e-11'),('MAX-Density change','-1e-12'),('Orbital Rotation','1.1e-8')]:
            bad=copy.deepcopy(cur);bad['summary'][key]['value_token']=value
            rejects(lambda b=bad:c.current_gate(b),'reject_residual_'+key)
        bad=copy.deepcopy(cur);bad['algorithm']='TRAH'
        rejects(lambda:c.current_gate(bad),'unknown_algorithm')
        bad=copy.deepcopy(cur);del bad['required_checks']['Orbital Rotation']
        rejects(lambda:c.current_gate(bad),'missing_current_required_field')
        p=out/'synthetic_energy.out';p.write_text('FINAL SINGLE POINT ENERGY -12.345600\n')
        printed=c.energy_printing(c.pin(p),-12.3456)
        check(printed['half_ulp_hartree']==0.0000005,'literal_energy_half_ULP')
        rejects(lambda:c.energy_printing(c.pin(p),-12.34),'wrong_energy_stage')
        p.write_text('FINAL SINGLE POINT ENERGY -12.345600\nFINAL SINGLE POINT ENERGY -12.345600\n')
        rejects(lambda:c.energy_printing(c.pin(p),-12.3456),'duplicate_print_token')
        atom={'charge':3,'symbols':['La'],'ghost_mask':[False],'expected_electrons':8}
        moment={'output':'O','property_file':'P','charge':3,'labels':['La'],'ghost_mask':[False],'native_electrons':8,
                'total_dipole_au':[1.,2.,-3.],'electronic_dipole_au':[-2.,-3.,-4.],'nuclear_dipole_au':[3.,5.,1.],
                'precise_origin_Bohr':[0.,0.,.123456789],'printed_origin_Bohr':[0.,0.,.123457],
                'print_bounds_au':{'total_dipole_au':1.5e-6}}
        c.moment_gate(moment,atom,'O','P');checks.append('precise_charged_origin_identity')
        bad=copy.deepcopy(moment);bad['charge']=0
        rejects(lambda:c.moment_gate(bad,atom,'O','P'),'changed_charged_origin_identity')
        # Deliberately asymmetric fabricated values exercise every signed component.
        terms={}; fake_inv={'rows':[{'cell_id':k} for k in ids],'cycles':{}}; models={'comparison':{'cycles':{}}}
        coeff=[1,-1,-1,-1,1,1,1]
        for metal,cycle in enumerate(c.CYCLES):
            subset=ids[7*metal:7*(metal+1)]; co=dict(zip(subset,coeff))
            for j,k in enumerate(subset):
                base=[float(j+1),float(2*j-5),float(3-j)]
                old={'electronic_energy_hartree':float(j*j+metal),'moment':{'charge':0,'total_dipole_au':base,
                     'print_bounds_au':{'total_dipole_au':1e-6}},'energy_printing':{'half_ulp_hartree':5e-13}}
                new=copy.deepcopy(old)
                if j==0:
                    new['electronic_energy_hartree']+=float(metal+1)
                    new['moment']['total_dipole_au']=[base[0]+metal+1,base[1]-2*(metal+1),base[2]+3*(metal+1)]
                terms[k]={'old':old,'new':new}
            # Closed-form expectations for the fabricated seven terms, independent of combine.
            E=helper.HA*(sum(coeff[j]*j*j for j in range(7))+metal)
            mu=[helper.CONVERSION*sum(coeff[j]*v for j,v in enumerate(vals))
                for vals in ([j+1 for j in range(7)],[2*j-5 for j in range(7)],[3-j for j in range(7)])]
            fake_inv['cycles'][cycle]={'reference_terms':co,'reference_electronic_E3_kcal_mol':E,
                'reference_mu3_Debye':mu,'reference_print_bound_Debye_per_component':helper.CONVERSION*7e-6}
            models['comparison']['cycles'][cycle]={'status':'complete','reference_electronic_E3_kcal_mol':E,
                'reference_mu3_Debye':mu,'models':{label:{'E3_kcal_mol':100.+i+metal,'mu3_Debye':[i+.5,2.,-4.-metal]}
                for i,label in enumerate(c.LABELS)}}
        original_models=copy.deepcopy(models)
        cycles,paired=c.arithmetic(terms,fake_inv,models,helper)
        check(abs(cycles[c.CYCLES[0]]['measured_protocol_shift']['E3_kcal_mol']-helper.HA)<1e-10,'seven_term_energy_shift')
        expected=[helper.CONVERSION,-2*helper.CONVERSION,3*helper.CONVERSION]
        check(all(abs(a-b)<1e-12 for a,b in zip(cycles[c.CYCLES[0]]['measured_protocol_shift']['mu3_Debye'],expected)), 'all_signed_vector_components')
        check(abs(paired['measured_protocol_shift']['E3_kcal_mol']-helper.HA)<1e-10,'Dy_minus_La_orientation')
        check(all(abs(a-b)<1e-12 for a,b in zip(paired['measured_protocol_shift']['mu3_Debye'],expected)), 'paired_signed_vector')
        check(models==original_models,'model_values_unchanged')
        check(set(paired['fixed_models'])==set(c.LABELS),'all_three_fixed_model_labels')
        check(abs(c.norm([1.,-2.,3.])-14**.5)<1e-14,'vector_norm_retains_three_components')
        rejects(lambda:c.vector([1.,2.]),'truncated_dipole')
        rejects(lambda:c.vector([1.,2.,float('nan')]),'nonfinite_dipole')
        wrong=copy.deepcopy(fake_inv);wrong['cycles'][c.CYCLES[0]]['reference_terms']['0']=-1
        rejects(lambda:c.arithmetic(terms,wrong,models,helper),'changed_seven_term_sign')
        # The actual source path is exercised only with null future results. It must
        # fail before scientific arithmetic and preserve all fourteen missing IDs.
        cfg=out/'unbound.json';cfg.write_text(json.dumps(template))
        rel=out/'synthetic_release.json';rel.write_text(json.dumps({'status':'approved_saved_reference_precision_comparison',
            'configuration':c.pin(cfg),'implementation':c.pin(c.__file__),'inputs':template['inputs']}))
        for k in c.THREADS:os.environ[k]='1'
        failed=c.run(cfg,rel,out/'incomplete_control')
        check(failed['status']=='failed_saved_reference_precision_comparison' and
              len(failed['admission_progress']['missing_states'])==14 and 'cycles' not in failed,'incomplete_receipt_no_cycles')
        result['status']='passed_reference_precision_comparison_controls'
    except Exception as exc:
        result['error']=repr(exc)
    result.update(checks=checks,check_count=len(checks),cost={'wall_seconds':time.monotonic()-tick,'CPU_seconds':time.process_time()-cpu,'threads':1})
    (out/'RESULT.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({'status':result['status'],'checks':len(checks),'error':result['error']}))
    return 0 if result['error'] is None else 1


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--out',required=True);args=ap.parse_args()
    raise SystemExit(main(args.out))
