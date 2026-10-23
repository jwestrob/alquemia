#!/usr/bin/env python3
"""Frozen raw admission plus independent actual seven-term and 25-gate arithmetic."""
import argparse
import importlib.util
import json
import math
import os
from pathlib import Path
import time

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[4]
W=ROOT/'workspaces/lanm_selectivity_free_energy_20261021/explicit_model_v1/metal_pair_response_v1/sensitivity_v1'
s=importlib.util.spec_from_file_location('sensitivity_review_admission',HERE/'RUN.py')
a=importlib.util.module_from_spec(s);s.loader.exec_module(a)


def cross(x,y):return [x[1]*y[2]-x[2]*y[1],x[2]*y[0]-x[0]*y[2],x[0]*y[1]-x[1]*y[0]]


def rigid_error(row,model,reference):
    """Independent water/ion category-balanced force and centroid-torque RMS."""
    raw=[];weighted=0.;water_count=sum(b['category']=='water' for b in row['bodies'])
    for kind in ('ion','water'):
        for body in row['bodies']:
            if body['category']!=kind:continue
            ids=body['indices'];positions=row['geometry']['positions_A']
            center=[sum(positions[i][j] for i in ids)/len(ids) for j in range(3)]
            difference={i:[model[i][j]-reference[i][j] for j in range(3)] for i in ids}
            net=[sum(difference[i][j] for i in ids) for j in range(3)]
            if kind=='ion':components=net;weight=1/6
            else:
                torques=[cross([positions[i][j]-center[j] for j in range(3)],difference[i]) for i in ids]
                components=net+[sum(t[j] for t in torques) for j in range(3)];weight=1/(12*water_count)
            raw.extend(components);weighted+=weight*sum(v*v for v in components)
    return math.sqrt(weighted),max(abs(v) for v in raw),raw


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--out',required=True,type=Path);args=parser.parse_args()
    args.out.mkdir(parents=True,exist_ok=False);tick=time.monotonic();cpu=time.process_time();checks=[]
    def check(name,ok):
        checks.append({'name':name,'passed':bool(ok)})
        a.require(ok,name)
    report={'status':'started','implementation':a.pin(__file__),'checks':checks,'new_native_calls':0,'error':None}
    try:
        check('one local thread',all(os.environ.get(k)=='1' for k in a.THREADS))
        cp=a.pin(W/'run_v1/CONFIG.json');rp=a.pin(W/'run_v1/execution/RESULT.json')
        c=a.read(cp);r=a.read(rp);m=a.read(c['manifest'])
        report.update(configuration=cp,result=rp,manifest=c['manifest'],adapter=c['implementation'])
        check('terminal fixed job',r['job_id']=='18413' and r['status']=='consumed_pair_response_sensitivity_complete' and r['error'] is None)
        a.release_gate(c,m)
        admitted=a.admit_result(c,m,r)
        check('complete saved raw admission',admitted['status']=='saved_sensitivity_admitted' and admitted['new_native_calls']==0 and len(admitted['rows'])==120 and r['native_program_invocations']==236 and r['reused_gradient_programs']==4 and r['internal_FD_energies']==0)
        check('four retained profiles',set(r['profiles'])==set(a.PROFILES) and all(not p['selected_candidate'] and not p['full_shell_B_computed'] for p in r['profiles'].values()))
        baseline=m['baseline_model'];references={k:v['reference'] for k,v in m['training_rows'].items()};summary={}
        # The frozen parser has admitted raw E/F/moments. Arithmetic below does
        # not call its summarizer or objective implementation.
        for profile in a.PROFILES:
            rows={gid:r['rows'][profile+'__'+gid] for gid in m['rows']}
            model={gid:{'energy':v['moments']['energy_kcal_mol'],'forces':v['gradient']['forces']} for gid,v in rows.items()}
            score=r['profiles'][profile]['original25_controls']
            def delta(values,g):return values[g+'_Dy']['energy']-values[g+'_La']['energy']
            energies={g:(delta(model,g)-delta(model,'aquo8_central'))-(delta(references,g)-delta(references,'aquo8_central')) for g in ('water_minus','water_central','water_plus','aquo9_central')}
            works={metal+'_'+sign:(model['water_'+sign+'_'+metal]['energy']-model['water_central_'+metal]['energy'])-(references['water_'+sign+'_'+metal]['energy']-references['water_central_'+metal]['energy']) for metal in ('La','Dy') for sign in ('minus','plus')}
            paired={sign:works['Dy_'+sign]-works['La_'+sign] for sign in ('minus','plus')}
            groups={'balanced_exchange':energies,'individual_radial_work':works,'paired_radial_work':paired};energy_pass=0;changes={'newly_passed':[],'newly_failed':[]}
            for group,values in groups.items():
                for key,error in values.items():
                    saved=score['energy_controls'][group][key];passed=abs(error)<=1
                    check(profile+' energy '+group+' '+key,abs(error-saved['error_kcal_mol'])<1e-10 and passed==saved['passed']);energy_pass+=passed
                    prior=m['baseline_score']['energy_controls'][group][key]['passed']
                    if passed!=prior:changes['newly_passed' if passed else 'newly_failed'].append('energy/'+group+'/'+key)
            force_pass=0;force_rms=[];force_checks={}
            for kind in ('La','Dy','paired'):
                for gid in ('water_minus','water_central','water_plus','aquo8_central','aquo9_central'):
                    if kind=='paired':
                        mf=[[dy-la for dy,la in zip(d,l)] for d,l in zip(model[gid+'_Dy']['forces'],model[gid+'_La']['forces'])]
                        rf=[[dy-la for dy,la in zip(d,l)] for d,l in zip(references[gid+'_Dy']['forces'],references[gid+'_La']['forces'])]
                    else:mf=model[gid+'_'+kind]['forces'];rf=references[gid+'_'+kind]['forces']
                    rms,maximum,raw=rigid_error(m['training_rows'][gid+'_La'],mf,rf);old=score['per_fixture'][kind][gid];passed=rms<=2 and maximum<=5
                    check(profile+' force '+kind+' '+gid,abs(rms-old['rigid_RMS_kcal_mol_A'])<1e-10 and abs(maximum-old['maximum_body_component_error_kcal_mol_A'])<1e-10 and max(abs(x-y) for x,y in zip(raw,old['raw_body_components_kcal_mol_A']))<1e-10 and passed==old['passed'])
                    force_rms.append(rms);force_pass+=passed;force_checks[kind+'__'+gid]={'RMS':rms,'maximum':maximum,'passed':passed}
                    prior=m['baseline_score']['per_fixture'][kind][gid]['passed']
                    if passed!=prior:changes['newly_passed' if passed else 'newly_failed'].append('force/'+kind+'/'+gid)
            objective=sum(v*v for v in list(energies.values())+list(works.values()))/32+sum(v*v for v in force_rms)/80
            check(profile+' independent objective and all25 gate',abs(objective-score['objective'])<1e-9 and (energy_pass==10 and force_pass==15)==score['development_checks_passed'])
            moments=dict(m['baseline_moments']);moments.update({gid:v['moments'] for gid,v in rows.items()});cycles={}
            for cid,cycle in m['trimer_cycles'].items():
                parent=cycle['parent_geometry_id'];metal=parent.rsplit('_',1)[1];kind=cid.rsplit('__',1)[1];i,j=cycle['selection']['water_indices']
                ordered=[(cid+'__metal_two_water',1),(parent+f'__dimer_{i:02d}',-1),(parent+f'__dimer_{j:02d}',-1),(cid+'__two_water',-1),('ion_'+metal,1),(parent+f'__water_{i:02d}',1),(parent+f'__water_{j:02d}',1)]
                check(profile+' '+cid+' literal seven-term signs',dict(ordered)==cycle['model_terms'])
                energy=sum(sign*moments[k]['energy_kcal_mol'] for k,sign in ordered)
                mu=[sum(sign*moments[k]['induced_sum_Debye'][axis] for k,sign in ordered) for axis in range(3)]
                low=[sum(sign*moments[k]['total_lab_Debye'][axis] for k,sign in ordered) for axis in range(3)]
                old=r['profiles'][profile]['trimer_cycles'][cid];ref=m['trimer_reference'][cid]
                tolerance=cycle['M_print_bound_Debye_per_component']+cycle['induced_print_bound_Debye_per_component']+14e-9
                check(profile+' '+cid+' independent E3/mu3/M',abs(energy-old['model_E3_kcal_mol'])<1e-10 and max(abs(x-y) for x,y in zip(mu,old['model_mu3_Debye']))<1e-12 and max(abs(x-y) for x,y in zip(low,old['corrected_M_mu3_Debye']))<1e-12 and max(abs(x-y) for x,y in zip(mu,low))<=tolerance)
                check(profile+' '+cid+' fixed reference and signs',old['reference_electronic_E3_kcal_mol']==ref['reference_electronic_E3_kcal_mol'] and abs((energy-ref['reference_electronic_E3_kcal_mol'])-old['model_minus_reference_E3_kcal_mol'])<1e-10 and old['reference_mu3_Debye']==ref['reference_mu3_Debye'])
                cycles[cid]={'E3':energy,'baseline_E3':ref['model_E3_kcal_mol'],'reference_E3':ref['reference_electronic_E3_kcal_mol'],'E3_error':energy-ref['reference_electronic_E3_kcal_mol'],'mu3_Debye':mu,'baseline_mu3_Debye':ref['model_mu3_Debye'],'reference_mu3_Debye':ref['reference_mu3_Debye']}
            centers={}
            for gid in a.CENTERS:
                metal=gid.rsplit('_',1)[1];terms=[(gid,1),('ion_'+metal,-1),(gid+'__shell',-1)]
                check(profile+' '+gid+' C cancellation',dict(terms)==m['center_cancellation'][gid]['terms'])
                energy=sum(sign*moments[k]['energy_kcal_mol'] for k,sign in terms);mu=[sum(sign*moments[k]['induced_sum_Debye'][axis] for k,sign in terms) for axis in range(3)];old=r['profiles'][profile]['centers'][gid]
                check(profile+' '+gid+' independent C',abs(energy-old['C_energy_kcal_mol'])<1e-10 and max(abs(x-y) for x,y in zip(mu,old['C_mu_Debye']))<1e-12)
                centers[gid]={'C_energy':energy,'C_mu_Debye':mu,'reference_C_mu_Debye':old['reference_C_mu_Debye'],'baseline_C_mu_Debye':old['baseline_C_mu_Debye']}
            summary[profile]={'objective':objective,'energy_gates_passed':energy_pass,'force_gates_passed':force_pass,'changes_from_baseline':changes,'energy_errors':groups,'force_checks':force_checks,'development_passed':energy_pass==10 and force_pass==15,'cycles':cycles,'centers':centers,'no_profile_selected':True}
        maximum=max(v['moments']['induction']['maximum_printed_RMS_residual_Debye'] for v in r['rows'].values())
        maximum=max(maximum,max(v.get('gradient_induction',{}).get('maximum_printed_RMS_residual_Debye',0) for v in r['rows'].values()))
        closure=max(abs(v['gradient']['energy']-v['moments']['energy_kcal_mol']) for v in r['rows'].values())
        check('all120 energy force closures and induction',closure<=1e-4 and maximum<=1e-8)
        report.update(status='actual_sensitivity_review_passed',raw_admission_status=admitted['status'],profiles=summary,raw_cells_admitted=120,new_programs_admitted=236,inherited_gradients_admitted=4,new_internal_FD_energies=0,maximum_new_induction_RMS_Debye=maximum,maximum_energy_force_program_energy_difference=closure,original_application_cost=r['cost'],limitations=['Fixed consumed sensitivity, no fit or selected candidate.','All four profiles retained, including worsening or failed physical checks.','Eight sampled E3 cycles do not reconstruct full-shell B.','Original provisional reference, SCF, ECP/frozen-f and printed-moment uncertainty limits remain.'])
    except Exception as exc:report.update(status='actual_sensitivity_review_failed',error=repr(exc))
    finally:
        report.update(check_count=len(checks),review_cost={'wall_seconds':time.monotonic()-tick,'CPU_seconds':time.process_time()-cpu,'numerical_threads':1})
        a.save(args.out/'REVIEW.json',report)
    print(json.dumps({k:report[k] for k in ('status','check_count','error')}))
    if report['error']:raise SystemExit(1)


if __name__=='__main__':main()
