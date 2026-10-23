#!/usr/bin/env python3
"""Saved-only EXCHPOL92-row replay and independent C/P/B summary; no native calls."""
import argparse
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import time
HERE=Path(__file__).resolve().parent;D=HERE.parent
SOURCE_SHA='9eb552e6ec65da01b3e62ac803740405ee7cc05ab45d1371294fdea5b3ca422d'
READY_SHA='5f4fa09cf51822faadc64b97d8c0ba8ed1f2fbf0bb909a772d6e2ce45325c7be'
SIGNS={'aquo8_central_La':1,'aquo8_central_Dy':-1,'aquo9_central_La':-1,'aquo9_central_Dy':1}
def pin(p):
    p=Path(p).resolve();return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def verify(p):assert pin(p['path'])==p,'Changed artifact';return Path(p['path'])
def read(p):return json.loads(verify(p).read_text())
def save(p,x):
    with Path(p).open('x') as f:json.dump(x,f,indent=2,allow_nan=False);f.write('\n')
def delta(v):return math.fsum(SIGNS[k]*v[k] for k in SIGNS)
def scalar_controls(score):
    energy={group+'/'+k:v for group,rows in score['energy_controls'].items() for k,v in rows.items()}
    forces={group+'/'+k:{f:v[f] for f in ('rigid_RMS_kcal_mol_A','maximum_body_component_error_kcal_mol_A','passed')} for group,rows in score['per_fixture'].items() for k,v in rows.items()}
    assert len(energy)==score['available_energy_controls']==score['declared_energy_controls']==10
    assert len(forces)==score['available_force_controls']==score['declared_force_controls']==15
    return dict(objective=score['objective'],energy_controls=energy,force_controls=forces,energy_failed=sum(not v['passed'] for v in energy.values()),force_failed=sum(not v['passed'] for v in forces.values()),all_development_controls_passed=score['development_checks_passed'])
def compare(a,b):
    assert len(a['forces'])==len(b['forces']) and all(len(v)==3 for v in a['forces']+b['forces'])
    de=a['energy']-b['energy'];df=max(abs(x-y) for u,v in zip(a['forces'],b['forces']) for x,y in zip(u,v))
    assert math.isfinite(de) and math.isfinite(df)
    return dict(energy_difference_kcal_mol=de,maximum_force_component_difference_kcal_mol_A=df,passed=abs(de)<=1e-4 and df<=1e-3)
def cycle(rows,manifest,label):
    C=delta({gid:rows[label+'__'+gid]['energy_kcal_mol'] for gid in SIGNS})
    P=math.fsum(SIGNS[meta['gid']]*rows[cid]['energy_kcal_mol'] for cid,meta in manifest['rows'].items() if meta['kind']=='pair' and cid.startswith(label+'__'))
    ions=manifest['analytic_bare_ions'];assert all(v['status']=='analytically_zero' and v['energy_kcal_mol']==0 and v['missing_data_substitution'] is False for v in ions.values())
    P-=ions['Dy']['energy_kcal_mol']-ions['La']['energy_kcal_mol']
    return dict(C=C,P=P,B=C-P)
def markdown(r):
    lines=['# Independent EXCHPOL actual review','',f"Status: **{r['status']}**. Native application: {r.get('application_status','unavailable')}.",'']
    if r['error']:lines+=['Review error: '+r['error'],'']
    if r.get('counts'):lines+=['Counts: '+str(r['counts']),'']
    if r.get('invariance'):lines+=['Zero-q and metal-absent energy/force checks all passed at the prospectively fixed tolerances.','']
    if r.get('fixed_q_results'):
        lines+=['| Fixed q | C | P | B | ΔC | ΔP | ΔB | Failed energy / force |','|---|---:|---:|---:|---:|---:|---:|---:|']
        for label,v in r['fixed_q_results'].items():lines.append('| '+label+' | '+' | '.join(f'{x:+.8f}' for x in [*[v['cycle'][k] for k in ('C','P','B')],*[v['change_from_baseline'][k] for k in ('C','P','B')]])+f" | {v['controls']['energy_failed']}/10; {v['controls']['force_failed']}/15 |")
        lines+=['','Energies are kcal/mol, with paired order Dy9 − La9 − Dy8 + La8. Both q values are reported; none is selected.']
    lines+=['','Numerical replay does not establish physical accuracy, a unique collective mechanism or a sampling model. Original targets and held validation gates are unchanged.', '', 'This review invoked no molecular program. Native costs are retained from the run; per-program child-CPU differences are not summed.']
    return '\n'.join(lines)+'\n'
def review(config,result_path,out):
    tick,cpu=time.monotonic(),time.process_time();out=Path(out);out.mkdir(parents=True,exist_ok=False)
    r=dict(status='unavailable',reviewer=pin(__file__),source=None,result=None,error=None,raw_replay_passed=False,counts=None,invariance=None,FD_gates=None,induction=None,fixed_q_results=None,new_native_calls=0,candidate_selected=False,fit_targets_changed=False,physical_accuracy_qualified=False)
    try:
        assert pin(D/'RUN.py')['sha256']==SOURCE_SHA and pin(D/'READY.json')['sha256']==READY_SHA
        ready=read(pin(D/'READY.json'));r['source']=ready['implementation'];r['source_review']=pin(D/'independent_execution_review_v1/REVIEW.json')
        cfg=json.loads(Path(config).read_text());r['configuration']=pin(config);assert cfg['implementation']==ready['implementation']
        s=importlib.util.spec_from_file_location('exch_saved_actual',verify(cfg['implementation']));m=importlib.util.module_from_spec(s);s.loader.exec_module(m)
        assert m.fingerprint(cfg)==ready['config_fingerprint'] and cfg['preparation_result']==ready['preparation_result']
        prep=read(cfg['preparation_result']);manifest=read(prep['manifest']);r['manifest']=prep['manifest']
        r['cell_statuses']={cid:'unavailable' for cid in manifest['rows']}
        if not Path(result_path).is_file():r['status']='unavailable_result_absent'
        else:
            r['result']=pin(result_path);actual=read(r['result']);r['application_status']=actual.get('status');r['native_cost']=actual.get('cost');r['job_id']=actual.get('job_id')
            r['cell_statuses']={cid:actual.get('rows',{}).get(cid,{}).get('status','missing') for cid in manifest['rows']}
            r['reported_error']=actual.get('error');r['reported_counts']={k:actual.get(k) for k in ('native_program_invocations','internal_FD_energies_requested','internal_FD_energies_completed')}
            if actual.get('status')=='running':r['status']='unavailable_nonterminal'
            elif actual.get('status')!='collective_response_diagnostic_complete':r['status']='unavailable_failed_native_diagnostic'
            else:
                admission=m.admit_result(cfg,manifest,actual);assert admission['status']=='saved_collective_response_replay_passed'
                rows=actual['rows'];assert len(rows)==92 and all(x['status']=='complete' for x in rows.values())
                receipts=[p for row in rows.values() for p in row['receipts']];assert len(receipts)==116
                r['raw_replay_passed']=True;r['counts']=dict(cells=92,program_receipts=len(receipts),requested_FD=actual['internal_FD_energies_requested'],completed_FD=actual['internal_FD_energies_completed'])
                assert r['counts']['requested_FD']==r['counts']['completed_FD']==636
                data=read(cfg['data']);invariance={f'zero_q_{metal}':compare(rows['zero_aquo9_'+metal]['model'],data['baseline_model']['aquo9_central_'+metal]) for metal in ('La','Dy')}
                invariance['metal_absent_water']=compare(rows['shell_zero']['model'],rows['shell_q05']['model'])
                for k,v in invariance.items():assert v['passed'] and all(v[f]==actual['prechecks']['comparisons'][k][f] for f in v)
                r['invariance']=invariance;r['FD_gates']={cid:rows[cid]['numerical_check'] for cid in manifest['stages']['derivatives']}
                assert len(r['FD_gates'])==4 and all(v['passed'] and v['step_A']==1e-5 and v['absolute_tolerance_kcal_mol_A']==1e-3 and v['max_absolute_error_kcal_mol_A']<=1e-3 for v in r['FD_gates'].values())
                induction=[x for row in rows.values() for x in row['induction'].values()];assert len(induction)==116
                r['induction']=dict(programs=116,total_printed_solutions=sum(v['solutions'] for v in induction),maximum_printed_residual_Debye=max(v['maximum_printed_RMS_residual_Debye'] for v in induction),maximum_iterations=max(v['maximum_iterations'] for v in induction),tolerance_Debye=1e-8)
                assert r['induction']['maximum_printed_residual_Debye']<=1e-8
                r['fixed_q_results']={}
                for label in ('q01','q05'):
                    values=cycle(rows,manifest,label);saved=actual['summaries'][label];assert all(abs(values[k]-saved['cycle'][k])<1e-10 for k in values)
                    changes={k:values[k]-manifest['baseline_cycle'][k] for k in values};assert all(abs(changes[k]-saved['cycle_change'][k])<1e-10 for k in values)
                    r['fixed_q_results'][label]=dict(q=saved['q'],cycle=values,change_from_baseline=changes,components=saved['components'],controls=scalar_controls(saved['score']),selected=False)
                r['scale_dependence']={k:r['fixed_q_results']['q05']['change_from_baseline'][k]/.05-r['fixed_q_results']['q01']['change_from_baseline'][k]/.01 for k in ('C','P','B')}
                assert all(abs(r['scale_dependence'][k]-actual['summaries']['scale_dependence'][k])<1e-9 for k in r['scale_dependence'])
                r['status']='passed_independent_actual_review'
    except Exception as exc:r.update(status='failed_or_unavailable_review',error=repr(exc))
    r['review_cost']=dict(wall_seconds=time.monotonic()-tick,CPU_seconds=time.process_time()-cpu,threads=1,new_native_calls=0)
    save(out/'REVIEW.json',r);(out/'REPORT.md').write_text(markdown(r));print(json.dumps({k:r[k] for k in ('status','raw_replay_passed','counts','error')}));return r
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--config',required=True);p.add_argument('--result',required=True);p.add_argument('--out',required=True);v=p.parse_args();raise SystemExit(0 if review(v.config,v.result,v.out)['status']=='passed_independent_actual_review' else 2)
