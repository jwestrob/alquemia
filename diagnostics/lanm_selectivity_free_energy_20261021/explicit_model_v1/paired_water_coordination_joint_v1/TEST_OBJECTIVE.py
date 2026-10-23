#!/usr/bin/env python3
"""Saved and analytic controls only; never invokes a molecular engine."""
import copy
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import time

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[3]
W=ROOT/'workspaces/lanm_selectivity_free_energy_20261021'
D=ROOT/'diagnostics/lanm_selectivity_free_energy_20261021'

def pin(p):
    p=Path(p).resolve();return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def load(p):
    assert pin(p['path'])==p
    return json.loads(Path(p['path']).read_text())
def module(p,name):
    s=importlib.util.spec_from_file_location(name,p);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m

def main():
    start,cpu=time.monotonic(),time.process_time()
    assert all(os.environ.get(k)=='1' for k in ('OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'))
    objective=module(HERE/'OBJECTIVE.py','joint_objective')
    helperpin=pin(D/'explicit_model_v1/matched_dry_v1/OBJECTIVE.py')
    assert helperpin['sha256']=='bf65a20a379811f8d35e2cb5ae78365a5d34022691b332bae5df223196078819'
    body=module(helperpin['path'],'original_body')
    actualpin=pin(W/'aqueous_qualification_design_v1/coordination_screen_v1/run_v1/execution/RESULT.json');actual=load(actualpin)
    assert actual['job_id']=='17670' and actual['status']=='complete'
    bound=load(actual['bound']);ccfg=load(bound['candidate_config']);transfer=load(ccfg['transfer_result']);tb=load(transfer['bound']);olddata=load(tb['data'])
    data={'rows':{g:copy.deepcopy(olddata['rows'][g] if not g.startswith('aquo9_') else bound['rows'][g]) for g in objective.CELLS}}
    baseline={g:copy.deepcopy(transfer['rows'][g]['model'] if not g.startswith('aquo9_') else actual['rows'][g]['model']) for g in objective.CELLS}
    pins={'actual17670':actualpin,'coordination_bound':actual['bound'],'transfer17612':ccfg['transfer_result'],'transfer_data':tb['data'],'body_helper':helperpin,'root_declaration':pin(HERE/'ROOT_DECLARATION.md')}
    checks=[]
    def check(ok,name):
        checks.append({'name':name,'passed':bool(ok)})
        if not ok:raise AssertionError(name)
    def close(a,b):return abs(a-b)<=1e-11*max(1,abs(a),abs(b))
    def reject(d,v,name):
        try:objective.evaluate(d,v,body)
        except (ValueError,KeyError,TypeError,IndexError):ok=True
        else:ok=False
        check(ok,name)
    score=objective.evaluate(data,baseline,body)
    check(len(score['weighted_residuals'])==413,'Actual exact413 residual denominator')
    check(score['declared_energy_controls']==score['available_energy_controls']==10 and score['declared_force_controls']==score['available_force_controls']==15,'Actual ten energy and fifteen force controls')
    check(score['development_partitions']['old_eight']['passed'] and not score['development_partitions']['added_W9']['passed'] and not score['development_checks_passed'],'Actual old8 pass and W9 fail remain separate')
    check(close(score['energy_errors_kcal_mol']['balanced_exchange']['aquo9_central'],actual['comparison']['energy']['error_kcal_mol']),'Actual W9 energy error exactly retained')
    for kind in ('La','Dy','paired'):
        z=score['per_fixture'][kind]['aquo9_central'];a=actual['comparison']['rigid'][kind]
        check(z['passed']==a['passed'] and close(z['rigid_RMS_kcal_mol_A'],a['RMS_kcal_mol_A']) and close(z['maximum_body_component_error_kcal_mol_A'],a['maximum_component_kcal_mol_A']),'Actual W9 '+kind+' control replay')
    references={g:{'energy':r['reference']['energy'],'forces':copy.deepcopy(r['reference']['forces'])} for g,r in data['rows'].items()}
    selfscore=objective.evaluate(data,references,body)
    check(selfscore['objective']==0 and selfscore['development_checks_passed'],'Actual reference self-comparison gives zero loss and all controls pass')
    zero=copy.deepcopy(data)
    for r in zero['rows'].values():r['reference']={'energy':0.,'forces':[[0.,0.,0.] for _ in range(r['geometry']['atoms'])]}
    zeros={g:copy.deepcopy(r['reference']) for g,r in zero['rows'].items()}
    synthetic=copy.deepcopy(zeros)
    for qi,q in enumerate(objective.GEOMETRIES):
        for metal in ('La','Dy'):
            gid=q+'_'+metal
            synthetic[gid]['energy']=(qi+1)*.13 if metal=='La' else (qi+1)**2*.17
            scale=1 if metal=='La' else -.7
            synthetic[gid]['forces']=[[scale*.003*(qi+1)*(i+1)*(k+1) for k in range(3)] for i in range(zero['rows'][gid]['geometry']['atoms'])]
    # Independent manual reconstruction: direct atom sums, cross products and
    # explicit scalar weights. No body response/projection/objective helpers.
    expected=[]
    def e(q,metal):return synthetic[q+'_'+metal]['energy']
    for q in ('water_minus','water_central','water_plus','aquo9_central'):
        expected.append(((e(q,'Dy')-e(q,'La'))-(e('aquo8_central','Dy')-e('aquo8_central','La')))*math.sqrt(.25*.5/4))
    for metal in ('La','Dy'):
        for sign in ('minus','plus'):expected.append((e('water_'+sign,metal)-e('water_central',metal))*math.sqrt(.25*.5/4))
    for kind in ('La','Dy','paired'):
        for q in objective.GEOMETRIES:
            row=zero['rows'][q+'_La'];xyz=row['geometry']['positions_A'];count=len(xyz);nw=(count-1)//3
            assert row['geometry']['metal_index']==0
            f=synthetic[q+'_'+kind]['forces'] if kind!='paired' else [[d-l for d,l in zip(a,b)] for a,b in zip(synthetic[q+'_Dy']['forces'],synthetic[q+'_La']['forces'])]
            expected.extend(v*math.sqrt(.25/(5*2*3))/2 for v in f[0])
            for first in range(1,count,3):
                ids=range(first,first+3);center=[sum(xyz[i][k] for i in ids)/3 for k in range(3)]
                net=[sum(f[i][k] for i in ids) for k in range(3)];torque=[0.,0.,0.]
                for i in ids:
                    x,y,z=[xyz[i][k]-center[k] for k in range(3)];fx,fy,fz=f[i]
                    t=[y*fz-z*fy,z*fx-x*fz,x*fy-y*fx]
                    torque=[a+b for a,b in zip(torque,t)]
                expected.extend(v*math.sqrt(.25/(5*2*nw*6))/2 for v in net+torque)
    test=objective.evaluate(zero,synthetic,body)
    check(len(expected)==len(test['weighted_residuals'])==413,'Manual projection has413 entries')
    check(all(close(a,b) for a,b in zip(expected,test['weighted_residuals'])),'All413 residual signs, body order, fixture order and numerical weights match independent manual construction')
    check(close(sum(v*v for v in expected),test['objective']),'Independent sum of squared residuals equals objective')
    for name,lo,hi in (('energy',0,8),('La',8,143),('Dy',143,278),('paired',278,413)):
        check(close(sum(v*v for v in expected[lo:hi]),.25*test['blocks'][name]),'Independent '+name+' block scaling and boundaries')
    shifted=copy.deepcopy(zeros)
    for q in objective.GEOMETRIES:
        for metal,offset in (('La',123.),('Dy',-71.)):shifted[q+'_'+metal]['energy']+=offset
    check(objective.evaluate(zero,shifted,body)['objective']==0,'Per-element constant offsets cancel in every balanced energy and radial work')
    forcebad=copy.deepcopy(zeros)
    for metal in ('La','Dy'):forcebad['aquo9_central_'+metal]['forces'][0][0]=5.2
    badscore=objective.evaluate(zero,forcebad,body)
    check(all(badscore['controls'][k]['category_balanced_rigid_RMS_kcal_mol_A']<2 for k in ('La','Dy')),'Synthetic W9 failure is below averaged RMS gate')
    check(not badscore['development_checks_passed'] and not badscore['per_fixture']['La']['aquo9_central']['passed'] and not badscore['per_fixture']['Dy']['aquo9_central']['passed'],'W9 individual failures cannot hide behind means')
    check(badscore['per_fixture']['paired']['aquo9_central']['passed'] and badscore['development_partitions']['old_eight']['passed'],'Shared force error cancels only in paired control; old8 unchanged')
    energybad=copy.deepcopy(zeros);energybad['aquo9_central_Dy']['energy']=1.1;badscore=objective.evaluate(zero,energybad,body)
    check(not badscore['energy_checks_passed'] and badscore['development_partitions']['old_eight']['passed'] and not badscore['development_checks_passed'],'W9 energy failure cannot hide in energy block mean')
    pairedonly=copy.deepcopy(zeros)
    for q in ('water_central','water_plus'):pairedonly[q+'_La']['energy']=.3;pairedonly[q+'_Dy']['energy']=-.3
    pairedonly['water_minus_La']['energy']=-.3;pairedonly['water_minus_Dy']['energy']=.3
    badscore=objective.evaluate(zero,pairedonly,body)
    check(all(abs(v)<=1 for name,vals in badscore['energy_errors_kcal_mol'].items() if name!='paired_radial_work' for v in vals.values()) and not badscore['energy_checks_passed'],'Paired radial acceptance check independently prevents pass')
    check(close(badscore['energy_errors_kcal_mol']['paired_radial_work']['minus'],1.2),'Paired radial work sign retained outside loss')
    for gid in objective.CELLS:
        missing=copy.deepcopy(zeros);missing.pop(gid);reject(zero,missing,'Missing model row rejected: '+gid)
    missing=copy.deepcopy(zero);missing['rows'].pop('aquo9_central_Dy');reject(missing,zeros,'Missing reference metadata rejected')
    for origin in ('model','reference'):
        for fault in ('energy_nan','energy_inf','energy_bool','force_nan','force_inf','short_force','missing_value'):
            d=copy.deepcopy(zero);v=copy.deepcopy(zeros);target=v['aquo9_central_La'] if origin=='model' else d['rows']['aquo9_central_La']['reference']
            if fault.startswith('energy_'):target['energy']={'energy_nan':float('nan'),'energy_inf':float('inf'),'energy_bool':True}[fault]
            elif fault=='force_nan':target['forces'][0][0]=float('nan')
            elif fault=='force_inf':target['forces'][0][0]=float('inf')
            elif fault=='short_force':target['forces'].pop()
            elif origin=='model':v['aquo9_central_La']=None
            else:d['rows']['aquo9_central_La']['reference']=None
            reject(d,v,origin+' '+fault+' fails closed')
    bad=copy.deepcopy(zero);bad['rows']['aquo9_central_Dy']['geometry']['positions_A'][1][0]+=.01;reject(bad,zeros,'Unmatched paired coordinates rejected')
    bad=copy.deepcopy(zero);bad['rows']['aquo9_central_La']['bodies'].append(copy.deepcopy(bad['rows']['aquo9_central_La']['bodies'][-1]));reject(bad,zeros,'Duplicated body cannot alter denominator')
    bad=copy.deepcopy(zero);bad['rows']['aquo9_central_La']['geometry']['charge']=2;reject(bad,zeros,'Changed water charge rejected')
    check(all(pin(p['path'])==p for p in pins.values()),'All saved input pins unchanged')
    out=W/'explicit_model_v1/paired_water_coordination_joint_v1/objective_controls_v1';out.mkdir(parents=True,exist_ok=False)
    result={'status':'passed','objective_source':pin(HERE/'OBJECTIVE.py'),'test_source':pin(__file__),'input_pins':pins,'checks':checks,'checks_passed':len(checks),'actual_baseline_score':score,'manual_residual_maximum_difference':max(abs(a-b) for a,b in zip(expected,test['weighted_residuals'])),'new_native_calls':0,'reserved_outcomes_opened':False,'local_threads':1,'wall_seconds':time.monotonic()-start,'CPU_seconds':time.process_time()-cpu}
    (out/'RESULT.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({'result':pin(out/'RESULT.json'),'checks':len(checks),'baseline_objective':score['objective'],'manual_maximum_difference':result['manual_residual_maximum_difference'],'wall_seconds':result['wall_seconds'],'CPU_seconds':result['CPU_seconds']}))

if __name__=='__main__':main()
