"""Finite source/artificial controls; no campaign/native/numerical helper import."""
import ast
import copy
import importlib.util
import json
import os
from pathlib import Path
import time
from types import SimpleNamespace

HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('new_consumed_regression',HERE/'RUN.py')
R=importlib.util.module_from_spec(spec);spec.loader.exec_module(R)

def run(out):
    tick,cpu=time.perf_counter(),time.process_time();checks=[]
    def check(name,x):
        assert x,name;checks.append(name)
    def reject(name,f):
        try:f()
        except ValueError:checks.append(name);return
        raise AssertionError(name)
    check('one numerical thread',all(os.environ.get(k)=='1' for k in (*R.THREADS,'BLIS_NUM_THREADS')))
    inputs=R.read(R.pin(HERE/'SOURCE_INPUTS.json'));parent=R.verify(inputs['parent_adapter'])
    def node(path,name):return next(n for n in ast.walk(ast.parse(path.read_text())) if isinstance(n,ast.FunctionDef) and n.name==name)
    for name in ('cycles','login_contract'):
        check('unchanged '+name,ast.dump(node(HERE/'RUN.py',name))==ast.dump(node(parent,name)))
    def native_body(path):
        n=node(path,'run');loop=next(x for x in ast.walk(n) if isinstance(x,ast.For) and ast.unparse(x.target)=='cid' and ast.unparse(x.iter)=='IDS')
        start=next(i for i,x in enumerate(loop.body) if ast.unparse(x).startswith("folder = out / 'cells' / cid"))
        return [ast.dump(x) for x in loop.body[start:]]
    check('unchanged render/invoke/parser/cap/failure body',native_body(HERE/'RUN.py')==native_body(parent))
    for status in ('started','failed_physical_gate_component_fit','bounded_component_fit_complete'):
        reject('reject '+status,lambda status=status:R.candidate_contract({}, {},dict(status=status,mode='run',candidate_frozen=True,error=None),{},{},{}))
    ev=SimpleNamespace(effective=lambda x,m,d:x[0 if m=='La' else 1])
    hist=[({'result':'old_component'},dict(parameters=[1.,2.])),({'result':'old_joint'},dict(parameters=[3.,4.]))]
    check('old component exact metal reuse',R.matching_history([1.,9.],'La',{},ev,hist)[0]['result']=='old_component')
    check('old joint exact metal reuse',R.matching_history([9.,4.],'Dy',{},ev,hist)[0]['result']=='old_joint')
    check('nearby state cannot reuse',R.matching_history([1.+1e-12,9.],'La',{},ev,hist) is None)
    rows={};previous={}
    for n in (8,9):
        for which in ('near','far'):
            for m,delta in [('La',0),('Dy',3)]:
                cid=f'aquo{n}_central_{m}__{which}'
                rows[cid]=dict(status='complete',reference_electronic_E3_kcal_mol=10+delta,reference_mu3_Debye=[delta,0,0],baseline_E3_kcal_mol=5+delta,baseline_mu3_Debye=[delta,1,0],candidate_E3_kcal_mol=11+2*delta,candidate_mu3_Debye=[delta,2,0])
                previous[cid]=dict(candidate_E3_kcal_mol=8+delta,candidate_mu3_Debye=[delta,3,0])
    comp=R.comparisons(rows,previous);pair=comp['paired_Dy_minus_La']['aquo8__near']
    check('fixed eight/four denominator',len(comp['cycles'])==8 and len(comp['paired_Dy_minus_La'])==4)
    check('Dy-minus-La energy error sign',pair['models']['new_terminal']['Dy_minus_La_E3_kcal_mol']==6 and pair['models']['new_terminal']['energy_error_kcal_mol']==3)
    check('paired moment is vector difference',pair['models']['new_terminal']['dipole_error_norm_Debye']==0 and comp['cycles']['aquo8_central_La__near']['models']['new_terminal']['dipole_error_norm_Debye']==2)
    check('old terminal and baseline stay separate',comp['cycles']['aquo8_central_La__near']['models']['old_component_terminal']['energy_error_kcal_mol']==-2 and comp['cycles']['aquo8_central_La__near']['models']['original_pre_response_baseline']['energy_error_kcal_mol']==-5)
    missing=copy.deepcopy(rows);missing['aquo9_central_Dy__far']['status']='unavailable';comp=R.comparisons(missing,previous)
    check('failure retains eight/four and explicit missing pair',len(comp['cycles'])==8 and len(comp['paired_Dy_minus_La'])==4 and comp['paired_Dy_minus_La']['aquo9__far']['status']=='unavailable' and 'new_terminal' not in comp['paired_Dy_minus_La']['aquo9__far']['models'])
    out=Path(out);out.mkdir(parents=True,exist_ok=False)
    result=dict(status='passed_physical_gate_postselection_source_controls',implementation=R.pin(HERE/'RUN.py'),test_source=R.pin(__file__),parent_source=inputs['parent_adapter'],checks=checks,checks_passed=len(checks),native_calls=0,solver_calls=0,real_candidate_bound=False,cost=dict(wall_seconds=time.perf_counter()-tick,CPU_seconds=time.process_time()-cpu,numerical_threads=1))
    R.save(out/'RESULT.json',result);print(json.dumps(result))

if __name__=='__main__':
    import sys
    run(sys.argv[1])
