#!/usr/bin/env python3
"""Artificial mapping/source checks only; solver construction is forbidden."""
import ast
import copy
import importlib.util
import json
from pathlib import Path
import sys
import time
import numpy as np

HERE=Path(__file__).resolve().parent


def load(path,name):
    spec=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module


def main(out):
    wall=time.monotonic();cpu=time.process_time();checks=[]
    new=load(HERE/'DIAGNOSE.py','joint_ceiling_diagnostic')
    old=load(HERE.parent/'conic_saved_comparison_v1/TEST_SOURCE.py','parent_artificial_controls')
    def forbidden(*args,**kwargs):raise AssertionError('No solver construction in source controls')
    new.clarabel.DefaultSolver=forbidden;old.C.clarabel.DefaultSolver=forbidden
    def check(name,ok):
        if not ok:raise AssertionError(name)
        checks.append(name)
    def rejects(name,call):
        try:call()
        except (ValueError,TypeError):checks.append(name);return
        raise AssertionError('Unexpected pass: '+name)
    p=old.artificial();full=old.C.assemble(p);model,removed,keep=new.remove_joint_ceilings(full)
    check('only_exact535_rows_removed',keep==list(range(535,1817)) and len(removed)==4)
    check('retained_A_b_literal_slices',np.array_equal(model['A'].toarray(),full['A'].toarray()[535:]) and np.array_equal(model['b'],full['b'][535:]))
    check('objective_and_P_unchanged',np.array_equal(model['c'],full['c']) and (model['P']!=full['P']).nnz==0)
    check('17SOC_plus856_sides',len(model['cones'])==18 and model['A'].shape==(1282,13))
    check('both_component_cones_retained',all(any(r['name']=='block/'+name for r in model['layout']) for name in ('P_energy','B_energy')))
    check('15_force_RMS_cones_retained',sum(r['name'].startswith('force/RMS/') for r in model['layout'])==15)
    check('all_original_nonnegative_layout_preserved',
          [r['name'] for r in model['layout'] if r['kind']=='nonnegative']==
          [r['name'] for r in full['layout'] if r['kind']=='nonnegative'])
    for i,(step,t) in enumerate(((0.,1.),(.01,.75),(-.02,.5),(.025,.25))):
        x=np.r_[[step*(-1 if j%2 else 1) for j in range(12)],t]
        got=new.correspondence(p,model,x,old.core)
        check('847_retained_constraint_probe_'+str(i),got['retained_source_inequalities']==847)
    x=np.r_[np.zeros(12),1.]
    rows=new.reports(p,x,'Solved',old.core)
    check('only_inherited4fractions',tuple(r['fraction'] for r in rows)==old.core.FRACTIONS)
    check('all25_and_all6blocks_reported',all(len(r['all25_original_gates'])==25 and set(r['all_six_blocks'])==set(old.C.BLOCKS) for r in rows))
    check('no_candidate_selection_authorized',all(r['selection_authorized'] is False and r['native_evaluated'] is False for r in rows))
    for name,b in p['blocks'].items():
        factor=1 if name in ('P_energy','B_energy') else 4
        check('score_scale_'+name,abs(rows[0]['all_six_blocks'][name]['origin_score']-factor*float(b['r']@b['r']))<=1e-14)
    q=copy.deepcopy(p);q['blocks']['old_EF']['r']*=1000
    rr=new.reports(q,x,'Solved',old.core)[0]
    check('removed_block_violation_visible_not_selectable',
          rr['retained_predicted_feasible'] and not rr['source_certificate_under_original_policy']['eligible']
          and rr['removed_source_margins']['block/old_EF']<0 and not rr['selection_authorized'])
    q=copy.deepcopy(p);q['energy']['r'][0]=2
    check('energy_gate_still_required',not new.reports(q,x,'Solved',old.core)[0]['retained_predicted_feasible'])
    q=copy.deepcopy(p);next(iter(q['forces'].values()))['r'][0]=10
    check('force_gate_still_required',not new.reports(q,x,'Solved',old.core)[0]['retained_predicted_feasible'])
    xx=x.copy();xx[0]=.03
    check('trust_box_still_required',not new.reports(p,xx,'Solved',old.core)[0]['retained_predicted_feasible'])
    q=dict(full);q['layout']=copy.deepcopy(full['layout']);q['layout'][0]['stop']-=1
    rejects('changed_removed_cone_rejected',lambda:new.remove_joint_ceilings(q))
    q=dict(full);q['layout']=copy.deepcopy(full['layout']);q['layout'][0]['name']='block/unapproved'
    rejects('missing_named_ceiling_rejected',lambda:new.remove_joint_ceilings(q))
    rejects('negative_t_rejected',lambda:new.reports(p,np.r_[np.zeros(12),-.1],'Solved',old.core))
    tree=ast.parse((HERE/'DIAGNOSE.py').read_text())
    check('one_solve_site_only',sum(isinstance(n,ast.Call) and isinstance(n.func,ast.Attribute) and n.func.attr=='solve' for n in ast.walk(tree))==1)
    source=(HERE/'DIAGNOSE.py').read_text()
    check('baseline_success_actual_review_gates',all(t in source for t in (
        "review['status'] == 'passed_saved_conic_actual_review'",
        "review['certificate_checks_passed'] is True",
        "review['solver_success'] is True", "review['solver_status'] == 'Solved'")))
    check('persistent_exclusive_attempt_marker',"open(cfg['one_attempt_marker'],'x')" in source)
    report={'status':'passed_joint_ceiling_removal_artificial_source_controls',
            'implementation':new.pin(HERE/'DIAGNOSE.py'),'test_implementation':new.pin(__file__),
            'parent_artificial_source':new.pin(HERE.parent/'conic_saved_comparison_v1/TEST_SOURCE.py'),
            'checks':checks,'check_count':len(checks),'solver_calls':0,'native_calls':0,
            'scientific_vectors_read':0,'campaign_outcomes_read':False,
            'cost':{'wall_seconds':time.monotonic()-wall,'CPU_seconds':time.process_time()-cpu,'threads':1}}
    out.mkdir(parents=True,exist_ok=False);(out/'RESULT.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'result':new.pin(out/'RESULT.json'),'checks':len(checks),'cost':report['cost']}))


if __name__=='__main__':main(Path(sys.argv[1]))
