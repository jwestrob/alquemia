#!/usr/bin/env python3
"""Narrow saved-format, literal inventory and finite-release controls; zero native calls."""
import argparse
import copy
import importlib.util
import json
from pathlib import Path
import time
HERE=Path(__file__).resolve().parent
s=importlib.util.spec_from_file_location('common_controls',HERE/'REFERENCE.py');a=importlib.util.module_from_spec(s);s.loader.exec_module(a)
def main(panel,out):
    tick,cpu=time.monotonic(),time.process_time();p=a.load(panel);checks={}
    def check(name,ok):a.b.check(ok,name);checks[name]=True
    def reject(name,fn):
        try:fn()
        except (AssertionError,ValueError,KeyError,FileNotFoundError):checks[name]=True;return
        raise AssertionError('Expected rejection: '+name)
    check('68_literal_full_parent_fragments',len(p['cells'])==68)
    check('34_each_real_role',all(sum(c['role']==r for c in p['cells'])==34 for r in a.ROLES))
    check('12_reused_raw_authorities',len(p['saved_reuse'])==12)
    check('all_parent_auxiliary_counts_preserved',all(v['native_aux_J_basis_function_counts']['basis_header']==[v['expected_aux_J_basis_functions']] and v['native_aux_J_basis_function_counts']['integral_header']==[v['expected_aux_J_basis_functions']] for v in p['saved_reuse'].values()))
    old=a.helper.load(a.b.verify(p['cp_panel']))
    for role in a.helper.ROLES:
        c=copy.deepcopy(next(c for c in old['cells'] if c['cell_id']==f'pbe0__aquo9__Dy__{role}'));c['expected_aux_J_basis_functions']=1008
        raw=a.b.verify(p['saved_reuse'][c['cell_id']]['output']).read_text();parsed=a.parse(raw,c,0)
        check('actual_'+role+'_format',parsed['status']=='complete' and parsed['checks']['full_parent_aux_J_basis'])
        check('actual_'+role+'_raw_false_flags_preserved',parsed['printed_scf_criteria']==a.helper.parse(raw,c,0)['printed_scf_criteria'])
        check('reject_wrong_aux_'+role,a.parse(raw,dict(c,expected_aux_J_basis_functions=1007),0)['status']=='failed_admission')
    text='# of basis functions in Aux-J ... 1008\nNumber of AuxJ basis functions .... 1008\n'
    check('both_aux_labels',a.auxiliary_admission(text,1008)[1])
    check('reject_missing_aux_label',not a.auxiliary_admission(text.splitlines()[0],1008)[1])
    check('reject_mixed_aux_counts',not a.auxiliary_admission(text+'Number of AuxJ basis functions .... 1007\n',1008)[1])
    pre=a.b.pin(Path(panel).parent/'PREFLIGHT.json')
    release=dict(status='approved_common_basis_decomposition_cells',implementation=p['adapter'],panel=a.b.pin(panel),cell_ids=[c['cell_id'] for c in p['cells']],maximum_native_calls=68,ranks=4,preparation=pre,independent_review=p['declaration'],stage='scouts',released_cell_ids=a.SCOUTS,execution_root='/unexecuted/saved_control/run_v1')
    check('exact_two_scout_subset',a.release_contract(p,panel,release)==release)
    reject('reject_larger_scout_subset',lambda:a.release_contract(p,panel,dict(release,released_cell_ids=release['cell_ids'])))
    reject('reject_repeat_scout',lambda:a.release_contract(p,panel,dict(release,released_cell_ids=[a.SCOUTS[0]]*2)))
    reject('reject_altered_total_cap',lambda:a.release_contract(p,panel,dict(release,maximum_native_calls=69)))
    reject('reject_altered_ranks',lambda:a.release_contract(p,panel,dict(release,ranks=8)))
    reject('remaining_requires_actual_scout_review',lambda:a.release_contract(p,panel,dict(release,stage='remaining',released_cell_ids=[x for x in release['cell_ids'] if x not in a.SCOUTS])))
    out.mkdir(parents=True,exist_ok=False)
    result=dict(status='saved_controls_passed',implementation=a.b.pin(__file__),adapter=a.b.pin(HERE/'REFERENCE.py'),panel=a.b.pin(panel),checks=checks,passed=len(checks),new_native_calls=0,local_threads=1,wall_seconds=time.monotonic()-tick,parent_cpu_seconds=time.process_time()-cpu,synthetic_release_never_executed=True)
    a.b.write(out/'RESULT.json',result);print(json.dumps({'status':result['status'],'passed':len(checks),'new_native_calls':0}))
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--panel',required=True);ap.add_argument('--out',required=True,type=Path);v=ap.parse_args();main(v.panel,v.out)
