#!/usr/bin/env python3
"""Pure arithmetic and missing-slot controls; no future result or native call."""
import argparse
import copy
import importlib.util
import json
from pathlib import Path
import time
HERE=Path(__file__).resolve().parent
s=importlib.util.spec_from_file_location('common_arithmetic_controls',HERE/'COLLECT.py');m=importlib.util.module_from_spec(s);s.loader.exec_module(m)
a=m.a;b=m.b

def main(panel,out):
    tick,cpu=time.monotonic(),time.process_time();p=a.load(panel);checks={}
    def check(name,ok):b.check(ok,name);checks[name]=True
    def reject(name,fn):
        try:fn()
        except (AssertionError,ValueError,KeyError):checks[name]=True;return
        raise AssertionError('Expected rejection: '+name)
    row=lambda x,d4=0.:dict(status='complete',electronic_energy_hartree=x,D4_energy_hartree=d4)
    new={};reused={};manual={};complex_d4={}
    for j,g in enumerate(a.GIDS):
        n=p['geometries'][g]['water_count'];ion=-1.-j*.13;shell=-10.-j*.17;E=-21.-j*.31;d4=.00001*(j+1)
        reused[g]=row(E,d4);reused[m.reused_id(g,'ion_ghost_shell')]=row(ion);reused[m.reused_id(g,'shell_ghost_ion')]=row(shell);complex_d4[g]=d4
        pairs=[]
        for i in range(1,n+1):
            d=-3.-j*.21-i*.007;w=-1.5-j*.11-i*.009
            new[m.cid(g,i,'metal_water')]=row(d);new[m.cid(g,i,'water_only')]=row(w);pairs.append(d-ion-w)
        C=E-ion-shell;P=sum(pairs);manual[g]={'C':C,'P':P,'B':C-P}
    signs=[1,-1,-1,1];manual_delta={q:sum(sign*manual[g][q] for sign,g in zip(signs,a.GIDS)) for q in ('C','P','B')}
    d4C=sum(sign*complex_d4[g] for sign,g in zip(signs,a.GIDS));d4P=.003;d4={'cluster_contrast':d4C,'pair_contrast':d4P,'cooperative_remainder':d4C-d4P}
    cp={'hartree':{'CP_electronic':manual_delta['C'],'original_D4_unchanged':d4C}}
    pair={'reference_decomposition':{'D4_energy_hartree':{'hartree':d4},'final_energy_hartree':{'kcal_mol':{k:0. for k in m.QUANTITIES}}},'model_decomposition':{k:0. for k in m.QUANTITIES}}
    result=m.calculate(p,new,reused,cp,pair)
    for k,q in zip(m.QUANTITIES,('C','P','B')):
        check('manual_'+q,abs(result['electronic_hartree'][k]-manual_delta[q])<1e-12)
        check('D4_once_'+q,abs(result['electronic_plus_original_D4_hartree'][k]-manual_delta[q]-d4[k])<1e-12)
    g=a.GIDS[3];n=9;changed=copy.deepcopy(reused);changed[m.reused_id(g,'ion_ghost_shell')]['electronic_energy_hartree']+=.125
    changed_cp=copy.deepcopy(cp);changed_cp['hartree']['CP_electronic']-=.125
    r=m.calculate(p,new,changed,changed_cp,pair)
    check('ion_N_fold_in_pair',abs(r['electronic_hartree']['pair_contrast']-result['electronic_hartree']['pair_contrast']+n*.125)<1e-12)
    check('ion_N_minus1_in_B',abs(r['electronic_hartree']['cooperative_remainder']-result['electronic_hartree']['cooperative_remainder']-(n-1)*.125)<1e-12)
    changed=copy.deepcopy(new);changed[m.cid(g,1,'water_only')]['electronic_energy_hartree']+=.25;r=m.calculate(p,changed,reused,cp,pair)
    check('individual_ghost_water_not_cancelled',abs(r['electronic_hartree']['cooperative_remainder']-result['electronic_hartree']['cooperative_remainder']-.25)<1e-12)
    missing=dict(new);missing.pop(next(iter(missing)));reject('missing_one_of68',lambda:m.calculate(p,missing,reused,cp,pair))
    missing=dict(reused);missing.pop(next(iter(missing)));reject('missing_one_of12',lambda:m.calculate(p,new,missing,cp,pair))
    for label,field,value in [('failed','status','failed_admission'),('nonfinite','electronic_energy_hartree',float('nan')),('ghost_D4','D4_energy_hartree',.1)]:
        bad=copy.deepcopy(new);bad[next(iter(bad))][field]=value;reject('reject_'+label,lambda:m.calculate(p,bad,reused,cp,pair))
    bad=copy.deepcopy(cp);bad['hartree']['CP_electronic']+=.01;reject('wrong_prior_CP8_C',lambda:m.calculate(p,new,reused,bad,pair))
    bad=copy.deepcopy(pair);bad['reference_decomposition']['D4_energy_hartree']['hartree']['pair_contrast']+=.01;reject('D4_closure_failure',lambda:m.calculate(p,new,reused,cp,bad))
    out.mkdir(parents=True,exist_ok=False)
    review=b.write(out/'SYNTHETIC_REVIEW.json',dict(status='passed_collection_review',implementation=b.pin(HERE/'COLLECT.py'),panel=b.pin(panel),synthetic_test_only=True))
    paths={c['cell_id']:str(out/'deliberately_missing'/c['cell_id']/'RESULT.json') for c in p['cells']}
    bad_id=p['cells'][0]['cell_id'];badpath=out/'truncated_RESULT.json';badpath.write_text('{');paths[bad_id]=str(badpath)
    cfg=dict(implementation=b.pin(HERE/'COLLECT.py'),adapter=p['adapter'],panel=b.pin(panel),independent_review=review,result_paths=paths)
    config=b.write(out/'SYNTHETIC_CONFIG.json',cfg)
    saved=m.admit_reused
    m.admit_reused=lambda key,*args:dict(p['saved_reuse'][key],synthetic_saved_control_passthrough=True)
    try:missing_result=m.collect(config['path'],out/'missing_and_malformed')
    finally:m.admit_reused=saved
    check('missing_and_malformed_keeps68slots',len(missing_result['new_fragments'])==68 and missing_result['new_complete']==0)
    check('reused12_denominator_preserved',len(missing_result['reused_terms'])==missing_result['reused_complete']==12)
    check('no_partial_decomposition',missing_result['decomposition'] is None and missing_result['status']=='common_basis_collection_incomplete')
    check('malformed_row_retains_error',missing_result['new_fragments'][bad_id]['status']=='failed_or_unavailable' and 'JSONDecodeError' in missing_result['new_fragments'][bad_id]['error'])
    check('other67rows_remain_missing',sum(v['status']=='missing' for v in missing_result['new_fragments'].values())==67)
    r=dict(status='saved_collection_controls_passed',implementation=b.pin(__file__),collector=b.pin(HERE/'COLLECT.py'),panel=b.pin(panel),checks=checks,passed=len(checks),wall_seconds=time.monotonic()-tick,parent_cpu_seconds=time.process_time()-cpu,local_threads=1,new_native_calls=0,synthetic_only=True,future_native_outputs_read=False)
    b.write(out/'RESULT.json',r);print(json.dumps({'status':r['status'],'passed':len(checks),'new_native_calls':0}))
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--panel',required=True);ap.add_argument('--out',type=Path,required=True);args=ap.parse_args();main(args.panel,args.out)
