#!/usr/bin/env python3
"""Bounded saved controls; no subprocess/native call or new molecular model."""
import argparse
import copy
import importlib.util
import json
import os
from pathlib import Path
import time

HERE=Path(__file__).resolve().parent
s=importlib.util.spec_from_file_location('sensitivity_author_control',HERE/'RUN.py')
r=importlib.util.module_from_spec(s);s.loader.exec_module(r)

def main(out,prepared):
    tick=time.monotonic();cpu=time.process_time();checks=[]
    def check(ok,name):
        r.require(ok,name);checks.append({'name':name,'passed':True})
    def reject(fn,name):
        try:fn()
        except (ValueError,KeyError,TypeError):check(True,name)
        else:raise AssertionError(name)
    check(all(os.environ.get(k)=='1' for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS')),'One local numerical thread')
    cfg=json.loads((prepared/'CONFIG_RUN_TEMPLATE.json').read_text());m,h=r.context(cfg)
    check(len(m['rows'])==30 and len(m['jobs'])==120,'Exact30×4 cells')
    check(len(m['reused_gradients'])==4 and m['new_analyze_programs']==120 and m['new_gradient_programs']==116,'Four reused gradients and236 new programs')
    check(sum(v is not None for v in m['baseline_forces'].values())==16,'16 known baseline force arrays')
    check(sum(v is None for v in m['baseline_forces'].values())==14,'14 unavailable baseline arrays remain null')
    check(len(m['trimer_cycles'])==8 and len(m['center_reference_C'])==4 and len(m['unchanged_moment_ids'])==26,'Full8/4 outcomes and26 unchanged fragments')
    for label,job in m['jobs'].items():
        base=m['rows'][job['cell']]
        text=r.verify(job['input']['key']).read_text()
        check(text==r.profile_key(r.verify(base['input']['key']).read_text(),job['profile']) and text.lower().count('polar-pair-response')==2,'Exact two-record-only key edit '+label)
    records={}
    for label,job in m['jobs'].items():
        cid=job['cell'];moment=m['baseline_moments'].get(cid)
        if moment is None:moment={'energy_kcal_mol':m['baseline_model'][cid]['energy']}
        force=m['baseline_forces'][cid]
        records[label]={'status':'complete','moments':copy.deepcopy(moment),'gradient':{'forces':force['forces'] if force else None}}
    replay=r.summarize(m,records,h)
    for profile,value in replay.items():
        check(value['original25_controls']==m['baseline_score'],'Unchanged original25 controls '+profile)
        for cid,x in value['trimer_cycles'].items():
            old=m['trimer_reference'][cid]
            check(abs(x['model_E3_kcal_mol']-old['model_E3_kcal_mol'])<1e-11 and max(abs(a-b) for a,b in zip(x['model_mu3_Debye'],old['model_mu3_Debye']))<1e-12,'Exact saved seven-term energy/moment '+profile+' '+cid)
        check(not value['full_shell_B_computed'] and not value['selected_candidate'],'No B reconstruction or profile selection '+profile)
    missing=dict(records);missing.pop(next(iter(missing)));reject(lambda:r.summarize(m,missing,h),'Missing profile cell rejects')
    bad=copy.deepcopy(records);bad[next(iter(bad))]['gradient']['forces'][0][0]=float('nan');reject(lambda:r.summarize(m,bad,h),'Nonfinite original force control rejects')
    label,old=next(iter(m['reused_gradients'].items()));row=m['rows'][m['jobs'][label]['cell']]
    fake={'status':'complete','input':old['input'],'receipts':[copy.deepcopy(old['receipt'])],'value':old['value'],'job':{'id':old['producer_record']}}
    fake['receipts'][0]['argv'][0]='/wrong/binary'
    reject(lambda:r.reuse_gradient(m['sources'],fake,m['jobs'][label]['input'],row,m['jobs'][label]['profile'],h),'Wrong producer executable rejects before raw parse')
    changed=copy.deepcopy(m['jobs'][label]['input']);changed['key']=m['rows'][m['jobs'][label]['cell']]['input']['key']
    reject(lambda:r.reuse_gradient(m['sources'],fake,changed,row,m['jobs'][label]['profile'],h),'Different effective key rejects reuse')
    count=r.Counter()
    for _ in range(236):count.reserve()
    reject(count.reserve,'237th program reservation rejects')
    count=r.Counter();count.stop.set();reject(count.reserve,'Global stop prevents a new program')
    release={'status':'approved_consumed_pair_response_sensitivity','implementation':cfg['implementation'],'manifest':cfg['manifest'],'plan':cfg['plan'],'maximum_native_programs':236,'workers':4,'internal_FD_energies':0,'fit_authorized':False,'qualification_review':m['sources']['qualification_review']}
    r.save(out/'FAKE_RELEASE.json',release);bound=dict(cfg,root_release=r.pin(out/'FAKE_RELEASE.json'))
    check(r.release_gate(bound,m)==release,'Concrete one-way release/config binding passes')
    release['maximum_native_programs']=237;r.save(out/'BAD_RELEASE.json',release)
    reject(lambda:r.release_gate(dict(cfg,root_release=r.pin(out/'BAD_RELEASE.json')),m),'Wrong released ceiling rejects')
    result={'status':'passed_saved_sensitivity_controls','implementation':r.pin(__file__),'adapter':r.pin(HERE/'RUN.py'),'manifest':r.pin(prepared/'MANIFEST.json'),'checks':checks,'check_count':len(checks),'new_native_calls':0,'compiler_calls':0,'error':None,'cost':{'wall_seconds':time.monotonic()-tick,'CPU_seconds':time.process_time()-cpu,'local_threads':1}}
    r.save(out/'RESULT.json',result);print(json.dumps({'checks':len(checks),'result':r.pin(out/'RESULT.json')}))

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,required=True);ap.add_argument('--prepared',type=Path,required=True);a=ap.parse_args();a.out.mkdir(parents=True,exist_ok=False);main(a.out,a.prepared)
