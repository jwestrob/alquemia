#!/usr/bin/env python3
"""One finite arithmetic/missing/provenance test set; no native output fabricated."""
import argparse
import copy
import importlib.util
import json
import math
from pathlib import Path
import time

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--out',required=True);a=ap.parse_args();out=Path(a.out);out.mkdir(parents=True,exist_ok=False)
    t=time.monotonic();cpu=time.process_time();source=Path(__file__).with_name('COLLECT.py');s=importlib.util.spec_from_file_location('cp_controls',source);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);checks=[]
    def check(ok,name):m.require(ok,name);checks.append(name)
    def rejects(fn,name):
        try:fn()
        except (ValueError,KeyError,TypeError):checks.append(name);return
        raise AssertionError('Expected rejection: '+name)
    complexes={g:{'electronic_energy_hartree':e,'D4_energy_hartree':d,'final_energy_hartree':e+d} for g,e,d in zip(m.GIDS,[-100.,-110.,-120.,-131.],[-.1,-.2,-.3,-.5])}
    G=dict(zip(m.GIDS,[-99.,-108.,-118.,-129.]));fragments={}
    for g in m.GIDS:
        for j,role in enumerate(m.ROLES):fragments[f'pbe0__aquo{g[4]}__{g.rsplit("_",1)[1]}__{role}']={'status':'complete','electronic_energy_hartree':-5. if j==0 else G[g]+5.,'D4_energy_hartree':0.}
    model=dict(zip(m.GIDS,[10.,9.,8.,6.]));ok={'own_ion':True,'own_shell':True};r=m.calculate(complexes,fragments,model,ok)
    expected={'original_electronic':-1.,'original_D4_unchanged':-.1,'original_total':-1.1,'delta_G_electronic':-2.,'CP_shift':2.,'CP_electronic':1.,'CP_total':.9}
    for k,v in expected.items():check(math.isclose(r['hartree'][k],v,abs_tol=1e-12),'cycle sign '+k)
    check(math.isclose(r['kcal_mol']['model_minus_CP'],-1.-.9*m.HA,abs_tol=1e-9),'model minus CP sign and units')
    check(math.isclose(r['kcal_mol']['model_minus_original'],-1.+1.1*m.HA,abs_tol=1e-9),'model minus original sign and units')
    check(r['kcal_mol']['original_D4_unchanged']==r['hartree']['original_D4_unchanged']*m.HA and not r['fragment_D4_added'],'original D4 kept separately')
    shifted=copy.deepcopy(complexes)
    for g in m.GIDS:
        offset=(7. if g.endswith('_La') else -13.)+(11. if g.startswith('aquo8') else 29.)
        for k in ('electronic_energy_hartree','final_energy_hartree'):shifted[g][k]+=offset
    other=m.calculate(shifted,fragments,model,ok);check(all(math.isclose(other['hartree'][k],r['hartree'][k],abs_tol=1e-12) for k in expected),'own ion and shell constants cancel with correct signs')
    for kind in ('complex','fragment','model'):
        cs,fs,ms=copy.deepcopy(complexes),copy.deepcopy(fragments),copy.deepcopy(model)
        if kind=='complex':cs.pop(m.GIDS[0])
        elif kind=='fragment':fs.pop(m.IDS[0])
        else:ms.pop(m.GIDS[0])
        rejects(lambda cs=cs,fs=fs,ms=ms:m.calculate(cs,fs,ms,ok),'missing '+kind+' never produces cycle')
    for value in (float('nan'),float('inf'),None):
        fs=copy.deepcopy(fragments);fs[m.IDS[0]]['electronic_energy_hartree']=value
        rejects(lambda fs=fs:m.calculate(complexes,fs,model,ok),'invalid fragment energy '+repr(value))
    fs=copy.deepcopy(fragments);fs[m.IDS[0]]['D4_energy_hartree']=-.1;rejects(lambda:m.calculate(complexes,fs,model,ok),'ghost D4 forbidden')
    fs=copy.deepcopy(fragments);fs[m.IDS[0]]['status']='failed_admission';rejects(lambda:m.calculate(complexes,fs,model,ok),'failed fragment forbidden')
    cs=copy.deepcopy(complexes);cs[m.GIDS[0]]['final_energy_hartree']+=.01;rejects(lambda:m.calculate(cs,fragments,model,ok),'original electronic D4 closure required')
    rejects(lambda:m.calculate(complexes,fragments,model,{'own_ion':False}),'cancellation failure blocks cycle')
    rejects(lambda:m.calculate(complexes,fragments,model,{}),'missing cancellation evidence blocks cycle')
    check(m.unavailable(out/'not_created.json')['status']=='missing','missing result explicit')
    badfile=out/'malformed.json';badfile.write_text('{');bad=m.unavailable(badfile,ValueError('malformed'))
    check(bad['status']=='failed_or_unavailable' and bad['result']==m.pin(badfile),'malformed result retains pin')
    pinfile=out/'pin_fixture.txt';pinfile.write_text('original');original=m.pin(pinfile);pinfile.write_text('changed');rejects(lambda:m.verify(original),'changed source pin rejected')
    fixture={'status':'complete','returncode':0,'checks':{'fixture':True},'cell_id':m.IDS[0],'cell':{'cell_id':m.IDS[0]},'panel':{'path':'fixture','sha256':'a'},'adapter':{'path':'bad','sha256':'b'}}
    target=out/'wrong_provenance.json';m.save(target,fixture)
    rejects(lambda:m.admit_fragment(target,fixture['cell'],{'adapter':{'path':'correct','sha256':'c'}},fixture['panel'],None),'wrong native producer rejected before raw admission')
    m.save(out/'RESULT.json',{'status':'saved_controls_passed','implementation':m.pin(__file__),'collector':m.pin(source),'checks':checks,'checks_passed':len(checks),'new_native_calls':0,'local_threads':1,'wall_seconds':time.monotonic()-t,'CPU_seconds':time.process_time()-cpu})
    print(json.dumps({'status':'saved_controls_passed','checks':len(checks),'out':str(out)}))

if __name__=='__main__':main()
