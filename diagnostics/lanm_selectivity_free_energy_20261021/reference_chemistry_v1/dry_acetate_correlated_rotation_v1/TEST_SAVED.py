#!/usr/bin/env python3
"""Bounded saved-parser, numerical-denominator and release controls; zero native."""
import copy
from decimal import Decimal as D
import importlib.util
import json
from pathlib import Path
import sys
import time

HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('rotation_adapter',HERE/'REFERENCE.py')
r=importlib.util.module_from_spec(spec);spec.loader.exec_module(r)

def main():
    tick=time.monotonic();cpu=time.process_time();panel=Path(sys.argv[1]);out=Path(sys.argv[2]);out.mkdir(parents=True,exist_ok=False)
    p=r.load(panel);checks=[]
    def check(name,ok):
        checks.append({'name':name,'passed':bool(ok)})
        if not ok:raise AssertionError(name)
    def rejects(name,fn):
        try:value=fn();rejected=isinstance(value,dict) and value.get('status')=='failed_admission'
        except Exception:rejected=True
        check(name,rejected)
    old=r.read(r.CANONICAL)
    for record in p['cost_controls']:
        saved=r.read(r.b.verify(record['result']));cell=next(c for c in old['cells'] if c['cell_id']==record['cell_id']);text=r.b.verify(saved['output']).read_text()
        replay=r.parse(text,cell,saved['returncode'])
        check(record['cell_id']+' saved canonical parser',replay['status']=='complete' and replay['energy_hartree']==saved['energy_hartree'] and all(replay['checks'].values()))
        check(record['cell_id']+' precise final token',abs(float(replay['energy_decimal_tokens_hartree']['CCSD_T'])-saved['energy_hartree'])<1e-8)
    bad=copy.deepcopy(cell);bad['expected_basis_functions']=190
    rejects('wrong orbital count rejected',lambda:r.parse(text,bad,0))
    rejects('nonzero native status rejected',lambda:r.parse(text,cell,1))
    token=replay['energy_decimal_tokens_hartree']['CCSD_T']
    rejects('final energy closure rejected',lambda:r.parse(text.replace(token,'-1.000000000000'),cell,0))
    rejects('nonfinite final energy rejected',lambda:r.parse(text.replace(token,'NaN'),cell,0))
    def rows(cubic=None):
        cubic=cubic or {};result={};ha=D(str(r.b.HA))
        for c in p['cells']:
            q=D(str(c['q_A']));s=D(2 if c['metal']=='La' else 5);e=D('-260')+(s*q+D(str(cubic.get(c['metal'],0)))*q**3)/ha
            result[c['cell_id']]={'energy_decimal_tokens_hartree':{'CCSD_T':str(e),'CCSD':str(e+D('.04')),'triples':'-.040000000'}}
        return result
    scored=r.score(rows());check('eight exact linear slopes pass three checks',scored['all_checks_passed'] and all(abs(scored['groups'][m]['steps']['0.005']['CCSD_T']['derivative_kcal_mol_A']-s)<1e-10 for m,s in [('La',2),('Dy',5),('paired',3)]))
    scored=r.score(rows({'La':1000,'Dy':-1000}));check('paired two-step failure not hidden by individual passes',scored['checks']['La']['passed'] and scored['checks']['Dy']['passed'] and not scored['checks']['paired']['passed'] and not scored['all_checks_passed'])
    missing=rows();missing.pop(r.ALL[-1]);scored=r.score(missing)
    check('missing endpoint keeps all three checks unavailable',not scored['all_checks_passed'] and all(v['passed'] is None for v in scored['checks'].values()))
    bad=rows();bad[r.ALL[0]]['energy_decimal_tokens_hartree']['CCSD_T']='NaN';rejects('nonfinite score rejected',lambda:r.score(bad))
    check('decimal energy bound follows token precision',r.half_token('-1.000000000000')==D('0.0000000000005'))
    inv={cid:None for cid in r.ALL};malformed=out/'truncated.json';malformed.write_text('{');inv[r.ALL[0]]=str(malformed);ip=out/'missing_inventory.json';r.b.write(ip,inv)
    r.collect(panel,ip,out/'missing_collection');result=r.read(out/'missing_collection/RESULT.json')
    check('malformed native record retains exact eight slots',len(result['slots'])==8 and result['admitted_cells']==0 and result['status']=='canonical_rotation_incomplete')
    check('missing collection cannot pass numerical denominator',result['scoring']['declared_checks']==3 and all(v['passed'] is None for v in result['scoring']['checks'].values()) and not result['scoring']['all_checks_passed'])
    review=out/'mock_source_review.json';r.b.write(review,{'status':'passed_dry_acetate_correlated_rotation_source_review','implementation':r.b.pin(r.__file__),'panel':r.b.pin(panel),'synthetic_control_only':True})
    root=(out/'mock_execution_root').resolve();release={'status':'released_dry_acetate_correlated_rotation_scout','stage':'scout','implementation':r.b.pin(r.__file__),'panel':r.b.pin(panel),'cell_ids':[r.SCOUT],'maximum_total_native_calls':8,'maximum_stage_calls':1,'ranks':4,'nodes':1,'execution_root':str(root),'task_output_directories':{r.SCOUT:str(root/'tasks'/r.SCOUT/'execution')},'source_review':r.b.pin(review)}
    rp=out/'mock_release.json';r.b.write(rp,release);check('source-only scout release joins',r.release_gate(p,r.b.pin(rp),r.SCOUT)==release)
    rejects('scout cannot authorize another endpoint',lambda:r.release_gate(p,r.b.pin(rp),r.REMAINING[0]))
    bad=copy.deepcopy(release);bad['nodes']=2;bp=out/'bad_nodes.json';r.b.write(bp,bad);rejects('two-node release rejected',lambda:r.release_gate(p,r.b.pin(bp),r.SCOUT))
    result={'status':'passed_saved_controls','implementation':r.b.pin(r.__file__),'panel':r.b.pin(panel),'test_source':r.b.pin(__file__),'checks':checks,'passed':len(checks),'new_native_calls':0,'cost':{'wall_seconds':time.monotonic()-tick,'CPU_seconds':time.process_time()-cpu,'threads':1}}
    print(json.dumps({'result':r.b.write(out/'RESULT.json',result),'checks':len(checks),'cost':result['cost']},indent=2))

if __name__=='__main__':main()
