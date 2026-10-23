#!/usr/bin/env python3
"""One bounded saved control set for literal cycles and reused rendering."""
import argparse
import copy
import importlib.util
import json
from pathlib import Path

p=argparse.ArgumentParser();p.add_argument('--manifest',required=True);p.add_argument('--out',required=True);args=p.parse_args()
here=Path(__file__).resolve().parent;s=importlib.util.spec_from_file_location('regression_controls',here/'RUN_v2.py');r=importlib.util.module_from_spec(s);s.loader.exec_module(r)
m=r.read(r.pin(args.manifest));ev=r.module(m['inputs']['evaluator']);cfg=r.read(m['evaluator_template']);em=r.read(cfg['manifest']);domain=r.read(cfg['domain'])
resolver=ev.module(cfg['resolver'],'regression_control_resolver');checks={}
def ck(ok,name):
    checks[name]=bool(ok)
    if not ok:raise ValueError(name)
ck(set(m['rows'])==set(r.IDS) and len(r.IDS)==6 and len(m['cycles'])==8,'six additional rows/eight cycles')
values={cid:dict(status='complete',energy=v['energy_kcal_mol'],mu_Debye=v['mu_Debye']) for cid,v in em['static'].items()}
values.update(em['baseline_cells']);values.update(m['baseline_cells'])
baseline=r.cycles(m['cycles'],values)
def printing_bound(k):
    if k in ('ion_La','ion_Dy'):
        assert em['static'][k]['energy_kcal_mol']==0 and em['static'][k]['mu_Debye']==[0.,0.,0.]
        return .5e-6+2e-9  # One explicitly admitted native ion site.
    if k in em['static']:return em['static'][k]['mu_print_bound']
    return (em['rows'].get(k) or m['rows'][k])['atoms']*.5e-6+2e-9
for cid,row in baseline.items():
    terms=m['cycles'][cid]['model_terms'];cabs=sum(abs(c) for c in terms.values())
    ck(row['status']=='complete' and abs(row['candidate_E3_kcal_mol']-row['baseline_E3_kcal_mol'])<=2*cabs*5e-9,'baseline energy cycle '+cid)
    ours=sum(abs(c)*printing_bound(k) for k,c in terms.items())
    bound=ours+row['original_model_print_bound_Debye_per_component']
    ck(max(abs(a-b) for a,b in zip(row['candidate_mu3_Debye'],row['baseline_mu3_Debye']))<=bound,'baseline common-origin moment cycle '+cid)
for cid,row in m['rows'].items():
    key,meta=ev.render(m['initial_vector'],row,domain,resolver,exports=False)
    ck('RESPONSE-MATRIX-EXPORT' not in key and meta['POLPAIR_unchanged'] and meta['export_enabled'] is False,'exact reused rendering '+cid)
one=next(iter(m['cycles']));missing=next(k for k in m['cycles'][one]['model_terms'] if k.endswith('__metal_two_water'))
bad=copy.deepcopy(values);del bad[missing];unavailable=r.cycles(m['cycles'],bad)
ck(unavailable[one]['status']=='unavailable' and 'candidate_E3_kcal_mol' not in unavailable[one],'missing term stays unavailable')
bad=copy.deepcopy(values);bad[missing]['energy']=float('nan')
try:r.cycles(m['cycles'],bad)
except ValueError:ck(True,'nonfinite cycle rejects')
else:ck(False,'nonfinite cycle rejects')
out=Path(args.out).resolve();out.mkdir(parents=True,exist_ok=False)
result=dict(status='passed_saved_regression_controls',implementation=r.pin(here/'RUN_v2.py'),test_implementation=r.pin(__file__),manifest=r.pin(args.manifest),checks=checks,count=len(checks),baseline_cycles=baseline,new_native_calls=0,new_parameter_molecular_points=0)
(out/'RESULT.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n');print(json.dumps({'result':r.pin(out/'RESULT.json'),'checks':len(checks)}))
