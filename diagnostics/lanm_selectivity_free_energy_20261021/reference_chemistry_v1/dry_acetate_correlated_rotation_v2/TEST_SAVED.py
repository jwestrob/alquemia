#!/usr/bin/env python3
"""Narrow actual-scout equality and malformed-metadata controls, zero native."""
import ast,copy,importlib.util,json,time
from pathlib import Path
H=Path(__file__).resolve().parent;s=importlib.util.spec_from_file_location('v2',H/'REFERENCE.py');m=importlib.util.module_from_spec(s);s.loader.exec_module(m)
checks={}
def check(k,v):checks[k]=bool(v);m.b.check(v,k)
def rejects(k,fn):
 try:fn()
 except Exception as e:checks[k]=dict(passed=True,error=repr(e));return
 raise ValueError('Expected rejection '+k)
def main():
 tick,cpu=time.monotonic(),time.process_time();x=m.read(H/'INPUTS.json');p=m.load(x['panel']['path']);r=m.read(m.b.verify(x['legacy_scout_result']));c=next(c for c in p['cells'] if c['cell_id']==m.SCOUT);raw=m.parse(m.b.verify(r['output']).read_text(),c,0)
 different=[k for k,v in raw.items() if r[k]!=v];check('only_basis_groups_representation',different==['basis_groups']);m.raw_equal(r,raw);check('complete_exact_JSON_replay',m.canonical({k:r[k] for k in raw})==m.canonical(raw))
 for name,edit in [('energy',lambda z:z.update(energy_hartree=z['energy_hartree']+1e-4)),('basis_count',lambda z:z['basis_groups'][next(iter(z['basis_groups']))].update(nbf=999) if isinstance(z['basis_groups'][next(iter(z['basis_groups']))],dict) else z['basis_groups'].update({next(iter(z['basis_groups'])):'wrong'})),('missing',lambda z:z.pop('energy_hartree')),('nonfinite',lambda z:z.update(energy_hartree=float('nan'))),('native_status',lambda z:z.update(status='failed_admission'))]:
  bad=copy.deepcopy(r);edit(bad);rejects(name+'_still_rejected',lambda:m.raw_equal(bad,raw))
 oldtree=ast.parse(Path(m.old.__file__).read_text());newtree=ast.parse(Path(m.__file__).read_text())
 for name in ('run','collect'):
  aa=next(n for n in oldtree.body if isinstance(n,ast.FunctionDef) and n.name==name);bb=next(n for n in newtree.body if isinstance(n,ast.FunctionDef) and n.name==name);check(name+'_unchanged_function_AST',ast.dump(aa)==ast.dump(bb))
 check('deck_parser_score_unchanged_helpers',m.deck is m.old.deck and m.parse is m.old.parse and m.score is m.old.score)
 check('exact8_and_remaining7',m.ALL==m.old.ALL and len(m.ALL)==8 and len(m.REMAINING)==7 and m.SCOUT not in m.REMAINING)
 a=m.read(H/'independent_actual_review_v1/CONFIG.json');rv=m.read(m.WS/'x') if False else m.read(Path.cwd()/'workspaces/lanm_selectivity_free_energy_20261021/reference_chemistry_v1/dry_acetate_correlated_rotation_v2/independent_actual_review_v1/actual_v1/REVIEW.json')
 check('actual_recovered_scout_full_review',rv['status']=='passed_dry_acetate_correlated_rotation_scout_review' and rv['result']==x['legacy_scout_result'] and rv['implementation']==m.b.pin(m.__file__) and all(q['passed'] for q in rv['checks']))
 folder=Path.cwd()/'workspaces/lanm_selectivity_free_energy_20261021/reference_chemistry_v1/dry_acetate_correlated_rotation_v2/saved_controls_v1';folder.mkdir(parents=True,exist_ok=False)
 out=dict(status='passed_saved_representation_controls',implementation=m.b.pin(m.__file__),test_implementation=m.b.pin(__file__),legacy_scout=x['legacy_scout_result'],checks=checks,check_count=len(checks),new_native_calls=0,cost=dict(wall_seconds=time.monotonic()-tick,CPU_seconds=time.process_time()-cpu,threads=1));m.b.write(folder/'RESULT.json',out);print(json.dumps(dict(result=m.b.pin(folder/'RESULT.json'),checks=len(checks))))
if __name__=='__main__':main()
