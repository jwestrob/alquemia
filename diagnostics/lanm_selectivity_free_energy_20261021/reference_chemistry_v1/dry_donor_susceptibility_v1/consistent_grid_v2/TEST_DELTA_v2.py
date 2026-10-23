#!/usr/bin/env python3
"""Focused saved/synthetic grid-version controls; no molecular execution."""
import ast, copy, importlib.util, json, re, sys, time
from pathlib import Path
HERE=Path(__file__).resolve().parent
s=importlib.util.spec_from_file_location('grid_candidate',HERE/'REFERENCE.py');m=importlib.util.module_from_spec(s);s.loader.exec_module(m)
tick,cpu=time.monotonic(),time.process_time()
panel=Path(sys.argv[1]);out=Path(sys.argv[2]);out.mkdir(parents=True,exist_ok=False)
p=m.load(panel);checks={}
def check(name,ok):
    checks[name]=bool(ok)
    if not ok:raise AssertionError(name)
def rejects(name,f):
    try:f()
    except Exception:check(name,True)
    else:check(name,False)

prior=m.read(p['prior_panel']);evidence=m.read(p['grid_evidence']);saved=m.read(p['failed_collection'])
for c,old in zip(p['cells'],prior['cells']):
    text=m.b.verify(c['template_input']).read_text()
    check(c['cell_id']+'__literal_only_grid_delta',text.replace(m.METHOD_BLOCK,'',1)==m.b.verify(old['template_input']).read_text() and text.count(m.METHOD_BLOCK)==1 and 'DoEQ' not in text)
check('Dy5_then_La5',all('__Dy__' in k for k in p['scout_cell_ids']) and all('__La__' in k for k in p['remaining_cell_ids']) and len(set(p['scout_cell_ids']+p['remaining_cell_ids']))==10)

# Exact unchanged function ASTs, including the energy/property and invocation paths.
before=ast.parse(m.b.verify(p['prior_adapter']).read_text());after=ast.parse(Path(m.__file__).read_text())
bd={f.name:f for f in before.body if isinstance(f,ast.FunctionDef)};ad={f.name:f for f in after.body if isinstance(f,ast.FunctionDef)}
names=('base_deck','old_center','one','array_tokens','field_properties','view','run','admit_result','compare_model')
check('unchanged_nine_core_function_ASTs',all(ast.dump(bd[n])==ast.dump(ad[n]) for n in names))
bs=copy.deepcopy(bd['score']);new=copy.deepcopy(ad['score'])
new.body[0].value.body.elts[0].value='La'
check('unchanged_score_except_scout_metal',ast.dump(bs)==ast.dump(new))

for metal in ('La','Dy'):
    c=next(c for c in p['cells'] if c['metal']==metal and c['field_tag']=='zero')
    raw=m.b.verify(evidence['current_failed_native_rows'][c['cell_id']]['output']).read_text()
    q=m.parse(raw,c,0)
    check(metal+'__old_recipe_cannot_admit',q['status']=='failed_admission' and not q['checks']['matching_physical_grids'] and not q['checks']['declared_grid_echo'])
    # Explicitly synthetic: change only the input echo and grid2 printed block.
    # This proves parser behavior, never a new scientific result.
    a=raw.index('\nGRIDX 2\n');z=raw.index('\nGRIDX 3\n');end=raw.index('Grids setup in',z)
    fixture=raw[:a]+raw[z:end].replace('GRIDX 3','GRIDX 2',1)+raw[z:]
    fixture=fixture.replace('****END OF INPUT****','\n'+''.join('| 900> '+line+'\n' for line in m.METHOD_BLOCK.rstrip().splitlines())+'****END OF INPUT****',1)
    q=m.parse(fixture,c,0)
    check(metal+'__synthetic_matching_grid_admitted',q['status']=='complete' and q['physical_grid_fingerprints']==c['expected_grid_fingerprints'])
    if metal=='Dy':
        bad=fixture.replace(' AngularGridX 2,4,4',' AngularGridX 2,3,4',1)
        check('wrong_declared_grid_rejected',m.parse(bad,c,0)['status']=='failed_admission')
        bad=fixture.replace('| 900> %method','| 900> %method\n| 900>  DoEQ false',1)
        check('unrelated_DoEQ_rejected',m.parse(bad,c,0)['status']=='failed_admission')
        bad=re.sub(r'(Maximum # iterations\s+MaxIter\s+\.{2,}\s*)300',r'\g<1>125',fixture)
        check('old_Mex_iteration_limit_rejected',m.parse(bad,c,0)['status']=='failed_admission')
        bad=fixture.replace('Recomputing exchange energy using gridx3','SYNTHETIC missing final integration')
        check('removed_final_integration_rejected',m.parse(bad,c,0)['status']=='failed_admission')
        bad=fixture.replace('\nGRIDX 3\n','\nSYNTHETIC MISSING GRID HEADER\n')
        q=m.parse(bad,c,0)
        check('missing_grid_grammar_unavailable',q['status']=='failed_admission' and q['grid_parse_error'] is not None)

check('original_14_of_16_failure_retained',m.score(saved['rows'],p,'full')==saved['response'] and sum(saved['response']['checks'].values())==14)
scout=m.score(saved['rows'],p,'scout')
check('Dy_actual_old_scout_remains_unqualified',len(scout['checks'])==6 and sum(scout['checks'].values())==4 and not scout['passed'])
for stage,ids,denom in [('scout',p['scout_cell_ids'],6),('full',[c['cell_id'] for c in p['cells']],16)]:
    inv={cid:str(out/(cid+'.missing.json')) for cid in ids}
    # One malformed file also stays unavailable, preserving all other slots.
    Path(inv[ids[0]]).write_text('{truncated')
    ip=out/(stage+'_inventory.json');m.save(ip,inv)
    code=m.collect(panel,ip,out/(stage+'_collection'),stage)
    result=json.loads((out/(stage+'_collection/RESULT.json')).read_text())
    check(stage+'__missing_malformed_full_denominator',code==1 and result['complete']==0 and result['declared']==len(ids) and len(result['rows'])==len(ids) and len(result['response']['checks'])==denom and not any(result['response']['checks'].values()) and result['new_native_calls']==0)

result=dict(status='passed',implementation=m.b.pin(__file__),adapter=p['adapter'],panel=m.b.pin(panel),checks=checks,passed=sum(checks.values()),declared=len(checks),
    fixture_notice='Modified raw text above is synthetic parser evidence only; no native calculation or new scientific result.',
    new_native_calls=0,cost=dict(wall_seconds=time.monotonic()-tick,CPU_seconds=time.process_time()-cpu,threads=1))
print(json.dumps(m.save(out/'RESULT.json',result)))
