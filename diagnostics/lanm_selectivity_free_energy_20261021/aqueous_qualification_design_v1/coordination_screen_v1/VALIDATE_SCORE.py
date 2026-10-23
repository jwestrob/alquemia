#!/usr/bin/env python3
"""Saved/synthetic scoring controls only; never invokes a molecular program."""
import copy
import importlib.util
import json
import math
from pathlib import Path
import time

HERE = Path(__file__).resolve().parent
s = importlib.util.spec_from_file_location('screen', HERE / 'RUN.py')
m = importlib.util.module_from_spec(s)
s.loader.exec_module(m)

def main():
    start, cpu = time.monotonic(), time.process_time()
    path = m.W / 'aqueous_qualification_design_v1/coordination_screen_v1/preparation_v2/BOUND.json'
    bp = m.pin(path)
    bound = m.read(bp)
    body = m.module(bound['body_helper'], 'screen_body')
    baseline = copy.deepcopy(bound)
    baseline['anchor_model'] = copy.deepcopy(baseline['anchor_reference'])
    values = {g:copy.deepcopy(baseline['rows'][g]['reference']) for g in m.CELLS}
    checks = []
    def check(ok, name):
        checks.append({'name':name, 'passed':bool(ok)})
        if not ok: raise AssertionError(name)
    def near(a,b): return abs(a-b) < 1e-7
    exact = m.score(baseline, values, body)
    check(exact['passed'] and exact['available_comparisons']==exact['declared_comparisons']==4, 'Self-reference has all four comparisons')
    check(exact['energy']['error_kcal_mol']==0, 'Self-reference energy is zero')
    check(all(r['RMS_kcal_mol_A']==r['maximum_component_kcal_mol_A']==0 for r in exact['rigid'].values()), 'Self-reference rigid responses are zero')
    signed = copy.deepcopy(values)
    signed[m.CELLS[0]]['energy'] += .4
    signed[m.CELLS[1]]['energy'] += .9
    b = copy.deepcopy(baseline)
    b['anchor_model']['La']['energy'] += .3
    b['anchor_model']['Dy']['energy'] += .1
    check(near(m.score(b,signed,body)['energy']['error_kcal_mol'], .7), 'Dy-minus-La and W9-minus-W8 sign')
    shifted = copy.deepcopy(values)
    b = copy.deepcopy(baseline)
    for metal, offset in (('La',123.0),('Dy',-71.0)):
        shifted['aquo9_central_'+metal]['energy'] += offset
        b['anchor_model'][metal]['energy'] += offset
    check(m.score(b,shifted,body)['passed'], 'Element-dependent free-fragment offsets cancel')
    high = copy.deepcopy(values)
    high[m.CELLS[1]]['energy'] += 2
    scored = m.score(baseline,high,body)
    check(near(scored['energy']['error_kcal_mol'],2) and not scored['passed'] and not scored['energy']['passed'], 'Two-kcal paired error fails without filtering')
    check(all(v['passed'] for v in scored['rigid'].values()), 'Energy failure leaves force checks independent')
    shared = copy.deepcopy(values)
    for g in m.CELLS: shared[g]['forces'][0][0] += 10
    scored = m.score(baseline,shared,body)
    check(not scored['rigid']['La']['passed'] and not scored['rigid']['Dy']['passed'] and not scored['passed'], 'Shared large force error cannot pass individual gates')
    check(scored['rigid']['paired']['passed'] and near(scored['rigid']['paired']['maximum_component_kcal_mol_A'],0), 'Shared force error cancels only in paired response')
    check(all(near(scored['rigid'][metal]['RMS_kcal_mol_A'],10/math.sqrt(6)) for metal in ('La','Dy')), 'Ion force RMS has fixed equal-category normalization')
    check(all(near(scored['rigid'][metal]['maximum_component_kcal_mol_A'],10) for metal in ('La','Dy')), 'Maximum is the unweighted body component')
    for retained in (None, m.CELLS[0], m.CELLS[1]):
        partial = {g:(values[g] if g==retained else None) for g in m.CELLS}
        scored = m.score(baseline,partial,body)
        expected = 0 if retained is None else 1
        check(scored['declared_comparisons']==4 and scored['available_comparisons']==expected and not scored['passed'], 'Missing slots retain denominator: '+str(retained))
        check(scored['energy']['passed'] is None and scored['rigid']['paired']['passed'] is None, 'Missing paired quantities remain unknown: '+str(retained))
    for kind in ('nonfinite_energy','nonfinite_force','short_force_array'):
        bad = copy.deepcopy(values)
        if kind=='nonfinite_energy': bad[m.CELLS[0]]['energy']=float('nan')
        elif kind=='nonfinite_force': bad[m.CELLS[0]]['forces'][0][0]=float('inf')
        else: bad[m.CELLS[0]]['forces'].pop()
        try: m.score(baseline,bad,body)
        except (ValueError,TypeError): rejected=True
        else: rejected=False
        check(rejected, 'Reject malformed candidate: '+kind)
    check(m.pin(path)==bp, 'Original binding remains unchanged')
    out=m.W/'aqueous_qualification_design_v1/coordination_screen_v1/scoring_controls_v1'
    out.mkdir(parents=True,exist_ok=False)
    result={'status':'passed','implementation':m.pin(HERE/'RUN.py'),'control_source':m.pin(__file__),'bound':bp,'checks':checks,'checks_passed':len(checks),'new_native_calls':0,'wall_seconds':time.monotonic()-start,'CPU_seconds':time.process_time()-cpu}
    m.save(out/'RESULT.json',result)
    print(json.dumps({'result':m.pin(out/'RESULT.json'),'checks':len(checks),'wall_seconds':result['wall_seconds']}))

if __name__=='__main__': main()
