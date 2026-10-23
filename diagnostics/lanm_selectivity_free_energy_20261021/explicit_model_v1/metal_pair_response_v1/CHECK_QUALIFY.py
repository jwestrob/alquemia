#!/usr/bin/env python3
"""Small saved-only controls for inventory, rigid paths, and native readback."""
import argparse
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import time

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('qualification_controls_target', HERE/'QUALIFY.py')
q = importlib.util.module_from_spec(spec)
spec.loader.exec_module(q)


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--manifest', required=True)
    p.add_argument('--out', required=True)
    args = p.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=False)
    start, cpu = time.monotonic(), time.process_time()
    checks = []
    def check(name, ok):
        checks.append({'name':name,'passed':bool(ok)})
        q.require(ok, name)
    def reject(name, fn):
        try:
            fn()
        except ValueError:
            check(name, True)
        else:
            check(name, False)
    result = {'status':'started','implementation':q.pin(__file__),'target':q.pin(HERE/'QUALIFY.py'),'manifest':q.pin(args.manifest),'checks':checks,'native_calls':0,'fixtures_are_saved_metadata_or_explicit_synthetic_readback':True,'error':None}
    try:
        check('one numerical thread', all(os.environ.get(k) == '1' for k in q.THREADS))
        b = q.read(result['manifest'])
        check('prepared source joins', b['implementation'] == result['target'] and b['plan'] == q.pin(HERE/'QUALIFICATION_PLAN.md'))
        jobs = q.inventory()
        check('56 unique jobs and 9 source fixtures', len(jobs) == len({j['id'] for j in jobs}) == 56 and len(b['rows']) == 9)
        check('40 testgrad 16 analyze 10 Cartesian', sum(j['program']=='testgrad' for j in jobs)==40 and sum(j['program']=='analyze' for j in jobs)==16 and sum(j['fd'] for j in jobs)==10)
        check('384 internal energies', sum(6*b['rows'][j['cell']]['atoms'] for j in jobs if j['fd']) == 384)
        for g, row in b['rows'].items():
            check(g+' immutable input bytes', all(q.pin(p['path']) == p and p['sha256'] == row['source_input'][k]['sha256'] for k,p in row['input'].items()))
        for g in q.NEAR:
            coords = b['rows'][g]['positions_A']
            for kind in ('radial','rotation'):
                plus,tangent = q.path_geometry(coords,kind,1e-6)
                minus,_ = q.path_geometry(coords,kind,-1e-6)
                derivative = [[(a-b)/2e-6 for a,b in zip(r,s)] for r,s in zip(plus,minus)]
                check(g+' '+kind+' unit tangent', abs(q.norm([x for r in tangent for x in r])-1)<1e-12)
                check(g+' '+kind+' signed geometry derivative', max(abs(a-b) for r,s in zip(derivative,tangent) for a,b in zip(r,s))<1e-9)
                for h in (1e-4,-1e-4,5e-5,-5e-5):
                    changed,_ = q.path_geometry(coords,kind,h)
                    distance = lambda xyz,a,b:q.norm([xyz[a][j]-xyz[b][j] for j in range(3)])
                    check(g+' '+kind+f' rigid water {h:g}', max(abs(distance(changed,a,z)-distance(coords,a,z)) for a,z in ((1,2),(1,3),(2,3)))<1e-12)
        # Synthetic native text deliberately tests admission, not molecular behavior.
        for g in (q.NEAR[0],q.NEAR[1],q.PURE,q.IONS[0]):
            row = b['rows'][g]
            key = q.verify(row['input']['key']).read_text()
            for profile in ('zero','direct','mutual','combined'):
                lines=[]
                for metal in (500,501):
                    if metal not in row['types'] or 36 not in row['types']:
                        lines.append(f'PAIR RESPONSE ABSENT TYPES: {metal} 36')
                    else:
                        gamma = [float(t[3]) for t in (line.split() for line in key.splitlines()) if len(t)>=4 and t[0].lower()=='polpair' and {int(t[1]),int(t[2])}=={metal,36}][-1]
                        dd,dm=q.PROFILES[profile]
                        lines.append(f'PAIR RESPONSE TYPES/GAMMA0/D/M: {metal} 36 {gamma:.17g} {gamma*(1+dd):.17g} {gamma*(1+dm):.17g}')
                text='\n'.join(lines)
                check(g+' '+profile+' synthetic native readback', q.response_readback(text,profile,row['types'],key)['passed'])
                reject(g+' '+profile+' missing native readback', lambda:q.response_readback('',profile,row['types'],key))
        row=b['rows'][q.NEAR[0]];key=q.verify(row['input']['key']).read_text()
        reject('no-op gamma readback', lambda:q.response_readback('PAIR RESPONSE TYPES/GAMMA0/D/M: 500 36 .17662704152555941 .17662704152555941 .17662704152555941\nPAIR RESPONSE ABSENT TYPES: 501 36','combined',row['types'],key))
        reject('nonfinite native readback', lambda:q.response_readback('PAIR RESPONSE TYPES/GAMMA0/D/M: 500 36 NaN NaN NaN\nPAIR RESPONSE ABSENT TYPES: 501 36','combined',row['types'],key))
        reject('absent keyword with active native report', lambda:q.response_readback('PAIR RESPONSE ABSENT TYPES: 500 36','absent',row['types'],key))
        check('legacy identity thresholds and exact observation', q.identity({'energy':0.,'forces':[[0.,0.,0.]]},{'energy':0.,'forces':[[0.,0.,0.]]})['printed_exact'])
        check('identity force failure retained', not q.identity({'energy':0.,'forces':[[0.,0.,0.]]},{'energy':0.,'forces':[[.002,0.,0.]]})['passed'])
        result['status']='saved_qualification_controls_passed'
    except Exception as exc:
        result.update(status='saved_qualification_controls_failed',error=repr(exc))
    finally:
        result.update(check_count=len(checks),cost={'wall_seconds':time.monotonic()-start,'CPU_seconds':time.process_time()-cpu,'numerical_threads':1})
        q.save(out/'RESULT.json',result)
    print(json.dumps({'status':result['status'],'checks':len(checks)}))
    if result['error']:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
