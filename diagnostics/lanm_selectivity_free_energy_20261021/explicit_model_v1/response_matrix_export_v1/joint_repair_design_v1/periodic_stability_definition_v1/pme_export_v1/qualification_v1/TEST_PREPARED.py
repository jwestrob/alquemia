#!/usr/bin/env python3
"""Finite source/input controls only; no compiler, FFT or native call."""
import ast
import hashlib
import json
from pathlib import Path
import re
import time


def pin(p):
    p=Path(p).resolve()
    return dict(path=str(p),sha256=hashlib.sha256(p.read_bytes()).hexdigest())


def main():
    start,cpu=time.monotonic(),time.process_time()
    d=Path(__file__).resolve().parent
    w=Path(str(d).replace('/diagnostics/','/workspaces/',1))
    out=w/'controls_v1';out.mkdir(parents=True,exist_ok=False)
    result=dict(status='failed',controls=[],native_calls=0,builds=0,FFT_calls=0)
    def check(name,value):
        result['controls'].append(dict(name=name,passed=bool(value)))
        if not value:raise AssertionError(name)
    try:
        m=json.loads((w/'preparation_v1/MANIFEST.json').read_text())
        for item in m['source_input'].values():check('source fixture pin',pin(item['path'])==item)
        check('neutral fixture',m['source_fixture_charge']==0)
        for mode,row in m['cells'].items():
            for item in row.values():check('prepared pin '+mode,pin(item['path'])==item)
            check('literal XYZ '+mode,row['xyz']['sha256']==m['source_input']['xyz']['sha256'])
            check('literal parameters '+mode,row['parameters']['sha256']==m['source_input']['parameters']['sha256'])
        off=Path(m['cells']['OFF']['key']['path']).read_text()
        on=Path(m['cells']['ON']['key']['path']).read_text()
        check('only off/on delta',on==off+'PME-NORM-EXPORT 1\n')
        check('original key prefix',off.startswith(Path(m['source_input']['key']['path']).read_text()))
        records=[v.split()[0].lower() for v in off.splitlines() if v.split() and not v.lstrip().startswith('#')]
        check('response/external/vacuum boundary absent',not set(records)&{'polar-local-response','polar-donor-response','external-field','ewald-boundary','response-matrix-export','dewald'})
        source=(d/'PME_HARNESS.f90').read_text()
        parent=json.loads((d.parent/'PATCH_INPUTS.json').read_text())
        native=Path(parent['parent_source']['path']).parent
        modules=re.findall(r'^\s*use\s+(\w+)',source,re.M|re.I)
        check('all imported modules exist',all(v=='ieee_arithmetic' or (native/(v+'.f')).is_file() for v in modules))
        check('one gradient call',source.lower().count('call gradient(')==1)
        check('all native reciprocal stages once',all(source.lower().count('call '+v+'(')==1 if v in ['grid_uind','fphi_uind'] else source.lower().count('call '+v+'\n')==1 for v in ['grid_uind','fftfront','fftback','fphi_uind']))
        check('no global uind/uinp assignment',not re.search(r'^\s*(uind|uinp)\([^\n]*?\)\s*=',source,re.M|re.I))
        check('bounded three modes and eighteen columns','do mode=1,3' in source and 'do col=1,18' in source)
        check('fresh files and finite readbacks',source.count("status='new'")==3 and 'ieee_is_finite' in source)
        check('free form line bound',max(map(len,source.splitlines()))<=132)
        for name in ['PREPARE.py','REPLAY.py']:
            tree=ast.parse((d/name).read_text())
            check(name+' has no process invocation',not any(isinstance(n,ast.Name) and n.id in ['subprocess','exec','eval'] for n in ast.walk(tree)))
        check('release/build/link unbound',all(m[k] is None for k in ['root_release','build_result','harness_executable']))
        result['status']='passed_source_and_prepared_fixture_controls'
    except Exception as error:result['error']=repr(error)
    result.update(implementation=pin(__file__),harness=pin(d/'PME_HARNESS.f90'),replay=pin(d/'REPLAY.py'),
                  preparation=pin(w/'preparation_v1/MANIFEST.json'),wall_seconds=time.monotonic()-start,CPU_seconds=time.process_time()-cpu,threads=1)
    (out/'RESULT.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'status':result['status'],'checks':len(result['controls']),'result':pin(out/'RESULT.json')}))
    if result['status'].startswith('failed'):raise SystemExit(1)


if __name__=='__main__':main()
