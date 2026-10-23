#!/usr/bin/env python3
"""Manufactured key/source checks only; no actual candidate or geometry processing."""
import hashlib
import importlib.util
import json
from pathlib import Path
import re
import time


def pin(path):
    path=Path(path).resolve()
    return dict(path=str(path),sha256=hashlib.sha256(path.read_bytes()).hexdigest())


def main():
    start,cpu=time.monotonic(),time.process_time()
    d=Path(__file__).resolve().parent
    w=Path(str(d).replace('/diagnostics/','/workspaces/',1))/'controls_v1'
    w.mkdir(parents=True,exist_ok=False)
    spec=importlib.util.spec_from_file_location('derive',d/'DERIVE.py')
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    checks=[]
    def check(name,value):
        checks.append(dict(name=name,passed=bool(value)))
        if not value:raise AssertionError(name)
    result=dict(status='failed',new_native_calls=0,builds=0,actual_candidate_reads=0,actual_box_evaluations=0)
    try:
        prefix='parameters manufactured.prm\npolar-eps 1e-8\npolarize 900 1.6597 .25\n'
        la='POLAR-LOCAL-RESPONSE 900 402 .1 -.1 .1 -.05\n'
        dy='POLAR-LOCAL-RESPONSE 901 402 -.1 .1 -.1 .05\n'
        source=prefix+la+'# intervening literal line\n'+dy
        text,receipt=module.derive(source)
        check('exact two removed lines',receipt['removed']==[dict(line_index=3,literal_line=la),dict(line_index=5,literal_line=dy)])
        check('literal rest and settings preserved',text==prefix+'# intervening literal line\n'+module.SETTINGS)
        check('reconstruction and no energy admission',receipt['exact_reconstruction'] and not receipt['energy_admitted'])
        for name,bad in [('missing Dy',source.replace(dy,'')),('duplicate',source+la),('wrong oxygen',source.replace('900 402','900 403')),
                         ('wrong namespace',source.replace('900 402','500 402')),('nonfinite',source.replace('.1 -.1 .1','nan -.1 .1')),
                         ('outside domain',source.replace('.1 -.1 .1','2 -.1 .1')),('preexisting PME',source+'pme-grid 64 64 64\n'),
                         ('external field',source+'external-field 0 0 0\n'),('no final newline',source[:-1])]:
            try:module.derive(bad)
            except (ValueError,OverflowError):rejected=True
            else:rejected=False
            check(name+' rejected',rejected)
        native=json.loads((d.parent/'PATCH_INPUTS.json').read_text())
        base=Path(native['parent_source']['path']).parent
        harness=(d/'INIT_PME.f90').read_text()
        modules=re.findall(r'^\s*use\s+(\w+)',harness,re.M|re.I)
        check('actual imported modules',all(v=='ieee_arithmetic' or (base/(v+'.f')).is_file() for v in modules))
        calls=re.findall(r'^\s*call\s+(\w+)',harness,re.M|re.I)
        check('exact bounded native call chain',calls==['initial','getxyz','mechanic','nextarg','gettext','upcase','chkpole','rotpole','localresponsesetup','udirect1','final'])
        check('no solver/energy/gradient/real field',not set(calls)&{'induce','induce0a','energy','gradient','dfield0','dfield0c','udirect2a','udirect2b'})
        check('normal response setup retained',calls.count('localresponsesetup')==1 and 'if (lractive.or.lrperiodic)' in harness)
        check('single reciprocal initialization',calls.count('udirect1')==1)
        check('dfield0c initialization assignments',all(s in harness for s in ['nfft1=nefft1','nfft2=nefft2','nfft3=nefft3','bsorder=bsporder','aewald=apewald']))
        check('free form line limit',max(map(len,harness.splitlines()))<=132)
        result['status']='passed_manufactured_key_and_source_controls'
    except Exception as error:result['error']=repr(error)
    result.update(checks=checks,implementation=pin(__file__),derive=pin(d/'DERIVE.py'),harness=pin(d/'INIT_PME.f90'),
                  wall_seconds=time.monotonic()-start,CPU_seconds=time.process_time()-cpu,threads=1)
    (w/'RESULT.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'status':result['status'],'checks':len(checks),'result':pin(w/'RESULT.json')}))
    if result['status'].startswith('failed'):raise SystemExit(1)


if __name__=='__main__':main()
