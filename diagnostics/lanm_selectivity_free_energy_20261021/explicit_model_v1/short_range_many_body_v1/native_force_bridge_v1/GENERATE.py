"""Generate private source clones only. Never compile, launch or edit parents."""
import argparse
import hashlib
import json
from pathlib import Path
import re
import time

PINS = {'epolar1.f':'b558f1b1b459da25c2639b41a4f2c1a095936f94356ac083a7cfb095b9301b76',
        'gradient.f':'1fed31b16d1fafb29a0bdc1f7a469d153a9850c69f1723f543bb06d5a50d592f',
        'induce.f':'35bdd1c9182ceb5994770e52221f5cb664b13cd1ddd42df1632d8953bc03d7ba'}


def pin(path):
    return {'path':str(path.resolve()),'sha256':hashlib.sha256(path.read_bytes()).hexdigest()}


def routine(source, name):
    start=re.search(r'^      subroutine '+re.escape(name)+r'(?:\s|\().*$',source,re.M)
    if not start: raise ValueError('missing '+name)
    end=re.search(r'^      end\s*$',source[start.end():],re.M)
    if not end:raise ValueError('missing end '+name)
    return source[start.start():start.end()+end.end()]+'\n'


def one(text, old, new):
    if text.count(old)!=1:raise ValueError('expected exactly one source anchor '+repr(old))
    return text.replace(old,new)


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--source',required=True);ap.add_argument('--out',required=True);args=ap.parse_args()
    start=time.monotonic();cpu=time.process_time();source=Path(args.source);out=Path(args.out)
    parents={name:pin(source/name) for name in PINS}
    for name in PINS:
        if parents[name]['sha256']!=PINS[name]:raise ValueError('parent pin mismatch '+name)
    rawp=routine((source/'epolar1.f').read_text(),'epolar1a')
    rawg=routine((source/'gradient.f').read_text(),'gradient')
    p=one(rawp,'      subroutine epolar1a\n','      subroutine srmb_epolar1a\n')
    p=one(p,'      call induce\n','c     Supplied fixed uind/uinp; initialized udir/udirp. No SCF.\n')
    g=one(rawg,'      subroutine gradient (energy,derivs)\n','      subroutine srmb_native_gradient (energy,derivs)\n')
    g=one(g,'      if (use_polar)  call epolar1\n','      if (use_polar)  call srmb_epolar1a\n')
    # Exact inverse joins guarantee preserved loops/torques and total assembly.
    inverse=p.replace('      subroutine srmb_epolar1a\n','      subroutine epolar1a\n').replace(
        'c     Supplied fixed uind/uinp; initialized udir/udirp. No SCF.\n','      call induce\n')
    if inverse!=rawp:raise AssertionError('polar clone changed inherited math')
    inverse=g.replace('      subroutine srmb_native_gradient (energy,derivs)\n','      subroutine gradient (energy,derivs)\n').replace(
        '      if (use_polar)  call srmb_epolar1a\n','      if (use_polar)  call epolar1\n')
    if inverse!=rawg:raise AssertionError('gradient clone changed inherited math')
    if re.search(r'^\s*call induce',p+g,re.M):raise AssertionError('unremoved induction call')
    out.mkdir(parents=True,exist_ok=False)
    (out/'srmb_epolar1a.f').write_text(p);(out/'srmb_native_gradient.f').write_text(g)
    result={'status':'generated_private_source_only','parents':parents,'generator':pin(Path(__file__)),
            'outputs':{name:pin(out/name) for name in ('srmb_epolar1a.f','srmb_native_gradient.f')},
            'joins':{'exact_inverse_polar_clone':True,'exact_inverse_gradient_clone':True,'direct_induce_calls':0},
            'compiler_calls':0,'native_calls':0,'numerical_tests':0,
            'cost':{'wall_seconds':time.monotonic()-start,'cpu_seconds':time.process_time()-cpu}}
    (out/'MANIFEST.json').write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'status':result['status'],'outputs':result['outputs'],'cost':result['cost']}))


if __name__=='__main__':main()
