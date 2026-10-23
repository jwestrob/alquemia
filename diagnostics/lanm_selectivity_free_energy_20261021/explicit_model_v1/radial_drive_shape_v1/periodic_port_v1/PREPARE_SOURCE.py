#!/usr/bin/env python3
"""Source-only, exclusive two-file overlay on the pinned periodic parent."""
import difflib
import hashlib
import json
from pathlib import Path
import shutil
import time

HERE=Path(__file__).resolve().parent

def pin(path):
    path=Path(path).resolve()
    return {'path':str(path),'sha256':hashlib.sha256(path.read_bytes()).hexdigest()}

def verify(ref):
    path=Path(ref['path'])
    if pin(path)!=ref:raise ValueError('Changed source authority '+str(path))
    return path

def save(path,obj):
    with Path(path).open('x') as stream:
        json.dump(obj,stream,indent=2,allow_nan=False);stream.write('\n')

def routine(text,name):
    start=text.index('      subroutine '+name+'\n')
    return text[start:text.index('      end\n',start)+len('      end\n')]

def transform(original,helper):
    changed=dict(original);edits=[]
    def replace(name,old,new):
        if changed[name].count(old)!=1:raise ValueError('Nonunique join '+old[:80])
        changed[name]=changed[name].replace(old,new,1)
        edits.append({'file':name,'before':old,'after':new})
    replace('polar.f','      real*8 :: lrcoef(4,2) = 0.0d0\n',
            '      real*8 :: lrcoef(4,2) = 0.0d0\n'
            'c     BIO-only independent radial drive, absent/zero inert.\n'
            '      real*8 :: lrpshape(2) = 0.0d0\n')
    old=routine(changed['induce.f'],'periodicresponsesetup');new=old
    new=new.replace('      handled=.false.\n',
                    '      lrpshape=0.0d0\n      handled=.false.\n',1)
    # Every new keyword enters the BIO parser, including malformed/mixed ones.
    marker="         if (trim(keyword).ne.'POLAR-LOCAL-RESPONSE') cycle\n"
    new=new.replace(marker,
        "         if (trim(keyword).eq.'POLAR-RADIAL-DRIVE')\n"
        '     &      handled=.true.\n'+marker,1)
    new=new.replace('c     Zero BIO coefficients add no restriction or arithmetic.\n',
        'c     Parse e before inactive bypass; e-only uses all normal guards.\n'
        '      call periodicradialsetup\n'
        'c     Zero BIO coefficients add no restriction or arithmetic.\n',1)
    replace('induce.f',old,new)
    old=routine(changed['induce.f'],'periodicresponsegeometry');new=old
    marker='      if (.not.all(ieee_is_finite(lrh)) .or.\n'
    if old.count(marker)!=1:raise ValueError('Geometry finite-check join')
    new=new.replace(marker,
        '      if (lrpshape(slot).ne.0.0d0)\n'
        '     &   call periodicradialfield(nw,ow,rr,vec,ff,ph,slot)\n'+marker)
    replace('induce.f',old,new)
    old=routine(changed['induce.f'],'periodicresponsegradient');new=old
    marker='         gg=fac*gg\n'
    if old.count(marker)!=1:raise ValueError('Gradient unit-conversion join')
    new=new.replace(marker,
        '         if (lrpshape(slot).ne.0.0d0)\n'
        '     &      call periodicradialgradient(i,lrpshape(slot),w,gg)\n'+marker)
    replace('induce.f',old,new)
    addition='\nc\n'+helper
    changed['induce.f']+=addition
    inverse=dict(changed);inverse['induce.f']=inverse['induce.f'][:-len(addition)]
    for edit in reversed(edits):
        if inverse[edit['file']].count(edit['after'])!=1:raise ValueError('Inverse join')
        inverse[edit['file']]=inverse[edit['file']].replace(edit['after'],edit['before'],1)
    if inverse!=original:raise ValueError('Parent byte restoration failed')
    return changed,edits

def main():
    tick,cpu=time.monotonic(),time.process_time()
    inp=json.loads((HERE/'INPUTS.json').read_text())
    for ref in inp['authorities'].values():verify(ref)
    parent=json.loads(verify(inp['authorities']['parent_manifest']).read_text())
    parent_root=Path(parent['source_root']);out=Path(inp['source_output'])
    out.mkdir(parents=True,exist_ok=False)
    originals={name:verify(inp['authorities']['parent_'+name[:-2]]).read_text()
               for name in ('polar.f','induce.f')}
    changed,edits=transform(originals,(HERE/'PERIODIC_RADIAL.f').read_text())
    copied=[]
    # Copy only the complete frozen source inventory, never a live build tree.
    for ref in parent['source_outputs']:
        source=verify(ref);relative=source.relative_to(parent_root)
        dest=out/relative;dest.parent.mkdir(parents=True,exist_ok=True)
        shutil.copyfile(source,dest);copied.append(relative)
    for name,text in changed.items():
        path=out/'source'/name
        if path.read_text()!=originals[name]:raise ValueError('Copied parent mismatch')
        path.write_text(text)
        for line in text.splitlines():
            if line and line[0].lower()!='c' and len(line)>72:
                raise ValueError('Fixed-form width '+name+': '+line)
        (out/(name+'.diff')).write_text(''.join(difflib.unified_diff(
            originals[name].splitlines(True),text.splitlines(True),
            fromfile=str(parent_root/'source'/name),tofile=str(path))))
    save(out/'EDITS.json',edits)
    outputs=[pin(out/p) for p in copied]
    differences=[str(p) for p,ref in zip(copied,parent['source_outputs'])
                 if pin(out/p)['sha256']!=ref['sha256']]
    if sorted(differences)!=['source/induce.f','source/polar.f']:
        raise ValueError('Unexpected source delta '+str(differences))
    result={'status':'prepared_periodic_radial_source_only_unbuilt',
        'producer':pin(__file__),'inputs':pin(HERE/'INPUTS.json'),
        'root_declaration':inp['authorities']['root_declaration'],
        'parent_manifest':inp['authorities']['parent_manifest'],
        'source_root':str(out),'source_outputs':outputs,'changed_files':differences,
        'inverse_patch_restores_exact_parent':True,'edits':pin(out/'EDITS.json'),
        'diffs':[pin(out/(name+'.diff')) for name in changed],
        'normal_periodic_guards_remain_fatal':True,
        'pending':['future e-aware state schema','candidate composition',
                   'cutoff12 exact-state/operator admission','build and finite qualification'],
        'actual':{'compilations':0,'native_calls':0,'solvers':0,'factorizations':0,
                  'wall_seconds':time.monotonic()-tick,'CPU_seconds':time.process_time()-cpu,
                  'threads':1}}
    save(out/'MANIFEST.json',result)
    print(json.dumps({'manifest':pin(out/'MANIFEST.json'),'changed_files':differences,
                      'actual':result['actual']}))

if __name__=='__main__':main()
