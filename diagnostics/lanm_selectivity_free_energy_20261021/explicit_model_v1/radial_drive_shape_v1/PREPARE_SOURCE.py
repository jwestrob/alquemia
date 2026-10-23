#!/usr/bin/env python3
"""Exact copied parent + two-file source delta. No compiler or native calls."""
import argparse, difflib, hashlib, json, shutil, time
from pathlib import Path

def pin(p):
    p=Path(p).resolve(); return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def checked(p):
    assert pin(p['path'])==p, p
    return Path(p['path'])
def dump(p,x): p.write_text(json.dumps(x,indent=2,allow_nan=False)+'\n')
def once(s,old,new):
    assert s.count(old)==1,(old,s.count(old));return s.replace(old,new)

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--inputs',required=True);ap.add_argument('--out',required=True);a=ap.parse_args()
    start=time.monotonic();cpu=time.process_time(); ip=Path(a.inputs); inp=json.loads(ip.read_text())
    for p in inp['sources'].values(): checked(p)
    src=inp['sources']['parent_source_manifest']; sm=json.loads(checked(src).read_text()); parent=Path(sm['source_root'])
    out=Path(a.out).resolve();out.mkdir(parents=True,exist_ok=False)
    for p in sm['source_outputs']+[sm['cmake']]:
        f=checked(p); dest=out/f.relative_to(parent);dest.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(f,dest)
    helper=checked(inp['sources']['radial_helper']).read_text()
    old=(parent/'source/induce.f').read_text();s=old
    # Deleting the bounded marked insertions restores the parent byte for byte.
    s=once(s,'c     Absent and identically zero records use literal original paths.',
        'c     RADIAL_SHAPE_BEGIN_SETUP\n      call radialdrivesetup\nc     RADIAL_SHAPE_END_SETUP\n'
        'c     Absent and identically zero records use literal original paths.')
    s=once(s,'      if (.not.all(ieee_is_finite(lrh)) .or.',
        'c     RADIAL_SHAPE_BEGIN_GEOMETRY\n      if (lrshape(slot).ne.0.0d0) then\n'
        '         call radialdrivegeometry(nw,ow,rr,vec,ff,fr,ph,pr,\n'
        '     &                            m,slot)\n      end if\nc     RADIAL_SHAPE_END_GEOMETRY\n'
        '      if (.not.all(ieee_is_finite(lrh)) .or.')
    s+='\nc     RADIAL_SHAPE_BEGIN_HELPER\n'+helper+'c     RADIAL_SHAPE_END_HELPER\n'
    (out/'source/induce.f').write_text(s)
    p=parent/'source/polar.f'; s=p.read_text();s=once(s,'      real*8 :: lrcoef(4,2) = 0.0d0',
        '      real*8 :: lrcoef(4,2) = 0.0d0\nc     RADIAL_SHAPE_BEGIN_MODULE\n'
        '      real*8 :: lrshape(2) = 0.0d0\nc     RADIAL_SHAPE_END_MODULE')
    (out/'source/polar.f').write_text(s)
    changed=[];patch=[]
    for p in sorted((out/'source').iterdir()):
        b=parent/'source'/p.name
        if p.read_bytes()!=b.read_bytes():
            changed.append(p.name);patch.extend(difflib.unified_diff(b.read_text().splitlines(True),p.read_text().splitlines(True),fromfile='parent/source/'+p.name,tofile='radial/source/'+p.name))
    assert changed==['induce.f','polar.f'],changed
    pp=out/'RADIAL_DRIVE.patch';pp.write_text(''.join(patch))
    result={'status':'prepared_radial_drive_source_only','implementation':pin(__file__),'inputs':pin(ip),'parent_source_manifest':src,'source_root':str(out),'changed_files':changed,'patch':pin(pp),'source_outputs':[pin(p) for p in sorted((out/'source').iterdir()) if p.is_file()],'cmake':pin(out/'cmake/CMakeLists.txt'),'build':None,'executables':None,'native_calls':0,'builds':0,'cost':{'wall_seconds':time.monotonic()-start,'CPU_seconds':time.process_time()-cpu,'numerical_threads':1}}
    dump(out/'MANIFEST.json',result);print(json.dumps({'manifest':pin(out/'MANIFEST.json'),'changed_files':changed,'cost':result['cost']}))
if __name__=='__main__':main()
