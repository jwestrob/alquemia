#!/usr/bin/env python3
"""Copy the pinned water source and add the isolated donor channel. No build/run."""
import difflib, hashlib, json, shutil, time
from pathlib import Path
HERE=Path(__file__).resolve().parent
ROOT=Path.cwd()
WS=ROOT/'workspaces/lanm_selectivity_free_energy_20261021/explicit_model_v1/response_matrix_export_v1/joint_repair_design_v1'
SOURCES=HERE.parent/'donor_extension_design_v1/SOURCES.json'
def pin(path):
    p=Path(path).resolve();return dict(path=str(p),sha256=hashlib.sha256(p.read_bytes()).hexdigest())
def verify(p):
    if pin(p['path'])!=p:raise ValueError('Changed pin: '+p['path'])
    return Path(p['path'])
def save(p,x):
    with p.open('x') as f:f.write(json.dumps(x,indent=2,allow_nan=False)+'\n')
    return pin(p)
def replace(text,old,new):
    if text.count(old)!=1:raise ValueError('Expected one source anchor: '+old[:70])
    return text.replace(old,new)

tick,cpu=time.monotonic(),time.process_time()
assert pin(SOURCES)['sha256']=='cf76658930a7aa8ab2d440efa7f7dc6d2bb4de34e4f37c3536896c53e7a35e5c'
sources=json.loads(SOURCES.read_text())
verify(sources['note']);manifest=json.loads(verify(sources['native_manifest']).read_text())
for item in manifest['source_outputs']:verify(item)
verify(manifest['cmake'])
for item in sources['native_sources'].values():verify(item)
verify(sources['native_helper']);verify(sources['donor_inventory'])
prior=Path(manifest['source_root']);out=WS/'donor_extension_native_v1/source_v1'
out.mkdir(parents=True,exist_ok=False)
for name in ('source','cmake'):shutil.copytree(prior/name,out/name)

f=out/'source/polar.f';s=f.read_text()
s=replace(s,'      logical :: lractive = .false.','      logical :: lractive = .false.\n      logical :: lrdonor = .false.\n      real*8 :: drcoef(2,2) = 0.0d0\n      real*8, allocatable :: drdb(:,:,:)')
f.write_text(s)

f=out/'source/induce.f';s=f.read_text()
s=replace(s,'      if (matrix_export)  call exportresponsematrix','      if (matrix_export) then\n         if (.not.lrdonor) call donorresponsepermanent\n         call exportresponsematrix\n      end if')
s=replace(s,'      logical requested,bad\n      logical seen(28)','      logical requested,bad,donor_fixture\n      logical donorresponsefixture\n      logical seen(28)')
s=replace(s,'      seen = .false.\n      nmetal = 0','      seen = .false.\n      nmetal = 0\n      donor_fixture=donorresponsefixture()')
s=replace(s,'         if ((type(i).ne.36 .and. type(i).ne.37 .and.\n     &        type(i).ne.500 .and. type(i).ne.501) .or.',
    '         if (((type(i).ne.36 .and. type(i).ne.37 .and.\n     &        type(i).ne.500 .and. type(i).ne.501)\n     &        .and. .not.donor_fixture) .or.')
s=replace(s,'      lrm = 0\n      seen = .false.','      lrm = 0\n      call donorresponsesetup\n      if (lrdonor) return\n      seen = .false.')
s=replace(s,'      dual=0.0d0\n      rhsdual=0.0d0','      if (lrdonor) call donorresponsepermanent\n      dual=0.0d0\n      rhsdual=0.0d0')
s=replace(s,'      if (.not.lractive .or. lrm.eq.0) return\n      fac=0.5d0*electric/dielec',
    '      if (.not.lractive .or. lrm.eq.0) return\n      if (lrdonor) then\n         call donorresponsegradient\n         return\n      end if\n      fac=0.5d0*electric/dielec')
s+='\nc\n'+(HERE/'DONOR_RESPONSE.f').read_text()
f.write_text(s)

f=out/'source/testgrad.f';s=f.read_text()
s=replace(s,'      use iounit\n      use solpot','      use iounit\n      use polar, only: lrdonor\n      use solpot')
s=replace(s,'            call gradient (eval,derivs)','            call gradient (eval,derivs)\n            if (lrdonor)\n     &         call donorresponsetotalgradient (eval,derivs)')
f.write_text(s)

changed=[];patch=[]
for f in sorted((out/'source').iterdir()):
    if not f.is_file():continue
    orig=prior/'source'/f.name
    if f.read_bytes()!=orig.read_bytes():
        changed.append(f.name)
        patch.extend(difflib.unified_diff(orig.read_text().splitlines(True),f.read_text().splitlines(True),fromfile='a/source/'+f.name,tofile='b/source/'+f.name))
assert changed==['induce.f','polar.f','testgrad.f'],changed
(out/'DONOR_RESPONSE.patch').write_text(''.join(patch))
# Native fixed-format statements must fit the same 72-column convention.
for name in changed:
    for lineno,line in enumerate((out/'source'/name).read_text().splitlines(),1):
        if line and line[0] not in 'cC*!' and len(line)>72:
            raise ValueError(f'Overlong executable fixed-format line {name}:{lineno}')
result=dict(status='prepared_donor_response_source_uncompiled',source_root=str(out),baseline_source_manifest=sources['native_manifest'],
    design_sources=pin(SOURCES),design=sources['note'],declaration=pin(HERE/'DECLARATION.md'),helper=pin(HERE/'DONOR_RESPONSE.f'),
    implementation=pin(__file__),patch=pin(out/'DONOR_RESPONSE.patch'),changed_files=changed,
    source_outputs=[pin(f) for sub in ('source','cmake') for f in sorted((out/sub).rglob('*')) if f.is_file()],cmake=pin(out/'cmake/CMakeLists.txt'),
    donor_inventory=sources['donor_inventory'],water_native_actual_qualification=sources['water_native_actual_qualification'],
    native_keyword='POLAR-DONOR-RESPONSE metaltype 403 a c',allowed_pairs=[[500,403],[501,403]],
    coefficient_scope='Finite native inputs with arithmetic/SPD guards; execution limited by separate fixed six-state qualification manifest.',
    new_native_calls=0,builds=0,cost=dict(wall_seconds=time.monotonic()-tick,CPU_seconds=time.process_time()-cpu,threads=1))
print(json.dumps(save(out/'MANIFEST.json',result)))
