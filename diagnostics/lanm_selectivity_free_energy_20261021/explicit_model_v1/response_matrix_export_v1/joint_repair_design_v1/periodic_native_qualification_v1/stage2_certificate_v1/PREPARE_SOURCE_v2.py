#!/usr/bin/env python3
"""Create an isolated source clone and retain its exact diff; never compile/run."""
import ast
import datetime
import difflib
import hashlib
import json
from pathlib import Path
import time

HERE=Path(__file__).resolve().parent

def pin(path):
    path=Path(path).resolve()
    with path.open('rb') as handle: digest=hashlib.file_digest(handle,'sha256').hexdigest()
    return {'path':str(path),'sha256':digest}

def checked(ref):
    if pin(ref['path'])!=ref: raise ValueError('changed source '+ref['path'])
    return Path(ref['path'])

def write(path,text):
    with path.open('x') as handle: handle.write(text)

def replace_once(text,old,new,label):
    if text.count(old)!=1: raise ValueError('nonunique source join '+label)
    return text.replace(old,new,1)

def run():
    wall,cpu=time.monotonic(),time.process_time()
    inp=json.loads((HERE/'SOURCE_INPUTS_v2.json').read_text())
    original=checked(inp['induce']).read_text()
    text=original
    text=replace_once(text,'      call localresponsesetup\n','      call localresponsesetup\n      call stage2_induce_enter\n','one induce entry')
    text=replace_once(text,'      call periodicresponsebackendguard\n','      call stage2_backend_check\n','setup guard callsite')
    text=replace_once(text,'      if (lrperiodic) call periodicresponsestability\n','      if (lrperiodic) call stage2_certify\n','initialized pre-PCG guard')
    begin=text.index('      subroutine dfield0 (field,fieldp)\n')
    end=text.index('      subroutine ufield0 (field,fieldp)\n',begin)
    body=text[begin:end]
    body=replace_once(body,'      if (lractive) then\n','      call stage2_before_direct(field,fieldp)\n      if (lractive) then\n','raw RHS observation')
    body=replace_once(body,'      return\n','      call stage2_after_direct(field,fieldp)\n      return\n','corrected RHS observation')
    text=text[:begin]+body+text[end:]
    begin=text.index('      subroutine ufield0 (field,fieldp)\n')
    end=text.index('      subroutine localresponsesetup\n',begin)
    body=text[begin:end]
    body=replace_once(body,'      if (use_ewald) then\n','      call stage2_mutual_enter\n      if (use_ewald) then\n','mutual counter')
    text=text[:begin]+body+text[end:]
    begin=text.index('      subroutine induce0a\n')
    end=text.index('      subroutine dfield0a (field,fieldp)\n',begin)
    body=text[begin:end]
    marker='c     print the results from the conjugate gradient iteration\nc\n'
    body=replace_once(body,marker,marker+'         call stage2_solver_readback(iter,eps,epsd,epsp,epsold)\nc\n','existing residual readback')
    text=text[:begin]+body+text[end:]
    # Removing observations and restoring both callsites must reconstruct every
    # byte of the scientific parent object, including the original fatal bodies.
    undo=text
    for line in ('      call stage2_induce_enter\n','      call stage2_before_direct(field,fieldp)\n','      call stage2_after_direct(field,fieldp)\n','      call stage2_mutual_enter\n','         call stage2_solver_readback(iter,eps,epsd,epsp,epsold)\nc\n'):
        undo=replace_once(undo,line,'','inverse observation')
    undo=replace_once(undo,'      call stage2_backend_check\n','      call periodicresponsebackendguard\n','inverse backend')
    undo=replace_once(undo,'      if (lrperiodic) call stage2_certify\n','      if (lrperiodic) call periodicresponsestability\n','inverse stability')
    if undo!=original: raise ValueError('scientific parent did not reconstruct exactly')
    parity=checked(inp['stage1_harness']).read_text()
    start=parity.index('   subroutine write_state\n')
    end=parity.index('\n   subroutine reciprocal_apply',start)
    state=parity[start:end]
    state=replace_once(state,'   subroutine write_state\n','subroutine stage2_write_state(direct1,direct2)\n   use atoms\n   use boxes\n   use ewald\n   use fft\n   use mpole\n   use pme\n   use polar\n   implicit none\n   integer :: ios,i\n   real*8 :: transform(3,3),direct1(3,n),direct2(3,n)\n','separate serializer')
    state=replace_once(state,"      open(newunit=u,file='parity_state.dat'","      do i=1,3\n         transform(1,i)=dble(nfft1)*recip(i,1)\n         transform(2,i)=dble(nfft2)*recip(i,2)\n         transform(3,i)=dble(nfft3)*recip(i,3)\n      end do\n      open(newunit=u,file='certificate_actual_state.dat'",'output name only')
    output=Path(inp['source_output'])
    output.mkdir(parents=True,exist_ok=False)
    write(output/'induce_qualification.f',text)
    write(output/'STATE_WRITER.f90',state+'\n')
    write(output/'induce_qualification.diff',''.join(difflib.unified_diff(original.splitlines(True),text.splitlines(True),fromfile=inp['induce']['path'],tofile=str(output/'induce_qualification.f'))))
    for name in ('GATE.f90','GRADIENT.f90','INITIALIZE.f90'):
        write(output/name,(HERE/name).read_text())
    checks=[]
    for path in output.glob('*.f*'):
        if path.suffix not in ('.f','.f90'): continue
        width=72 if path.suffix=='.f' else 132
        if any(len(x)>width for x in path.read_text().splitlines() if not(path.suffix=='.f' and x[:1].lower()=='c')):
            raise ValueError('source line width '+path.name)
        checks.append('line width '+path.name)
    checks.append('inverse clone edits reproduce every original induce byte')
    for name in ('periodicresponsebackendguard','periodicresponsestability'):
        old=original.split('      subroutine '+name+'\n')[1].split('      end\n')[0]
        new=text.split('      subroutine '+name+'\n')[1].split('      end\n')[0]
        if new!=old or '      call fatal\n' not in new: raise ValueError('original guard changed')
        checks.append('literal original fatal body '+name)
    g=(output/'GRADIENT.f90').read_text().lower()
    init=(output/'INITIALIZE.f90').read_text().lower()
    if g.count('   call gradient(')!=1 or any('call '+x in init for x in ('gradient','induce','energy','ufield')):
        raise ValueError('native call inventory')
    checks.append('one gradient reporter; zero initializer solve/E/G/mutual calls')
    ast.parse(Path(__file__).read_text())
    result={'status':'prepared_stage2_source_only_unbuilt_bindings_missing','source_inputs':pin(HERE/'SOURCE_INPUTS_v2.json'),'producer':pin(__file__),'parent_induce':inp['induce'],'files':{p.name:pin(p) for p in sorted(output.iterdir()) if p.is_file()},'checks':checks,'scientific_parent_byte_reconstruction':True,'original_fatal_guard_bodies_unchanged':True,'programs_planned':{'initialization':32,'active_gradient':33,'bypass_gradient':3,'total_stage2':68,'new_extra_residual_operator_applications':0},'pending':{'stage1_actual_admission':None,'wrapper':None,'collector':None,'build':None,'root_execution_release':None},'actual':{'builds':0,'native_programs':0,'assemblies':0,'factorizations':0,'solves':0,'wall_seconds':time.monotonic()-wall,'CPU_seconds':time.process_time()-cpu,'threads':1},'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()}
    write(output/'MANIFEST.json',json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({'manifest':pin(output/'MANIFEST.json'),'checks':len(checks),'actual':result['actual']}))

if __name__=='__main__': run()
