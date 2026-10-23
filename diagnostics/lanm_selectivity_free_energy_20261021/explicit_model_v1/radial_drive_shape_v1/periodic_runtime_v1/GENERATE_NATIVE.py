"""Isolated qualification sources; no build or runtime invocation."""
import difflib,hashlib,json,time
from pathlib import Path
HERE=Path(__file__).resolve().parent

def pin(p):
 p=Path(p).resolve();return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def read(p):
 if pin(p['path'])!=p:raise ValueError('Changed source '+p['path'])
 return Path(p['path']).read_text()
def replace(s,a,b):
 if s.count(a)!=1:raise ValueError('Nonunique source anchor '+a[:80])
 return s.replace(a,b,1)
def changed(s,edits):
 original=s
 for a,b in edits:s=replace(s,a,b)
 undo=s
 for a,b in reversed(edits):undo=replace(undo,b,a)
 if undo!=original:raise ValueError('Inverse mismatch')
 return s

def derive(parents):
 manifest=json.loads(read(parents['radial_native_manifest']))
 # The radial manifest is an immutable complete source inventory.
 candidates=[]
 def visit(x):
  if isinstance(x,dict):
   if set(x)=={'path','sha256'} and x['path'].endswith('/source/induce.f'):candidates.append(x)
   else:
    for v in x.values():visit(v)
  elif isinstance(x,list):
   for v in x:visit(v)
 visit(manifest)
 candidates=[p for p in candidates if '/periodic_port_v1/source_v2/' in p['path']]
 if len({p['path'] for p in candidates})!=1:raise ValueError('Exact radial induce parent')
 induce_pin=candidates[0];original=read(induce_pin)
 start=original.index('      subroutine periodicresponsesetup (handled)\n');stop=original.index('      subroutine periodicresponsebackendguard\n',start)
 setup=original[start:stop]
 setup_edits=[('subroutine periodicresponsesetup (handled)','subroutine parityresponsesetup (handled)'),
 ('c     This source has no admission of a real periodic execution.\n      call periodicresponsebackendguard\n','c     Standalone explicit-operator initialization only.\nc     Normal setup and both fatal guards remain unchanged.\n')]
 parity_setup=changed(setup,setup_edits)
 edits=[('      call localresponsesetup\n','      call localresponsesetup\n      call stage2_induce_enter\n'),
 ('      call periodicresponsebackendguard\n','      call stage2_backend_check\n'),
 ('      if (lrperiodic) call periodicresponsestability\n','      if (lrperiodic) call stage2_certify\n')]
 qualified=changed(original,edits)
 # Keep each insertion local to its exact old routine; invert the whole routine.
 for begin,end,changes in [
  ('      subroutine dfield0 (field,fieldp)\n','      subroutine ufield0 (field,fieldp)\n',[
   ('      if (lractive) then\n','      call stage2_before_direct(field,fieldp)\n      if (lractive) then\n'),
   ('      return\n','      call stage2_after_direct(field,fieldp)\n      return\n')]),
  ('      subroutine ufield0 (field,fieldp)\n','      subroutine localresponsesetup\n',[
   ('      if (use_ewald) then\n','      call stage2_mutual_enter\n      if (use_ewald) then\n')]),
  ('      subroutine induce0a\n','      subroutine dfield0a (field,fieldp)\n',[
   ('c     print the results from the conjugate gradient iteration\nc\n','c     print the results from the conjugate gradient iteration\nc\n         call stage2_solver_readback(iter,eps,epsd,epsp,epsold)\nc\n')])]:
  a=qualified.index(begin);b=qualified.index(end,a);old=qualified[a:b];new=changed(old,changes)
  qualified=replace(qualified,old,new);edits.append((old,new))
 undo=qualified
 for a,b in reversed(edits):undo=replace(undo,b,a)
 if undo!=original:raise ValueError('Whole radial source inverse')
 for name in ('periodicresponsebackendguard','periodicresponsestability'):
  a=original.split('      subroutine '+name+'\n')[1].split('      end\n')[0]
  b=qualified.split('      subroutine '+name+'\n')[1].split('      end\n')[0]
  if a!=b or '      call fatal\n' not in b:raise ValueError('Normal fatal body changed')
 old=read(parents['stage1_harness']);parity_edits=[]
 def pe(a,b):parity_edits.append((a,b))
 pe('   logical :: handled,exists\n','   logical :: handled,exists\n   real*8 :: old_coef(4,2),old_e(2)\n')
 pe("   if (xbox/=49.324d0.or.ybox/=49.324d0.or.zbox/=49.324d0) error stop 'Frozen cell changed'", "   if (min(xbox,ybox,zbox)<44d0) error stop 'Declared orthogonal cell domain'")
 pe('if (ewaldcut/=9d0.or.vdwcut/=9d0.or.vdwtaper/=8.1d0)','if (ewaldcut/=12d0.or.vdwcut/=9d0.or.vdwtaper/=8.1d0)')
 pe("   if (apewald/=0.4d0.or.aeewald/=0.4d0) error stop 'Ewald exponent changed'\n", "   if (apewald/=0.4d0.or.aeewald/=0.4d0) error stop 'Ewald exponent changed'\n   if (poleps/=1d-8.or.politer/=100) error stop 'Induction controls changed'\n")
 a=old.index('   ! Exact input identity');b=old.index('   allocate(probes',a)
 pe(old[a:b],"""   ! Full V2 state is checked after native direct initialization, before probes.
   if (any(.not.ieee_is_finite(lrcoef)).or.any(.not.ieee_is_finite(lrpshape))) &
       error stop 'Nonfinite response coefficients'
   if (any(abs(lrpshape)>1d0)) error stop 'Radial coefficient domain'
   if (any(.not.ieee_is_finite(lrh)).or.any(.not.ieee_is_finite(lrbeta))) error stop 'Response arrays'
   if (any(lrbeta(2:n)/=0d0)) error stop 'H must occur once on primary metal'
   do i=1,n
      if (type(i)/=402.and.any(lrh(:,i)/=0d0)) error stop 'Unexpected h support'
   end do

""")
 pe('   call write_state\n',"   old_coef=lrcoef\n   old_e=lrpshape\n   call stage1_write_state_v2(direct1,direct2)\n   call compare_expected_v2\n")
 a=old.index('   subroutine write_state\n');b=old.index('\n   subroutine reciprocal_apply',a)
 pe(old[a:b],"""   subroutine compare_expected_v2
      integer :: actual,expected,ia,ib,lines
      character*512 :: a,b
      open(newunit=actual,file='parity_state.dat',status='old',action='read',iostat=ia)
      open(newunit=expected,file='expected_state.dat',status='old',action='read',iostat=ib)
      if (ia/=0.or.ib/=0) error stop 'Fresh expected state absent'
      lines=0
      do
         read(actual,'(a)',iostat=ia) a
         read(expected,'(a)',iostat=ib) b
         if (ia<0.and.ib<0) exit
         if (ia/=0.or.ib/=0.or.a/=b) error stop 'Exact V2 initialized state mismatch'
         lines=lines+1
         if (lines==1.and.trim(a)/='PERIODIC_PARITY_STATE_V2') error stop 'V2 state required'
      end do
      close(actual)
      close(expected)
      if (lines<262144) error stop 'Incomplete V2 state'
   end subroutine
""")
 pe('   subroutine check_state\n','   subroutine check_state\n      if (any(lrcoef/=old_coef).or.any(lrpshape/=old_e)) error stop \'Model identity changed\'\n')
 parity=changed(old,parity_edits)
 identity=json.loads(read(parents['identity_source_manifest']))
 files={'PARITY_V2.f90':parity,'PARITY_SETUP.f':parity_setup,'induce_qualification.f':qualified,
        'INITIALIZE.f90':read(parents['stage2_initializer_source'])}
 for name,p in identity['outputs'].items():files[name]=read(p)
 # The existing analytic reporter itself is not changed by V2 metadata.
 s2=Path(parents['stage2_initializer_source']['path']).parent
 reporter=s2/'GRADIENT.f90';files['GRADIENT.f90']=reporter.read_text()
 inverse={'parent_induce':induce_pin,'parent_reporter':pin(reporter),'setup':setup_edits,'qualification':edits,'parity':parity_edits}
 for name,text in files.items():
  width=72 if name.endswith('.f') else 132
  if any(len(line)>width for line in text.splitlines() if not(name.endswith('.f') and line[:1].lower()=='c')):raise ValueError('Fortran line width '+name)
 return files,inverse

def main(output):
 wall=time.monotonic();cpu=time.process_time();decl=json.loads((HERE/'ROOT_DECLARATION.json').read_text());parents=decl['parents']
 for p in parents.values():read(p)
 files,inverse=derive(parents);out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False)
 for n,t in files.items():(out/n).write_text(t)
 (out/'INVERSE.json').write_text(json.dumps(inverse,indent=2)+'\n')
 result={'status':'prepared_periodic_v2_runtime_native_source_only','implementation':pin(__file__),'declaration':pin(HERE/'ROOT_DECLARATION.json'),
 'files':{n:pin(out/n) for n in files},'inverse':pin(out/'INVERSE.json'),'normal_fatal_guards_unchanged':True,
 'counts':{'compilations':0,'native_calls':0,'FFT_replays':0,'assemblies':0,'factorizations':0,'solves':0},
 'cost':{'wall_seconds':time.monotonic()-wall,'CPU_seconds':time.process_time()-cpu,'threads':1}}
 (out/'MANIFEST.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({'manifest':pin(out/'MANIFEST.json'),'cost':result['cost']}))
if __name__=='__main__':
 import sys
 main(sys.argv[1])
