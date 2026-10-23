"""Generate isolated identity sources only; no compiler/native/runtime entry."""
import difflib, hashlib, json, time
from pathlib import Path

HERE=Path(__file__).resolve().parent

def pin(path):
    path=Path(path).resolve()
    return {'path':str(path),'sha256':hashlib.sha256(path.read_bytes()).hexdigest()}

def checked(p):
    if pin(p['path'])!=p:raise ValueError('Source pin mismatch '+p['path'])
    return Path(p['path']).read_text()

def patch(source, edits):
    out=source
    for old,new in edits:
        if out.count(old)!=1:raise ValueError('Nonunique source anchor '+repr(old))
        out=out.replace(old,new,1)
    restored=out
    for old,new in reversed(edits):
        if restored.count(new)!=1:raise ValueError('Nonunique inverse anchor')
        restored=restored.replace(new,old,1)
    if restored!=source:raise ValueError('Inverse source mismatch')
    return out

def sources(parents):
    writer=checked(parents['state_writer_actual_parent'])
    gate=checked(parents['cutoff12_gate_actual_parent'])
    writer_edits=[
        ('subroutine stage2_write_state(direct1,direct2)', 'subroutine periodic_write_state_v2(filename,direct1,direct2)'),
        ('   use polar\n','   use polar\n   use limits\n   use polpot\n'),
        ('   integer :: ios,i\n','   integer :: ios,i\n   character(*),intent(in) :: filename\n'),
        ("file='certificate_actual_state.dat'",'file=filename'),
        ('PERIODIC_PARITY_STATE_V1','PERIODIC_PARITY_STATE_V2'),
        ("      write(u,'(a)') 'MODE FFTW TINFOIL MUTUAL ALL_ACTIVE SIMPLE_LOOP'\n", "      write(u,'(a)') 'MODE FFTW TINFOIL MUTUAL ALL_ACTIVE SIMPLE_LOOP'\n"
         "      write(u,'(a,3(1x,es26.17e3))') 'CELL ',xbox,ybox,zbox\n"
         "      write(u,'(a,3(1x,es26.17e3))') 'CUTOFFS ',ewaldcut,vdwcut,vdwtaper\n"
         "      write(u,'(a,2(1x,es26.17e3))') 'EWALD_SPLIT ',aeewald,apewald\n"
         "      write(u,'(a,1x,es26.17e3,1x,i12)') 'INDUCTION ',poleps,politer\n"),
        ('      do ia=1,n\n', "      do slot=1,2\n"
         "         write(u,'(a,i8,1x,es26.17e3)') 'RADIAL_COEFFICIENT ',slot+899,lrpshape(slot)\n"
         "      end do\n      do ia=1,n\n"),
        ("write(u,'(a,3i12)') 'END ',n,n*15,nfft1*nfft2*nfft3", "write(u,'(a,4i12)') 'END ',n,n*15,nfft1*nfft2*nfft3,2")]
    radial_check="""   subroutine check_radial_identity
      if (any(.not.ieee_is_finite(lrpshape)).or.any(abs(lrpshape)>1d0)) &
         error stop 'Stage2 radial coefficient domain'
      if (any(lrpshape/=0d0).and..not.lractive) error stop 'Stage2 radial bypass forbidden'
   end subroutine

"""
    gate_edits=[
        ('STAGE2_BOUND_INPUT_V1','STAGE2_BOUND_INPUT_V2'),
        ("      read(u,'(a)',iostat=ios) role\n", "      read(u,'(a)',iostat=ios) line\n"
         "      if (ios/=0.or.trim(line)/='PERIODIC_PARITY_STATE_V2') error stop 'Stage2 state schema'\n"
         "      read(u,'(a)',iostat=ios) role\n"),
        ('   subroutine backend_check\n',radial_check+'   subroutine backend_check\n      call check_radial_identity\n'),
        ('   subroutine induce_enter\n','   subroutine induce_enter\n      call check_radial_identity\n'),
        ("      write(u,'(a,1x,a)') 'ROLE',trim(role)\n", "      call check_radial_identity\n"
         "      write(u,'(a)') 'STATE_SCHEMA PERIODIC_PARITY_STATE_V2'\n"
         "      do i=1,2\n"
         "         write(u,'(a,i8,1x,es26.17e3)') 'RADIAL_COEFFICIENT ',i+899,lrpshape(i)\n"
         "      end do\n      write(u,'(a,1x,a)') 'ROLE',trim(role)\n")]
    wrappers="""! Both isolated callers use the exact same serializer, raw Ed and raw Ep.
subroutine stage1_write_state_v2(direct1,direct2)
   use atoms, only: n
   implicit none
   real*8,intent(in) :: direct1(3,n),direct2(3,n)
   call periodic_write_state_v2('parity_state.dat',direct1,direct2)
end subroutine
subroutine stage2_write_state(direct1,direct2)
   use atoms, only: n
   implicit none
   real*8,intent(in) :: direct1(3,n),direct2(3,n)
   call periodic_write_state_v2('certificate_actual_state.dat',direct1,direct2)
end subroutine
"""
    outputs={'STATE_WRITER_V2.f90':patch(writer,writer_edits),'STATE_ENTRY_POINTS.f90':wrappers,
             'GATE_V2.f90':patch(gate,gate_edits)}
    edits={'writer':writer_edits,'gate':gate_edits}
    diffs=''.join(difflib.unified_diff(writer.splitlines(True),outputs['STATE_WRITER_V2.f90'].splitlines(True),fromfile='parent/STATE_WRITER.f90',tofile='STATE_WRITER_V2.f90'))
    diffs+=''.join(difflib.unified_diff(gate.splitlines(True),outputs['GATE_V2.f90'].splitlines(True),fromfile='parent/GATE.f90',tofile='GATE_V2.f90'))
    return outputs,edits,diffs

def main(output):
    wall=time.monotonic();cpu=time.process_time()
    declaration=json.loads((HERE/'ROOT_DECLARATION.json').read_text())
    for p in declaration['parents'].values():checked(p)
    outputs,edits,diffs=sources(declaration['parents'])
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False)
    for name,body in outputs.items():(out/name).write_text(body)
    (out/'PATCH.diff').write_text(diffs)
    (out/'INVERSE_EDITS.json').write_text(json.dumps(edits,indent=2)+'\n')
    manifest={'status':'prepared_periodic_identity_v2_source_only','implementation':pin(__file__),
        'declaration':pin(HERE/'ROOT_DECLARATION.json'),'parents':declaration['parents'],
        'outputs':{n:pin(out/n) for n in outputs},'patch':pin(out/'PATCH.diff'),'inverse_edits':pin(out/'INVERSE_EDITS.json'),
        'inverse_restores_both_parents':True,'normal_engine_files_changed':0,
        'cost':{'wall_seconds':time.monotonic()-wall,'CPU_seconds':time.process_time()-cpu,'threads':1},
        'counts':{'compilations':0,'native_calls':0,'FFT_replays':0,'solver_calls':0,'factorizations':0},
        'runnable':False}
    (out/'MANIFEST.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(json.dumps({'manifest':pin(out/'MANIFEST.json'),'cost':manifest['cost']}))
if __name__=='__main__':
    import sys
    main(sys.argv[1])
