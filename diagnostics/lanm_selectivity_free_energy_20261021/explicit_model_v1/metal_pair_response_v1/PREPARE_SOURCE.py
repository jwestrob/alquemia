#!/usr/bin/env python3
"""Make an isolated source-only prototype. No compiler or molecule is run."""
from pathlib import Path
import difflib
import hashlib
import json
import shutil

ROOT = Path(__file__).resolve().parents[4]
D = Path(__file__).resolve().parent
BASE = ROOT / 'workspaces/lanm_selectivity_free_energy_20261021/explicit_model_v1/tinker_source'
OUT = ROOT / 'workspaces/lanm_selectivity_free_energy_20261021/explicit_model_v1/metal_pair_response_v1/source_v3'

def replace(s, a, b, count=1):
    assert s.count(a) == count, (a[:100], s.count(a), count)
    return s.replace(a, b)

def pin(p):
    return {'path': str(p.resolve()), 'sha256': hashlib.sha256(p.read_bytes()).hexdigest()}

def main():
    assert not OUT.exists(), 'Immutable source version already exists'
    OUT.mkdir(parents=True)
    shutil.copytree(BASE / 'source', OUT / 'source')
    shutil.copytree(BASE / 'cmake', OUT / 'cmake')
    original = {n: (BASE / 'source' / n).read_text() for n in
                ['polar.f', 'kpolar.f', 'damping.f', 'induce.f', 'epolar1.f']}
    changed = dict(original)

    s = changed['polar.f']
    s = replace(s, '      integer npolar\n', '''      integer npolar
c     Pair-scoped CPU prototype; zero/absent preserves original paths.
      integer npairdm = 0
      logical, allocatable :: pairdm(:,:)
      real*8, allocatable :: pairdmd(:,:),pairdmm(:,:)
''')
    changed['polar.f'] = s

    s = changed['kpolar.f']
    s = replace(s, '      subroutine kpolar\n', '''      subroutine kpolar
      use, intrinsic :: ieee_arithmetic, only: ieee_is_finite,
     &                                        ieee_value,ieee_quiet_nan
''')
    s = replace(s, '      integer nlist,npg\n', '''      integer nlist,npg
      integer dmios,dmextra,dmslot,dmtype
      logical dmseen(2)
      real*8 dmdirect,dmmutual,dmgamma
      character*240 dmtail
''')
    s = replace(s, '      polprt = .false.\n', '''      polprt = .false.
      npairdm = 0
      dmseen = .false.
''')
    s = replace(s, '      allocate (thdval(nlist,nlist))\n', '''      allocate (thdval(nlist,nlist))
      if (allocated(pairdm))  deallocate (pairdm)
      if (allocated(pairdmd))  deallocate (pairdmd)
      if (allocated(pairdmm))  deallocate (pairdmm)
      allocate (pairdm(nlist,nlist))
      allocate (pairdmd(nlist,nlist),pairdmm(nlist,nlist))
      pairdm = .false.
      pairdmd = 0.0d0
      pairdmm = 0.0d0
''')
    s = replace(s, '  210 continue\nc\nc     perform deallocation of some local arrays\n', '''  210 continue
c
c     Opt-in fractional changes to original cubic Thole coefficients.
c     This prototype admits only the exact SI metal/water-O type pairs.
c     Missing types are allowed, so identical keys work on fragments.
c
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         if (trim(keyword) .eq. 'POLAR-PAIR-RESPONSE') then
            string = record(next:240)
            ia = 0
            ib = 0
            dmdirect = ieee_value(0.0d0,ieee_quiet_nan)
            dmmutual = ieee_value(0.0d0,ieee_quiet_nan)
            read (string,*,iostat=dmios) ia,ib,dmdirect,dmmutual
            if (dmios .ne. 0) then
               write (iout,*) 'PAIR RESPONSE: expected 4 arguments'
               call fatal
            end if
            dmtail = ' '
            read (string,*,iostat=dmextra)
     &                       ia,ib,dmdirect,dmmutual,dmtail
            if (dmextra .eq. 0) then
               write (iout,*) 'PAIR RESPONSE: extra argument'
               call fatal
            end if
            dmtype = max(ia,ib)
            if (min(ia,ib).ne.36 .or.
     &          (dmtype.ne.500 .and. dmtype.ne.501)) then
               write (iout,*) 'PAIR RESPONSE: unsupported type pair'
               call fatal
            end if
            dmslot = dmtype - 499
            if (dmseen(dmslot)) then
               write (iout,*) 'PAIR RESPONSE: duplicate type pair'
               call fatal
            end if
            dmseen(dmslot) = .true.
            if (.not.ieee_is_finite(dmdirect) .or.
     &          .not.ieee_is_finite(dmmutual)) then
               write (iout,*) 'PAIR RESPONSE: nonfinite increment'
               call fatal
            end if
            if (dmdirect.le.-1.0d0 .or. dmmutual.le.-1.0d0) then
               write (iout,*) 'PAIR RESPONSE: increment must exceed -1'
               call fatal
            end if
            if (dmdirect.ne.0.0d0 .or. dmmutual.ne.0.0d0)
     &         npairdm = npairdm + 1
            ii = rlist(ia)
            kk = rlist(ib)
            if (ii.ne.0 .and. kk.ne.0) then
               dmgamma = thlval(ii,kk)
               if (.not.ieee_is_finite(dmgamma) .or.
     &             dmgamma.le.0.0d0) then
                  write (iout,*) 'PAIR RESPONSE: invalid original gamma'
                  call fatal
               end if
               pairdmd(ii,kk) = dmgamma * (1.0d0+dmdirect)
               pairdmm(ii,kk) = dmgamma * (1.0d0+dmmutual)
               if (.not.ieee_is_finite(pairdmd(ii,kk)) .or.
     &             .not.ieee_is_finite(pairdmm(ii,kk)) .or.
     &             pairdmd(ii,kk).le.0.0d0 .or.
     &             pairdmm(ii,kk).le.0.0d0) then
                  write (iout,*) 'PAIR RESPONSE: invalid changed gamma'
                  call fatal
               end if
               pairdmd(kk,ii) = pairdmd(ii,kk)
               pairdmm(kk,ii) = pairdmm(ii,kk)
               pairdm(ii,kk) = dmdirect.ne.0.0d0 .or.
     &                         dmmutual.ne.0.0d0
               pairdm(kk,ii) = pairdm(ii,kk)
               write (iout,*) 'PAIR RESPONSE TYPES/GAMMA0/D/M:',
     &                        ia,ib,dmgamma,
     &                        pairdmd(ii,kk),pairdmm(ii,kk)
            else
               write (iout,*) 'PAIR RESPONSE ABSENT TYPES:',ia,ib
            end if
         end if
      end do
c
c     perform deallocation of some local arrays
''')
    changed['kpolar.f'] = s

    s = changed['damping.f']
    for name, mode in [('dampthole', 0), ('damptholed', 1)]:
        start = s.index('      subroutine ' + name + ' (')
        end = s.index('      end\n', start) + len('      end\n')
        chunk = s[start:end]
        chunk = replace(chunk, '      real*8 dmpik(*)\n', '''      real*8 dmpik(*)
      logical dmselected
''')
        chunk = replace(chunk, 'c\nc\nc     initialize the Thole damping factors to a value of one\n', f'''c
c
      call dampresponse (i,k,{mode},rorder,r,dmpik,dmselected)
      if (dmselected)  return
c
c     initialize the Thole damping factors to a value of one
''')
        s = s[:start] + chunk + s[end:]
    s += (D / 'PAIR_RESPONSE_HELPERS.f').read_text()
    changed['damping.f'] = s

    s = changed['induce.f']
    s = replace(s, 'c     choose the method for computation of induced dipoles\nc\n', '''c     choose the method for computation of induced dipoles
c
      if (npairdm .ne. 0)  call checkpairresponse
c
''')
    changed['induce.f'] = s

    s = changed['epolar1.f']
    start = s.index('      subroutine epolar1a\n')
    end = s.index('      subroutine epolar1b\n')
    chunk = s[start:end]
    chunk = replace(chunk, '      character*6 mode\n', '''      character*6 mode
      logical dmselected
''')
    # Only the nonperiodic double-loop block is admitted. Replica/list/
    # Ewald routes remain original and the runtime guard forbids them.
    marker = '''                  if (use_tholed) then
                     pgamma = thdval(it,kt)'''
    assert chunk.count(marker) == 2
    chunk = chunk.replace(marker, '''                  call pairresponse (i,k,1,pgamma,dmselected)
                  if (dmselected) then
                     call derivresponse (i,k,1,r,xr,yr,zr,
     &                      rr3,rr5,rr7,sc3,sc5,sc7,rc3,rc5,rc7)
                  else if (use_tholed) then
                     pgamma = thdval(it,kt)''', 1)
    marker = '''               if (use_tholed) then
                  sc3 = 1.0d0'''
    assert chunk.count(marker) == 2
    chunk = chunk.replace(marker, '''               call pairresponse (i,k,0,pgamma,dmselected)
               if (dmselected) then
                  call derivresponse (i,k,0,r,xr,yr,zr,
     &                   rr3,rr5,rr7,sc3,sc5,sc7,rc3,rc5,rc7)
               else if (use_tholed) then
                  sc3 = 1.0d0''', 1)
    s = s[:start] + chunk + s[end:]
    changed['epolar1.f'] = s

    diff = []
    for name, text in changed.items():
        (OUT / 'source' / name).write_text(text)
        diff.extend(difflib.unified_diff(original[name].splitlines(True),
                    text.splitlines(True), 'a/source/'+name,
                    'b/source/'+name))
    (OUT / 'PROTOTYPE.patch').write_text(''.join(diff))
    inputs = [pin(BASE/'source'/p.name) for p in sorted((OUT/'source').glob('*')) if p.is_file()]
    outputs = [pin(p) for p in sorted((OUT/'source').glob('*')) if p.is_file()]
    all_changes = [p.name for p in (OUT/'source').glob('*') if p.is_file()
                   and p.read_bytes() != (BASE/'source'/p.name).read_bytes()]
    assert sorted(all_changes) == sorted(changed)
    long_added = [line[1:] for line in diff if line.startswith('+') and
                  not line.startswith('+++') and len(line.rstrip('\n')) > 73
                  and line[1:2] not in ('c', 'C', '*')]
    assert not long_added, long_added
    result = {'status': 'source_prepared_uncompiled', 'implementation': pin(Path(__file__)),
              'helpers': pin(D/'PAIR_RESPONSE_HELPERS.f'), 'source_root': str(OUT),
              'baseline_root': str(BASE), 'changed_files': sorted(all_changes),
              'baseline_inputs': inputs, 'source_outputs': outputs,
              'cmake': pin(OUT/'cmake/CMakeLists.txt'), 'patch': pin(OUT/'PROTOTYPE.patch'),
              'new_native_calls': 0, 'compiler_calls': 0,
              'scientific_parameters_selected': False,
              'validation': {'only_declared_five_files_changed': True,
                             'added_fixed_form_code_within_72_columns': True},
              'error': None}
    (OUT/'MANIFEST.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({'manifest': pin(OUT/'MANIFEST.json'), 'patch': result['patch'],
                      'changed_files': result['changed_files']}))

if __name__ == '__main__':
    main()
