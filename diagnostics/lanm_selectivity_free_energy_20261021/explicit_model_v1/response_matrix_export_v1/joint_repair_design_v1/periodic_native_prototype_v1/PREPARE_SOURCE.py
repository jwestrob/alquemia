#!/usr/bin/env python3
"""Isolated source preparation; no compiler, native, box, or candidate API."""
from pathlib import Path
import difflib
import hashlib
import json
import shutil
import time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[5]
OUT = ROOT/'workspaces/lanm_selectivity_free_energy_20261021/explicit_model_v1/response_matrix_export_v1/joint_repair_design_v1/periodic_native_prototype_v1/source_v1'


def pin(p):
    p = Path(p).resolve()
    return {'path': str(p), 'sha256': hashlib.sha256(p.read_bytes()).hexdigest()}


def replace(text, old, new, count=1):
    assert text.count(old) == count, (old[:120], text.count(old), count)
    return text.replace(old, new)


def routine(text, name):
    start = text.index('      subroutine '+name+'\n')
    end = text.index('      end\n', start)+len('      end\n')
    return text[start:end]


def transform(baseline):
    changed = {name: (baseline/'source'/name).read_text()
               for name in ('polar.f', 'induce.f', 'epolar1.f', 'epolar3.f')}
    changed['polar.f'] = replace(changed['polar.f'],
        '      logical :: lractive = .false.\n',
        '      logical :: lractive = .false.\n'
        '      logical :: lrperiodic = .false.\n'
        '      integer :: lrpnwater = 0\n')
    original = changed['induce.f']
    text = original
    # BIO predispatch is separate; original parser/body is otherwise literal.
    text = replace(text, '      logical bad,seen(2),sites(28)\n',
        '      logical bad,seen(2),sites(28),periodic_handled\n')
    text = replace(text, '      lractive = .false.\n',
        '      lrperiodic = .false.\n'
        '      call periodicresponsesetup(periodic_handled)\n'
        '      if (periodic_handled) return\n'
        '      lractive = .false.\n')
    # Native direct PME preparation precedes this future certificate point.
    anchor = '''      call dfield0 (field,fieldp)
c
c     add external electric field to the direct field values
'''
    assert text.count(anchor) == 1
    text = replace(text, anchor, '''      call dfield0 (field,fieldp)
c     The source-only periodic path has no implemented stability backend.
      if (lrperiodic) call periodicresponsestability
c
c     add external electric field to the direct field values
''')
    # Separate dynamic periodic contraction; unchanged vacuum body follows.
    oldgradient = routine(original, 'localresponsegradient')
    newgradient = replace(oldgradient,
        '      if (.not.lractive .or. lrm.eq.0 .or. lrnw.eq.0) return\n',
        '      if (lrperiodic) then\n'
        '         call periodicresponsegradient\n'
        '         return\n'
        '      end if\n'
        '      if (.not.lractive .or. lrm.eq.0 .or. lrnw.eq.0) return\n')
    text = replace(text, oldgradient, newgradient)
    # Generate the separate geometry from the admitted algebra. Only storage,
    # namespace, compact radial kernel and nearest-image displacement differ.
    geo = routine(original, 'localresponsegeometry')
    geo = replace(geo, 'subroutine localresponsegeometry',
                       'subroutine periodicresponsegeometry')
    geo = replace(geo, '      use atoms\n',
                      '      use atoms\n      use boxes, only: xbox,ybox,zbox\n')
    geo = replace(geo, '      integer i,j,m,slot,nw,oi,oj,ow(28)\n',
        '      integer i,j,m,slot,nw,oi,oj,nalloc\n'
        '      integer, allocatable :: ow(:)\n')
    geo = replace(geo, '      real*8 a,b,c,d,r,t,p,fzero,s,q,eta,g,sigma,other\n',
        '      real*8 a,b,c,d,r,s,q,eta,g,sigma,other\n'
        '      real*8 sides(3)\n')
    geo = replace(geo, '''      real*8 ff(28),fr(28),ph(28),pr(28),et(28),os(28)
      real*8 rr(28),vec(3,28),df(3,28),ds(28)
      real*8 pref,prefix(0:28),suffix(1:29)
''', '''      real*8, allocatable :: ff(:),fr(:),ph(:),pr(:),et(:),os(:)
      real*8, allocatable :: rr(:),vec(:,:),df(:,:),ds(:)
      real*8, allocatable :: prefix(:),suffix(:)
      real*8 pref
''')
    geo = replace(geo, '      if (lrm.eq.0) return\n', '''      if (lrm.eq.0) return
      nalloc=lrpnwater
      allocate(ow(nalloc),ff(nalloc),fr(nalloc),ph(nalloc))
      allocate(pr(nalloc),et(nalloc),os(nalloc),rr(nalloc))
      allocate(vec(3,nalloc),df(3,nalloc),ds(nalloc))
      allocate(prefix(0:nalloc),suffix(1:nalloc+1))
      sides=(/xbox,ybox,zbox/)
''')
    geo = replace(geo, '      slot=type(m)-499\n',
                      '      slot=type(m)-899\n')
    geo = replace(geo, '''      fzero=1.0d0/(((1.0d0+5.0d0+25.0d0/3.0d0)
     &              *exp(-5.0d0))**2)
''', '')
    geo = replace(geo, '         if (type(oi).ne.36) cycle\n',
                      '         if (type(oi).ne.402) cycle\n')
    geo = replace(geo, '         nw=nw+1\n', '''         nw=nw+1
         if (nw.gt.nalloc) then
            write(iout,*) 'PERIODIC RESPONSE: water count changed'
            call fatal
         end if
''')
    geo = replace(geo,
        '         vec(:,i)=(/x(oi)-x(m),y(oi)-y(m),z(oi)-z(m)/)\n',
        '         vec(:,i)=(/x(oi)-x(m),y(oi)-y(m),z(oi)-z(m)/)\n'
        'c        Nearest orthogonal image; ties are outside support.\n'
        '         vec(:,i)=vec(:,i)-sides*anint(vec(:,i)/sides)\n')
    geo = replace(geo, '''         t=2.0d0*r
         p=1.0d0+t+t*t/3.0d0
         ff(i)=fzero*(p*exp(-t))**2
         fr(i)=-(4.0d0/3.0d0)*ff(i)*t*(1.0d0+t)/p
''', '         call periodicresponsekernel(r,ff(i),fr(i))\n')
    geo = replace(geo, '''         kappa=gamma/prod**3
         zeta=kappa*r**3
''', '''         kappa=gamma/prod**3
         if (.not.ieee_is_finite(kappa) .or. kappa.le.0.0d0) then
            write(iout,*) 'PERIODIC RESPONSE: invalid damping length'
            call fatal
         end if
c        Outside compact support all phi contractions vanish exactly.
         if (r.ge.21.0d0) then
            ph(i)=0.0d0
            pr(i)=0.0d0
            cycle
         end if
         zeta=kappa*r**3
         if (.not.ieee_is_finite(zeta)) then
            write(iout,*) 'PERIODIC RESPONSE: damping overflow'
            call fatal
         end if
''')
    # Automatic deallocation of local allocatables also covers nw=0 return.
    grad = replace(oldgradient, 'subroutine localresponsegradient',
                               'subroutine periodicresponsegradient')
    grad = replace(grad,
        '      real*8 zi(28),term(28),prefix(0:28),suffix(1:29)\n',
        '      real*8, allocatable :: zi(:),term(:),prefix(:),suffix(:)\n')
    grad = replace(grad,
        '      if (.not.lractive .or. lrm.eq.0 .or. lrnw.eq.0) return\n',
        '      if (.not.lractive .or. lrm.eq.0 .or. lrnw.eq.0) return\n'
        '      allocate(zi(lrnw),term(lrnw))\n'
        '      allocate(prefix(0:lrnw),suffix(1:lrnw+1))\n')
    grad = replace(grad, '      slot=type(lrm)-499\n',
                        '      slot=type(lrm)-899\n')
    text += '\nc\n'+(HERE/'PERIODIC_SETUP.f').read_text()+'c\n'+geo+'c\n'+grad
    changed['induce.f'] = text
    # One common gradient hook after all native derivatives and frame torques.
    text = changed['epolar1.f']
    text = replace(text, '      if (lractive) call localresponsegradient\n', '')
    text = replace(text, '      subroutine epolar1\n',
                         '      subroutine epolar1\n      use polar\n')
    text = replace(text, '''      if (use_expol) then
         call dexpol
      end if
      return
''', '''      if (use_expol) then
         call dexpol
      end if
      if (lractive) call localresponsegradient
      return
''')
    changed['epolar1.f'] = text
    # Pairwise E analysis uses only the added p-field; nonpairwise already has it.
    text = changed['epolar3.f']
    text = replace(text, '      subroutine epolar3\n',
                         '      subroutine epolar3\n      use polar\n')
    correction = '''c     Reduced stationary energy needs the added p-field once.
      if (lractive) then
         do i = 1, n
            e = -f*dot_product(uind(:,i),lrh(:,i))
            ep = ep + e
            einter = einter + e
            aep(i) = aep(i) + e
            if (e.ne.0.0d0) nep = nep + 1
         end do
      end if
'''
    text = replace(text, correction, '')
    text = replace(text, '''      else
         call epolar3e
      end if
      return
''', '''      else
         call epolar3e
      end if
      if (lractive.and.pairwise) call localresponseanalysis
      return
''')
    text += '''c
      subroutine localresponseanalysis
      use action
      use analyz
      use atoms
      use chgpot
      use energi
      use inter
      use polar
      use units
      implicit none
      integer i
      real*8 e,f
      f=0.5d0*electric/dielec
'''+correction+'''      return
      end
'''
    changed['epolar3.f'] = text
    return changed


def main():
    start = time.monotonic(); cpu = time.process_time()
    inputs = json.loads((HERE/'SOURCE_INPUTS.json').read_text())
    for p in inputs.values():
        assert pin(p['path']) == p
    parent = json.loads(Path(inputs['parent_manifest']['path']).read_text())
    baseline = Path(parent['source_root'])
    for p in parent['source_outputs']:
        assert pin(p['path']) == p
    changed = transform(baseline)
    assert not OUT.exists(), 'Source version already exists'
    OUT.mkdir(parents=True)
    for name in ('source', 'cmake'):
        shutil.copytree(baseline/name, OUT/name)
    patch = []
    for name, text in sorted(changed.items()):
        old = (baseline/'source'/name).read_text()
        (OUT/'source'/name).write_text(text)
        patch.extend(difflib.unified_diff(old.splitlines(True), text.splitlines(True),
                     fromfile='a/source/'+name, tofile='b/source/'+name))
    (OUT/'PERIODIC_NATIVE_PROTOTYPE.patch').write_text(''.join(patch))
    manifest = {
        'status': 'source_only_periodic_prototype_execution_disabled',
        'source_root': str(OUT), 'parent_manifest': inputs['parent_manifest'],
        'parent_actual_admission': inputs['parent_actual_admission'],
        'inputs': inputs, 'source_inputs': pin(HERE/'SOURCE_INPUTS.json'),
        'declaration': pin(HERE/'DECLARATION.md'),
        'generator': pin(__file__), 'helper': pin(HERE/'PERIODIC_SETUP.f'),
        'changed_files': sorted(changed),
        'source_outputs': [pin(p) for part in ('source', 'cmake')
                           for p in sorted((OUT/part).iterdir()) if p.is_file()],
        'patch': pin(OUT/'PERIODIC_NATIVE_PROTOTYPE.patch'),
        'new_native_calls': 0, 'builds': 0, 'periodic_execution_enabled': False,
        'unimplemented_fatal_hooks': ['periodicresponsebackendguard',
                                      'periodicresponsestability'],
        'candidate_bound': False, 'actual_box_evaluations': 0,
        'cost': {'wall_seconds': time.monotonic()-start,
                 'CPU_seconds': time.process_time()-cpu, 'threads': 1}}
    (OUT/'MANIFEST.json').write_text(json.dumps(manifest, indent=2)+'\n')
    print(json.dumps({'manifest': pin(OUT/'MANIFEST.json'),
                      'changed_files': sorted(changed)}))


if __name__ == '__main__':
    main()
