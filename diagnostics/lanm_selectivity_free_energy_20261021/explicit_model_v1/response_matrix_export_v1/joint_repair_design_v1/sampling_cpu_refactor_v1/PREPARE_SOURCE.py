#!/usr/bin/env python3
"""Generate an isolated unchanged-physics CPU refactor. No compiler/native calls."""
from pathlib import Path
import difflib
import hashlib
import json
import re
import shutil
import time

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[5]
OUT=ROOT/'workspaces/lanm_selectivity_free_energy_20261021/explicit_model_v1/response_matrix_export_v1/joint_repair_design_v1/sampling_cpu_refactor_v1/source_v1'
def pin(p):
 p=Path(p).resolve();return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def replace(text,old,new,count=1):
 assert text.count(old)==count,(old[:100],text.count(old),count)
 return text.replace(old,new)
def main():
 tick=time.monotonic();cpu=time.process_time();inputs=json.loads((HERE/'SOURCE_INPUTS.json').read_text())
 for p in inputs.values():assert pin(p['path'])==p
 parent=json.loads(Path(inputs['baseline_source_manifest']['path']).read_text());baseline=Path(parent['source_root'])
 for p in parent['source_outputs']:assert pin(p['path'])==p
 OUT.mkdir(parents=True,exist_ok=False)
 for name in ('source','cmake'):shutil.copytree(baseline/name,OUT/name)
 p=OUT/'source/polar.f';text=p.read_text()
 text=replace(text,'      real*8, allocatable :: lrdh(:,:,:,:)\n','''c     Selected-pair geometry for O(Q) conservative contraction.
      integer :: lrnw = 0
      integer, allocatable :: lrow(:)
      real*8, allocatable :: lrf(:),lros(:),lreta(:),lrr(:)
      real*8, allocatable :: lrphi(:),lrphir(:)
      real*8, allocatable :: lrvec(:,:),lrdf(:,:)
''');p.write_text(text)
 p=OUT/'source/induce.f';text=p.read_text()
 # One common post-dispatch hook; the original field kernels are unmodified
 # apart from removing their old opt-in additions. External field ordering stays.
 text=replace(text,'''      if (lractive) then
         field(:,1:n) = field(:,1:n) + lrh
         fieldp(:,1:n) = fieldp(:,1:n) + lrh
      end if
''','')
 text=replace(text,'''      if (lractive) then
         do i = 1, n
            field(:,i) = field(:,i) - lrbeta(i)*uind(:,i)
            fieldp(:,i) = fieldp(:,i) - lrbeta(i)*uinp(:,i)
         end do
      end if
''','')
 # The four solver dispatches retain their original branch ordering.
 pattern=r'(?m)^( +)if \(use_ewald\) then\n\1   call ([du])field0c \(field,fieldp\)\n\1else if \(use_mlist\) then\n\1   call \2field0b \(field,fieldp\)\n\1else\n\1   call \2field0a \(field,fieldp\)\n\1end if'
 text,n=re.subn(pattern,lambda m:m[1]+'call '+m[2]+'field0 (field,fieldp)',text);assert n==4,n
 # Exporter and SPD diagnostics must see the same post-dispatch operator.
 text=replace(text,'      call dfield0a (field,fieldp)','      call dfield0 (field,fieldp)')
 text=replace(text,'            call ufield0a (field,fieldp)','            call ufield0 (field,fieldp)')
 text=replace(text,'      call ufield0a (field,fieldp)','      call ufield0 (field,fieldp)')
 text=replace(text,'            call ufield0a(fld,fldp)','            call ufield0(fld,fldp)')
 wrappers='''c     Common complete native fields. Enabled mode guards remain in setup.
c     Static external fields are still added by the original caller after h.
      subroutine dfield0 (field,fieldp)
      use atoms
      use limits
      use polar
      implicit none
      real*8 field(3,n),fieldp(3,n)
      if (use_ewald) then
         call dfield0c (field,fieldp)
      else if (use_mlist) then
         call dfield0b (field,fieldp)
      else
         call dfield0a (field,fieldp)
      end if
      if (lractive) then
         field(:,1:n)=field(:,1:n)+lrh
         fieldp(:,1:n)=fieldp(:,1:n)+lrh
      end if
      return
      end
c
      subroutine ufield0 (field,fieldp)
      use atoms
      use limits
      use polar
      implicit none
      integer i
      real*8 field(3,n),fieldp(3,n)
      if (use_ewald) then
         call ufield0c (field,fieldp)
      else if (use_mlist) then
         call ufield0b (field,fieldp)
      else
         call ufield0a (field,fieldp)
      end if
      if (lractive) then
         do i=1,n
            field(:,i)=field(:,i)-lrbeta(i)*uind(:,i)
            fieldp(:,i)=fieldp(:,i)-lrbeta(i)*uinp(:,i)
         end do
      end if
      return
      end
c
'''
 marker='c     Regularized local response: isolated original CPU MUTUAL route.\n'
 text=replace(text,marker,wrappers+marker)
 text=replace(text,'      lrm = 0\n','      lrm = 0\n      lrnw = 0\n')
 text=replace(text,'''      if (allocated(lrh)) deallocate(lrh,lrbeta,lrdb,lrdh)
      allocate(lrh(3,n),lrbeta(n),lrdb(3,n),lrdh(3,3,n,n))
      lrh=0.0d0
      lrbeta=0.0d0
      lrdb=0.0d0
      lrdh=0.0d0
''','''      if (allocated(lrh)) then
         deallocate(lrh,lrbeta,lrdb,lrow,lrf,lros,lreta,lrr)
         deallocate(lrphi,lrphir,lrvec,lrdf)
      end if
      allocate(lrh(3,n),lrbeta(n),lrdb(3,n),lrow(n))
      allocate(lrf(n),lros(n),lreta(n),lrr(n),lrphi(n),lrphir(n))
      allocate(lrvec(3,n),lrdf(3,n))
      lrh=0.0d0
      lrbeta=0.0d0
      lrdb=0.0d0
      lrow=0
      lrf=0.0d0
      lros=0.0d0
      lreta=0.0d0
      lrr=0.0d0
      lrphi=0.0d0
      lrphir=0.0d0
      lrvec=0.0d0
      lrdf=0.0d0
''')
 text=replace(text,'      integer i,j,k,l,m,slot,nw,oi,oj,ow(28)\n','      integer i,j,m,slot,nw,oi,oj,ow(28)\n')
 text=replace(text,'      real*8 rr(28),vec(3,28),df(3,28),ds(28),jac(3,3)\n      real*8 pref,w,wp,bij\n','      real*8 rr(28),vec(3,28),df(3,28),ds(28)\n      real*8 pref,prefix(0:28),suffix(1:29)\n')
 text=replace(text,'''      s=sum(ff(1:nw))
      q=0.0d0
      do i=1,nw
         os(i)=0.0d0
         do j=1,nw
            if (j.ne.i) os(i)=os(i)+ff(j)
            if (j.gt.i) q=q+2.0d0*ff(i)*ff(j)
         end do
         et(i)=os(i)/(1.0d0+os(i))
      end do
''','''c     O(Q) nonnegative environment sums; no total-minus-self loss.
      prefix(0)=0.0d0
      q=0.0d0
      do i=1,nw
         q=q+2.0d0*ff(i)*prefix(i-1)
         prefix(i)=prefix(i-1)+ff(i)
      end do
      suffix(nw+1)=0.0d0
      do i=nw,1,-1
         suffix(i)=suffix(i+1)+ff(i)
      end do
      s=prefix(nw)
      do i=1,nw
         os(i)=prefix(i-1)+suffix(i+1)
         et(i)=os(i)/(1.0d0+os(i))
      end do
''')
 start=text.index('      do i=1,nw\n         oi=ow(i)\n         pref=a+b*et(i)')
 end=text.index('      if (.not.all(ieee_is_finite(lrh))',start)
 text=text[:start]+'''      lrnw=nw
      lrow(1:nw)=ow(1:nw)
      lrf(1:nw)=ff(1:nw)
      lros(1:nw)=os(1:nw)
      lreta(1:nw)=et(1:nw)
      lrr(1:nw)=rr(1:nw)
      lrphi(1:nw)=ph(1:nw)
      lrphir(1:nw)=pr(1:nw)
      lrvec(:,1:nw)=vec(:,1:nw)
      lrdf(:,1:nw)=df(:,1:nw)
      do i=1,nw
         oi=ow(i)
         pref=a+b*et(i)
         lrh(:,oi)=3.0d0*ff(i)*pref*ph(i)*vec(:,i)
      end do
'''+text[end:]
 text=replace(text,'     &    .not.all(ieee_is_finite(lrdh))) then\n','''     &    .not.all(ieee_is_finite(lrf)) .or.
     &    .not.all(ieee_is_finite(lros)) .or.
     &    .not.all(ieee_is_finite(lreta)) .or.
     &    .not.all(ieee_is_finite(lrr)) .or.
     &    .not.all(ieee_is_finite(lrphi)) .or.
     &    .not.all(ieee_is_finite(lrphir)) .or.
     &    .not.all(ieee_is_finite(lrvec)) .or.
     &    .not.all(ieee_is_finite(lrdf))) then
''')
 start=text.index('      subroutine localresponsegradient\n')
 assert text[start:].count('      end\n')==1
 text=text[:start]+'''      subroutine localresponsegradient
      use atoms
      use chgpot
      use deriv
      use polar
      use virial
      use units
      implicit none
      integer i,oi,k,l,slot
      real*8 fac,cross,gg(3),w(3),a,b,pref,dotw,ci,phioverr
      real*8 zi(28),term(28),prefix(0:28),suffix(1:29)
      if (.not.lractive .or. lrm.eq.0 .or. lrnw.eq.0) return
      fac=0.5d0*electric/dielec
      cross=dot_product(uind(:,lrm),uinp(:,lrm))
      slot=type(lrm)-499
      a=lrcoef(1,slot)
      b=lrcoef(2,slot)
c     First pass contracts each site's induced half-sum with h's
c     radial vector. Prefix/suffix excludes self without subtraction.
      prefix(0)=0.0d0
      do i=1,lrnw
         oi=lrow(i)
         w=uind(:,oi)+uinp(:,oi)
         zi(i)=3.0d0*lrphi(i)*dot_product(lrvec(:,i),w)
         term(i)=b*lrf(i)*zi(i)/(1.0d0+lros(i))**2
         prefix(i)=prefix(i-1)+term(i)
      end do
      suffix(lrnw+1)=0.0d0
      do i=lrnw,1,-1
         suffix(i)=suffix(i+1)+term(i)
      end do
c     Second pass includes every other-water chain derivative.
c     Add its exact metal reaction and symmetric native virial.
      do i=1,lrnw
         oi=lrow(i)
         w=uind(:,oi)+uinp(:,oi)
         pref=a+b*lreta(i)
         dotw=dot_product(lrvec(:,i),w)
         ci=pref*zi(i)+prefix(i-1)+suffix(i+1)
         phioverr=0.0d0
         if (lrr(i).gt.0.0d0) phioverr=lrphir(i)/lrr(i)
         gg=cross*lrdb(:,oi)-ci*lrdf(:,i)
     &      -3.0d0*lrf(i)*pref*(lrphi(i)*w
     &                         +phioverr*dotw*lrvec(:,i))
         gg=fac*gg
         dep(:,oi)=dep(:,oi)+gg
         dep(:,lrm)=dep(:,lrm)-gg
         do k=1,3
            do l=1,3
               vir(k,l)=vir(k,l)+0.5d0*(lrvec(k,i)*gg(l)
     &                                  +lrvec(l,i)*gg(k))
            end do
         end do
      end do
      return
      end
'''
 assert 'lrdh' not in text
 p.write_text(text)
 changed=[];patch=[]
 for old in sorted((baseline/'source').iterdir()):
  new=OUT/'source'/old.name
  if old.read_bytes()!=new.read_bytes():
   changed.append(old.name);patch.extend(difflib.unified_diff(old.read_text().splitlines(True),new.read_text().splitlines(True),fromfile='a/source/'+old.name,tofile='b/source/'+old.name))
 assert changed==['induce.f','polar.f'],changed
 (OUT/'SAMPLING_CPU_REFACTOR.patch').write_text(''.join(patch))
 outputs=[pin(p) for name in ('source','cmake') for p in sorted((OUT/name).iterdir()) if p.is_file()]
 manifest={'status':'prepared_vacuum_cpu_refactor_source_uncompiled','source_root':str(OUT),'baseline_source_manifest':inputs['baseline_source_manifest'],'baseline_actual_qualification':inputs['baseline_actual_qualification'],'source_inputs':pin(HERE/'SOURCE_INPUTS.json'),'declaration':pin(HERE/'DECLARATION.md'),'generator':pin(__file__),'changed_files':changed,'source_outputs':outputs,'cmake':pin(OUT/'cmake/CMakeLists.txt'),'patch':pin(OUT/'SAMPLING_CPU_REFACTOR.patch'),'new_native_calls':0,'builds':0,'enabled_scope_unchanged':True,'periodic_execution_enabled':False,'dense_SPD_certificate_retained':True,'cost':{'wall_seconds':time.monotonic()-tick,'CPU_seconds':time.process_time()-cpu,'threads':1}}
 (OUT/'MANIFEST.json').write_text(json.dumps(manifest,indent=2)+'\n');print(json.dumps({'manifest':pin(OUT/'MANIFEST.json'),'changed_files':changed,'source_outputs':len(outputs)}))
if __name__=='__main__':main()
