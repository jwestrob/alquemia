c     BIO-only radial drive; the vacuum namespace is unchanged.
      subroutine periodicradialsetup
      use, intrinsic :: ieee_arithmetic
      use iounit
      use keys
      use polar
      implicit none
      integer i,next,ia,ib,ios,slot
      real*8 e,nanval
      logical seen(2),bad
      character*240 record,keyword,string,extra
      seen=.false.
      nanval=ieee_value(0.0d0,ieee_quiet_nan)
      do i=1,nkey
         next=1
         record=keyline(i)
         call gettext(record,keyword,next)
         call upcase(keyword)
         if (trim(keyword).ne.'POLAR-RADIAL-DRIVE') cycle
         ia=-1
         ib=-1
         e=nanval
         string=record(next:240)
         read(string,*,iostat=ios) ia,ib,e
         bad=ios.ne.0 .or. (ia.ne.900.and.ia.ne.901)
     &       .or. ib.ne.402 .or. .not.ieee_is_finite(e)
         if (bad) then
            write(iout,*) 'PERIODIC RADIAL: invalid/mixed record'
            call fatal
         end if
         read(string,*,iostat=ios) ia,ib,e,extra
         if (ios.ge.0) then
            write(iout,*) 'PERIODIC RADIAL: extra record tokens'
            call fatal
         end if
         slot=ia-899
         if (seen(slot) .or. abs(e).gt.1.0d0) then
            write(iout,*) 'PERIODIC RADIAL: duplicate/domain error'
            call fatal
         end if
         seen(slot)=.true.
         lrpshape(slot)=e
         if (e.ne.0.0d0) then
            lractive=.true.
            write(iout,'(a,2i8,es26.17e3)')
     &       'PERIODIC_RADIAL_DRIVE_PARAMETERS ',ia,ib,e
         end if
      end do
      return
      end
c
c     Reuse the current compact geometry; no dense derivative storage.
      subroutine periodicradialfield(nw,ow,rr,vec,ff,ph,slot)
      use polar
      implicit none
      integer nw,ow(nw),slot,i,oi
      real*8 rr(nw),vec(3,nw),ff(nw),ph(nw),e,psi
      e=lrpshape(slot)
      if (e.eq.0.0d0) return
      do i=1,nw
c        Skip outside support before forming the polynomial shape.
         if (rr(i).ge.21.0d0) cycle
         oi=ow(i)
         psi=(rr(i)*rr(i)-6.25d0)/6.25d0
         lrh(:,oi)=lrh(:,oi)
     &      +3.0d0*e*ff(i)*ph(i)*psi*vec(:,i)
      end do
      return
      end
c
c     Add -J_h^T(u+v) before the parent's electric/(2*dielec).
c     The caller adds opposite metal reaction and image virial once.
      subroutine periodicradialgradient(i,e,w,gg)
      use polar
      implicit none
      integer i
      real*8 e,w(3),gg(3),psi,dotw,phioverr
      if (e.eq.0.0d0 .or. lrr(i).ge.21.0d0) return
      psi=(lrr(i)*lrr(i)-6.25d0)/6.25d0
      dotw=dot_product(lrvec(:,i),w)
      phioverr=0.0d0
      if (lrr(i).gt.0.0d0) phioverr=lrphir(i)/lrr(i)
      gg=gg-3.0d0*e*(lrphi(i)*psi*dotw*lrdf(:,i)
     &   +lrf(i)*(lrphi(i)*psi*w
     &   +(psi*phioverr+2.0d0*lrphi(i)/6.25d0)
     &   *dotw*lrvec(:,i)))
      return
      end
