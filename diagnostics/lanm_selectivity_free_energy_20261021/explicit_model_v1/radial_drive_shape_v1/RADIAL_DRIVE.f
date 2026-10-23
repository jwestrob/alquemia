c     Disjoint radial shape; e=0 calls no added geometry arithmetic.
      subroutine radialdrivesetup
      use, intrinsic :: ieee_arithmetic
      use iounit
      use keys
      use polar
      implicit none
      integer i,next,ia,ib,ios,slot
      real*8 e,nanval
      logical seen(2),bad
      character*240 record,keyword,string,extra
      lrshape=0.0d0
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
         bad=ios.ne.0 .or. (ia.ne.500.and.ia.ne.501)
     &       .or. ib.ne.36 .or. .not.ieee_is_finite(e)
         if (bad) then
            write(iout,*) 'RADIAL DRIVE: invalid record'
            call fatal
         end if
         read(string,*,iostat=ios) ia,ib,e,extra
         if (ios.ge.0) then
            write(iout,*) 'RADIAL DRIVE: extra record tokens'
            call fatal
         end if
         slot=ia-499
         if (seen(slot) .or. abs(e).gt.1.0d0) then
            write(iout,*) 'RADIAL DRIVE: duplicate/domain error'
            call fatal
         end if
         seen(slot)=.true.
         lrshape(slot)=e
         if (e.ne.0.0d0) then
            lractive=.true.
            write(iout,'(a,2i8,es26.17e3)')
     &       'RADIAL_DRIVE_PARAMETERS ',ia,ib,e
         end if
      end do
      return
      end
c
      subroutine radialdrivegeometry(nw,ow,rr,vec,ff,fr,ph,pr,
     &                               m,slot)
      use iounit
      use polar
      implicit none
      integer nw,ow(28),m,slot,i,k,l,oi
      real*8 rr(28),vec(3,28),ff(28),fr(28),ph(28),pr(28)
      real*8 e,psi,b,bp,jac(3,3),dh(3)
      real*8,parameter :: r0sq=6.25d0
      e=lrshape(slot)
      if (e.eq.0.0d0) return
      do i=1,nw
         oi=ow(i)
         psi=(rr(i)*rr(i)-r0sq)/r0sq
         b=3.0d0*ff(i)*ph(i)*psi
         bp=3.0d0*((fr(i)*ph(i)+ff(i)*pr(i))*psi
     &      +2.0d0*rr(i)*ff(i)*ph(i)/r0sq)
         dh=e*b*vec(:,i)
         lrh(:,oi)=lrh(:,oi)+dh
         jac=0.0d0
         do k=1,3
            jac(k,k)=e*b
            if (rr(i).gt.0.0d0) then
               do l=1,3
                  jac(k,l)=jac(k,l)+e*bp*vec(k,i)*vec(l,i)
     &                                        /rr(i)
               end do
            end if
         end do
         lrdh(:,:,oi,oi)=lrdh(:,:,oi,oi)+jac
         lrdh(:,:,oi,m)=lrdh(:,:,oi,m)-jac
         write(iout,'(a,2i8,5es26.17e3)')
     &    'RADIAL_DRIVE_PAIR ',m,oi,e,psi,dh
      end do
      return
      end
