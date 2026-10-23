c
c     Isolated CPU prototype: pair-scoped original cubic Thole kernels.
c     mode=1 selects direct driving fields; mode=0 selects mutual T.
c     The same symmetric selector is used by fields and derivatives.
c
      subroutine pairresponse (i,k,mode,gamma,selected)
      use polar
      implicit none
      integer i,k,mode,it,kt
      real*8 gamma
      logical selected
      selected = .false.
      gamma = 0.0d0
      if (npairdm .eq. 0)  return
      if (.not.allocated(pairdm))  return
      it = jpolar(i)
      kt = jpolar(k)
      if (.not.pairdm(it,kt))  return
      selected = .true.
      if (mode .eq. 1) then
         gamma = pairdmd(it,kt)
      else
         gamma = pairdmm(it,kt)
      end if
      return
      end
c
      subroutine dampresponse (i,k,mode,rorder,r,dm,selected)
      use polar
      implicit none
      integer i,k,j,mode,rorder
      real*8 r,dm(*),gamma,damp,t,t2,t3,expt
      logical selected
      call pairresponse (i,k,mode,gamma,selected)
      if (.not.selected)  return
      do j = 1, rorder
         dm(j) = 1.0d0
      end do
      damp = pdamp(i) * pdamp(k)
      if (damp.ne.0.0d0 .and. gamma.ne.0.0d0) then
         t = gamma * (r/damp)**3
         if (t .lt. 50.0d0) then
            expt = exp(-t)
            dm(3) = 1.0d0 - expt
            dm(5) = 1.0d0 - expt*(1.0d0+t)
            if (rorder .ge. 7) then
               t2 = t * t
               dm(7) = 1.0d0 - expt*(1.0d0+t+0.6d0*t2)
               if (rorder .ge. 9) then
                  t3 = t * t2
                  dm(9) = 1.0d0 - expt*(1.0d0+t
     &                    +(18.0d0/35.0d0)*t2+(9.0d0/35.0d0)*t3)
               end if
            end if
         end if
      end if
      return
      end
c
c     rcN is grad(fN)*rrN, in the existing epolar1a convention.
c     Return all factors separately for each selected D or M kernel.
c
      subroutine derivresponse (i,k,mode,r,xr,yr,zr,
     &                  rr3,rr5,rr7,sc3,sc5,sc7,rc3,rc5,rc7)
      use polar
      implicit none
      integer i,k,j,mode
      real*8 r,xr,yr,zr,rr3,rr5,rr7,sc3,sc5,sc7
      real*8 rc3(3),rc5(3),rc7(3),dm(9)
      real*8 gamma,damp,t,expt,temp3,temp5,temp7,r2
      logical selected
      call dampresponse (i,k,mode,7,r,dm,selected)
      if (.not.selected)  return
      sc3 = dm(3)
      sc5 = dm(5)
      sc7 = dm(7)
      do j = 1, 3
         rc3(j) = 0.0d0
         rc5(j) = 0.0d0
         rc7(j) = 0.0d0
      end do
      call pairresponse (i,k,mode,gamma,selected)
      damp = pdamp(i) * pdamp(k)
      if (damp.ne.0.0d0 .and. gamma.ne.0.0d0) then
         t = gamma * (r/damp)**3
         if (t .lt. 50.0d0) then
            expt = exp(-t)
            r2 = r * r
            temp3 = t * expt * rr5
            temp5 = 3.0d0 * t / r2
            temp7 = (-1.0d0+3.0d0*t) / r2
            rc3(1) = xr * temp3
            rc3(2) = yr * temp3
            rc3(3) = zr * temp3
            do j = 1, 3
               rc5(j) = rc3(j) * temp5
               rc7(j) = rc5(j) * temp7
            end do
         end if
      end if
      return
      end
c
c     Called at induction entry, after all mechanic flags are assigned.
c     This does not activate a global damping mode. Unsupported modes
c     are rejected only when the new keyword is present.
c
      subroutine checkpairresponse
      use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
      use bound
      use extfld
      use iounit
      use limits
      use mplpot
      use mpole
      use polar
      use polpot
      use potent
      use solpot
      use uprior
      implicit none
      integer ii,i
      logical bad
      bad = use_bounds .or. use_replica .or. use_polymer
     &      .or. use_ewald .or. use_mlist .or. use_ulist
     &      .or. use_solv .or. use_chgpen .or. use_expol
     &      .or. use_tholed .or. use_chgflx .or. use_exfld
     &      .or. use_pred .or. poltyp.ne.'MUTUAL'
     &      .or. .not.use_thole .or. mpolecut.lt.1.0d12
      if (bad) then
         write (iout,*) 'PAIR RESPONSE: unsupported execution mode'
         call fatal
      end if
      do ii = 1, npole
         i = ipole(ii)
         if (.not.ieee_is_finite(polarity(i)) .or.
     &       polarity(i).le.0.0d0) then
            write (iout,*) 'PAIR RESPONSE: invalid site polarizability'
            call fatal
         end if
      end do
      return
      end
