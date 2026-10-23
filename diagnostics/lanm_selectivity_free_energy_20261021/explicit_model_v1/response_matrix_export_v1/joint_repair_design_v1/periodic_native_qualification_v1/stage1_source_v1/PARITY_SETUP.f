      subroutine parityresponsesetup (handled)
      use, intrinsic :: ieee_arithmetic
      use atomid
      use atoms
      use bound
      use boxes
      use couple
      use ewald
      use extfld
      use ielscf
      use iounit
      use keys
      use limits
      use mplpot
      use mpole
      use polar
      use polpot
      use potent
      use uprior
      use usage
      implicit none
      logical handled,bad,seen(2)
      logical, allocatable :: sites(:)
      integer i,j,k,next,ia,ib,slot,ios,extra,nmet,nwat
      integer nh,nna,ncl,oj,hj
      real*8 a,b,c,d,nanval,dom
      character*240 record,keyword,string
c     Discover only the new namespace; leave old parsing untouched.
      handled=.false.
      do i=1,nkey
         next=1
         record=keyline(i)
         call gettext(record,keyword,next)
         call upcase(keyword)
         if (trim(keyword).ne.'POLAR-LOCAL-RESPONSE') cycle
         ia=-1
         read(record(next:240),*,iostat=ios) ia
         if (ios.eq.0 .and. (ia.eq.900.or.ia.eq.901)) then
            handled=.true.
         end if
      end do
      if (.not.handled) return
      lractive=.false.
      lrperiodic=.false.
      lrcoef=0.0d0
      lrm=0
      lrnw=0
      lrpnwater=0
      seen=.false.
      nanval=ieee_value(0.0d0,ieee_quiet_nan)
      dom=log(4.0d0)
      do i=1,nkey
         next=1
         record=keyline(i)
         call gettext(record,keyword,next)
         call upcase(keyword)
         if (trim(keyword).ne.'POLAR-LOCAL-RESPONSE') cycle
         ia=-1
         ib=-1
         a=nanval
         b=nanval
         c=nanval
         d=nanval
         string=record(next:240)
         read(string,*,iostat=ios) ia,ib,a,b,c,d
         bad=ios.ne.0 .or. (ia.ne.900.and.ia.ne.901)
     &       .or. ib.ne.402
         bad=bad.or..not.all(ieee_is_finite((/a,b,c,d/)))
         if (bad) then
            write(iout,*) 'PERIODIC RESPONSE: invalid/mixed record'
            call fatal
         end if
         read(string,*,iostat=ios) ia,ib,a,b,c,d,extra
         if (ios.ge.0) then
            write(iout,*) 'PERIODIC RESPONSE: extra record tokens'
            call fatal
         end if
         slot=ia-899
         if (seen(slot) .or. abs(a).gt.1.0d0 .or.
     &       abs(b).gt.1.0d0 .or.
     &       abs(c)+abs(d).gt.0.5d0*dom+1.0d-15) then
            write(iout,*) 'PERIODIC RESPONSE: duplicate/domain error'
            call fatal
         end if
         seen(slot)=.true.
         lrcoef(:,slot)=(/a,b,c,d/)
         if (any(lrcoef(:,slot).ne.0.0d0)) lractive=.true.
      end do
c     Zero BIO coefficients add no restriction or arithmetic.
      if (.not.lractive) return
      do i=1,nkey
         next=1
         record=keyline(i)
         call gettext(record,keyword,next)
         call upcase(keyword)
         if (trim(keyword).eq.'POLAR-PAIR-RESPONSE' .or.
     &       trim(keyword).eq.'POLAR-DONOR-RESPONSE' .or.
     &       trim(keyword).eq.'RESPONSE-MATRIX-EXPORT') then
            write(iout,*) 'PERIODIC RESPONSE: unsupported overlay'
            call fatal
         end if
      end do
      bad=.not.use_bounds .or. .not.use_ewald
     &    .or. .not.orthogonal .or. monoclinic .or. triclinic
     &    .or. octahedron .or. dodecadron .or. nonprism
     &    .or. use_replica .or. use_polymer .or. use_ulist
     &    .or. use_solv .or. use_chgpen .or. use_expol
     &    .or. use_tholed .or. use_chgflx .or. use_exfreq
     &    .or. use_pred .or. use_ielscf
     &    .or. poltyp.ne.'MUTUAL' .or. .not.use_thole
     &    .or. .not.use_polar .or. boundary.ne.'TINFOIL'
     &    .or. n.lt.1 .or. n.gt.11958 .or. npole.ne.n
     &    .or. nuse.ne.n
      if (.not.all(ieee_is_finite((/xbox,ybox,zbox/)))) bad=.true.
      if (min(xbox,ybox,zbox).lt.44.0d0) bad=.true.
      if (bad) then
         write(iout,*) 'PERIODIC RESPONSE: unsupported native mode'
         call fatal
      end if
      allocate(sites(n))
      sites=.false.
      nmet=0
      nwat=0
      nh=0
      nna=0
      ncl=0
      do j=1,npole
         i=ipole(j)
         if (i.lt.1 .or. i.gt.n) then
            write(iout,*) 'PERIODIC RESPONSE: invalid site index'
            call fatal
         end if
         if (sites(i) .or. .not.douind(i)) then
            write(iout,*) 'PERIODIC RESPONSE: inactive/duplicate site'
            call fatal
         end if
         sites(i)=.true.
         if (.not.ieee_is_finite(polarity(i)) .or.
     &       polarity(i).lt.1.0d-8 .or.
     &       .not.all(ieee_is_finite((/x(i),y(i),z(i)/)))) then
            write(iout,*) 'PERIODIC RESPONSE: invalid alpha/geometry'
            call fatal
         end if
         if (type(i).eq.900 .or. type(i).eq.901) then
            nmet=nmet+1
            lrm=i
            if (n12(i).ne.0 .or. rpole(1,i).ne.3.0d0 .or.
     &          (type(i).eq.900.and.atomic(i).ne.57) .or.
     &          (type(i).eq.901.and.atomic(i).ne.66)) then
               write(iout,*) 'PERIODIC RESPONSE: metal identity'
               call fatal
            end if
         else if (type(i).eq.402) then
            nwat=nwat+1
            if (n12(i).ne.2 .or. atomic(i).ne.8) then
               write(iout,*) 'PERIODIC RESPONSE: incomplete water O'
               call fatal
            end if
            if (i12(1,i).eq.i12(2,i)) then
               write(iout,*) 'PERIODIC RESPONSE: duplicate water H'
               call fatal
            end if
            do k=1,2
               hj=i12(k,i)
               if (hj.lt.1 .or. hj.gt.n) then
                  write(iout,*) 'PERIODIC RESPONSE: bad O bond index'
                  call fatal
               end if
               if (type(hj).ne.403 .or. n12(hj).ne.1) then
                  write(iout,*) 'PERIODIC RESPONSE: water O graph'
                  call fatal
               end if
               if (i12(1,hj).ne.i) then
                  write(iout,*) 'PERIODIC RESPONSE: O/H reciprocity'
                  call fatal
               end if
            end do
c           Exact paired charges preserve literal neutral water.
            if (rpole(1,i)+rpole(1,i12(1,i))+
     &          rpole(1,i12(2,i)).ne.0.0d0) then
               write(iout,*) 'PERIODIC RESPONSE: nonneutral water'
               call fatal
            end if
         else if (type(i).eq.403) then
            nh=nh+1
            if (n12(i).ne.1 .or. atomic(i).ne.1) then
               write(iout,*) 'PERIODIC RESPONSE: incomplete water H'
               call fatal
            end if
            oj=i12(1,i)
            if (oj.lt.1 .or. oj.gt.n) then
               write(iout,*) 'PERIODIC RESPONSE: bad H bond index'
               call fatal
            end if
            if (type(oj).ne.402 .or. n12(oj).ne.2) then
               write(iout,*) 'PERIODIC RESPONSE: water H graph'
               call fatal
            end if
            if (all(i12(1:2,oj).ne.i)) then
               write(iout,*) 'PERIODIC RESPONSE: H/O reciprocity'
               call fatal
            end if
         else if (type(i).eq.405 .or. type(i).eq.414) then
            if (n12(i).ne.0) then
               write(iout,*) 'PERIODIC RESPONSE: bonded salt'
               call fatal
            end if
            if (type(i).eq.405) then
               nna=nna+1
               if (atomic(i).ne.11 .or. rpole(1,i).ne.1.0d0)
     &            call fatal
            else
               ncl=ncl+1
               if (atomic(i).ne.17 .or. rpole(1,i).ne.-1.0d0)
     &            call fatal
            end if
         else
            write(iout,*) 'PERIODIC RESPONSE: unsupported atom type'
            call fatal
         end if
      end do
      deallocate(sites)
      if (nmet.ne.1 .or. nwat.gt.3980 .or. nh.ne.2*nwat .or.
     &    nna.ne.7 .or. ncl.ne.10) then
         write(iout,*) 'PERIODIC RESPONSE: population/neutrality scope'
         call fatal
      end if
      if (use_exfld) then
         if (.not.all(ieee_is_finite(texfld))) then
            write(iout,*) 'PERIODIC RESPONSE: nonfinite static field'
            call fatal
         end if
      end if
      lrperiodic=.true.
      lrpnwater=nwat
      if (allocated(lrh)) then
         deallocate(lrh,lrbeta,lrdb,lrow,lrf,lros,lreta,lrr)
         deallocate(lrphi,lrphir,lrvec,lrdf)
      end if
      allocate(lrh(3,n),lrbeta(n),lrdb(3,n),lrow(nwat))
      allocate(lrf(nwat),lros(nwat),lreta(nwat),lrr(nwat))
      allocate(lrphi(nwat),lrphir(nwat))
      allocate(lrvec(3,nwat),lrdf(3,nwat))
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
      do slot=1,2
         if (seen(slot)) then
            write(iout,'(a,2i8,4es26.17e3)')
     &       'LOCAL_RESPONSE_PARAMETERS ',slot+899,402,lrcoef(:,slot)
         end if
      end do
      call periodicresponsegeometry
c     Standalone explicit-operator initialization only.
c     Normal setup and both fatal guards remain unchanged.
      return
      end
c
