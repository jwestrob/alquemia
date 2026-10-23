c
c     Opt-in reporting only; original baseline kernels are unchanged.
c
      subroutine responsematrixflag (requested)
      use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
      use atoms
      use bound
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
      implicit none
      integer i,j,next,count,nmetal
      logical requested,bad
      logical seen(28)
      character*240 record,keyword
      requested = .false.
      count = 0
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         if (trim(keyword) .eq. 'RESPONSE-MATRIX-EXPORT') then
            count = count + 1
            if (len_trim(record(next:240)) .ne. 0) then
               write (iout,*) 'MATRIX EXPORT: no arguments allowed'
               call fatal
            end if
         end if
      end do
      if (count .eq. 0)  return
      if (count .ne. 1) then
         write (iout,*) 'MATRIX EXPORT: duplicate reporting keyword'
         call fatal
      end if
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         if (trim(keyword) .eq. 'POLAR-PAIR-RESPONSE') then
            write (iout,*) 'MATRIX EXPORT: original baseline only'
            call fatal
         end if
      end do
      bad = use_bounds .or. use_replica .or. use_polymer
     &      .or. use_ewald .or. use_mlist .or. use_ulist
     &      .or. use_solv .or. use_chgpen .or. use_expol
     &      .or. use_tholed .or. use_chgflx .or. use_exfreq
     &      .or. use_pred .or. use_ielscf
     &      .or. poltyp.ne.'MUTUAL' .or. .not.use_thole
     &      .or. .not.use_polar .or. mpolecut.lt.1.0d12
     &      .or. n.lt.1 .or. n.gt.28 .or. npole.ne.n
      if (bad) then
         write (iout,*) 'MATRIX EXPORT: unsupported native mode'
         call fatal
      end if
      seen = .false.
      nmetal = 0
      do j = 1, npole
         i = ipole(j)
         if (i.lt.1 .or. i.gt.n) then
            write (iout,*) 'MATRIX EXPORT: invalid site order'
            call fatal
         end if
         if (seen(i) .or. .not.douind(i)) then
            write (iout,*) 'MATRIX EXPORT: inactive/duplicate site'
            call fatal
         end if
         seen(i) = .true.
         if (type(i).eq.500 .or. type(i).eq.501)
     &      nmetal = nmetal + 1
         if ((type(i).ne.36 .and. type(i).ne.37 .and.
     &        type(i).ne.500 .and. type(i).ne.501) .or.
     &       .not.ieee_is_finite(polarity(i)) .or.
     &       polarity(i).lt.1.0d-8 .or.
     &       .not.ieee_is_finite(x(i)) .or.
     &       .not.ieee_is_finite(y(i)) .or.
     &       .not.ieee_is_finite(z(i))) then
            write (iout,*) 'MATRIX EXPORT: unsupported type/alpha'
            call fatal
         end if
      end do
      if (nmetal .gt. 1) then
         write (iout,*) 'MATRIX EXPORT: at most one metal'
         call fatal
      end if
      if (use_exfld) then
         if (.not.all(ieee_is_finite(texfld))) then
            write (iout,*) 'MATRIX EXPORT: nonfinite static field'
            call fatal
         end if
      end if
      requested = .true.
      return
      end
c
      subroutine exportresponsematrix
      use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
      use atomid
      use atoms
      use cell
      use chgpot
      use extfld
      use iounit
      use mpole
      use polar
      use polpot
      use shunt
      use units
      implicit none
      integer i,j,k,ii,jj,ic,ir,dim,iu,ios,oldcell
      integer freeunit
      integer, save :: sequence = 0
      real*8 oldshunt(18),applied(3)
      real*8 dualerr,directerr,symerr,operatorerr
      real*8 residuald,residualp,value,valuep
      real*8, allocatable :: savedu(:,:),savedv(:,:)
      real*8, allocatable :: savedd(:,:),savedp(:,:)
      real*8, allocatable :: field(:,:),fieldp(:,:)
      real*8, allocatable :: am(:,:),ed(:),ep(:),uv(:),vv(:)
      logical finite,restored
      character*40 filename
      sequence = sequence + 1
      if (sequence .gt. 2) then
         write (iout,*) 'MATRIX EXPORT: only two solves per program'
         call fatal
      end if
      dim = 3*npole
      allocate (savedu(3,n),savedv(3,n),savedd(3,n),savedp(3,n))
      allocate (field(3,n),fieldp(3,n),am(dim,dim))
      allocate (ed(dim),ep(dim),uv(dim),vv(dim))
      savedu = uind(:,1:n)
      savedv = uinp(:,1:n)
      savedd = udir(:,1:n)
      savedp = udirp(:,1:n)
      oldcell = ncell
      oldshunt = (/off,off2,cut,cut2,c0,c1,c2,c3,c4,c5,
     &             f0,f1,f2,f3,f4,f5,f6,f7/)
      applied = 0.0d0
      if (use_exfld)  applied = texfld
c     Direct fields are the actual solver RHS, not a new model.
      call dfield0a (field,fieldp)
      directerr = 0.0d0
      do ii = 1, npole
         i = ipole(ii)
         do j = 1, 3
            ir = 3*(ii-1)+j
            ed(ir) = savedd(j,i)/polarity(i)
            ep(ir) = savedp(j,i)/polarity(i)
            uv(ir) = savedu(j,i)
            vv(ir) = savedv(j,i)
            directerr = max(directerr,
     &                    abs(ed(ir)-field(j,i)-applied(j)),
     &                    abs(ep(ir)-fieldp(j,i)-applied(j)))
         end do
      end do
c     Exactly the PCG mutual operator, on unit basis dipoles.
      dualerr = 0.0d0
      do jj = 1, npole
         k = ipole(jj)
         do j = 1, 3
            ic = 3*(jj-1)+j
            uind(:,1:n) = 0.0d0
            uinp(:,1:n) = 0.0d0
            uind(j,k) = 1.0d0
            uinp(j,k) = 1.0d0
            call ufield0a (field,fieldp)
            do ii = 1, npole
               i = ipole(ii)
               do iu = 1, 3
                  ir = 3*(ii-1)+iu
                  am(ir,ic) = -field(iu,i)
                  if (ir .eq. ic)
     &               am(ir,ic) = am(ir,ic)+1.0d0/polarity(i)
                  dualerr = max(dualerr,
     &                          abs(field(iu,i)-fieldp(iu,i)))
               end do
            end do
         end do
      end do
c     Restore native dipoles before the independent operator closure.
      uind(:,1:n) = savedu
      uinp(:,1:n) = savedv
      call ufield0a (field,fieldp)
      symerr = 0.0d0
      do ir = 1, dim
         do ic = 1, dim
            symerr = max(symerr,abs(am(ir,ic)-am(ic,ir)))
         end do
      end do
      operatorerr = 0.0d0
      residuald = 0.0d0
      residualp = 0.0d0
      do ii = 1, npole
         i = ipole(ii)
         do j = 1, 3
            ir = 3*(ii-1)+j
            value = dot_product(am(ir,:),uv)
            valuep = dot_product(am(ir,:),vv)
            operatorerr = max(operatorerr,
     &         abs(value-(savedu(j,i)/polarity(i)-field(j,i))),
     &         abs(valuep-(savedv(j,i)/polarity(i)-fieldp(j,i))))
            residuald = max(residuald,abs(value-ed(ir)))
            residualp = max(residualp,abs(valuep-ep(ir)))
         end do
      end do
c     switch(MPOLE) also changes shunt and ncell via replica(off).
c     No other global state is written by these guarded kernels.
      uind(:,1:n) = savedu
      uinp(:,1:n) = savedv
      udir(:,1:n) = savedd
      udirp(:,1:n) = savedp
      off=oldshunt(1)
      off2=oldshunt(2)
      cut=oldshunt(3)
      cut2=oldshunt(4)
      c0=oldshunt(5)
      c1=oldshunt(6)
      c2=oldshunt(7)
      c3=oldshunt(8)
      c4=oldshunt(9)
      c5=oldshunt(10)
      f0=oldshunt(11)
      f1=oldshunt(12)
      f2=oldshunt(13)
      f3=oldshunt(14)
      f4=oldshunt(15)
      f5=oldshunt(16)
      f6=oldshunt(17)
      f7=oldshunt(18)
      ncell=oldcell
      restored = all(uind(:,1:n).eq.savedu) .and.
     &           all(uinp(:,1:n).eq.savedv) .and.
     &           all(udir(:,1:n).eq.savedd) .and.
     &           all(udirp(:,1:n).eq.savedp) .and.
     &           ncell.eq.oldcell .and.
     &           all((/off,off2,cut,cut2,c0,c1,c2,c3,c4,c5,
     &                 f0,f1,f2,f3,f4,f5,f6,f7/).eq.oldshunt)
      finite = all(ieee_is_finite(am)) .and.
     &         all(ieee_is_finite(ed)) .and.
     &         all(ieee_is_finite(ep)) .and.
     &         all(ieee_is_finite(uv)) .and.
     &         all(ieee_is_finite(vv))
      if (.not.finite .or. .not.restored) then
         write (iout,*) 'MATRIX EXPORT: invalid/restoration failure'
         call fatal
      end if
      write (filename,'(a,i3.3,a)')
     &                    'response_matrix_',sequence,'.dat'
      iu = freeunit()
      open (unit=iu,file=trim(filename),status='new',
     &      action='write',iostat=ios)
      if (ios .ne. 0) then
         write (iout,*) 'MATRIX EXPORT: refuse output overwrite'
         call fatal
      end if
      write (iu,'(a)') 'TINKER_RESPONSE_MATRIX_V1'
      write (iu,'(a,4i8)') 'COUNTS ',sequence,n,npole,dim
      write (iu,'(a)') 'UNITS A=Angstrom^-3 RHS=e/Angstrom^2 U=eA'
      write (iu,'(a,3es26.17e3)') 'CONSTANTS ',electric,dielec,debye
      write (iu,'(a,l2,3es26.17e3)') 'EXTERNAL ',use_exfld,applied
      write (iu,'(a,6es26.17e3)') 'CHECKS ',dualerr,directerr,
     &                    symerr,operatorerr,residuald,residualp
      write (iu,'(a,l2)') 'RESTORED ',restored
      do ii = 1, npole
         i = ipole(ii)
         write (iu,'(a,3i8,1x,a3,4es26.17e3)') 'SITE ',
     &                 ii,i,type(i),name(i),polarity(i),x(i),y(i),z(i)
      end do
      do ir = 1, dim
         write (iu,'(a,i8,4es26.17e3)') 'VECTOR ',
     &                                    ir,ed(ir),ep(ir),uv(ir),vv(ir)
         write (iu,'(a,i8,100es26.17e3)') 'A ',
     &                                    ir,(am(ir,j),j=1,dim)
      end do
      write (iu,'(a)') 'END_TINKER_RESPONSE_MATRIX_V1'
      close (iu)
      write (iout,*) 'Response Matrix Written To : ',trim(filename)
      deallocate (savedu,savedv,savedd,savedp,field,fieldp)
      deallocate (am,ed,ep,uv,vv)
      return
      end
