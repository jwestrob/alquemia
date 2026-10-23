! Isolated qualification program. It does not change the library or Hamiltonian.
program pme_harness
   use atoms
   use bound
   use boxes
   use ewald
   use fft
   use mpole
   use pme
   use polar
   use polpot
   use ieee_arithmetic
   implicit none
   integer :: i,j,k,b,slot,unit,ios,mode,col,other,countapp
   real*8 :: energy,a(3,3),out1(3,6),out2(3,6),in1(3,6),in2(3,6)
   real*8 :: grad(3,6),old_u(3,6),old_v(3,6)
   real*8, allocatable :: old_q(:,:,:),old_t1(:,:,:),old_t2(:,:,:),old_t3(:,:,:)
   logical :: exists,on
   character*16 :: arg

   call initial
   call getxyz
   call mechanic
   call nextarg(arg,exists)
   if (.not.exists) error stop 'Expected OFF or ON'
   on = trim(arg)=='ON'
   if (.not.on .and. trim(arg)/='OFF') error stop 'Unknown harness mode'
   if (n/=6 .or. npole/=6) error stop 'Exact consumed six-site fixture only'
   if (lractive .or. lrperiodic) error stop 'Response must be inactive'
   if (.not.use_bounds .or. .not.orthogonal .or. nonprism) error stop 'Orthorhombic only'
   if (boundary/='TINFOIL' .or. poltyp/='MUTUAL') error stop 'Boundary/polarization'
   if (xbox/=44d0 .or. ybox/=44d0 .or. zbox/=44d0) error stop 'Cell changed'
   if (any(.not.douind(1:6)) .or. any(polarity(1:6)<1d-8)) error stop 'Active sites'
   if (any(ipole(1:6)/=[1,2,3,4,5,6])) error stop 'Site order changed'
   call gradient(energy,grad)
   if (.not.ieee_is_finite(energy) .or. any(.not.ieee_is_finite(grad))) error stop 'Nonfinite EF'
   open(newunit=unit,file='harness_ef.dat',status='new',action='write',iostat=ios)
   if (ios/=0) error stop 'EF file collision/open failure'
   write(unit,'(a)') 'PME_HARNESS_EF_V1'
   write(unit,'(a,es26.17e3)') 'ENERGY ',energy
   do i=1,6
      write(unit,'(a,i3,3(1x,es26.17e3))') 'GRADIENT ',i,grad(:,i)
   end do
   write(unit,'(a)') 'END 6'
   close(unit,iostat=ios)
   if (ios/=0) error stop 'EF close failure'
   inquire(file='pme_norm_export_001.dat',exist=exists)
   if (exists.neqv.on) error stop 'Exporter off/on state mismatch'
   if (.not.on) then
      call final
      stop
   end if
   ! These are the current initialized arrays used below; the collector must
   ! compare every token to the earlier export from the same gradient call.
   if (nfft1/=32.or.nfft2/=32.or.nfft3/=32.or.bsorder/=5) error stop 'Grid/order'
   if (aewald/=0.3d0 .or. ffttyp/='FFTW') error stop 'PME settings changed'
   do i=1,3
      a(1,i)=dble(nfft1)*recip(i,1)
      a(2,i)=dble(nfft2)*recip(i,2)
      a(3,i)=dble(nfft3)*recip(i,3)
   end do
   allocate(old_q,source=qfac)
   allocate(old_t1,source=thetai1)
   allocate(old_t2,source=thetai2)
   allocate(old_t3,source=thetai3)
   old_u=uind(:,1:6)
   old_v=uinp(:,1:6)
   open(newunit=unit,file='harness_state.dat',status='new',action='write',iostat=ios)
   if (ios/=0) error stop 'State file collision/open failure'
   write(unit,'(a)') 'PME_HARNESS_STATE_V1'
   write(unit,'(a,4i6,2(1x,es26.17e3))') 'STATE ',nfft1,nfft2,nfft3,bsorder,aewald,volbox
   do i=1,3
      write(unit,'(a,i3,3(1x,es26.17e3))') 'TRANS ',i,a(i,:)
   end do
   do i=1,6
      write(unit,'(a,5i6,4(1x,es26.17e3))') 'SITE ',i,type(i),igrid(:,i),polarity(i),x(i),y(i),z(i)
      do slot=1,5
         write(unit,'(a,3i4,2(1x,es26.17e3))') 'THETA ',i,1,slot,thetai1(1:2,slot,i)
         write(unit,'(a,3i4,2(1x,es26.17e3))') 'THETA ',i,2,slot,thetai2(1:2,slot,i)
         write(unit,'(a,3i4,2(1x,es26.17e3))') 'THETA ',i,3,slot,thetai3(1:2,slot,i)
      end do
   end do
   do k=1,nfft3
      do j=1,nfft2
         do i=1,nfft1
            write(unit,'(a,3i5,1x,es26.17e3)') 'QFAC ',i,j,k,qfac(i,j,k)
         end do
      end do
   end do
   write(unit,'(a)') 'END 6 90 32768'
   close(unit,iostat=ios)
   if (ios/=0) error stop 'State close failure'
   open(newunit=unit,file='harness_apply.dat',status='new',action='write',iostat=ios)
   if (ios/=0) error stop 'Application file collision/open failure'
   write(unit,'(a)') 'PME_HARNESS_APPLY_V1'
   countapp=0
   ! Mode 1: u basis, v zero. Mode 2: u zero, v basis.
   ! Mode 3: u basis column j, v cyclic next column. Both are unequal.
   do mode=1,3
      do col=1,18
         in1=0d0
         in2=0d0
         i=(col-1)/3+1
         b=mod(col-1,3)+1
         if (mode/=2) in1(b,i)=sqrt(polarity(i))
         if (mode==2) in2(b,i)=sqrt(polarity(i))
         if (mode==3) then
            other=mod(col,18)+1
            i=(other-1)/3+1
            b=mod(other-1,3)+1
            in2(b,i)=sqrt(polarity(i))
         end if
         call reciprocal_apply(in1,in2,out1,out2,a)
         countapp=countapp+1
         if (any(.not.ieee_is_finite(out1)).or.any(.not.ieee_is_finite(out2))) error stop 'Nonfinite apply'
         do i=1,6
            write(unit,'(a,3i4,6(1x,es26.17e3))') 'APPLY ',mode,col,i,out1(:,i),out2(:,i)
         end do
      end do
   end do
   if (any(qfac/=old_q).or.any(thetai1/=old_t1).or.any(thetai2/=old_t2).or. &
       any(thetai3/=old_t3).or.any(uind(:,1:6)/=old_u).or.any(uinp(:,1:6)/=old_v)) &
       error stop 'Persistent initialized state changed'
   write(unit,'(a,i4)') 'END ',countapp
   close(unit,iostat=ios)
   if (ios/=0) error stop 'Application close failure'
   call final
contains
   subroutine reciprocal_apply(v1,v2,r1,r2,transform)
      real*8,intent(in) :: v1(3,6),v2(3,6),transform(3,3)
      real*8,intent(out) :: r1(3,6),r2(3,6)
      real*8 :: fu1(3,6),fu2(3,6),phi1(10,6),phi2(10,6),phis(20,6)
      integer :: ia,ix,iy,iz
      ! Same native chain and Cartesian/fractional maps as umutual1.
      fu1=matmul(transform,v1)
      fu2=matmul(transform,v2)
      call grid_uind(fu1,fu2)
      call fftfront
      do iz=1,nfft3
         do iy=1,nfft2
            do ix=1,nfft1
               qgrid(:,ix,iy,iz)=qfac(ix,iy,iz)*qgrid(:,ix,iy,iz)
            end do
         end do
      end do
      call fftback
      call fphi_uind(phi1,phi2,phis)
      do ia=1,6
         r1(:,ia)=sqrt(polarity(ia))*matmul(transpose(transform),phi1(2:4,ia))
         r2(:,ia)=sqrt(polarity(ia))*matmul(transpose(transform),phi2(2:4,ia))
      end do
   end subroutine
end program
