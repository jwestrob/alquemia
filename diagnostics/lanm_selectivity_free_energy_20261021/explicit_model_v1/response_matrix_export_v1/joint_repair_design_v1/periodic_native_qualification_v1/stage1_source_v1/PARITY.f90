! Standalone explicit-operator probe. No induced solve, energy or gradient.
! Both normal periodic fatal guards remain in the immutable linked library.
program periodic_parity
   use atoms
   use bound
   use boxes
   use ewald
   use extfld
   use fft
   use limits
   use mplpot
   use mpole
   use pme
   use polar
   use polpot
   use ieee_arithmetic
   implicit none
   integer :: i,j,k,unit,ios,mode,col,other,atom,typ,ns,ng,np,complete_count,reference_count
   real*8 :: transform(3,3),expected_xyz(3),expected_alpha,expected_beta,expected_coef(4,2)
   real*8,allocatable :: probes(:,:,:),v1(:,:),v2(:,:),field1(:,:),field2(:,:)
   real*8,allocatable :: full1(:,:),full2(:,:),r1(:,:),r2(:,:),direct1(:,:),direct2(:,:)
   real*8,allocatable :: old_q(:,:,:),old_grid(:,:,:,:),old_t1(:,:,:),old_t2(:,:,:),old_t3(:,:,:)
   real*8,allocatable :: old_h(:,:),old_beta(:),old_alpha(:),fu1(:,:),fu2(:,:)
   real*8,allocatable :: phi1(:,:),phi2(:,:),phis(:,:)
   integer,allocatable :: old_igrid(:,:)
   logical :: handled,exists
   character*80 :: header,label

   call initial
   call getxyz
   call mechanic
   if (n/=11958.or.npole/=n) error stop 'Exact La BIO site count required'
   if (count(type(1:n)==900)/=1.or.type(1)/=900.or.any(type(1:n)==901)) error stop 'Exact La site map'
   if (count(type(1:n)==402)/=3980.or.count(type(1:n)==403)/=7960) error stop 'Water inventory'
   if (count(type(1:n)==405)/=7.or.count(type(1:n)==414)/=10) error stop 'Salt inventory'
   if (any(ipole(1:n)/=[(i,i=1,n)])) error stop 'Atom-major active map required'
   if (.not.use_bounds.or..not.orthogonal.or.nonprism.or..not.use_ewald) error stop 'Orthogonal periodic Ewald only'
   if (use_mlist.or.use_ulist.or.use_replica.or.use_exfld) error stop 'Exact simple-loop no-field mode required'
   if (boundary/='TINFOIL'.or.ffttyp/='FFTW'.or.poltyp/='MUTUAL') error stop 'Native mode changed'
   if (xbox/=49.324d0.or.ybox/=49.324d0.or.zbox/=49.324d0) error stop 'Frozen cell changed'
   if (ewaldcut/=9d0.or.vdwcut/=9d0.or.vdwtaper/=8.1d0) error stop 'Frozen cutoff changed'
   if (nefft1/=64.or.nefft2/=64.or.nefft3/=64.or.bsporder/=5.or.bseorder/=5) error stop 'PME grid/order'
   if (apewald/=0.4d0.or.aeewald/=0.4d0) error stop 'Ewald exponent changed'
   if (any(.not.douind(1:n)).or.any(polarity(1:n)<1d-8)) error stop 'Active alpha domain'
   if (u1scale/=1d0.or.u2scale/=1d0.or.u3scale/=1d0.or.u4scale/=1d0) error stop 'Mutual scales changed'
   call chkpole
   call rotpole('MPOLE')
   call parityresponsesetup(handled)
   if (.not.handled.or..not.lractive.or..not.lrperiodic.or.lrm/=1.or.lrpnwater/=3980) error stop 'Setup scope'
   ! Exact input identity is separately SHA-bound by the wrapper; compare the
   ! actual in-memory atom/coordinate/alpha map before any operator application.
   open(newunit=unit,file='expected_state.dat',status='old',action='read',iostat=ios)
   if (ios/=0) error stop 'Missing expected state'
   read(unit,'(a)',iostat=ios) header
   if (ios/=0.or.trim(header)/='PERIODIC_PARITY_EXPECTED_V1') error stop 'Expected header'
   read(unit,*,iostat=ios) ns,expected_beta
   if (ios/=0.or.ns/=n.or..not.ieee_is_finite(expected_beta)) error stop 'Expected dimensions/beta'
   do j=1,2
      read(unit,*,iostat=ios) expected_coef(:,j)
      if (ios/=0) error stop 'Expected coefficients missing'
   end do
   if (any(lrcoef/=expected_coef).or.any(.not.ieee_is_finite(lrh)).or. &
       any(.not.ieee_is_finite(lrbeta))) error stop 'Native response coefficients/arrays'
   if (abs(lrbeta(1)*polarity(1)-expected_beta)>1d-11*max(1d0,abs(expected_beta))) &
       error stop 'Native H differs from admitted static matrix'
   if (any(lrbeta(2:n)/=0d0)) error stop 'H must occur once on primary metal'
   do i=1,n
      read(unit,*,iostat=ios) atom,typ,expected_xyz,expected_alpha
      if (ios/=0.or.atom/=i.or.typ/=type(i)) error stop 'Expected site map'
      if (any(expected_xyz/=[x(i),y(i),z(i)]).or.expected_alpha/=polarity(i)) error stop 'Expected xyz/alpha'
      if (type(i)/=402.and.any(lrh(:,i)/=0d0)) error stop 'Unexpected h support'
   end do
   read(unit,'(a)',iostat=ios) header
   if (ios/=0.or.trim(header)/='END') error stop 'Expected tail'
   read(unit,'(a)',iostat=ios) header
   if (ios>=0) error stop 'Unexpected trailing state records'
   close(unit)

   allocate(probes(3,n,10),v1(3,n),v2(3,n),field1(3,n),field2(3,n),full1(3,n),full2(3,n))
   allocate(r1(3,n),r2(3,n),direct1(3,n),direct2(3,n),fu1(3,n),fu2(3,n))
   allocate(phi1(10,n),phi2(10,n),phis(20,n))
   open(newunit=unit,file='probes.dat',status='old',action='read',iostat=ios)
   if (ios/=0) error stop 'Missing frozen probes'
   read(unit,'(a)',iostat=ios) header
   if (ios/=0.or.trim(header)/='PERIODIC_PARITY_PROBES_V1') error stop 'Probe header'
   read(unit,*,iostat=ios) ns,np
   if (ios/=0.or.ns/=n.or.np/=10) error stop 'Probe inventory'
   do col=1,10
      read(unit,*,iostat=ios) j,label
      if (ios/=0.or.j/=col) error stop 'Probe label/order'
      do i=1,n
         read(unit,*,iostat=ios) atom,probes(:,i,col)
         if (ios/=0.or.atom/=i.or.any(.not.ieee_is_finite(probes(:,i,col)))) error stop 'Probe vector'
      end do
   end do
   read(unit,'(a)',iostat=ios) header
   if (ios/=0.or.trim(header)/='END') error stop 'Probe tail'
   read(unit,'(a)',iostat=ios) header
   if (ios>=0) error stop 'Unexpected trailing probe records'
   close(unit)

   ! dfield0c initializes the same native qfac/splines and unequal direct RHS.
   ! No induced solver is called. The normal dfield0 wrapper adds h to both RHS.
   call dfield0c(direct1,direct2)
   if (any(.not.ieee_is_finite(direct1)).or.any(.not.ieee_is_finite(direct2))) error stop 'Direct fields'
   inquire(file='pme_norm_export_001.dat',exist=exists)
   if (.not.exists) error stop 'Initialized native snapshot missing'
   do i=1,3
      transform(1,i)=dble(nfft1)*recip(i,1)
      transform(2,i)=dble(nfft2)*recip(i,2)
      transform(3,i)=dble(nfft3)*recip(i,3)
   end do
   allocate(old_q,source=qfac)
   allocate(old_grid,source=qgrid)
   allocate(old_t1,source=thetai1)
   allocate(old_t2,source=thetai2)
   allocate(old_t3,source=thetai3)
   allocate(old_igrid,source=igrid)
   allocate(old_h,source=lrh)
   allocate(old_beta,source=lrbeta)
   allocate(old_alpha,source=polarity)
   ! kpolar allocates u/v but does not define their values. Establish a zero
   ! explicit-operator scratch state rather than reading undefined dipoles.
   uind=0d0
   uinp=0d0
   call write_state
   open(newunit=unit,file='parity_apply.dat',status='new',action='write',iostat=ios)
   if (ios/=0) error stop 'Application output collision'
   write(unit,'(a)') 'PERIODIC_PARITY_APPLY_V1'
   write(unit,'(a)') 'CHANNELS real=uind imaginary=uinp NORMALIZED_M=I-C*(T-H)*C RECIPROCAL=positive_R'
   complete_count=0
   reference_count=0
   do mode=1,3
      do col=1,10
         v1=0d0
         v2=0d0
         if (mode/=2) v1=probes(:,:,col)
         if (mode==2) v2=probes(:,:,col)
         other=mod(col,10)+1
         if (mode==3) v2=probes(:,:,other)
         do i=1,n
            uind(:,i)=sqrt(polarity(i))*v1(:,i)
            uinp(:,i)=sqrt(polarity(i))*v2(:,i)
         end do
         qgrid=old_grid
         call ufield0(field1,field2)
         complete_count=complete_count+1
         do i=1,n
            full1(:,i)=v1(:,i)-sqrt(polarity(i))*field1(:,i)
            full2(:,i)=v2(:,i)-sqrt(polarity(i))*field2(:,i)
         end do
         call check_state
         do i=1,n
            if (any(uind(:,i)/=sqrt(polarity(i))*v1(:,i)).or. &
                any(uinp(:,i)/=sqrt(polarity(i))*v2(:,i))) error stop 'Probe inputs mutated'
         end do
         qgrid=old_grid
         call reciprocal_apply(uind,uinp,r1,r2)
         reference_count=reference_count+1
         call check_state
         if (any(.not.ieee_is_finite(full1)).or.any(.not.ieee_is_finite(full2)).or. &
             any(.not.ieee_is_finite(r1)).or.any(.not.ieee_is_finite(r2))) error stop 'Nonfinite action'
         do i=1,n
            write(unit,'(a,3i8,12(1x,es26.17e3))') 'APPLY ',mode,col,i,full1(:,i),full2(:,i),r1(:,i),r2(:,i)
         end do
         uind=0d0
         uinp=0d0
         qgrid=old_grid
      end do
   end do
   if (complete_count/=30.or.reference_count/=30) error stop 'Exact operator denominator'
   write(unit,'(a,2i8)') 'END ',complete_count,reference_count
   close(unit,iostat=ios)
   if (ios/=0) error stop 'Application output close'
   call check_state
   if (any(uind/=0d0).or.any(uinp/=0d0).or.any(qgrid/=old_grid)) error stop 'Scratch restore failed'
   open(newunit=unit,file='parity_receipt.dat',status='new',action='write',iostat=ios)
   if (ios/=0) error stop 'Receipt collision'
   write(unit,'(a)') 'PERIODIC_PARITY_RECEIPT_V1'
   write(unit,'(a)') 'COUNTS full=30 separate_reciprocal=30 total_reciprocal=60 direct_initializations=1'
   write(unit,'(a)') 'COUNTS induced_solves=0 energies=0 gradients=0 factorizations=0'
   write(unit,'(a)') 'NORMAL_GUARDS_UNCHANGED TRUE'
   write(unit,'(a)') 'STATE_RESTORED TRUE'
   write(unit,'(a)') 'END'
   close(unit,iostat=ios)
   if (ios/=0) error stop 'Receipt close'
   call final
contains
   subroutine check_state
      if (any(qfac/=old_q).or.any(thetai1/=old_t1).or.any(thetai2/=old_t2).or. &
          any(thetai3/=old_t3).or.any(igrid/=old_igrid).or.any(lrh/=old_h).or. &
          any(lrbeta/=old_beta).or.any(polarity/=old_alpha)) error stop 'Persistent native state changed'
   end subroutine

   subroutine write_state
      integer :: u,ia,axis,slot,ix,iy,iz
      open(newunit=u,file='parity_state.dat',status='new',action='write',iostat=ios)
      if (ios/=0) error stop 'State output collision'
      write(u,'(a)') 'PERIODIC_PARITY_STATE_V1'
      write(u,'(a)') 'UNITS xyz=angstrom alpha=angstrom^3 beta=angstrom^-3 h=electron/angstrom^2'
      write(u,'(a)') 'RHS LABELS raw_Ed raw_Ep; add_shared_h_to_both; no_solved_dipoles'
      write(u,'(a,4i8,2(1x,es26.17e3))') 'GRID ',nfft1,nfft2,nfft3,bsorder,aewald,volbox
      write(u,'(a)') 'MODE FFTW TINFOIL MUTUAL ALL_ACTIVE SIMPLE_LOOP'
      do axis=1,3
         write(u,'(a,i8,3(1x,es26.17e3))') 'TRANS ',axis,transform(axis,:)
      end do
      do slot=1,2
         write(u,'(a,i8,4(1x,es26.17e3))') 'COEFFICIENT ',slot+899,lrcoef(:,slot)
      end do
      do ia=1,n
         write(u,'(a,3i8,6(1x,es26.17e3))') 'SITE ',ia,type(ia),ipole(ia), &
              x(ia),y(ia),z(ia),polarity(ia),pdamp(ia),thole(ia)
         write(u,'(a,i8,4(1x,es26.17e3))') 'RESPONSE ',ia,lrbeta(ia),lrh(:,ia)
         write(u,'(a,i8,6(1x,es26.17e3))') 'DIRECT ',ia,direct1(:,ia),direct2(:,ia)
         if (type(ia)==402) write(u,'(a,i8,1x,es26.17e3)') 'PAIR_GAMMA ',ia,thlval(jpolar(1),jpolar(ia))
         write(u,'(a,4i8)') 'IGRID ',ia,igrid(:,ia)
         do slot=1,5
            write(u,'(a,3i8,2(1x,es26.17e3))') 'THETA ',ia,1,slot,thetai1(1:2,slot,ia)
            write(u,'(a,3i8,2(1x,es26.17e3))') 'THETA ',ia,2,slot,thetai2(1:2,slot,ia)
            write(u,'(a,3i8,2(1x,es26.17e3))') 'THETA ',ia,3,slot,thetai3(1:2,slot,ia)
         end do
      end do
      do iz=1,nfft3
         do iy=1,nfft2
            do ix=1,nfft1
               write(u,'(a,3i8,1x,es26.17e3)') 'QFAC ',ix,iy,iz,qfac(ix,iy,iz)
            end do
         end do
      end do
      write(u,'(a,3i12)') 'END ',n,n*15,nfft1*nfft2*nfft3
      close(u,iostat=ios)
      if (ios/=0) error stop 'State close'
   end subroutine

   subroutine reciprocal_apply(input1,input2,output1,output2)
      real*8,intent(in) :: input1(3,n),input2(3,n)
      real*8,intent(out) :: output1(3,n),output2(3,n)
      integer :: ia,ix,iy,iz
      ! Exact admitted harness chain; qfac remains unsymmetrized and two packed
      ! channels remain unequal. This is positive normalized reciprocal R.
      fu1=matmul(transform,input1)
      fu2=matmul(transform,input2)
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
      do ia=1,n
         output1(:,ia)=sqrt(polarity(ia))*matmul(transpose(transform),phi1(2:4,ia))
         output2(:,ia)=sqrt(polarity(ia))*matmul(transpose(transform),phi2(2:4,ia))
      end do
   end subroutine
end program
