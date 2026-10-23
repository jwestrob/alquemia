! Operator-only initialization; no induced solve, energy, gradient or dynamics.
program init_pme
   use atoms
   use bound
   use boxes
   use ewald
   use fft
   use keys
   use limits
   use mplpot
   use mpole
   use pme
   use polar
   use polpot
   use ieee_arithmetic
   implicit none
   integer :: i,next,unit,ios
   real*8,allocatable :: field(:,:)
   logical :: exists
   character*240 :: record,keyword
   character*16 :: metal
   call initial
   call getxyz
   call mechanic
   call nextarg(metal,exists)
   if (.not.exists.or.(trim(metal)/='La'.and.trim(metal)/='Dy')) error stop 'Exact metal selector required'
   if (n/=11958.or.npole/=11958) error stop 'Bound BIO fixture count mismatch'
   if (count(type(1:n)==402)/=3980.or.count(type(1:n)==403)/=7960) error stop 'Water counts'
   if (count(type(1:n)==405)/=7.or.count(type(1:n)==414)/=10) error stop 'Salt counts'
   if (trim(metal)=='La') then
      if (count(type(1:n)==900)/=1) error stop 'La selector mismatch'
   else
      if (count(type(1:n)==901)/=1) error stop 'Dy selector mismatch'
   end if
   do i=1,nkey
      record=keyline(i)
      next=1
      call gettext(record,keyword,next)
      call upcase(keyword)
      if (trim(keyword)=='POLAR-LOCAL-RESPONSE'.or.trim(keyword)=='POLAR-DONOR-RESPONSE'.or. &
          trim(keyword)=='RESPONSE-MATRIX-EXPORT'.or.trim(keyword)=='EXTERNAL-FIELD') &
          error stop 'Only explicitly derived operator input is accepted'
   end do
   if (.not.use_bounds.or..not.orthogonal.or.nonprism.or..not.use_ewald) error stop 'Periodic orthogonal Ewald only'
   if (boundary/='TINFOIL'.or.ffttyp/='FFTW'.or.poltyp/='MUTUAL') error stop 'Native convention mismatch'
   if (use_expol.or.use_chgpen.or.use_dewald) error stop 'Unqualified response/penetration/dispersion route'
   if (xbox/=49.324d0.or.ybox/=49.324d0.or.zbox/=49.324d0) error stop 'Cell changed'
   if (ewaldcut/=9d0.or.vdwcut/=9d0.or.vdwtaper/=8.1d0) error stop 'Prospective cutoff mismatch'
   if (nefft1/=64.or.nefft2/=64.or.nefft3/=64.or.bsporder/=5.or.bseorder/=5) error stop 'PME grid/order mismatch'
   if (apewald/=0.4d0.or.aeewald/=0.4d0) error stop 'Ewald alpha mismatch'
   if (any(.not.douind(1:n)).or.any(polarity(1:n)<1d-8)) error stop 'Active polarizability domain'
   ! Use the same native frame initialization as the multipole path.
   call chkpole
   call rotpole('MPOLE')
   ! Explicit absence is checked above, and the normal setup remains called.
   ! A forgotten response record is not silently skipped to evade its guards.
   call localresponsesetup
   if (lractive.or.lrperiodic) error stop 'Response must remain disabled'
   ! Exact assignments from dfield0c; only its reciprocal initialization is needed.
   nfft1=nefft1
   nfft2=nefft2
   nfft3=nefft3
   bsorder=bsporder
   aewald=apewald
   allocate(field(3,n))
   field=0d0
   call udirect1(field)
   if (any(.not.ieee_is_finite(field))) error stop 'Nonfinite initialized reciprocal permanent field'
   inquire(file='pme_norm_export_001.dat',exist=exists)
   if (.not.exists) error stop 'Required first initialized export missing'
   open(newunit=unit,file='operator_initialization.dat',status='new',action='write',iostat=ios)
   if (ios/=0) error stop 'Initialization receipt collision/open failure'
   write(unit,'(a)') 'PME_OPERATOR_ONLY_INITIALIZATION_V1'
   write(unit,'(a,1x,a)') 'METAL',trim(metal)
   write(unit,'(a,2i8)') 'SITES ',n,npole
   write(unit,'(a)') 'COUNTS reciprocal_initializations=1 induced_solves=0 energies=0 gradients=0'
   write(unit,'(a)') 'ENERGY_MODEL_ADMISSION FALSE'
   write(unit,'(a)') 'FULL_STABILITY_ADMISSION FALSE'
   write(unit,'(a)') 'END'
   close(unit,iostat=ios)
   if (ios/=0) error stop 'Initialization receipt close failure'
   call final
end program
