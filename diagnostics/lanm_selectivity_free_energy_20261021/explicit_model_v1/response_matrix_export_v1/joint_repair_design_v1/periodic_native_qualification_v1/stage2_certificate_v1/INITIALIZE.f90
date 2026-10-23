program stage2_certificate_initialization
   use atoms
   use polar
   use stage2_gate_state, only: begin_gate,certify_state,role,direct_calls,mutual_calls,solve_entries
   implicit none
   real*8,allocatable :: raw1(:,:),raw2(:,:)
   integer :: u,ios
   call initial
   call getxyz
   call mechanic
   call begin_gate
   if (role/='INITIALIZE') error stop 'Initialization requires explicit initializer ticket'
   call chkpole
   call rotpole('MPOLE')
   call localresponsesetup
   uind=0d0
   uinp=0d0
   allocate(raw1(3,n),raw2(3,n))
   call dfield0(raw1,raw2)
   call certify_state
   if (direct_calls/=1.or.mutual_calls/=0.or.solve_entries/=0) error stop 'Initialization count'
   open(newunit=u,file='stage2_initialization_receipt.dat',status='new',action='write',iostat=ios)
   if (ios/=0) error stop 'Initialization receipt collision'
   write(u,'(a)') 'STAGE2_INITIALIZATION_V1'
   write(u,'(a)') 'COUNTS direct=1 mutual=0 solve=0 energy=0 gradient=0 factor=0'
   write(u,'(a)') 'STATE_EXPORTED_NOT_CERTIFIED'
   write(u,'(a)') 'END'
   close(u,iostat=ios)
   if (ios/=0) error stop 'Initialization receipt close'
   call final
end program
