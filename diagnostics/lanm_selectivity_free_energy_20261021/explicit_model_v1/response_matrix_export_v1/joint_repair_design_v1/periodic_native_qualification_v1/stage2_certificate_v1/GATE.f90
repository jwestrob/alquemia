! Qualification-only observations and exact initialized-state gate.
! No certificate is created here: the external wrapper must admit and pin it.
module stage2_gate_state
   use atoms
   use boxes
   use bound
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
   logical :: started=.false.,backend_seen=.false.,certified=.false.,residual_seen=.false.
   integer :: direct_calls=0,mutual_calls=0,solve_entries=0,iterations=-1
   real*8 :: native_eps=0d0,native_epsd=0d0,native_epsp=0d0,native_epsold=0d0,hook_error=0d0
   character*40 :: role=''
   character*128 :: state_id=''
   character*64 :: binding=''
   real*8,allocatable :: raw1(:,:),raw2(:,:),rhs1(:,:),rhs2(:,:)
contains
   subroutine begin_gate
      integer :: u,ios,j
      character*240 :: line
      if (started) error stop 'Stage2 repeated entry'
      open(newunit=u,file='stage2_ticket.dat',status='old',action='read',iostat=ios)
      if (ios/=0) error stop 'Stage2 ticket absent'
      read(u,'(a)',iostat=ios) line
      if (ios/=0.or.trim(line)/='STAGE2_BOUND_INPUT_V1') error stop 'Stage2 ticket header'
      read(u,'(a)',iostat=ios) role
      if (ios/=0) error stop 'Stage2 ticket role'
      select case(trim(role))
      case('INITIALIZE','ACTIVE','BYPASS_OLD_ABSENT','BYPASS_CLONE_ABSENT','BYPASS_CLONE_ZERO')
      case default
         error stop 'Stage2 unsupported role'
      end select
      read(u,'(a)',iostat=ios) state_id
      if (ios/=0.or.len_trim(state_id)==0) error stop 'Stage2 state id'
      read(u,'(a)',iostat=ios) binding
      if (ios/=0.or.len_trim(binding)/=64) error stop 'Stage2 binding digest'
      do j=1,64
         if (index('0123456789abcdef',binding(j:j))==0) error stop 'Stage2 digest syntax'
      end do
      read(u,'(a)',iostat=ios) line
      if (ios/=0.or.trim(line)/='END') error stop 'Stage2 ticket tail'
      read(u,'(a)',iostat=ios) line
      if (ios>=0) error stop 'Stage2 trailing ticket'
      close(u)
      if (n/=11958.or.npole/=n) error stop 'Stage2 exact site denominator'
      if (any(ipole(1:n)/=[(j,j=1,n)])) error stop 'Stage2 site map'
      if (count(type(1:n)==900)/=1.or.type(1)/=900.or.any(type(1:n)==901)) error stop 'Stage2 La identity'
      if (count(type(1:n)==402)/=3980.or.count(type(1:n)==403)/=7960) error stop 'Stage2 water inventory'
      if (count(type(1:n)==405)/=7.or.count(type(1:n)==414)/=10) error stop 'Stage2 salt inventory'
      if (.not.use_bounds.or..not.orthogonal.or.nonprism.or..not.use_ewald) error stop 'Stage2 cell mode'
      if (use_mlist.or.use_ulist.or.use_replica.or.use_exfld) error stop 'Stage2 fixed field/list mode'
      if (boundary/='TINFOIL'.or.ffttyp/='FFTW'.or.poltyp/='MUTUAL') error stop 'Stage2 native mode'
      if (ewaldcut/=9d0.or.vdwcut/=9d0.or.vdwtaper/=8.1d0) error stop 'Stage2 cutoff changed'
      if (nefft1/=64.or.nefft2/=64.or.nefft3/=64.or.bsporder/=5.or.bseorder/=5) error stop 'Stage2 PME grid'
      if (apewald/=0.4d0.or.aeewald/=0.4d0) error stop 'Stage2 Ewald exponent'
      if (poleps/=1d-8.or.politer/=100) error stop 'Stage2 induction controls'
      if (any(.not.douind(1:n)).or.any(polarity(1:n)<1d-8)) error stop 'Stage2 alpha domain'
      if (any(.not.ieee_is_finite(polarity(1:n)))) error stop 'Stage2 nonfinite alpha'
      if (u1scale/=1d0.or.u2scale/=1d0.or.u3scale/=1d0.or.u4scale/=1d0) error stop 'Stage2 mutual scales'
      started=.true.
   end subroutine

   subroutine backend_check
      if (.not.started.or.backend_seen) error stop 'Stage2 backend sequence'
      if (role/='ACTIVE'.and.role/='INITIALIZE') error stop 'Stage2 active response in bypass'
      if (.not.lractive.or..not.lrperiodic.or.lrm/=1.or.lrpnwater/=3980) error stop 'Stage2 response state'
      backend_seen=.true.
   end subroutine

   subroutine induce_enter
      if (.not.started.or.role=='INITIALIZE') error stop 'Stage2 unauthorized solve entry'
      solve_entries=solve_entries+1
      if (solve_entries/=1) error stop 'Stage2 repeated induced solve'
      if (role=='ACTIVE') then
         if (.not.backend_seen.or.certified) error stop 'Stage2 pre-solve gate order'
      else
         if (lractive.or.lrperiodic.or.backend_seen) error stop 'Stage2 bypass is active'
      end if
   end subroutine

   subroutine before_direct(a,b)
      real*8,intent(in) :: a(3,n),b(3,n)
      if (.not.started) error stop 'Stage2 direct without ticket'
      direct_calls=direct_calls+1
      if (direct_calls/=1) error stop 'Stage2 repeated direct initialization'
      if (any(.not.ieee_is_finite(a)).or.any(.not.ieee_is_finite(b))) error stop 'Stage2 raw direct nonfinite'
      allocate(raw1,source=a)
      allocate(raw2,source=b)
   end subroutine

   subroutine after_direct(a,b)
      real*8,intent(in) :: a(3,n),b(3,n)
      real*8 :: tolerance,scale,error1,error2
      if (direct_calls/=1.or..not.allocated(raw1)) error stop 'Stage2 direct observation order'
      if (any(.not.ieee_is_finite(a)).or.any(.not.ieee_is_finite(b))) error stop 'Stage2 RHS nonfinite'
      allocate(rhs1,source=a)
      allocate(rhs2,source=b)
      if (lractive) then
         if (any(.not.ieee_is_finite(lrh))) error stop 'Stage2 h nonfinite'
         scale=max(1d0,maxval(abs(a)),maxval(abs(b)),maxval(abs(raw1)),maxval(abs(raw2)),maxval(abs(lrh)))
         tolerance=128d0*dble(3*n)*epsilon(1d0)*scale
         error1=maxval(abs((a-raw1)-lrh))
         error2=maxval(abs((b-raw2)-lrh))
         hook_error=max(error1,error2)
         if (hook_error>tolerance) error stop 'Stage2 shared direct h mismatch'
      else
         if (any(a/=raw1).or.any(b/=raw2)) error stop 'Stage2 bypass direct arithmetic'
      end if
   end subroutine

   subroutine certify_state
      integer :: actual,expected,ia,ib,lines
      character*512 :: a,b
      if (.not.started.or..not.backend_seen.or.certified.or.direct_calls/=1) error stop 'Stage2 certificate order'
      if (.not.allocated(rhs1)) error stop 'Stage2 RHS not observed'
      call stage2_write_state(raw1,raw2)
      if (role=='INITIALIZE') return
      if (role/='ACTIVE'.or.solve_entries/=1) error stop 'Stage2 certificate role'
      open(newunit=actual,file='certificate_actual_state.dat',status='old',action='read',iostat=ia)
      open(newunit=expected,file='expected_state.dat',status='old',action='read',iostat=ib)
      if (ia/=0.or.ib/=0) error stop 'Stage2 expected certificate state absent'
      lines=0
      do
         read(actual,'(a)',iostat=ia) a
         read(expected,'(a)',iostat=ib) b
         if (ia<0.and.ib<0) exit
         if (ia/=0.or.ib/=0.or.a/=b) error stop 'Stage2 exact initialized certificate mismatch'
         lines=lines+1
      end do
      close(actual)
      close(expected)
      if (lines<262144) error stop 'Stage2 incomplete initialized certificate'
      certified=.true.
   end subroutine

   subroutine mutual_enter
      if (.not.started.or.role=='INITIALIZE') error stop 'Stage2 unauthorized mutual application'
      if (solve_entries/=1) error stop 'Stage2 mutual outside sole solve'
      if (role=='ACTIVE'.and..not.certified) error stop 'Stage2 mutual before state certificate'
      mutual_calls=mutual_calls+1
   end subroutine

   subroutine solver_readback(it,eps,epsd,epsp,epsold)
      integer,intent(in) :: it
      real*8,intent(in) :: eps,epsd,epsp,epsold
      if (residual_seen.or.solve_entries/=1) error stop 'Stage2 solver readback order'
      if (.not.all(ieee_is_finite([eps,epsd,epsp,epsold]))) error stop 'Stage2 solver residual nonfinite'
      iterations=it
      native_eps=eps
      native_epsd=epsd
      native_epsp=epsp
      native_epsold=epsold
      residual_seen=.true.
   end subroutine

   subroutine finish_gate(u)
      integer,intent(in) :: u
      integer :: i
      write(u,'(a,1x,a)') 'ROLE',trim(role)
      write(u,'(a,1x,a)') 'STATE_ID',trim(state_id)
      write(u,'(a,1x,a)') 'INPUT_BINDING',binding
      if (role=='BYPASS_OLD_ABSENT') then
         if (lractive.or.lrperiodic) error stop 'Stage2 old library bypass active'
         write(u,'(a)') 'INTERNAL_COUNTERS unavailable_unchanged_old_archive'
         return
      end if
      if (solve_entries/=1.or.direct_calls/=1.or..not.residual_seen) error stop 'Stage2 solve accounting'
      if (role=='ACTIVE'.and..not.certified) error stop 'Stage2 missing certificate admission'
      if (role/='ACTIVE'.and.(lractive.or.lrperiodic)) error stop 'Stage2 bypass changed'
      write(u,'(a,3i12)') 'COUNTS induced_solve direct mutual ',solve_entries,direct_calls,mutual_calls
      write(u,'(a)') 'COUNTS extra_diagnostic_operator_applications=0 internal_FD=0'
      write(u,'(a,i8,4(1x,es26.17e3))') 'PCG_PRE_PEEK ',iterations,native_eps,native_epsd,native_epsp,native_epsold
      if (role=='ACTIVE') then
         write(u,'(a,1x,es26.17e3)') 'SHARED_H_READBACK_MAX_ERROR ',hook_error
         do i=1,n
            write(u,'(a,i8,6(1x,es26.17e3))') 'RAW_DIRECT ',i,raw1(:,i),raw2(:,i)
            write(u,'(a,i8,6(1x,es26.17e3))') 'CORRECTED_RHS ',i,rhs1(:,i),rhs2(:,i)
            write(u,'(a,i8,7(1x,es26.17e3))') 'FINAL_POST_PEEK ',i,polarity(i),uind(:,i),uinp(:,i)
            write(u,'(a,i8,4(1x,es26.17e3))') 'CORRECTION ',i,lrbeta(i),lrh(:,i)
         end do
      end if
   end subroutine
end module

! Explicit-shape external callbacks keep original fixed-form solver arithmetic intact.
subroutine stage2_backend_check
   use stage2_gate_state, only: backend_check
   call backend_check
end subroutine
subroutine stage2_induce_enter
   use stage2_gate_state, only: induce_enter
   call induce_enter
end subroutine
subroutine stage2_before_direct(a,b)
   use atoms, only: n
   use stage2_gate_state, only: before_direct
   real*8 a(3,n),b(3,n)
   call before_direct(a,b)
end subroutine
subroutine stage2_after_direct(a,b)
   use atoms, only: n
   use stage2_gate_state, only: after_direct
   real*8 a(3,n),b(3,n)
   call after_direct(a,b)
end subroutine
subroutine stage2_certify
   use stage2_gate_state, only: certify_state
   call certify_state
end subroutine
subroutine stage2_mutual_enter
   use stage2_gate_state, only: mutual_enter
   call mutual_enter
end subroutine
subroutine stage2_solver_readback(it,eps,epsd,epsp,epsold)
   use stage2_gate_state, only: solver_readback
   integer it
   real*8 eps,epsd,epsp,epsold
   call solver_readback(it,eps,epsd,epsp,epsold)
end subroutine
