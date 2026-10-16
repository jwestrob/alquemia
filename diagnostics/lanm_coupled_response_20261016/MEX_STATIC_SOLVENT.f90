! Fixed-charge energy and analytic gradient; derived from the preserved native replay.
! No SCC calculator. q is held fixed; gradient is dG/dR in Hartree/Bohr.
program replay
   use, intrinsic :: iso_fortran_env, only: output_unit
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use xtb_mctc_accuracy, only: wp
   use mctc_io_convert, only: aatoau
   use xtb_mctc_symbols, only: toNumber
   use xtb_type_environment, only: TEnvironment, init_env => init
   use xtb_solv_input, only: TSolvInput
   use xtb_solv_model, only: TSolvModel, init_model => init, newBornModel, info
   use xtb_solv_gbsa, only: TBorn
   use xtb_solv_kernel, only: gbKernel
   use xtb_solv_state, only: solutionState
   implicit none
   type(TEnvironment) :: env
   type(TSolvInput) :: input
   type(TSolvModel) :: model
   type(TBorn) :: born
   integer :: n, i, j, u, stat
   integer, allocatable :: num(:)
   real(wp), allocatable :: xyz(:,:), q(:), qsh(:), be(:,:), hb(:), grad(:,:)
   real(wp) :: ge, gh, gs, gshift, bound_e, bound_h
   real(wp), parameter :: dq = 0.0_wp ! exact declared FF model charges, serialized with round-trip precision
   logical :: failed
   character(len=4096) :: xyzfile, qfile, line
   character(len=4) :: symbol

   if (command_argument_count() /= 2) error stop 'expected exact XYZ and charges paths'
   call get_command_argument(1, xyzfile)
   call get_command_argument(2, qfile)
   open(newunit=u,file=trim(xyzfile),status='old',action='read')
   read(u,*) n
   if (n /= 2049) error stop 'declared Mex physical scope requires 2049 atoms'
   read(u,'(a)') line
   allocate(num(n),xyz(3,n),q(n),qsh(0),be(n,n),hb(n),grad(3,n))
   do i=1,n
      read(u,*) symbol,xyz(:,i)
      num(i)=toNumber(symbol)
      if (num(i)<1) error stop 'unknown atomic identity'
   end do
   close(u)
   xyz=xyz*aatoau ! actual xtb_io_reader -> mctc_io_read_xyz conversion
   open(newunit=u,file=trim(qfile),status='old',action='read')
   do i=1,n
      read(u,*) q(i)
   end do
   read(u,*,iostat=stat) line
   if (stat >= 0) error stop 'unexpected extra charge records'
   close(u)
   if (.not.all(ieee_is_finite(q))) error stop 'nonfinite charges'
   call init_env(env)
   input%solvent='water'
   input%alpb=.true.
   input%cosmo=.false.
   input%tmcosmo=.false.
   input%kernel=gbKernel%p16
   input%state=solutionState%gsolv
   input%temperature=298.15_wp
   input%nAng=230
   input%ionStrength=0.0_wp
   input%ionRad=0.0_wp
   call init_model(model,env,input,2)
   call env%check(failed)
   if (failed) then
      call env%show('native model setup failed')
      error stop 'native model setup failed'
   end if
   if (model%paramFile /= 'internal GFN2-xTB/ALPB') error stop 'unexpected external solvent parameters'
   call info(model,output_unit)
   call newBornModel(model,env,born,num)
   call born%update(env,num,xyz)
   call env%check(failed)
   if (failed) error stop 'native cavity update failed'
   ! Exact pinned getEnergyParts never reads qsh. A length-zero unused argument
   ! supplies no invented physical shell populations; no SCC is invoked.
   call born%getEnergyParts(env,q,qsh,ge,gh,gs,gshift)
   call env%check(failed)
   if (failed) error stop 'native component evaluation failed'
   grad=0.0_wp
   call born%addGradient(env,num,xyz,q,qsh,grad)
   call env%check(failed)
   if (failed .or. .not.all(ieee_is_finite(grad))) error stop 'native gradient failed'
   write(*,'(a,4(1x,es25.17))') 'COMPONENTS',ge,gh,gs,gshift
   open(newunit=u,file='gradient.dat',status='new',action='write')
   do i=1,n
      write(u,'(3(es25.17,1x))') grad(:,i)
   end do
   close(u)
   ! Bounds only: energy values above come entirely from the native routine.
   ! For q*=q+e, |e_i|<=dq, the quadratic error is bounded by
   ! dq*||Bq||_1 + dq^2/2*sum_ij|B_ij|. Native B includes 2*hbw
   ! on its diagonal; remove that only for the electrostatic error bound.
   be=born%bornMat
   hb=0.0_wp
   if (born%lhb) hb=born%hbw
   do i=1,n
      be(i,i)=be(i,i)-2.0_wp*hb(i)
   end do
   bound_e=dq*sum(abs(matmul(be,q)))+0.5_wp*dq*dq*sum(abs(be))
   bound_h=2.0_wp*dq*sum(abs(hb*q))+dq*dq*sum(abs(hb))
   write(*,'(a,2(1x,es25.17))') 'ROUNDING_BOUNDS',bound_e,bound_h
   write(*,'(a,1x,es25.17)') 'CHARGE_SUM',sum(q)
   write(*,'(a,4(1x,es25.17))') 'CAVITY_GLOBAL',born%aDet,born%keps,born%alpbet,born%dielectricConst
end program replay
