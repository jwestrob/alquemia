! Pinned native Tinker electrostatic/GK frontend. One energy() call in run mode.
! preflight mode performs parameter initialization only. No gradients or motion.
program alquemia_tinker_framework_solver
  use, intrinsic :: ieee_arithmetic
  use atoms
  use atomid
  use mpole
  use polar
  use polpot
  use potent
  use solpot
  use solute
  use gkstuf
  use nonpol
  use energi
  use chgpot
  use limits
  use bound
  use inform
  use omp_lib, only: omp_get_max_threads
  implicit none
  integer :: i, j, frozen_count, mask_unit, ios, index
  integer(kind=8) :: start_clock, end_clock, clock_rate
  real*8 :: energy, total_energy, begin_cpu, end_cpu
  character(len=2048) :: mask_path
  character(len=32) :: action
  logical, allocatable :: frozen(:)

  if (command_argument_count() /= 3) error stop 'expected XYZ MASK preflight|run'
  call get_command_argument(2, mask_path)
  call get_command_argument(3, action)
  if (trim(action)/='preflight' .and. trim(action)/='run') error stop 'unknown mode'
  call initial
  call getxyz
  call mechanic
  if (solvtyp/='GK' .or. poltyp/='MUTUAL') error stop 'wrong solvent/polarization'
  if (dielec/=1d0 .or. gkc/=2.455d0) error stop 'wrong dielectric/GKC'
  if (use_bounds .or. use_ewald .or. mpolecut/=1d12) error stop 'wrong boundary/cutoff'
  if (.not.use_mpole .or. .not.use_polar .or. .not.use_solv .or. .not.use_born) &
    error stop 'missing electrostatic/GK term'
  if (use_bond.or.use_angle.or.use_strbnd.or.use_urey.or.use_angang.or.use_opbend.or. &
      use_opdist.or.use_improp.or.use_imptor.or.use_tors.or.use_pitors.or.use_strtor.or. &
      use_angtor.or.use_tortor.or.use_vdw.or.use_repel.or.use_disp.or.use_charge.or. &
      use_chgdpl.or.use_dipole.or.use_chgtrn.or.use_chgflx.or.use_rxnfld.or.use_metal.or. &
      use_geom.or.use_extra.or.use_orbit.or.use_mutate.or.use_expol) error stop 'extra active term'
  allocate(frozen(n)); frozen=.false.
  open(newunit=mask_unit,file=trim(mask_path),status='old',action='read',iostat=ios)
  if (ios/=0) error stop 'missing source mask'
  read(mask_unit,*,iostat=ios) frozen_count
  if (ios/=0) error stop 'missing mask count'
  if (frozen_count<0 .or. frozen_count>=n) error stop 'invalid source count'
  do j=1,frozen_count
    read(mask_unit,*,iostat=ios) index
    if (ios/=0) error stop 'missing source index'
    if (index<1.or.index>n) error stop 'invalid source index'
    if (frozen(index)) error stop 'duplicate source index'
    frozen(index)=.true.
    douind(index)=.false.
  end do
  read(mask_unit,*,iostat=ios) index
  if (ios>=0) error stop 'extra/malformed source data'
  close(mask_unit)
  write(*,*) 'ALQUEMIA_PARAMETERS',n,npole,npolar,frozen_count,omp_get_max_threads()
  ! The native GK kernel hardcodes bulk dielectric78.3; it is source-pinned.
  write(*,*) 'GLOBAL',dielec,78.3d0,gkc,descoff,doffset,poleps,politer,mpolecut
  write(*,*) 'GK_FLAGS',trim(borntyp),useneck,usetanh
  write(*,*) 'NONPOLAR',cavprb,solvprs,surften,epso,epsh,rmino,rminh,awater,slevy,shctd,dspoff
  do i=1,n
    if (.not.all(ieee_is_finite([rsolv(i),rdescr(i),shct(i),sneck(i),pole(1,i), &
        polarity(i),pdamp(i)]))) error stop 'nonfinite input parameters'
    if (rsolv(i)<=0.or.rdescr(i)<=0.or.shct(i)<0) error stop 'invalid radius/scale'
    write(*,*) 'PREP',i,atomic(i),x(i),y(i),z(i),rsolv(i),rdescr(i),shct(i),sneck(i), &
      pole(1,i),polarity(i),pdamp(i),douind(i),radcav(i),raddsp(i),epsdsp(i),cdsp(i)
  end do
  if (trim(action)=='preflight') then
    write(*,*) 'ALQUEMIA_PREFLIGHT_COMPLETE',n
    call final
    stop
  end if
  debug=.true.
  write(*,*) 'ALQUEMIA_ENERGY_BEGIN'
  flush(6)
  call system_clock(start_clock,clock_rate)
  call cpu_time(begin_cpu)
  total_energy=energy()
  call cpu_time(end_cpu)
  call system_clock(end_clock)
  write(*,*) 'ENERGIES',total_energy,em,ep,es,eb,ea,eba,eub,eaa,eopb,eopd,eid,eit, &
    et,ept,ebt,eat,ett,ev,er,edsp,ec,ecd,ed,ect,erxf,elf,eg,ex
  do i=1,n
    write(*,*) 'RESPONSE',i,rborn(i),(uind(j,i),j=1,3),(uinp(j,i),j=1,3), &
      (uinds(j,i),j=1,3),(uinps(j,i),j=1,3)
  end do
  write(*,*) 'KERNEL_TIMING',real(end_clock-start_clock,8)/real(clock_rate,8),end_cpu-begin_cpu
  write(*,*) 'ALQUEMIA_ENERGY_COMPLETE',n,frozen_count
  call final
end program alquemia_tinker_framework_solver
