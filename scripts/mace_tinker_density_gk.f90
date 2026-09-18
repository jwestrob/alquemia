! Isolated native static GK, supplied-field response, and vacuum identity frontend.
program alquemia_tinker_density_gk
  use, intrinsic :: ieee_arithmetic
  use atoms
  use couple
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
  real*8 :: frame(3,3)
  logical :: planar
  integer :: k, count_input
  real*8 :: new_charge,energy,total_energy,begin_cpu,end_cpu,eps_final
  real*8, allocatable :: supplied(:,:,:),expected_born(:),fd(:,:),fp(:,:),fs(:,:),fps(:,:)
  character(len=2048) :: field_path
  character(len=32) :: operation
  character(len=6) :: switch_mode
  character(len=2048) :: mask_path
  character(len=2048) :: override_path
  logical, allocatable :: frozen(:)

  if (command_argument_count()<4.or.command_argument_count()>5) error stop 'expected XYZ MASK CHARGES OP [FIELDS]'
  call get_command_argument(2, mask_path)
  call get_command_argument(3, override_path)
  call get_command_argument(4, operation)
  if (trim(operation)/='static'.and.trim(operation)/='solve'.and.trim(operation)/='identity') error stop 'invalid mode'
  if (trim(operation)=='solve') then
    if (command_argument_count()/=5) error stop 'solve requires fields'
    call get_command_argument(5,field_path)
  end if
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
  open(newunit=mask_unit,file=trim(override_path),status='old',action='read',iostat=ios)
  if (ios/=0) error stop 'missing charge override file'
  read(mask_unit,*,iostat=ios) count_input
  if (ios/=0.or.count_input/=n) error stop 'charge inventory mismatch'
  do i=1,n
    write(*,*) 'ORIGINAL_POLE',i,pole(:,i)
    read(mask_unit,*,iostat=ios) index,new_charge
    if (ios/=0.or.index/=i) error stop 'charge atom order mismatch'
    if (.not.ieee_is_finite(new_charge)) error stop 'nonfinite source/environment charge'
    pole(1,i)=new_charge
    if (frozen(i)) pole(2:13,i)=0d0
    ! kpolar may have removed an initially neutral nonresponsive QM source.
    ipole(i)=i
    pollist(i)=i
    polsiz(i)=13
    mono0(i)=pole(1,i)
  end do
  npole=n
  read(mask_unit,*,iostat=ios) index
  if (ios>=0) error stop 'extra/malformed charge override data'
  close(mask_unit)
  if (trim(operation)=='identity') then
    do i=1,n
      if (.not.frozen(i)) pole(:,i)=0d0
    end do
    polarity=0d0;pdamp=0d0;douind=.false.;npolar=0
    use_polar=.false.;use_solv=.false.;use_born=.false.
  end if
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
  do i=1,n
    write(*,*) 'RAW_POLE',i,pole(:,i)
    write(*,*) 'IDENTITY',i,type(i),class(i),atomic(i),mass(i)
    write(*,*) 'BONDS',i,n12(i),(i12(j,i),j=1,n12(i))
    write(*,*) 'POLE_INDEX',i,ipole(i),pollist(i),polsiz(i)
  end do
  call chkpole
  call rotpole('MPOLE')
  do i=1,n
    call rotmat(i,frame,planar)
    write(*,*) 'AXIS',i,trim(polaxe(i)),zaxis(i),xaxis(i),yaxis(i)
    write(*,*) 'FRAME',i,planar,((frame(j,k),k=1,3),j=1,3)
    write(*,*) 'LOCAL_POLE',i,pole(:,i)
    write(*,*) 'GLOBAL_POLE',i,rpole(:,i)
  end do
  write(*,*) 'ALQUEMIA_MOMENTS_COMPLETE',n
  uind=0d0;uinp=0d0;uinds=0d0;uinps=0d0
  udir=0d0;udirp=0d0;udirs=0d0;udirps=0d0
  write(*,*) 'CONSTANTS',electric,dielec
  write(*,*) 'ACTIVE_TERMS ',trim(operation),use_polar,use_solv,use_born
  debug=.true.
  call system_clock(start_clock,clock_rate)
  call cpu_time(begin_cpu)
  if (trim(operation)=='solve') then
    allocate(supplied(3,n,4),expected_born(n))
    open(newunit=mask_unit,file=trim(field_path),status='old',action='read',iostat=ios)
    if (ios/=0) error stop 'missing supplied field file'
    read(mask_unit,*,iostat=ios) count_input
    if (ios/=0.or.count_input/=n) error stop 'field count mismatch'
    do i=1,n
      read(mask_unit,*,iostat=ios) index,expected_born(i),supplied(:,i,1), &
        supplied(:,i,2),supplied(:,i,3),supplied(:,i,4)
      if (ios/=0.or.index/=i) error stop 'field order mismatch'
      if (.not.all(ieee_is_finite(supplied(:,i,:)))) error stop 'nonfinite supplied field'
    end do
    read(mask_unit,*,iostat=ios) index
    if (ios>=0) error stop 'extra field rows'
    close(mask_unit)
    call replica(0d0)
    call born
    if (.not.all(ieee_is_finite(expected_born))) error stop 'nonfinite supplied cavity'
    if (any(abs(rborn-expected_born)>1d-12)) error stop 'supplied field cavity mismatch'
    write(*,*) 'ALQUEMIA_RESPONSE_BEGIN'
    flush(6)
    call alquemia_induce0c(supplied,eps_final)
    write(*,*) 'EXACT_FINAL_RMS_DEBYE',eps_final
  else
    use_polar=.false.
    write(*,*) 'ALQUEMIA_ENERGY_BEGIN'
    flush(6)
    total_energy=energy()
    write(*,*) 'ENERGIES',total_energy,em,ep,es,eb,ea,eba,eub,eaa,eopb,eopd,eid,eit, &
      et,ept,ebt,eat,ett,ev,er,edsp,ec,ecd,ed,ect,erxf,elf,eg,ex
    if (trim(operation)=='static') then
      allocate(fd(3,n),fp(3,n),fs(3,n),fps(3,n))
      switch_mode='MPOLE'
      call switch(switch_mode)
      write(*,*) 'ALQUEMIA_FIELDS_BEGIN'
      call dfield0d(fd,fp,fs,fps)
      do i=1,n
        write(*,*) 'FIELD',i,rborn(i),fd(:,i),fp(:,i),fs(:,i),fps(:,i)
      end do
    end if
  end if
  call cpu_time(end_cpu)
  call system_clock(end_clock)
  if (trim(operation)/='identity') then
    do i=1,n
      write(*,*) 'RESPONSE',i,rborn(i),uind(:,i),uinp(:,i),uinds(:,i),uinps(:,i)
    end do
  else
    do i=1,n
      write(*,*) 'VACUUM_RESPONSE',i,uind(:,i),uinp(:,i),uinds(:,i),uinps(:,i)
    end do
  end if
  write(*,*) 'KERNEL_TIMING',real(end_clock-start_clock,8)/real(clock_rate,8),end_cpu-begin_cpu
  write(*,*) 'ALQUEMIA_DENSITY_GK_COMPLETE ',trim(operation),n
  call final
end program alquemia_tinker_density_gk
