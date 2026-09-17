! Parameter-only adapter for a pinned Tinker library. No energies or forces.
! A requested frozen-source list is applied after mechanic because the inspected
! upstream kpolar resets douind for sites with nonzero polarizability.
program alquemia_tinker_parameter_probe
  use atoms
  use atomid
  use couple
  use mpole
  use polar
  implicit none
  integer :: i, j, frozen_count, mask_unit, ios, index
  character(len=2048) :: mask_path
  logical, allocatable :: seen(:)

  if (command_argument_count() /= 2) error stop 'expected XYZ and frozen-index file'
  call get_command_argument(2, mask_path)
  call initial
  call getxyz
  call mechanic
  write(*,*) 'ALQUEMIA_NATIVE_BEGIN', n, npole, npolar
  do i=1,n
    write(*,*) 'ATOM', i, type(i), atomic(i), x(i), y(i), z(i), &
      polarity(i), pdamp(i), thole(i), douind(i)
    write(*,*) 'AXIS', i, trim(polaxe(i)), zaxis(i), xaxis(i), yaxis(i)
    write(*,*) 'POLE', i, (pole(j,i), j=1,13)
    write(*,*) 'BONDS', i, n12(i), (i12(j,i), j=1,n12(i))
  end do
  allocate(seen(n)); seen=.false.
  open(newunit=mask_unit, file=trim(mask_path), status='old', action='read', iostat=ios)
  if (ios /= 0) error stop 'cannot open frozen-index file'
  read(mask_unit,*,iostat=ios) frozen_count
  if (ios /= 0) error stop 'missing frozen-index count'
  if (frozen_count < 1 .or. frozen_count >= n) error stop 'invalid frozen-index count'
  do j=1,frozen_count
    read(mask_unit,*,iostat=ios) index
    if (ios /= 0) error stop 'missing frozen index'
    if (index < 1 .or. index > n) error stop 'frozen index outside physical system'
    if (seen(index)) error stop 'duplicate frozen index'
    seen(index)=.true.
    write(*,*) 'FROZEN_BEFORE', index, douind(index), polarity(index), pdamp(index)
    douind(index)=.false.
    write(*,*) 'FROZEN_AFTER', index, douind(index), polarity(index), pdamp(index)
  end do
  read(mask_unit,*,iostat=ios) index
  if (ios >= 0) error stop 'extra or malformed frozen-index data'
  close(mask_unit)
  write(*,*) 'ALQUEMIA_PARAMETER_ONLY_COMPLETE', n, frozen_count
  call final
end program alquemia_tinker_parameter_probe
