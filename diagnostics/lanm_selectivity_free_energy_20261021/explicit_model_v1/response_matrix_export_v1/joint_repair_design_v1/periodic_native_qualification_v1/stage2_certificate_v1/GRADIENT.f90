program stage2_one_gradient
   use atoms
   use energi, only: ep
   use virial, only: vir
   use chgpot, only: dielec,electric
   use stage2_gate_state, only: begin_gate,finish_gate,role
   use ieee_arithmetic
   implicit none
   real*8 :: total_energy
   real*8,allocatable :: g(:,:)
   integer :: i,u,ios
   call initial
   call getxyz
   call mechanic
   call begin_gate
   if (role=='INITIALIZE') error stop 'Initialization ticket cannot release gradient'
   allocate(g(3,n))
   ! Exactly one native analytic gradient call; no testgrad/energy/FD loop.
   call gradient(total_energy,g)
   if (.not.ieee_is_finite(total_energy).or..not.ieee_is_finite(ep)) error stop 'Stage2 nonfinite energy'
   if (any(.not.ieee_is_finite(g)).or.any(.not.ieee_is_finite(vir))) error stop 'Stage2 nonfinite G/virial'
   open(newunit=u,file='stage2_gradient.dat',status='new',action='write',iostat=ios)
   if (ios/=0) error stop 'Stage2 gradient output collision'
   write(u,'(a)') 'STAGE2_ONE_ANALYTIC_GRADIENT_V1'
   write(u,'(a)') 'UNITS E=kcal/mol G=positive_dE/dr_kcal/mol/angstrom VIRIAL=kcal/mol'
   write(u,'(a,1x,es26.17e3)') 'TOTAL_ENERGY ',total_energy
   write(u,'(a,1x,es26.17e3)') 'POLARIZATION_ENERGY ',ep
   write(u,'(a,2(1x,es26.17e3))') 'CONTRACTION_CONSTANTS electric dielec ',electric,dielec
   do i=1,n
      write(u,'(a,i8,3(1x,es26.17e3))') 'TOTAL_GRADIENT ',i,g(:,i)
   end do
   do i=1,3
      write(u,'(a,i8,3(1x,es26.17e3))') 'TOTAL_VIRIAL ',i,vir(i,:)
   end do
   call finish_gate(u)
   write(u,'(a)') 'END gradient_calls=1 energy_only_calls=0 internal_FD=0'
   close(u,iostat=ios)
   if (ios/=0) error stop 'Stage2 gradient output close'
   call final
end program
