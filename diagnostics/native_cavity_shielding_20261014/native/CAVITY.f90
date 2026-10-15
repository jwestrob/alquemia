! Same native xTB6.7.1 factory; one cavity update, no SCC or dense wrapper copy.
program cavity
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
   integer :: n, nm, i, j, k, u, stat, c, v, nc
   integer, allocatable :: num(:), core(:)
   logical, allocatable :: cap(:), uncertain(:), real_uncertain(:)
   real(wp), allocatable :: xyz(:,:), q(:), qzero(:), qsh(:), hb(:)
   real(wp) :: full(4), repeated(4), zero(4), blocks(3), shape_blocks(3)
   real(wp) :: b, coef, qr, qc, ge_bound, hb_bound, zero_ge_bound, zero_hb_bound
   real(wp) :: shape_bound, zero_shape_bound
   real(wp), parameter :: dq=5.0e-9_wp
   logical :: failed
   character(len=4096) :: xyzfile,qfile,maskfile,line
   character(len=4) :: symbol

   if (command_argument_count()/=3) error stop 'expected XYZ charges masks'
   call get_command_argument(1,xyzfile)
   call get_command_argument(2,qfile)
   call get_command_argument(3,maskfile)
   open(newunit=u,file=trim(xyzfile),status='old',action='read')
   read(u,*) n
   if (n<1) error stop 'empty cavity'
   read(u,'(a)') line
   allocate(num(n),xyz(3,n),q(n),qzero(n),qsh(0),hb(n),cap(n),uncertain(n),real_uncertain(n))
   do i=1,n
      read(u,*) symbol,xyz(:,i)
      num(i)=toNumber(symbol)
      if (num(i)<1) error stop 'invalid element'
   end do
   close(u)
   if (.not.all(ieee_is_finite(xyz))) error stop 'nonfinite coordinates'
   xyz=xyz*aatoau
   open(newunit=u,file=trim(qfile),status='old',action='read')
   do i=1,n
      read(u,*) q(i)
   end do
   read(u,*,iostat=stat) line
   if (stat>=0) error stop 'extra charge records'
   close(u)
   if (.not.all(ieee_is_finite(q))) error stop 'nonfinite charge'
   open(newunit=u,file=trim(maskfile),status='old',action='read')
   read(u,*) nm
   if (nm/=n) error stop 'mask atom count mismatch'
   do i=1,n
      read(u,*) c,v
      if (min(c,v)<0.or.max(c,v)>1) error stop 'nonbinary mask'
      cap(i)=c==1
      uncertain(i)=v==1
      if (cap(i).and..not.uncertain(i)) error stop 'cap outside saved core'
      if (.not.uncertain(i).and.q(i)/=0.0_wp) error stop 'nonzero added environment charge'
   end do
   read(u,*,iostat=stat) line
   if (stat>=0) error stop 'extra mask records'
   close(u)
   nc=count(uncertain)
   allocate(core(nc))
   k=0
   do i=1,n
      if (.not.uncertain(i)) cycle
      k=k+1;core(k)=i
   end do
   qzero=q
   where(cap) qzero=0.0_wp
   real_uncertain=uncertain.and..not.cap
   call init_env(env)
   input%solvent='water';input%alpb=.true.;input%cosmo=.false.;input%tmcosmo=.false.
   input%kernel=gbKernel%p16;input%state=solutionState%gsolv
   input%temperature=298.15_wp;input%nAng=230;input%ionStrength=0.0_wp;input%ionRad=0.0_wp
   call init_model(model,env,input,2)
   call env%check(failed)
   if (failed) error stop 'native model initialization failed'
   if (model%paramFile/='internal GFN2-xTB/ALPB') error stop 'external parameter override'
   call info(model,output_unit)
   call newBornModel(model,env,born,num)
   call born%update(env,num,xyz)
   call env%check(failed)
   if (failed) error stop 'native cavity update failed'
   ! qsh is unused by the exact archived native getEnergyParts implementation.
   call born%getEnergyParts(env,q,qsh,full(1),full(2),full(3),full(4))
   call born%getEnergyParts(env,q,qsh,repeated(1),repeated(2),repeated(3),repeated(4))
   call born%getEnergyParts(env,qzero,qsh,zero(1),zero(2),zero(3),zero(4))
   hb=0.0_wp
   if (born%lhb) hb=born%hbw
   coef=born%keps*born%alpbet/born%aDet
   qr=sum(q,mask=.not.cap);qc=sum(q,mask=cap)
   shape_blocks=[0.5_wp*coef*qr**2,coef*qr*qc,0.5_wp*coef*qc**2]
   blocks=0.0_wp
   ! Only saved core charges can be nonzero. Never allocate or dump another N*N array.
   do k=1,nc
      i=core(k)
      do c=1,nc
         j=core(c);b=electrostatic(i,j)
         if (.not.cap(i).and..not.cap(j)) then
            v=1
         else if (cap(i).and.cap(j)) then
            v=3
         else
            v=2
         end if
         blocks(v)=blocks(v)+0.5_wp*q(i)*b*q(j)
      end do
   end do
   call rounding(q,uncertain,ge_bound,hb_bound,shape_bound)
   call rounding(qzero,real_uncertain,zero_ge_bound,zero_hb_bound,zero_shape_bound)
   write(*,'(a,4(1x,es25.17))') 'COMPONENTS',full
   write(*,'(a,4(1x,es25.17))') 'REPEAT_COMPONENTS',repeated
   write(*,'(a,4(1x,es25.17))') 'CAP_ZERO_COMPONENTS',zero
   write(*,'(a,3(1x,es25.17))') 'ELECTROSTATIC_BLOCKS',blocks
   write(*,'(a,3(1x,es25.17))') 'SHAPE_BLOCKS',shape_blocks
   write(*,'(a,2(1x,es25.17))') 'SHAPE',0.5_wp*coef*sum(q)**2,0.5_wp*coef*sum(qzero)**2
   write(*,'(a,3(1x,es25.17))') 'ROUNDING_BOUNDS',ge_bound,hb_bound,shape_bound
   write(*,'(a,3(1x,es25.17))') 'CAP_ZERO_ROUNDING_BOUNDS',zero_ge_bound,zero_hb_bound,zero_shape_bound
   write(*,'(a,3(1x,es25.17))') 'CHARGE_SUMS',sum(q),qr,qc
   write(*,'(a,4(1x,es25.17))') 'CAVITY_GLOBAL',born%aDet,born%keps,born%alpbet,born%dielectricConst
   write(*,'(a,3(1x,i0))') 'COUNTS',n,nc,count(cap)
   open(newunit=u,file='cavity.tsv',status='new',action='write')
   write(u,'(a)') '# index Z cap uncertain x_bohr y_bohr z_bohr Born_radius vdwr rho sasa hbw'
   do i=1,n
      write(u,'(4(i0,1x),8(es25.17,1x))') i,num(i),merge(1,0,cap(i)),merge(1,0,uncertain(i)), &
         xyz(:,i),born%brad(i),born%vdwr(i),born%rho(i),born%sasa(i),hb(i)
   end do
   close(u)
contains
   function electrostatic(i,j) result(value)
      integer,intent(in) :: i,j
      real(wp) :: value
      value=born%bornMat(i,j)
      if (i==j) value=value-2.0_wp*hb(i)
   end function
   subroutine rounding(charges,mask,be,bh,bs)
      real(wp),intent(in) :: charges(:)
      logical,intent(in) :: mask(:)
      real(wp),intent(out) :: be,bh,bs
      real(wp) :: linear,quadratic,product,epsq
      integer :: a,b,ii,jj
      linear=0.0_wp;quadratic=0.0_wp
      do a=1,nc
         ii=core(a)
         if (.not.mask(ii)) cycle
         product=0.0_wp
         do b=1,nc
            jj=core(b)
            product=product+electrostatic(ii,jj)*charges(jj)
            if (mask(jj)) quadratic=quadratic+abs(electrostatic(ii,jj))
         end do
         linear=linear+abs(product)
      end do
      be=dq*linear+0.5_wp*dq*dq*quadratic
      bh=2.0_wp*dq*sum(abs(hb*charges),mask=mask)+dq*dq*sum(abs(hb),mask=mask)
      epsq=dq*count(mask)
      bs=abs(coef)*(abs(sum(charges))*epsq+0.5_wp*epsq**2)
   end subroutine
end program
