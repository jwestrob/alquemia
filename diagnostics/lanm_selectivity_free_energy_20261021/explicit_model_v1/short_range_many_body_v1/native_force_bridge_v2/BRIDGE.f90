! Source-only finite vacuum bridge; not built or numerically qualified.
! Link the two private fixed-form clones before the exact original archive.
module srmb_force_state
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use atoms
  use atomid
  use bound
  use chgpot
  use deriv
  use energi
  use extfld
  use freeze
  use ielscf
  use keys
  use limits
  use mplpot
  use mpole
  use polar
  use polpcg
  use polpot
  use poltcg
  use potent
  use units
  use uprior
  use usage
  use virial
  implicit none
  real*8, allocatable :: savedxyz(:,:), savedalpha(:), sd(:,:), sp(:,:), a0(:,:)
  integer, allocatable :: savedsite(:), savedtype(:)
  logical :: prepared=.false., consumed=.false.
  integer :: direct_calls=0, mutual_calls=0, gradient_calls=0
contains
  subroutine insist(ok,message)
    logical, intent(in) :: ok
    character(*), intent(in) :: message
    if (.not.ok) then
      write(*,'(a)') 'SRMB_ERROR '//message
      error stop 1
    endif
  end subroutine

  subroutine guard()
    integer :: ii, atom
    logical, allocatable :: seen(:)
    ! Preserve all native switches; report them even if a later guard rejects.
    write(*,'(a,5(1x,l1))') 'SRMB_GUARD_CHGTRN_RXNFLD_ORBIT_EXTRA_METAL', &
                          use_chgtrn,use_rxnfld,use_orbit,use_extra,use_metal
    write(*,'(a)') 'SRMB_EXTRA_POLICY SOURCE_BOUND_ZERO_ONLY_EXTRA1_28ee5bf8'
    flush(6)
    call insist(n>=1.and.n<=28.and.npole==n,'finite 1..28 all-atom response scope')
    call insist(poltyp=='MUTUAL'.and.use_thole.and.use_polar.and.use_mpole,'MUTUAL Thole multipoles required')
    call insist(.not.(use_bounds.or.use_replica.or.use_polymer.or.use_ewald.or.use_list.or.use_mlist.or.use_ulist), &
                'vacuum direct-loop scope only')
    call insist(.not.(use_solv.or.use_chgpen.or.use_expol.or.use_tholed.or.use_chgflx), &
                'unsupported modified electrostatics')
    call insist(.not.(use_pred.or.use_ielscf.or.use_exfreq.or.use_exfld),'no history/dynamics/external field')
    ! field.f defaults use_extra=T. Exact linked extra1.f only zeros ex/dex;
    ! kextra.f is empty. The source/build contract pins this archive and noop.
    ! Keep use_extra and its original gradient dispatch, never turn physics off.
    call insist(.not.(use_chgtrn.or.use_rxnfld.or.use_orbit.or.use_metal), &
                'unsupported additional coupled terms')
    call insist(.not.(use_charge.or.use_chgdpl.or.use_dipole),'multipole permanent electrostatics only')
    call insist(mpolecut>=1d12.and.dielec==1d0,'baseline cutoff and unit dielectric')
    call insist(nwat4==0.and.nuse==n.and.all(use(1:n)),'real active atoms, no four-site redistribution')
    call insist(allocated(uind).and.allocated(uinp).and.allocated(udir).and.allocated(udirp), &
                'mechanic must allocate polarization arrays')
    call insist(all(ieee_is_finite(polarity(1:n))).and.all(polarity(1:n)>=1d-8),'positive finite intrinsic alpha')
    call insist(all(ieee_is_finite(x(1:n))).and.all(ieee_is_finite(y(1:n))).and.all(ieee_is_finite(z(1:n))), &
                'finite geometry')
    allocate(seen(n));seen=.false.
    do ii=1,npole
      atom=ipole(ii)
      call insist(atom>=1.and.atom<=n,'site index range')
      call insist(.not.seen(atom).and.douind(atom),'unique active response sites')
      seen(atom)=.true.
    enddo
    deallocate(seen)
  end subroutine

  subroutine prepare_baseline()
    intrinsic :: maxval
    integer :: ii,jj,atom,source,axis,row,col,j,iu,ios
    real*8, allocatable :: fd(:,:),fp(:,:)
    real*8 :: dualerr,symerr,permanent(3),charge
    call insist(.not.prepared,'prepare exactly once')
    call guard()
    ! MUTUAL does not use TCG arrays; epolar1a nevertheless loops to tcgnab.
    tcgnab=0
    call chkpole
    call rotpole('MPOLE')
    allocate(savedxyz(3,n),savedalpha(n),savedsite(npole),savedtype(n))
    savedxyz(1,:)=x(1:n);savedxyz(2,:)=y(1:n);savedxyz(3,:)=z(1:n)
    savedalpha=polarity(1:n);savedsite=ipole(1:npole);savedtype=type(1:n)
    allocate(sd(3,n),sp(3,n),fd(3,n),fp(3,n),a0(3*npole,3*npole))
    uind(:,1:n)=0d0;uinp(:,1:n)=0d0;udir(:,1:n)=0d0;udirp(:,1:n)=0d0
    call dfield0a(sd,sp)
    direct_calls=direct_calls+1
    call insist(all(ieee_is_finite(sd)).and.all(ieee_is_finite(sp)),'finite direct fields')
    do ii=1,npole
      atom=ipole(ii)
      udir(:,atom)=polarity(atom)*sd(:,atom)
      udirp(:,atom)=polarity(atom)*sp(:,atom)
    enddo
    ! Exactly the original exportresponsematrix unit-column mapping, no SCF.
    dualerr=0d0
    do jj=1,npole
      source=ipole(jj)
      do axis=1,3
        col=3*(jj-1)+axis
        uind(:,1:n)=0d0;uinp(:,1:n)=0d0
        uind(axis,source)=1d0;uinp(axis,source)=1d0
        call ufield0a(fd,fp)
        mutual_calls=mutual_calls+1
        do ii=1,npole
          atom=ipole(ii)
          do j=1,3
            row=3*(ii-1)+j
            a0(row,col)=-fd(j,atom)
            if(row==col) a0(row,col)=a0(row,col)+1d0/polarity(atom)
            dualerr=max(dualerr,abs(fd(j,atom)-fp(j,atom)))
          enddo
        enddo
      enddo
    enddo
    uind(:,1:n)=0d0;uinp(:,1:n)=0d0
    symerr=maxval(abs(a0-transpose(a0)))
    call insist(all(ieee_is_finite(a0)).and.ieee_is_finite(dualerr),'finite raw operator')
    permanent=0d0;charge=0d0
    do ii=1,npole
      atom=ipole(ii)
      charge=charge+rpole(1,atom)
      permanent=permanent+rpole(1,atom)*savedxyz(:,atom)+rpole(2:4,atom)
    enddo
    open(newunit=iu,file='baseline_state.dat',status='new',action='write',iostat=ios)
    call insist(ios==0,'baseline output collision/open error')
    write(iu,'(a)') 'SRMB_BASELINE_V1'
    write(iu,'(a,3(1x,i0))') 'COUNTS',n,npole,3*npole
    write(iu,'(a,4(1x,es26.17e3))') 'CONSTANTS',electric,dielec,debye,bohr
    write(iu,'(a,2(1x,es26.17e3))') 'OPERATOR_CHECKS',symerr,dualerr
    write(iu,'(a)') 'PERMANENT_ORIGIN_A 0 0 0'
    write(iu,'(a,3(1x,es26.17e3))') 'PERMANENT_MU_EA',permanent
    write(iu,'(a,1x,es26.17e3)') 'CHARGE',charge
    write(iu,'(a)') 'MODE VACUUM MUTUAL THOLE DIELECTRIC_1 ZERO_FIELD ALL_ACTIVE NO_VIRTUAL_SITES'
    write(iu,'(a,3(1x,l1))') 'FLAGS_DPEQUAL_PCGGUESS_TCGGUESS',dpequal,pcgguess,tcgguess
    write(iu,'(a,3(1x,i0))') 'CALLS',direct_calls,mutual_calls,0
    do ii=1,npole
      atom=ipole(ii)
      write(iu,'(a,3(1x,i0),4(1x,es26.17e3))') 'SITE',ii,atom,type(atom), &
                 polarity(atom),x(atom),y(atom),z(atom)
      do axis=1,3
        row=3*(ii-1)+axis
        write(iu,'(a,1x,i0,2(1x,es26.17e3))') 'VECTOR',row,sd(axis,atom),sp(axis,atom)
        write(iu,'(a,1x,i0,1x,*(es26.17e3,1x))') 'A',row,a0(row,:)
      enddo
    enddo
    write(iu,'(a)') 'END_SRMB_BASELINE_V1'
    close(iu)
    deallocate(fd,fp)
    prepared=.true.
    write(*,'(a)') 'SRMB_READY'
    flush(6)
  end subroutine

  subroutine fixed_gradient(filename)
    character(*), intent(in) :: filename
    character(80) :: tag
    integer :: iu,ios,nn,ns,ii,jj,atom,at,ty,out,j
    real*8 :: al,xyz(3),edread(3),epread(3),ur(3),vr(3),energy
    real*8, allocatable :: us(:,:),vs(:,:),derivs(:,:)
    call insist(prepared.and..not.consumed,'one contraction after prepare')
    call guard()
    call insist(tcgnab==0,'TCG scratch count must stay zero')
    call insist(all(x(1:n)==savedxyz(1,:)).and.all(y(1:n)==savedxyz(2,:)).and. &
                all(z(1:n)==savedxyz(3,:)).and.all(type(1:n)==savedtype).and. &
                all(polarity(1:n)==savedalpha).and.all(ipole(1:npole)==savedsite),'prepared state changed')
    allocate(us(3,n),vs(3,n),derivs(3,n));us=0d0;vs=0d0
    open(newunit=iu,file=trim(filename),status='old',action='read',iostat=ios)
    call insist(ios==0,'supplied uv file unavailable')
    read(iu,*,iostat=ios) tag
    call insist(ios==0.and.tag=='SRMB_UV_V1','uv header')
    read(iu,*,iostat=ios) tag,nn,ns
    call insist(ios==0.and.tag=='COUNTS'.and.nn==n.and.ns==npole,'uv counts')
    do ii=1,npole
      atom=ipole(ii)
      read(iu,*,iostat=ios) tag,jj,at,ty,al,xyz
      call insist(ios==0.and.tag=='SITE','uv site row')
      call insist(jj==ii.and.at==atom.and.ty==type(atom).and.al==polarity(atom).and. &
                  all(xyz==savedxyz(:,atom)),'uv exact site/geometry/alpha identity')
      read(iu,*,iostat=ios) tag,edread
      call insist(ios==0.and.tag=='ED'.and.all(edread==sd(:,atom)),'uv direct field identity')
      read(iu,*,iostat=ios) tag,epread
      call insist(ios==0.and.tag=='EP'.and.all(epread==sp(:,atom)),'uv energy field identity')
      read(iu,*,iostat=ios) tag,ur
      call insist(ios==0.and.tag=='U'.and.all(ieee_is_finite(ur)),'finite supplied u')
      read(iu,*,iostat=ios) tag,vr
      call insist(ios==0.and.tag=='V'.and.all(ieee_is_finite(vr)),'finite supplied v')
      us(:,atom)=ur;vs(:,atom)=vr
    enddo
    read(iu,*,iostat=ios) tag
    call insist(ios==0.and.tag=='END_SRMB_UV_V1','uv closing tag')
    read(iu,*,iostat=ios) tag
    call insist(ios<0,'uv trailing content')
    close(iu)
    consumed=.true.
    uind(:,1:n)=us;uinp(:,1:n)=vs
    call srmb_native_gradient(energy,derivs)
    gradient_calls=gradient_calls+1
    call insist(ex==0d0.and.all(dex(:,1:n)==0d0),'source-bound no-op extra1 changed energy or gradient')
    write(*,'(a)') 'SRMB_EXTRA_ZERO_ENERGY_GRADIENT T'
    flush(6)
    call insist(all(uind(:,1:n)==us).and.all(uinp(:,1:n)==vs),'gradient changed supplied dipoles')
    call insist(ieee_is_finite(energy).and.ieee_is_finite(ep).and.all(ieee_is_finite(derivs)).and. &
                all(ieee_is_finite(dep(:,1:n))).and.all(ieee_is_finite(vir)),'finite force readback')
    open(newunit=out,file='force_bridge_result.dat',status='new',action='write',iostat=ios)
    call insist(ios==0,'force output collision/open error')
    write(out,'(a)') 'SRMB_FIXED_UV_GRADIENT_V1'
    write(out,'(a,3(1x,i0))') 'COUNTS',n,npole,3*npole
    write(out,'(a,4(1x,i0))') 'CALLS',direct_calls,mutual_calls,gradient_calls,0
    write(out,'(a,1x,es26.17e3)') 'OTHER_ENERGY',energy-ep
    write(out,'(a,1x,es26.17e3)') 'RAW_REDUCED_POLARIZATION',ep
    write(out,'(a,1x,es26.17e3)') 'RAW_NATIVE_ENERGY',energy
    write(out,'(a)') 'SUPPLIED_UV_UNCHANGED T'
    do ii=1,npole
      atom=ipole(ii)
      write(out,'(a,2(1x,i0),6(1x,es26.17e3))') 'UV',ii,atom,us(:,atom),vs(:,atom)
    enddo
    do atom=1,n
      write(out,'(a,1x,i0,3(1x,es26.17e3))') 'BASELINE_GRADIENT',atom,derivs(:,atom)
      write(out,'(a,1x,i0,3(1x,es26.17e3))') 'POLAR_CONTRACTION_GRADIENT',atom,dep(:,atom)
    enddo
    do j=1,3
      write(out,'(a,1x,i0,3(1x,es26.17e3))') 'RAW_NATIVE_VIRIAL',j,vir(j,:)
    enddo
    write(out,'(a)') 'END_SRMB_FIXED_UV_GRADIENT_V1'
    close(out)
    deallocate(us,vs,derivs)
  end subroutine
end module

program srmb_bridge
  use srmb_force_state
  implicit none
  character(1024) :: filename
  integer :: ios
  call initial
  call getxyz
  call mechanic
  call prepare_baseline()
  read(*,'(a)',iostat=ios) filename
  call insist(ios==0.and.len_trim(filename)>0,'stdin supplied-uv pathname required')
  call fixed_gradient(trim(filename))
  call final
end program
