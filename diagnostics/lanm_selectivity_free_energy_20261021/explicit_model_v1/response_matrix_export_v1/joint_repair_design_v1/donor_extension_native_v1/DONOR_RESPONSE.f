c     Separate acetate O403 response; original water helper is retained.
      logical function donorresponsefixture ()
      use atoms
      use couple
      implicit none
      integer i,j,k,l,nmet,ct(4),neighbors(4)
      donorresponsefixture=.false.
      if (n.ne.7 .and. n.ne.8) return
      ct=0
      nmet=0
      do i=1,n
         if (type(i).eq.500 .or. type(i).eq.501) then
            nmet=nmet+1
            if (n12(i).ne.0) return
         else if (type(i).ge.401 .and. type(i).le.404) then
            ct(type(i)-400)=ct(type(i)-400)+1
         else
            return
         end if
         if (n12(i).lt.0 .or. n12(i).gt.4) return
         do j=1,n12(i)
            k=i12(j,i)
            if (k.lt.1 .or. k.gt.n .or. k.eq.i) return
            do l=1,j-1
               if (i12(l,i).eq.k) return
            end do
         end do
      end do
      if (nmet.gt.1 .or. n.ne.7+nmet) return
      if (any(ct.ne.(/1,1,2,3/))) return
      do i=1,n
         neighbors=0
         do j=1,n12(i)
            k=i12(j,i)
            if (.not.any(i12(1:n12(k),k).eq.i)) return
            if (type(k).lt.401 .or. type(k).gt.404) return
            neighbors(type(k)-400)=neighbors(type(k)-400)+1
         end do
         if (type(i).eq.401) then
            if (any(neighbors.ne.(/0,1,0,3/))) return
         else if (type(i).eq.402) then
            if (any(neighbors.ne.(/1,0,2,0/))) return
         else if (type(i).eq.403) then
            if (any(neighbors.ne.(/0,1,0,0/))) return
         else if (type(i).eq.404) then
            if (any(neighbors.ne.(/1,0,0,0/))) return
         end if
      end do
      donorresponsefixture=.true.
      return
      end
c
      subroutine donorresponsesetup
      use, intrinsic :: ieee_arithmetic
      use atoms
      use bound
      use couple
      use extfld
      use ielscf
      use iounit
      use keys
      use limits
      use mplpot
      use mpole
      use polar
      use polpot
      use potent
      use uprior
      use usage
      implicit none
      integer i,j,next,ia,ib,slot,ios,extra,nmet,no,m
      real*8 a,c,nanval
      logical bad,seen(2),sites(8),donorresponsefixture
      character*240 record,keyword,string
      lrdonor=.false.
      drcoef=0.0d0
      seen=.false.
      nanval=ieee_value(0.0d0,ieee_quiet_nan)
      do i=1,nkey
         next=1
         record=keyline(i)
         call gettext(record,keyword,next)
         call upcase(keyword)
         if (trim(keyword).ne.'POLAR-DONOR-RESPONSE') cycle
         ia=-1
         ib=-1
         a=nanval
         c=nanval
         string=record(next:240)
         read(string,*,iostat=ios) ia,ib,a,c
         bad=ios.ne.0 .or. (ia.ne.500.and.ia.ne.501)
     &       .or. ib.ne.403
         bad=bad.or..not.all(ieee_is_finite((/a,c/)))
         if (bad) then
            write(iout,*) 'DONOR RESPONSE: invalid record'
            call fatal
         end if
         read(string,*,iostat=ios) ia,ib,a,c,extra
         if (ios.ge.0) then
            write(iout,*) 'DONOR RESPONSE: extra record tokens'
            call fatal
         end if
         slot=ia-499
         if (seen(slot)) then
            write(iout,*) 'DONOR RESPONSE: duplicate record'
            call fatal
         end if
         seen(slot)=.true.
         drcoef(:,slot)=(/a,c/)
      end do
c     No new arithmetic/guards on absent or exactly zero channel.
      if (all(drcoef.eq.0.0d0)) return
      nmet=0
      no=0
      m=0
      do i=1,n
         if (type(i).eq.500 .or. type(i).eq.501) then
            nmet=nmet+1
            m=i
         end if
         if (type(i).eq.403) no=no+1
      end do
      if (nmet.eq.0 .or. no.eq.0) return
      if (nmet.ne.1) then
         write(iout,*) 'DONOR RESPONSE: single metal scope'
         call fatal
      end if
      slot=type(m)-499
      if (all(drcoef(:,slot).eq.0.0d0)) return
      if (n.ne.8 .or. .not.donorresponsefixture()) then
         write(iout,*) 'DONOR RESPONSE: exact acetate graph required'
         call fatal
      end if
      do i=1,nkey
         next=1
         record=keyline(i)
         call gettext(record,keyword,next)
         call upcase(keyword)
         if (trim(keyword).eq.'POLAR-PAIR-RESPONSE' .or.
     &       trim(keyword).eq.'POLAR-LOCAL-RESPONSE') then
            write(iout,*) 'DONOR RESPONSE: mixed channels unsupported'
            call fatal
         end if
      end do
      bad=use_bounds .or. use_replica .or. use_polymer
     &    .or. use_ewald .or. use_mlist .or. use_ulist
     &    .or. use_solv .or. use_chgpen .or. use_expol
     &    .or. use_tholed .or. use_chgflx .or. use_exfreq
     &    .or. use_pred .or. use_ielscf
     &    .or. poltyp.ne.'MUTUAL' .or. .not.use_thole
     &    .or. .not.use_polar .or. mpolecut.lt.1.0d12
     &    .or. npole.ne.n .or. nuse.ne.n
      if (bad) then
         write(iout,*) 'DONOR RESPONSE: unsupported native mode'
         call fatal
      end if
      sites=.false.
      do j=1,npole
         i=ipole(j)
         if (i.lt.1 .or. i.gt.n) then
            write(iout,*) 'DONOR RESPONSE: invalid site index'
            call fatal
         end if
         if (sites(i) .or. .not.douind(i)) then
            write(iout,*) 'DONOR RESPONSE: inactive/duplicate site'
            call fatal
         end if
         sites(i)=.true.
         if (.not.ieee_is_finite(polarity(i)) .or.
     &       polarity(i).lt.1.0d-8 .or.
     &       .not.all(ieee_is_finite((/x(i),y(i),z(i)/)))) then
            write(iout,*) 'DONOR RESPONSE: invalid alpha/geometry'
            call fatal
         end if
      end do
      if (rpole(1,m).ne.3.0d0) then
         write(iout,*) 'DONOR RESPONSE: metal charge must be3'
         call fatal
      end if
      if (use_exfld) then
         if (.not.all(ieee_is_finite(texfld))) then
            write(iout,*) 'DONOR RESPONSE: nonfinite static field'
            call fatal
         end if
      end if
      lrm=m
      lractive=.true.
      lrdonor=.true.
      if (allocated(lrh)) deallocate(lrh,lrbeta,lrdb,lrdh)
      if (allocated(drdb)) deallocate(drdb)
      allocate(lrh(3,n),lrbeta(n),lrdb(3,n),lrdh(3,3,n,n))
      allocate(drdb(3,n,n))
      lrh=0.0d0
      lrbeta=0.0d0
      lrdb=0.0d0
      lrdh=0.0d0
      drdb=0.0d0
      do slot=1,2
         if (seen(slot)) then
            write(iout,'(a,2i8,2es26.17e3)')
     &       'DONOR_RESPONSE_PARAMETERS ',slot+499,403,drcoef(:,slot)
         end if
      end do
      call donorresponsegeometry
      call localresponsespectrum
      return
      end
c
      subroutine donorresponsegeometry
      use, intrinsic :: ieee_arithmetic
      use atoms
      use iounit
      use mpole
      use polar
      implicit none
      integer i,k,l,m,slot
      real*8 a,c,r,t,p,fzero,f,fr,phi,phir,gamma,prod,kappa
      real*8 zeta,si,w,wr,sigma,qm,vec(3),df(3),jac(3,3)
      m=lrm
      slot=type(m)-499
      a=drcoef(1,slot)
      c=drcoef(2,slot)
      qm=rpole(1,m)
      fzero=1.0d0/(((1.0d0+5.0d0+25.0d0/3.0d0)
     &              *exp(-5.0d0))**2)
      do i=1,n
         if (type(i).ne.403) cycle
         vec=(/x(i)-x(m),y(i)-y(m),z(i)-z(m)/)
         r=sqrt(dot_product(vec,vec))
         t=2.0d0*r
         p=1.0d0+t+t*t/3.0d0
         f=fzero*(p*exp(-t))**2
         fr=-(4.0d0/3.0d0)*f*t*(1.0d0+t)/p
         df=0.0d0
         if (r.gt.0.0d0) df=fr*vec/r
         gamma=thlval(jpolar(m),jpolar(i))
         prod=pdamp(m)*pdamp(i)
         if (.not.ieee_is_finite(gamma) .or. gamma.le.0.0d0
     &       .or. .not.ieee_is_finite(prod) .or. prod.le.0.0d0)
     &       then
            write(iout,*) 'DONOR RESPONSE: invalid fixed damping'
            call fatal
         end if
         kappa=gamma/prod**3
         zeta=kappa*r**3
c        Identical stable contact series to the frozen water helper.
         if (zeta.lt.1.0d-4) then
            phi=kappa*(1.0d0-zeta/2.0d0+zeta*zeta/6.0d0
     &         -zeta**3/24.0d0+zeta**4/120.0d0-zeta**5/720.0d0)
            phir=kappa*kappa*r*r*(-1.5d0+zeta
     &         -3.0d0*zeta*zeta/8.0d0+zeta**3/10.0d0
     &         -zeta**4/48.0d0)
         else
            si=1.0d0-exp(-zeta)
            phi=si/r**3
            phir=3.0d0*(zeta*exp(-zeta)-si)/r**4
         end if
         w=2.0d0*f/(1.0d0+f)
         wr=2.0d0*fr/(1.0d0+f)**2
         sigma=c*w
         if (.not.ieee_is_finite(sigma) .or.
     &       sigma.gt.log(huge(1.0d0))) then
            write(iout,*) 'DONOR RESPONSE: nonfinite stiffness'
            call fatal
         end if
         lrbeta(i)=(exp(sigma)-1.0d0)/polarity(i)
         if (r.gt.0.0d0) then
            drdb(:,i,i)=exp(sigma)*c*wr*vec/(polarity(i)*r)
            drdb(:,i,m)=-drdb(:,i,i)
         end if
         lrh(:,i)=a*qm*f*phi*vec
         do k=1,3
            do l=1,3
               jac(k,l)=a*qm*phi*vec(k)*df(l)
               if (k.eq.l) jac(k,l)=jac(k,l)+a*qm*f*phi
               if (r.gt.0.0d0) jac(k,l)=jac(k,l)
     &                  +a*qm*f*phir*vec(k)*vec(l)/r
            end do
         end do
         lrdh(:,:,i,i)=jac
         lrdh(:,:,i,m)=-jac
         write(iout,'(a,2i8,6es26.17e3)')
     &    'DONOR_RESPONSE_PAIR ',m,i,gamma,kappa,r,f,w,lrbeta(i)
         write(iout,'(a,i8,3es26.17e3)')
     &    'DONOR_RESPONSE_FIELD ',i,lrh(:,i)
      end do
      if (.not.all(ieee_is_finite(lrh)) .or.
     &    .not.all(ieee_is_finite(lrbeta)) .or.
     &    .not.all(ieee_is_finite(drdb)) .or.
     &    .not.all(ieee_is_finite(lrdh))) then
         write(iout,*) 'DONOR RESPONSE: nonfinite geometry functions'
         call fatal
      end if
      return
      end
c
      subroutine donorresponsegradient
      use atoms
      use chgpot
      use deriv
      use polar
      use virial
      use units
      implicit none
      integer i,j,k,l
      real*8 fac,gg(3),vij(3),cross
      if (.not.lrdonor .or. lrm.eq.0) return
      fac=0.5d0*electric/dielec
      do j=1,n
         gg=0.0d0
         do i=1,n
            cross=dot_product(uind(:,i),uinp(:,i))
            gg=gg+cross*drdb(:,i,j)
            do k=1,3
               do l=1,3
                  gg(l)=gg(l)-(uind(k,i)+uinp(k,i))*lrdh(k,l,i,j)
               end do
            end do
         end do
         gg=fac*gg
         dep(:,j)=dep(:,j)+gg
         vij=(/x(j)-x(lrm),y(j)-y(lrm),z(j)-z(lrm)/)
         do k=1,3
            do l=1,3
               vir(k,l)=vir(k,l)+0.5d0*(vij(k)*gg(l)
     &                                  +vij(l)*gg(k))
            end do
         end do
      end do
      return
      end
c
c     Reporting only, called after the complete native analytic gradient.
      subroutine donorresponsetotalgradient (energy,derivs)
      use atoms
      use iounit
      use polar
      use virial
      implicit none
      integer i
      real*8 energy,derivs(3,n)
      if (.not.lrdonor) return
      write(iout,'(a,es26.17e3)')
     & 'DONOR_RESPONSE_TOTAL_ENERGY ',energy
      do i=1,n
         write(iout,'(a,i8,3es26.17e3)')
     &    'DONOR_RESPONSE_TOTAL_GRADIENT ',i,derivs(:,i)
      end do
      do i=1,3
         write(iout,'(a,i8,3es26.17e3)')
     &    'DONOR_RESPONSE_TOTAL_VIRIAL ',i,vir(i,:)
      end do
      return
      end
c
c     Current rpole is permanent and already rotated, before moments adds u.
c     Lab-origin permanent moment in e Angstrom, even on zero/off exports.
      subroutine donorresponsepermanent
      use atoms
      use iounit
      use mpole
      implicit none
      integer i,j
      real*8 mu(3),charge
      logical eligible,donorresponsefixture
      eligible=donorresponsefixture()
      if (n.eq.1) eligible=type(1).eq.500.or.type(1).eq.501
      if (.not.eligible) return
      mu=0.0d0
      charge=0.0d0
      do j=1,npole
         i=ipole(j)
         mu=mu+rpole(1,i)*(/x(i),y(i),z(i)/)+rpole(2:4,i)
         charge=charge+rpole(1,i)
      end do
      write(iout,'(a,3es26.17e3)')
     & 'DONOR_RESPONSE_PERMANENT_LAB_EA ',mu
      write(iout,'(a,es26.17e3)')
     & 'DONOR_RESPONSE_CHARGE ',charge
      return
      end
