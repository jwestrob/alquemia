c     Regularized local response: isolated original CPU MUTUAL route.
c     New h/H are separate from intrinsic alpha and native Thole T.
      subroutine localresponsesetup
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
      integer i,j,k,next,ia,ib,slot,ios,extra,nmet
      real*8 a,b,c,d,nanval,dom
      logical bad,seen(2),sites(28)
      character*240 record,keyword,string
      lractive = .false.
      lrcoef = 0.0d0
      lrm = 0
      seen = .false.
      nanval = ieee_value(0.0d0,ieee_quiet_nan)
      dom = log(4.0d0)
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         if (trim(keyword) .ne. 'POLAR-LOCAL-RESPONSE') cycle
         ia = -1
         ib = -1
         a = nanval
         b = nanval
         c = nanval
         d = nanval
         string = record(next:240)
         read (string,*,iostat=ios) ia,ib,a,b,c,d
         bad = ios.ne.0 .or. (ia.ne.500.and.ia.ne.501)
     &       .or. ib.ne.36
         bad = bad .or. .not.all(ieee_is_finite((/a,b,c,d/)))
         if (bad) then
            write (iout,*) 'LOCAL RESPONSE: invalid record'
            call fatal
         end if
         read (string,*,iostat=ios) ia,ib,a,b,c,d,extra
         if (ios .ge. 0) then
            write (iout,*) 'LOCAL RESPONSE: extra record tokens'
            call fatal
         end if
         slot = ia - 499
         if (seen(slot) .or. abs(a).gt.1.0d0 .or.
     &       abs(b).gt.1.0d0 .or.
     &       abs(c)+abs(d).gt.0.5d0*dom+1.0d-15) then
            write (iout,*) 'LOCAL RESPONSE: duplicate/domain error'
            call fatal
         end if
         seen(slot) = .true.
         lrcoef(:,slot) = (/a,b,c,d/)
         if (any(lrcoef(:,slot).ne.0.0d0)) lractive=.true.
      end do
c     Absent and identically zero records use literal original paths.
      if (.not.lractive) return
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         if (trim(keyword).eq.'POLAR-PAIR-RESPONSE') then
            write (iout,*) 'LOCAL RESPONSE: mixed models unsupported'
            call fatal
         end if
      end do
      bad = use_bounds .or. use_replica .or. use_polymer
     &      .or. use_ewald .or. use_mlist .or. use_ulist
     &      .or. use_solv .or. use_chgpen .or. use_expol
     &      .or. use_tholed .or. use_chgflx .or. use_exfreq
     &      .or. use_pred .or. use_ielscf
     &      .or. poltyp.ne.'MUTUAL' .or. .not.use_thole
     &      .or. .not.use_polar .or. mpolecut.lt.1.0d12
     &      .or. n.lt.1 .or. n.gt.28 .or. npole.ne.n
     &      .or. nuse.ne.n
      if (bad) then
         write (iout,*) 'LOCAL RESPONSE: unsupported native mode'
         call fatal
      end if
      sites = .false.
      nmet = 0
      do j = 1, npole
         i = ipole(j)
         if (i.lt.1 .or. i.gt.n) then
            write (iout,*) 'LOCAL RESPONSE: invalid site index'
            call fatal
         end if
         if (sites(i) .or. .not.douind(i)) then
            write (iout,*) 'LOCAL RESPONSE: inactive/duplicate site'
            call fatal
         end if
         sites(i) = .true.
         if (.not.ieee_is_finite(polarity(i)) .or.
     &       polarity(i).lt.1.0d-8 .or.
     &       .not.all(ieee_is_finite((/x(i),y(i),z(i)/)))) then
            write (iout,*) 'LOCAL RESPONSE: invalid alpha/geometry'
            call fatal
         end if
         if (type(i).eq.500 .or. type(i).eq.501) then
            nmet = nmet + 1
            lrm = i
            if (rpole(1,i).ne.3.0d0) then
               write (iout,*) 'LOCAL RESPONSE: metal charge must be3'
               call fatal
            end if
            if (n12(i).ne.0) then
               write (iout,*) 'LOCAL RESPONSE: bonded metal'
               call fatal
            end if
         else if (type(i).eq.36) then
            if (n12(i).ne.2) then
               write (iout,*) 'LOCAL RESPONSE: incomplete water O'
               call fatal
            end if
            do k = 1, 2
               if (type(i12(k,i)).ne.37) then
                  write (iout,*) 'LOCAL RESPONSE: water O graph'
                  call fatal
               end if
            end do
         else if (type(i).eq.37) then
            if (n12(i).ne.1) then
               write (iout,*) 'LOCAL RESPONSE: incomplete water H'
               call fatal
            end if
            if (type(i12(1,i)).ne.36) then
               write (iout,*) 'LOCAL RESPONSE: water H graph'
               call fatal
            end if
         else
            write (iout,*) 'LOCAL RESPONSE: unsupported atom type'
            call fatal
         end if
      end do
      if (nmet.gt.1) then
         write (iout,*) 'LOCAL RESPONSE: single metal scope'
         call fatal
      end if
      if (use_exfld) then
         if (.not.all(ieee_is_finite(texfld))) then
            write (iout,*) 'LOCAL RESPONSE: nonfinite static field'
            call fatal
         end if
      end if
      if (allocated(lrh)) deallocate(lrh,lrbeta,lrdb,lrdh)
      allocate(lrh(3,n),lrbeta(n),lrdb(3,n),lrdh(3,3,n,n))
      lrh=0.0d0
      lrbeta=0.0d0
      lrdb=0.0d0
      lrdh=0.0d0
      do slot = 1, 2
         if (seen(slot)) then
            write(iout,'(a,2i8,4es26.17e3)')
     &       'LOCAL_RESPONSE_PARAMETERS ',slot+499,36,lrcoef(:,slot)
         end if
      end do
      call localresponsegeometry
      call localresponsespectrum
      return
      end
c
      subroutine localresponsegeometry
      use, intrinsic :: ieee_arithmetic
      use atoms
      use iounit
      use polar
      implicit none
      integer i,j,k,l,m,slot,nw,oi,oj,ow(28)
      real*8 a,b,c,d,r,t,p,fzero,s,q,eta,g,sigma,other
      real*8 prod,kappa,gamma,zeta,phi,phir,si,dereta,dsig
      real*8 ff(28),fr(28),ph(28),pr(28),et(28),os(28)
      real*8 rr(28),vec(3,28),df(3,28),ds(28),jac(3,3)
      real*8 pref,w,wp,bij
      if (lrm.eq.0) return
      m=lrm
      slot=type(m)-499
      a=lrcoef(1,slot)
      b=lrcoef(2,slot)
      c=lrcoef(3,slot)
      d=lrcoef(4,slot)
      nw=0
      fzero=1.0d0/(((1.0d0+5.0d0+25.0d0/3.0d0)
     &              *exp(-5.0d0))**2)
      do oi=1,n
         if (type(oi).ne.36) cycle
         nw=nw+1
         i=nw
         ow(i)=oi
         vec(:,i)=(/x(oi)-x(m),y(oi)-y(m),z(oi)-z(m)/)
         r=sqrt(dot_product(vec(:,i),vec(:,i)))
         rr(i)=r
         t=2.0d0*r
         p=1.0d0+t+t*t/3.0d0
         ff(i)=fzero*(p*exp(-t))**2
         fr(i)=-(4.0d0/3.0d0)*ff(i)*t*(1.0d0+t)/p
         if (r.eq.0.0d0) then
            df(:,i)=0.0d0
         else
            df(:,i)=fr(i)*vec(:,i)/r
         end if
         gamma=thlval(jpolar(m),jpolar(oi))
         prod=pdamp(m)*pdamp(oi)
         if (.not.ieee_is_finite(gamma) .or. gamma.le.0.0d0
     &       .or. .not.ieee_is_finite(prod) .or. prod.le.0.0d0)
     &       then
            write(iout,*) 'LOCAL RESPONSE: invalid fixed damping'
            call fatal
         end if
         kappa=gamma/prod**3
         zeta=kappa*r**3
c        Small-z series evaluate the same phi, no physical clipping.
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
         ph(i)=phi
         pr(i)=phir
         write(iout,'(a,2i8,3es26.17e3)')
     &     'LOCAL_RESPONSE_PAIR ',m,oi,gamma,kappa,r
      end do
      if (nw.eq.0) then
         write(iout,'(a,i8,4es26.17e3)')
     &    'LOCAL_RESPONSE_ENVIRONMENT ',m,0.0d0,0.0d0,0.0d0,0.0d0
         return
      end if
      s=sum(ff(1:nw))
      q=0.0d0
      do i=1,nw
         os(i)=0.0d0
         do j=1,nw
            if (j.ne.i) os(i)=os(i)+ff(j)
            if (j.gt.i) q=q+2.0d0*ff(i)*ff(j)
         end do
         et(i)=os(i)/(1.0d0+os(i))
      end do
      eta=q/(1.0d0+s*s)
      g=2.0d0*s/(1.0d0+s)
      sigma=g*(c+d*eta)
      lrbeta(m)=(exp(sigma)-1.0d0)/polarity(m)
      do j=1,nw
         dereta=(2.0d0*os(j)*(1.0d0+s*s)-2.0d0*s*q)
     &                /(1.0d0+s*s)**2
         dsig=2.0d0*(c+d*eta)/(1.0d0+s)**2+g*d*dereta
         ds(j)=exp(sigma)*dsig/polarity(m)
         oj=ow(j)
         lrdb(:,oj)=ds(j)*df(:,j)
         lrdb(:,m)=lrdb(:,m)-lrdb(:,oj)
      end do
      do i=1,nw
         oi=ow(i)
         pref=a+b*et(i)
         lrh(:,oi)=3.0d0*ff(i)*pref*ph(i)*vec(:,i)
         do j=1,nw
            oj=ow(j)
            if (i.eq.j) then
               bij=pref
            else
               bij=ff(i)*b/(1.0d0+os(i))**2
            end if
            do k=1,3
               do l=1,3
                  jac(k,l)=3.0d0*ph(i)*vec(k,i)*bij*df(l,j)
                  if (i.eq.j) then
                     if (k.eq.l) jac(k,l)=jac(k,l)
     &                    +3.0d0*ff(i)*pref*ph(i)
                     if (rr(i).gt.0.0d0) jac(k,l)=jac(k,l)
     &                    +3.0d0*ff(i)*pref*pr(i)
     &                    *vec(k,i)*vec(l,i)/rr(i)
                  end if
               end do
            end do
            lrdh(:,:,oi,oj)=jac
            lrdh(:,:,oi,m)=lrdh(:,:,oi,m)-jac
         end do
      end do
      if (.not.all(ieee_is_finite(lrh)) .or.
     &    .not.all(ieee_is_finite(lrbeta)) .or.
     &    .not.all(ieee_is_finite(lrdb)) .or.
     &    .not.all(ieee_is_finite(lrdh))) then
         write(iout,*) 'LOCAL RESPONSE: nonfinite geometry functions'
         call fatal
      end if
      write(iout,'(a,i8,4es26.17e3)')
     & 'LOCAL_RESPONSE_ENVIRONMENT ',m,s,eta,sigma,lrbeta(m)
      return
      end
c
c     Native Jacobi on the actual normalized PCG operator, pre-solve.
c     induce0a clears both trial vectors immediately afterward.
      subroutine localresponsespectrum
      use, intrinsic :: ieee_arithmetic
      use atoms
      use cell
      use iounit
      use mpole
      use polar
      use shunt
      implicit none
      integer i,j,k,ii,jj,ir,ic,dim,oldcell
      real*8 oldshunt(18),lo,hi,sym,err,ortho,value
      real*8, allocatable :: am(:,:),original(:,:)
      real*8, allocatable :: ev(:),eigvec(:,:)
      real*8, allocatable :: fld(:,:),fldp(:,:)
      dim=3*npole
      allocate(am(dim,dim),original(dim,dim))
      allocate(ev(dim),eigvec(dim,dim))
      allocate(fld(3,n),fldp(3,n))
      oldcell=ncell
      oldshunt=(/off,off2,cut,cut2,c0,c1,c2,c3,c4,c5,
     &           f0,f1,f2,f3,f4,f5,f6,f7/)
      do jj=1,npole
         k=ipole(jj)
         do j=1,3
            ic=3*(jj-1)+j
            uind(:,1:n)=0.0d0
            uinp(:,1:n)=0.0d0
            uind(j,k)=1.0d0
            uinp(j,k)=1.0d0
            call ufield0a(fld,fldp)
            do ii=1,npole
               i=ipole(ii)
               do ir=1,3
                  am(3*(ii-1)+ir,ic)=-fld(ir,i)
                  if (3*(ii-1)+ir.eq.ic)
     &              am(ic,ic)=am(ic,ic)+1.0d0/polarity(i)
                  am(3*(ii-1)+ir,ic)=am(3*(ii-1)+ir,ic)
     &                *sqrt(polarity(i)*polarity(k))
               end do
            end do
         end do
      end do
      sym=0.0d0
      do i=1,dim
         do j=1,dim
            sym=max(sym,abs(am(i,j)-am(j,i)))
         end do
      end do
      if (.not.all(ieee_is_finite(am))) then
         write(iout,*) 'LOCAL RESPONSE: nonfinite response operator'
         call fatal
      end if
c     No symmetrization, eigenvalue clipping or alternate solver.
      original=am
      call jacobi(dim,am,ev,eigvec)
      lo=minval(ev)
      hi=-huge(hi)
      do i=1,dim
         hi=max(hi,ev(i))
      end do
      err=0.0d0
      ortho=0.0d0
      do i=1,dim
         do j=1,dim
            value=dot_product(original(i,:),eigvec(:,j))
     &              -ev(j)*eigvec(i,j)
            err=max(err,abs(value))
            value=dot_product(eigvec(:,i),eigvec(:,j))
            if (i.eq.j) value=value-1.0d0
            ortho=max(ortho,abs(value))
         end do
      end do
      write(iout,'(a,5es26.17e3)')
     &    'LOCAL_RESPONSE_SPECTRUM ',lo,hi,sym,err,ortho
      if (sym.gt.1.0d-10 .or. ortho.gt.1.0d-10 .or.
     &    err.gt.1.0d-10*max(1.0d0,hi) .or.
     &    .not.all(ieee_is_finite(eigvec))) then
         write(iout,*) 'LOCAL RESPONSE: eigensystem arithmetic error'
         call fatal
      end if
      if (.not.all(ieee_is_finite(ev)) .or. hi.le.0.0d0
     &    .or. lo.le.1.0d-6*hi) then
         write(iout,*) 'LOCAL RESPONSE: normalized SPD failure'
         call fatal
      end if
      uind(:,1:n)=0.0d0
      uinp(:,1:n)=0.0d0
      ncell=oldcell
      off=oldshunt(1)
      off2=oldshunt(2)
      cut=oldshunt(3)
      cut2=oldshunt(4)
      c0=oldshunt(5)
      c1=oldshunt(6)
      c2=oldshunt(7)
      c3=oldshunt(8)
      c4=oldshunt(9)
      c5=oldshunt(10)
      f0=oldshunt(11)
      f1=oldshunt(12)
      f2=oldshunt(13)
      f3=oldshunt(14)
      f4=oldshunt(15)
      f5=oldshunt(16)
      f6=oldshunt(17)
      f7=oldshunt(18)
      deallocate(am,original,ev,eigvec,fld,fldp)
      return
      end
c
      subroutine localresponsereport
      use atoms
      use extfld
      use iounit
      use mpole
      use polar
      implicit none
      integer i,j,ii
      real*8 dual,rhsdual
      dual=0.0d0
      rhsdual=0.0d0
      do ii=1,npole
         i=ipole(ii)
         do j=1,3
            dual=max(dual,abs(uind(j,i)-uinp(j,i)))
            rhsdual=max(rhsdual,abs(udir(j,i)-udirp(j,i))
     &                             /polarity(i))
         end do
         write(iout,'(a,i8,6es26.17e3)')
     &    'LOCAL_RESPONSE_INDUCED ',i,uind(:,i),uinp(:,i)
      end do
      write(iout,'(a,2es26.17e3)')
     &    'LOCAL_RESPONSE_DUAL ',rhsdual,dual
      return
      end
c
      subroutine localresponsegradient
      use atoms
      use chgpot
      use deriv
      use polar
      use virial
      use units
      implicit none
      integer i,j,k,l
      real*8 fac,cross,gg(3),vij(3),vg(3,3)
      if (.not.lractive .or. lrm.eq.0) return
      fac=0.5d0*electric/dielec
      cross=dot_product(uind(:,lrm),uinp(:,lrm))
      do j=1,n
         gg=cross*lrdb(:,j)
         do i=1,n
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
