      subroutine pot(symb,x,y,z,v,dvdx,dvdy,dvdz,natom,maxatom)

c NEW "AUTOFIT" PESs

      implicit real*8(a-h,o-z)
      dimension x(maxatom),y(maxatom),z(maxatom)
      dimension x2(maxatom),y2(maxatom),z2(maxatom)
      dimension x3(maxatom),y3(maxatom),z3(maxatom)
      dimension dvdx(maxatom),dvdy(maxatom),dvdz(maxatom)
      dimension dvdx2(maxatom),dvdy2(maxatom),dvdz2(maxatom)
      dimension dvdx3(maxatom),dvdy3(maxatom),dvdz3(maxatom)
      parameter(autoev=27.2113961d0)
      parameter(autocmi=219474.63067d0)
      parameter(autoang=0.529177249d0)
      character*2 symb(maxatom),symb2(maxatom),symb3(maxatom)
      integer at(maxatom),at2(maxatom),at3(maxatom)
      integer itinkxyz
      logical ltinkm

      dimension tmpprint(50)
      common/tmp/tmpprint

      v=0.d0
      v2=0.d0
      v3=0.d0
      do i=1,natom
      symb2(i)="xx"
      symb3(i)="xx"
      x2(i)=0.d0
      y2(i)=0.d0
      z2(i)=0.d0
      x3(i)=0.d0
      y3(i)=0.d0
      z3(i)=0.d0
      dvdz2(i)=0.d0
      dvdy2(i)=0.d0
      dvdx2(i)=0.d0
      dvdz3(i)=0.d0
      dvdy3(i)=0.d0
      dvdx3(i)=0.d0
      dvdz(i)=0.d0
      dvdy(i)=0.d0
      dvdx(i)=0.d0
      enddo

      do i=1,natom
      at(i)=0
      if ((symb(i).eq."H").or.
     &    (symb(i).eq."h"))  at(i)=1       ! hydrogen
      if ((symb(i).eq."C").or.
     &    (symb(i).eq."c"))  at(i)=2       ! carbon
      if ((symb(i).eq."N").or.
     &    (symb(i).eq."N"))  at(i)=2       ! carbon
      if ((symb(i).eq."O").or.
     &    (symb(i).eq."o"))  at(i)=4       ! oxygen
      if ((symb(i).eq."HO").or.
     &    (symb(i).eq."Ho"))  at(i)=5       ! OH H
      if ((symb(i).eq."He").or.
     &    (symb(i).eq."he").or.
     &    (symb(i).eq."HE")) at(i)=21      ! helium
      if ((symb(i).eq."Ar").or.
     &    (symb(i).eq."ar").or.
     &    (symb(i).eq."AR")) at(i)=23      ! argon
      if ((symb(i).eq."N2").or.
     &    (symb(i).eq."n2")) at(i)=26      ! N2  ! label your bath atoms N2, not N
      if ((symb(i).eq."O2").or.
     &    (symb(i).eq."o2")) at(i)=29      ! 
      if ((symb(i).eq."O3").or.
     &    (symb(i).eq."o3")) at(i)=31      ! 
      if ((symb(i).eq."O4").or.
     &    (symb(i).eq."o4")) at(i)=31      ! 
      if ((symb(i).eq."H2").or.
     &    (symb(i).eq."H2")) at(i)=32      ! 
      if ((symb(i).eq."C2").or.
     &    (symb(i).eq."C2")) at(i)=33      ! 
      if ((symb(i).eq."C3").or.
     &    (symb(i).eq."C3")) at(i)=33      ! 
      if ((symb(i).eq."C4").or.
     &    (symb(i).eq."C4")) at(i)=33      ! 
      if ((symb(i).eq."O1").or.
     &    (symb(i).eq."o1")) at(i)=33      ! 
      if ((symb(i).eq."H1").or.
     &    (symb(i).eq."h1")) at(i)=33      ! 
      if ((symb(i).eq."N3").or.
     &    (symb(i).eq."n3")) at(i)=33      ! 
      if ((symb(i).eq."H3").or.
     &    (symb(i).eq."h3")) at(i)=33      ! 
      if (at(i).eq.0) then ! atom not found
           write(6,*)"Atom # ",i," (",symb(i),") not found"
           stop
      endif
      enddo

      natom2=0
      natom3=0
      ltinkm=.false.
      do i=1,natom
        if (at(i).le.20) then ! collect target atoms
          natom2=natom2+1
          x2(natom2)=x(i)
          y2(natom2)=y(i)
          z2(natom2)=z(i)
          symb2(natom2)=symb(i)
          at2(natom2)=at(i)
        else ! collect bath atoms
          natom3=natom3+1
          x3(natom3)=x(i)
          y3(natom3)=y(i)
          z3(natom3)=z(i)
          symb3(natom3)=symb(i)
          at3(natom3)=at(i)
        endif
        if (at(i).ge.30) ltinkm=.true.
      enddo

      if (natom3.ne.0.and.natom2.ne.0) 
     &   call lsfit(at,x,y,z,v,dvdx,dvdy,dvdz,natom,maxatom) !  Interaction PES

      itinkxyz=0
      if (natom3.eq.0) itinkxyz=1   ! only targets
      if (natom2.eq.0) itinkxyz=2   ! only baths
      if (ltinkm)
     &    call tinkerpot(symb,x,y,z,v2,dvdx2,dvdy2,dvdz2,
     &          natom,maxatom,itinkxyz)
      if (.not.ltinkm.and.natom2.ne.0)
     &    call tinkerpot(symb2,x2,y2,z2,v2,dvdx2,dvdy2,dvdz2,
     &          natom2,maxatom,itinkxyz)
      if (.not.ltinkm.and.natom3.ne.0)
     &   call bath(at3,x3,y3,z3,v3,dvdx3,dvdy3,dvdz3,natom3,maxatom) !  1D diatomic baths

      tmpprint(1)=v    ! interaction
      tmpprint(3)=v2   ! target internal
      tmpprint(4)=v3   ! bath internal
      tmpprint(11)=v    ! interaction

      v=v+v2+v3

      do i=1,natom
        dvdx(i)=dvdx(i)+dvdx2(i)
        dvdy(i)=dvdy(i)+dvdy2(i)
        dvdz(i)=dvdz(i)+dvdz2(i)
      enddo
      do i=natom2+1,natom
        dvdx(i)=dvdx(i)+dvdx3(i-natom2)
        dvdy(i)=dvdy(i)+dvdy3(i-natom2)
        dvdz(i)=dvdz(i)+dvdz3(i-natom2)
      enddo

c hack int PES
      if (natom3.ne.0.and.natom2.ne.0) then
      do i=1,natom2
      x2(i)=x(i)+10.
      y2(i)=y(i)+10.
      z2(i)=z(i)+10.
      enddo
      do i=natom2+1,natom
      x2(i)=x(i)
      y2(i)=y(i)
      z2(i)=z(i)
      enddo
      call lsfit(at,x2,y2,z2,v2,dvdx2,dvdy2,dvdz2,natom,maxatom) !  Interaction PES
      tmpprint(1)=tmpprint(1)-v2    ! interaction
      endif

      return

      end


! ONE DIMENSIONAL DIATOMIC BATHS
      subroutine bath(at,x,y,z,v,dvdx,dvdy,dvdz,natom,maxatom)

      implicit real*8(a-h,o-z)
      dimension x(maxatom),y(maxatom),z(maxatom)
      dimension dvdx(maxatom),dvdy(maxatom),dvdz(maxatom)
      integer at(maxatom)
      parameter(autocmi=219474.63067d0)
      parameter(autoang=0.529177249d0)
      parameter(autoev=27.2113961d0)

      v=0.
      do i=1,natom
        dvdx(i)=0.d0
        dvdy(i)=0.d0
        dvdz(i)=0.d0
      enddo

      if (natom.eq.1) return

      if (natom.gt.2) then
         print *,"Can't handle more than 2 bath atoms"
         stop
      endif

      if (natom.eq.2) then
      dx=x(1)-x(2)
      dy=y(1)-y(2)
      dz=z(1)-z(2)
      rr=dsqrt(dx*dx+dy*dy+dz*dz)

        if (at(1).eq.25.and.at(2).eq.25) then
! H2 bath
!       From Hack's fit (eq 8 in Hack, Truhlar, JCP 110, 4315 (1999))
!                        to Kolos and Wolniewicz JCP 43, 2429 (1965)
!       Rmin = 1.40121 au, Vmin = -4.74772265 eV relative to H+H
        c1=139.7160d0        ! eV
        c2=-123.8978d0       ! eV / bohr
        c3=3.4031d0          ! 1 / bohr
        c4=-6.8725d0         ! eV / bohr**2
        c5=-23.0440d0        ! eV / bohr**3
        c6=2.032d0           ! 1 / bohr

        v=(c1+c2*rr)*dexp(-c3*rr)
     &   +(c4+c5*rr)*dexp(-c6*rr)*rr**2
c       move zero from asymptote to minimum
        v=v+4.74772265
        v=v/autoev

        dvdr=((c1+c2*rr)*(-c3)+c2)*dexp(-c3*rr)
     &      +((c4+c5*rr)*(-c6)+c5)*dexp(-c6*rr)*rr**2
     &       +(c4+c5*rr)*dexp(-c6*rr)*rr*2.d0
        dvdr=dvdr/autoev

        elseif (at(1).eq.29.and.at(2).eq.29) then
! O2 bath
!       fit to MRCI+Q/CBS(AQZ,A5Z) full valence
!       Jasper April 3, 2012
        de=42046.5d0 ! exp De in cm-1 
        re=1.2075d0 ! exp in A
        c1= 2.6938139d0 ! my fit
        c2= 0.384763939d0
        c3= 0.812506485d0

        yy=rr*autoang-re
        beta = c1+c2*yy+c3*yy**2
        v = de*(1.d0-dexp(-beta*yy))**2    ! A and cm-1
        v=v/autocmi  ! convert to au

c        print *,rr,yy,beta,v

        dvdr=c1+2.d0*c2*yy+3.d0*c3*yy**2
        dvdr=dvdr*2.d0*de*(1.d0-dexp(-beta*yy))*dexp(-beta*yy)
        dvdr=dvdr*autoang/autocmi  ! convert to au

c        print *,dvdr

        elseif (at(1).eq.26.and.at(2).eq.26) then
! N2 bath
!       fit to MRCI+Q/CBS(AQZ,A5Z) full valence
!       agrees reasonably well with more complicated form of LeRoy (JCP
!       125, 164310 (2006))
!       Jasper June 9, 2010
        de=79845.d0 ! exp De in cm-1 (Ronin, Luanay, Larzillier, 
!                                     PRL 53, 159 (1984), as quoted by
!                                     LeRoy)
        re=1.097679d0 ! exp in A
        c1=2.68872341 ! my fit
        c2=0.240070803
        c3=0.472261727

        yy=rr*autoang-re
        beta = c1+c2*yy+c3*yy**2
        v = de*(1.d0-dexp(-beta*yy))**2    ! A and cm-1
        v=v/autocmi  ! convert to au

c        print *,rr,yy,beta,v

        dvdr=c1+2.d0*c2*yy+3.d0*c3*yy**2
        dvdr=dvdr*2.d0*de*(1.d0-dexp(-beta*yy))*dexp(-beta*yy)
        dvdr=dvdr*autoang/autocmi  ! convert to au

c        print *,dvdr

        elseif ((at(1).eq.27.and.at(2).eq.28).or.
     &          (at(1).eq.28.and.at(2).eq.27)) then
! CO bath
!       Morse. Fit to RKR data of PAUL H. KRUPENIE and STANLEY WEISSMAN,
!       J. CHEM. PHYS. 43, 1529 (1965)
!       with De = 11.06 eV
        de=11.06d0    ! exp De in eV
        de=de/autoev
        re=1.128322d0 ! exp in A
        re=re/autoang
        beta=1.d0/0.428d0  ! my fit in 1/A
        beta=beta*autoang

        yy=rr-re
        v = de*(1.d0-dexp(-beta*yy))**2
        dvdr=2.d0*de*(1.d0-dexp(-beta*yy))*dexp(-beta*yy)*beta

!       elseif (at(1).eq.??.and.at(2).eq.??) then
! OTHER DIATOMIC BATHS HERE
        else
        print *,"Don't know this diatomic bath"
        stop
        endif

      dvdx(1) =  dvdr*dx/rr
      dvdx(2) = -dvdr*dx/rr
      dvdy(1) =  dvdr*dy/rr
      dvdy(2) = -dvdr*dy/rr
      dvdz(1) =  dvdr*dz/rr
      dvdz(2) = -dvdr*dz/rr

      endif

      return
      end




      subroutine rgexp(at,x,y,z,v,dvdx,dvdy,dvdz,natom,maxatom)
c Rare Gas exp6 potential subroutine
c loops over geometry and looks for Rg-X interactions
c returns the full Rg-target intermolecular potential and its derivatives

      implicit real*8(a-h,o-z)
      dimension x(maxatom),y(maxatom),z(maxatom)
      dimension dvdx(maxatom),dvdy(maxatom),dvdz(maxatom)
      parameter(autoev=27.2113961d0)
      parameter(autocmi=219474.63067d0)
      parameter(autokcal=627.509d0)
      parameter(autoang=0.529177249d0)
      integer at(maxatom)
      logical troya,cutoff

      integer nfitparams,mfitparams,ix
      parameter (mfitparams=50)
      double precision ccc(mfitparams)
      common/amfit/ccc,nfitparams

      save/amfit/

      v1=0.d0
      v=0.d0
      do i=1,natom
      dvdz(i)=0.d0
      dvdy(i)=0.d0
      dvdx(i)=0.d0
      enddo

      do 1 i=1,natom
      do 2 j=i+1,natom

      m1=min(at(i),at(j))
      m2=max(at(i),at(j))
      troya=.false.   ! do or don't use Troya's form
      cutoff=.false.   ! do or don't use cutoff

      if (m1.ge.21) then ! two rare gases, skip this pair
         go to 2
      endif
      if (m2.le.20) then ! no rare gas, skip this pair
         go to 2
      endif

      if (m2.eq.m2) then  ! ANYTHING
        ix=m1
        if (ix.gt.2) ix=ix-1
        aa = ccc((ix-1)*4+1)
        bb = ccc((ix-1)*4+2)
        cc = ccc((ix-1)*4+3)
        rrc = ccc((ix-1)*4+4)
        aa=(10.d0**aa)
        cutoff=.true.
      endif

      dx=x(i)-x(j)
      dy=y(i)-y(j)
      dz=z(i)-z(j)
      rr=dsqrt(dx*dx+dy*dy+dz*dz)
      rra=rr*autoang

! NOTE CANNOT HAVE BOTH TROYA FORM AND CUTOFF FORM

      if (troya) then   ! Troya uses different form & units
        v=aa*dexp(-rra*bb)+cc/rra**6
        v=v/autokcal
        dvdr = -aa*bb*dexp(-rra*bb)-6.d0*cc/rra**7
        dvdr=dvdr/autokcal*autoang
      elseif (cutoff) then  ! cutoff 1/R**-6 at short distances
        v=aa*dexp(-rra/bb)-(cc**6/(rra**6+rrc**6))
        v=v/autocmi
        dvdr = -aa/bb*dexp(-rra/bb)
     &      +6.d0*(cc**6)*(rra**5)/(rra**6+rrc**6)**2
        dvdr=dvdr/autocmi*autoang
      else
        v=aa*dexp(-rra/bb)-(cc/rra)**6
        v=v/autocmi
        dvdr = -aa/bb*dexp(-rra/bb)+(6.d0/rra)*(cc/rra)**6
        dvdr=dvdr/autocmi*autoang
      endif
        v1=v1+v

c      print *,m1,m2,rra,v*autocmi,v1*autocmi

c derivs = sum over all bonds (DV/DRij * DRij/DXi = DV/DRij * (Xi-Xj)/Rij)
      dvdx(i) = dvdx(i) + dvdr*dx/rr
      dvdx(j) = dvdx(j) - dvdr*dx/rr
      dvdy(i) = dvdy(i) + dvdr*dy/rr
      dvdy(j) = dvdy(j) - dvdr*dy/rr
      dvdz(i) = dvdz(i) + dvdr*dz/rr
      dvdz(j) = dvdz(j) - dvdr*dz/rr

    2 continue
    1 continue
      v=v1
c      print *,'rgexp',v1
      return
      end

c **********************************************************************
c **********************************************************************
      subroutine tinkerpot(symb,xx,yy,zz,pema,dxx,dyy,dzz,nat,mnat,idum)

      use sizes
      use atoms
      use files
      use inform
      use iounit

      implicit none
      double precision autokcal,autoang
      parameter(autokcal=627.509d0)
      parameter(autoang=0.52917706d0)

      integer iprepot,mnat,nat,i,idum,itinkxyz,itinkxyzlast
      double precision xx(mnat), yy(mnat), zz(mnat),
     &                dxx(mnat),dyy(mnat),dzz(mnat),
     &                pema,energy,derivs(3,maxatm)
      character*2 symb(mnat),dum

      integer nfitparams,mfitparams
      parameter (mfitparams=50)
      double precision ccc(mfitparams)
      common/amfit/ccc,nfitparams
      common/tinkerhack/itinkxyz

      save iprepot,itinkxyzlast
      entry prepot
      if (iprepot.eq.0) then
      call prepot2
      call prepot3
c     set up the structure and mechanics calculation
      call initial
      call getxyz
      call getkey
      call mechanic
      iprepot=1
      itinkxyz=0
      itinkxyzlast=0
      return
      endif

      itinkxyz=idum
      if (itinkxyz.ne.itinkxyzlast) then
      call initial
      call getxyz
      call getkey
      call mechanic
      endif

      do i=1,n
      x(i)=xx(i)*autoang
      y(i)=yy(i)*autoang
      z(i)=zz(i)*autoang
      enddo

      call gradient (energy,derivs)
      pema=energy/autokcal

      do i=1,nat
        dxx(i)=derivs(1,i)/autokcal*autoang
        dyy(i)=derivs(2,i)/autokcal*autoang
        dzz(i)=derivs(3,i)/autokcal*autoang
      enddo

      itinkxyzlast=itinkxyz
      return

      end
c **********************************************************************
c **********************************************************************




********************************************

      subroutine prepot2

      implicit double precision(a-h,o-z)
      include 'paramls.inc'
      dimension iagroup(maxatom),
     &  ind(maxterm,maxpair),nfound(maxatom),
     &  iatom(maxperm,maxatom),
     &  idum(maxatom),nngroup(maxatom),
     &  idum2(maxperm,maxatom),
     &  idum3(maxperm,maxatom),nperm0(maxperm),nperm1(maxperm),
     &  basis(maxterm),ibasis(maxterm),r(maxpair),
     &  dbasisdr(maxterm,maxpair),
     &  rrr(maxdata,maxpair),index(maxatom,maxatom),ix(maxperm,maxpair)
      character*2 symb(maxatom),dum
      logical lreadbasis
 
      common/foox/rrr,nncoef,natom1

      save npairs,nterms,ind,ibasis

c      if (.false.) then
c      else
      open(55,file="basis.dat")
      read(55,*)natom1,npairs,nncoef,nterms
      read(55,*)
      do i=1,nterms
      read(55,*)k,ibasis(k),dum,(ind(k,j),j=1,npairs)
      enddo
      close(55)
c      endif

      return

      entry funcs1(iii,basis,ncoef,dbasisdr)

      do i=1,ncoef
      basis(i)=0.d0
      do j=1,npairs
      dbasisdr(i,j)=0.d0
      enddo
      enddo

      do j=1,npairs
      r(j)=dexp(-rrr(iii,j)*autoang)
      enddo

      do i=1,nterms
      arg=1.d0
      do j=1,npairs
      if (ind(i,j).ne.0) arg=arg*(r(j)**ind(i,j))
      enddo
      basis(ibasis(i))=basis(ibasis(i))+arg
      do j=1,npairs
      dbasisdr(ibasis(i),j)=dbasisdr(ibasis(i),j)    ! dV/dy * dy/dr
     &      -arg*dble(ind(i,j))*autoang
      enddo
      enddo

      return 
      end

***************************************************

      recursive subroutine heapp(ia,size,n,iia,ii)

      include 'paramls.inc'
      integer i,n,size,ii
      integer ia(maxatom)
      integer iia(maxperm,maxatom)
      integer iagroup(maxatom)

      if (size.eq.1) then
         ii=ii+1
         do i=1,n
           iia(ii,i)=ia(i)
         enddo
        return
      endif

      do i=1,size
        call heapp(ia,size-1,n,iia,ii)
        if (mod(size,2).eq.1) then
          tmp=ia(1)
          ia(1)=ia(size)
          ia(size)=tmp
        else
          tmp=ia(i)
          ia(i)=ia(size)
          ia(size)=tmp
      endif

      enddo

      end subroutine

***************************************************
      subroutine lsfit(at,x,y,z,v,dvdx,dvdy,dvdz,natom,maxatomx)
c
      implicit double precision (a-h,o-z)
c
      include 'paramls.inc'
      dimension coef(maxterm),sig(maxdata)
      dimension basis(maxterm)
      dimension dbasisdr(maxterm,maxpair),dvdr(maxpair)

      dimension vv(maxdata),rrr(maxdata,maxpair)
      dimension vv2(maxdata),rrr2(maxdata,maxpair)
      dimension rcom2(maxdata),xprint(50,20)
      dimension x(maxatom),y(maxatom),z(maxatom)

      dimension dvdx(maxatom),dvdy(maxatom),dvdz(maxatom)
      parameter(autoev=27.2113961d0)
      parameter(autocmi=219474.63067d0)
      parameter(autokcal=627.509d0)
      integer at(maxatom)

      logical lnumerical


      character*2 dum

      common/foox/rrr,nncoef,natom1

      save iprepot3
      entry prepot3
      if (iprepot3.eq.0) then
      open(77,file="coef.dat")
      do k=1,nncoef
      read(77,*)i,coef(k)
      enddo
      iprepot3=1
      return
      endif

c      print *
c      do j=1,natom
c      print *,at(j),x(j)*autoang,y(j)*autoang,z(j)*autoang
c      enddo

      ii=0
      if (natom1.eq.natom) then
      do j=1,natom  ! molecule
      do k=j+1,natom  ! bath
      ii=ii+1
      rrr(1,ii)=dsqrt((x(j)-x(k))**2+(y(j)-y(k))**2+(z(j)-z(k))**2)
      enddo  
      enddo  
      else
      do j=1,natom1  ! molecule
      do k=natom1+1,natom  ! bath
      ii=ii+1
      rrr(1,ii)=dsqrt((x(j)-x(k))**2+(y(j)-y(k))**2+(z(j)-z(k))**2)
      enddo  
      enddo  
      endif

      ncoef=nncoef
      call funcs1(1,basis,ncoef,dbasisdr) 
      v=0.d0
      do j=1,ncoef
         v=v+coef(j)*basis(j)
      enddo
      v=v/autocmi

      do k=1,ii
      dvdr(k)=0.d0
      do j=1,ncoef
      dvdr(k)=dvdr(k)+coef(j)*dbasisdr(j,k)
      enddo
      dvdr(k)=dvdr(k)/autocmi
      enddo

      lnumerical=.true.
      lnumerical=.false.
      resp=0.0001
      ii=0
      do j=1,natom1  ! molecule
      ijk=natom1
      if (natom1.eq.natom) ijk=j
      do k=ijk+1,natom  ! bath
      dx=x(j)-x(k)
      dy=y(j)-y(k)
      dz=z(j)-z(k)
      ii=ii+1

      if (lnumerical) then
      rrr(1,ii)=rrr(1,ii)+resp
      call funcs1(1,basis,ncoef,dbasisdr) 
      vp=0.d0
      do l=1,ncoef
         vp=vp+coef(l)*basis(l)
      enddo
      vp=vp/autocmi
      rrr(1,ii)=rrr(1,ii)-2.d0*resp
      call funcs1(1,basis,ncoef,dbasisdr) 
      vm=0.d0
      do l=1,ncoef
         vm=vm+coef(l)*basis(l)
      enddo
      vm=vm/autocmi
      rrr(1,ii)=rrr(1,ii)+resp
      dtmpdrr=(vp-vm)/(2.d0*resp)
      else
      dtmpdrr=dvdr(ii)
      endif
      dvdx(j) = dvdx(j) + dtmpdrr*dx/rrr(1,ii)
      dvdx(k) = dvdx(k) - dtmpdrr*dx/rrr(1,ii)
      dvdy(j) = dvdy(j) + dtmpdrr*dy/rrr(1,ii)
      dvdy(k) = dvdy(k) - dtmpdrr*dy/rrr(1,ii)
      dvdz(j) = dvdz(j) + dtmpdrr*dz/rrr(1,ii)
      dvdz(k) = dvdz(k) - dtmpdrr*dz/rrr(1,ii)
      enddo
      enddo

      return
 
      end

