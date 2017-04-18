! ************************************************************************
! * module legsym - direct and indirect Legendre transformation routines *
! * using symmetric and antisymmetric fourier coefficients               *
! * E. Kirk 08-Sep-2010 tested for T21 - T682 resolutions 32 & 64 bit    *
! ************************************************************************

module legsym
! ************************
! * Legendre Polynomials *
! ************************

integer :: ntru
integer :: ntp1
integer :: ncsp
integer :: nesp
integer :: nlon
integer :: nlpp
integer :: nhpp
integer :: nlat
integer :: nlev
real    :: plavor

real   , allocatable :: qi(:,:) ! P(m,n) = Associated Legendre Polynomials, used in sp2fc
real   , allocatable :: qj(:,:) ! Q(m,n) = Used for d/d(mu)                 used in sp2fcdmu
real   , allocatable :: qc(:,:) ! P(m,n) * gwd                              used in fc2sp
real   , allocatable :: qe(:,:) ! Q(mn,) * gwd / cos2                       used in mktend
real   , allocatable :: qq(:,:) ! P(m,n) * gwd / cos2 * n * (n+1) / 2  "
real   , allocatable :: qu(:,:) ! P(m,n) / (n*(n+1)) * m                    used in dv2uv
real   , allocatable :: qv(:,:) ! Q(m,n) / (n*(n+1))                        used in dv2uv
complex, allocatable :: qx(:,:) ! P(m,n) * gwd / cos2 * m                   used in mktend
end module legsym

! =================
! SUBROUTINE LEGINI
! =================

subroutine legini(klat,klpp,kesp,klev,vorpla,sid,gwd)
use legsym
implicit none

integer :: klat
integer :: klpp
integer :: kesp
integer :: klev

real          :: vorpla   ! planetary vorticity
real (kind=8) :: sid(*)   ! sin(phi)
real (kind=8) :: gwd(*)   ! Gaussian weight (phi)

integer :: jlat ! Latitude
integer :: lm
integer :: m
integer :: n

real (kind=8) :: amsq
real (kind=8) :: z1
real (kind=8) :: z2
real (kind=8) :: z3
real (kind=8) :: f1m
real (kind=8) :: f2m
real (kind=8) :: znn1
real (kind=8) :: zsin    ! sin
real (kind=8) :: zcsq    ! cos2
real (kind=8) :: zcos    ! cos
real (kind=8) :: zgwd    ! gw
real (kind=8) :: zgwdcsq ! gw / cos2
real (kind=8) :: zpli(kesp)
real (kind=8) :: zpld(kesp)

nlat = klat
nlpp = klpp
nhpp = klpp / 2
nlon = klat + klat
ntru = (nlon - 1) /3
ntp1 = ntru + 1
ncsp = ((ntru + 1) * (ntru + 2)) / 2
nesp = kesp
nlev = klev

plavor = vorpla

allocate(qi(ncsp,nhpp))
allocate(qj(ncsp,nhpp))
allocate(qc(ncsp,nhpp))
allocate(qe(ncsp,nhpp))
allocate(qx(ncsp,nhpp))
allocate(qq(ncsp,nhpp))
allocate(qu(ncsp,nhpp))
allocate(qv(ncsp,nhpp))

do jlat = 1 , nhpp

! set p(0,0) and p(0,1)

   zgwd    = gwd(jlat)            ! gaussian weight - from inigau
   zsin    = sid(jlat)            ! sin(phi) - from inigau
   zcsq    = 1.0_8 - zsin * zsin  ! cos(phi) squared
   zgwdcsq = zgwd / zcsq          ! weight / cos squared
   zcos    = sqrt(zcsq)           ! cos(phi)
   f1m     = sqrt(1.5_8)
   zpli(1) = sqrt(0.5_8)
   zpli(2) = f1m * zsin
   zpld(1) = 0.0
   lm      = 2

! loop over wavenumbers

   do m = 0 , ntru
      if (m > 0) then
         lm  = lm + 1
         f2m = -f1m * sqrt(zcsq / (m+m))
         f1m =  f2m * sqrt(m+m + 3.0_8)
         zpli(lm) = f2m
         if (lm < ncsp) then
            lm = lm + 1
            zpli(lm  ) =       f1m * zsin
            zpld(lm-1) =  -m * f2m * zsin
         endif ! (lm < ncsp)
      endif ! (m > 0)

      amsq = m * m

      do n = m+2 , ntru
         lm = lm + 1
         z1 = sqrt(((n-1)*(n-1) - amsq) / (4*(n-1)*(n-1)-1))
         z2 = zsin * zpli(lm-1) - z1 * zpli(lm-2)
         zpli(lm  ) = z2 * sqrt((4*n*n-1) / (n*n-amsq))
         zpld(lm-1) = (1-n) * z2 + n * z1 * zpli(lm-2)
      enddo ! n

      if (lm < ncsp) then ! mode (m,ntru)
         z3 = sqrt((ntru*ntru-amsq) / (4*ntru*ntru-1))
         zpld(lm)=-ntru*zsin*zpli(lm) + (ntru+ntru+1)*zpli(lm-1)*z3
      else                ! mode (ntru,ntru)
         zpld(lm)=-ntru*zsin*zpli(lm)
      endif
   enddo ! m

   lm = 0
   do m = 0 , ntru
      do n = m , ntru
           lm = lm + 1
           znn1 = 0.0
           if (n > 0) znn1 = 1.0_8 / (n*(n+1))
           qi(lm,jlat) = zpli(lm)
           qj(lm,jlat) = zpld(lm)
           qc(lm,jlat) = zpli(lm) * zgwd
           qu(lm,jlat) = zpli(lm) * znn1 * m
           qv(lm,jlat) = zpld(lm) * znn1
           qe(lm,jlat) = zpld(lm) * zgwdcsq
           qq(lm,jlat) = zpli(lm) * zgwdcsq * n * (n+1) * 0.5_8
           qx(lm,jlat) = zpli(lm) * zgwdcsq * m * (0.0,1.0)
      enddo ! n
   enddo ! m
enddo ! jlat
return
end subroutine legini


! ================
! SUBROUTINE FC2SP
! ================

subroutine fc2sp(fc,sp)
use legsym
implicit none
complex, dimension(nlon,nhpp), intent(in ) :: fc
complex, dimension(nesp/2),    intent(out) :: sp
! local
integer :: l ! Index for latitude
integer :: m ! Index for zonal wavenumber
integer :: w ! Index for spherical harmonic
integer :: e ! Index for last wavenumber

sp(:) = (0.0,0.0)

do l = 1 , nhpp
   w = 1
   do m = 1 , ntp1
      e = w + ntp1 - m
      sp(w  :e:2) = sp(w  :e:2) + qc(w  :e:2,l) * (fc(m,l) + fc(m+nlat,l)) ! XW: qc derived from module "legsym"
      sp(w+1:e:2) = sp(w+1:e:2) + qc(w+1:e:2,l) * (fc(m,l) - fc(m+nlat,l))
      w = e + 1
   enddo ! m
enddo ! l
end subroutine fc2sp


! ================
! SUBROUTINE SP2FC
! ================

subroutine sp2fc(sp,fc) ! Spectral to Fourier
use legsym
implicit none
complex, dimension(ncsp),      intent(in ) :: sp   ! Coefficients of spherical harmonics
complex, dimension(nlon,nhpp), intent(out) :: fc   ! Fourier coefficients
! local
integer :: l ! Loop index for latitude
integer :: m ! Loop index for zonal wavenumber m
integer :: w ! Index for spectral mode
integer :: e ! Index for last wavenumber
complex :: fs,fa

fc(:,:) = (0.0,0.0)

do l = 1 , nhpp
  w = 1  
  do m = 1 ,ntp1
    e = w + ntp1 - m
    fs = dot_product(qi(w  :e:2,l),sp(w  :e:2)) ! XW: qi derived from module "legsym"
    fa = dot_product(qi(w+1:e:2,l),sp(w+1:e:2))
    fc(m     ,l) = fs + fa
    fc(m+nlat,l) = fs - fa
    w = e + 1
  enddo ! m
enddo ! l
end subroutine sp2fc

! XW(2017/4/14): apply sp2fc on all 2D grid (lon,lat)
! for outsp to produce new output in GrADS format
subroutine sp2fc_all(sp,fc) ! Spectral to Fourier
use legsym
implicit none
complex, dimension(ncsp),        intent(in ) :: sp   ! Coefficients of spherical harmonics
complex, dimension(nlon,nlat/2), intent(out) :: fc   ! Fourier coefficients
! local
integer :: l ! Loop index for latitude
integer :: m ! Loop index for zonal wavenumber m
integer :: w ! Index for spectral mode
integer :: e ! Index for last wavenumber
complex :: fs,fa

fc(:,:) = (0.0,0.0)

do l = 1 , nlat/2
  w = 1  
  do m = 1 ,ntp1
    e = w + ntp1 - m
    fs = dot_product(qi(w  :e:2,l),sp(w  :e:2)) ! XW: qi derived from module "legsym"
    fa = dot_product(qi(w+1:e:2,l),sp(w+1:e:2))
    fc(m     ,l) = fs + fa
    fc(m+nlat,l) = fs - fa
    w = e + 1
  enddo ! m
enddo ! l
end subroutine sp2fc_all



! ===================
! SUBROUTINE SP2FCDMU
! ===================

subroutine sp2fcdmu(sp,fc) ! Spectral to Fourier d/dmu
use legsym
implicit none
complex, dimension(ncsp),      intent(in)  :: sp   ! Coefficients of spherical harmonics
complex, dimension(nlon,nhpp), intent(out) :: fc   ! Fourier coefficients
! local
integer :: l ! Loop index for latitude
integer :: m ! Loop index for zonal wavenumber m
integer :: w ! Index for spectral mode
integer :: e ! Index for last wavenumber
complex :: fs,fa

fc(:,:) = (0.0,0.0)

do l = 1 , nhpp
  w = 1  
  do m = 1 , ntp1
    e = w + ntp1 - m
    fs = dot_product(qj(w  :e:2,l),sp(w  :e:2))
    fa = dot_product(qj(w+1:e:2,l),sp(w+1:e:2))
    fc(m     ,l) = fa + fs
    fc(m+nlat,l) = fa - fs
    w = e + 1
  enddo ! m
enddo ! l
end subroutine sp2fcdmu


! ================
! SUBROUTINE DV2UV
! ================

! This is an alternative subroutine for computing U and V from Div and Vor.
! It looks much prettier than the regular one, but unfortunately it is slower
! if compiling with "gfortran". 
! I leave it (unused) in this module for educational purposes.

subroutine dv2uv_alt(pd,pz,pu,pv)
use legsym
implicit none
complex, parameter :: i = (0.0,1.0)
complex :: pd(nesp/2)    ! Spherical harmonics  of divergence
complex :: pz(nesp/2)    ! Spherical harmonics  of vorticity
complex :: pu(nlon,nhpp) ! Fourier coefficients of u
complex :: pv(nlon,nhpp) ! Fourier coefficients of v
complex :: zsave,uds,vds,uzs,vzs,uda,vda,uza,vza

integer :: l ! Loop index for latitude
integer :: m ! Loop index for zonal wavenumber m
integer :: w ! Loop index for spectral mode
integer :: e ! End index

pu(:,:) = (0.0,0.0)
pv(:,:) = (0.0,0.0)

zsave = pz(2)                      ! Save mode(0,1) of vorticity
pz(2) = zsave - cmplx(plavor,0.0)  ! Convert pz from absolute to relative vorticity

do l = 1 , nhpp
  w = 1  
  do m = 1 ,ntp1
    e = w + ntp1 - m
    uds = i * dot_product(qu(w  :e:2,l),pd(w:  e:2))
    vds =     dot_product(qv(w  :e:2,l),pd(w:  e:2))
    uzs = i * dot_product(qu(w  :e:2,l),pz(w:  e:2))
    vzs =     dot_product(qv(w  :e:2,l),pz(w:  e:2))
    uda = i * dot_product(qu(w+1:e:2,l),pd(w+1:e:2))
    vda =     dot_product(qv(w+1:e:2,l),pd(w+1:e:2))
    uza = i * dot_product(qu(w+1:e:2,l),pz(w+1:e:2))
    vza =     dot_product(qv(w+1:e:2,l),pz(w+1:e:2))
    pu(m     ,l) =  vzs - uds + vza - uda
    pu(m+nlat,l) = -vzs - uds + vza + uda
    pv(m     ,l) = -vds - uzs - vda - uza
    pv(m+nlat,l) =  vds - uzs - vda + uza
    w = e + 1
  enddo ! m
enddo ! l
pz(2) = zsave
return
end subroutine dv2uv_alt


! ================
! SUBROUTINE DV2UV
! ================

subroutine dv2uv(pd,pz,pu,pv)
use legsym
implicit none

real :: pd(2,nesp/2)    ! Spherical harmonics  of divergence
real :: pz(2,nesp/2)    ! Spherical harmonics  of vorticity
real :: pu(2,nlon,nhpp) ! Fourier coefficients of u
real :: pv(2,nlon,nhpp) ! Fourier coefficients of v
real :: zsave
real :: unr,uni,usr,usi,vnr,vni,vsr,vsi
real :: zdr,zdi,zzr,zzi

integer :: l ! Loop index for latitude
integer :: m ! Loop index for zonal wavenumber m
integer :: n ! Loop index for total wavenumber n
integer :: w ! Loop index for spectral mode

pu(:,:,:) = 0.0
pv(:,:,:) = 0.0

zsave = pz(1,2)           ! Save mode(0,1) of vorticity
pz(1,2) = zsave - plavor  ! Convert pz from absolute to relative vorticity

do l = 1 , nhpp
   w = 1
   do m = 1 , ntp1
      unr = 0.0 ! u - north - real
      uni = 0.0 ! u - north - imag
      usr = 0.0 ! u - south - real
      usi = 0.0 ! u - south - imag
      vnr = 0.0 ! v - north - real
      vni = 0.0 ! v - north - imag
      vsr = 0.0 ! v - south - real
      vsi = 0.0 ! v - south - imag

!     process two modes per iteration, one symmetric (m+n = even) and one anti
!     we start the loop with (n=m), so the starting mode is always symmetric

      do n = m , ntru , 2
         zdr = qu(w,l) * pd(1,w)  ! symmetric mode
         zdi = qu(w,l) * pd(2,w)
         zzr = qv(w,l) * pz(1,w)
         zzi = qv(w,l) * pz(2,w)
         unr = unr + zzr + zdi
         uni = uni + zzi - zdr
         usr = usr - zzr + zdi
         usi = usi - zzi - zdr
         zzr = qu(w,l) * pz(1,w)
         zzi = qu(w,l) * pz(2,w)
         zdr = qv(w,l) * pd(1,w)
         zdi = qv(w,l) * pd(2,w)
         vnr = vnr + zzi - zdr
         vni = vni - zzr - zdi
         vsr = vsr + zzi + zdr
         vsi = vsi - zzr + zdi
         w = w + 1
         zdr = qu(w,l) * pd(1,w)  ! antisymmetric mode
         zdi = qu(w,l) * pd(2,w)
         zzr = qv(w,l) * pz(1,w)
         zzi = qv(w,l) * pz(2,w)
         unr = unr + zzr + zdi
         uni = uni + zzi - zdr
         usr = usr + zzr - zdi
         usi = usi + zzi + zdr
         zzr = qu(w,l) * pz(1,w)
         zzi = qu(w,l) * pz(2,w)
         zdr = qv(w,l) * pd(1,w)
         zdi = qv(w,l) * pd(2,w)
         vnr = vnr + zzi - zdr
         vni = vni - zzr - zdi
         vsr = vsr - zzi - zdr
         vsi = vsi + zzr - zdi
         w = w + 1
      enddo
      if (n == ntp1) then         ! additional symmetric mode
         zdr = qu(w,l) * pd(1,w)  ! if (ntp1-m) is even
         zdi = qu(w,l) * pd(2,w)
         zzr = qv(w,l) * pz(1,w)
         zzi = qv(w,l) * pz(2,w)
         unr = unr + zzr + zdi
         uni = uni + zzi - zdr
         usr = usr - zzr + zdi
         usi = usi - zzi - zdr
         zzr = qu(w,l) * pz(1,w)
         zzi = qu(w,l) * pz(2,w)
         zdr = qv(w,l) * pd(1,w)
         zdi = qv(w,l) * pd(2,w)
         vnr = vnr + zzi - zdr
         vni = vni - zzr - zdi
         vsr = vsr + zzi + zdr
         vsi = vsi - zzr + zdi
         w = w + 1
      endif
      pu(1,m     ,l) = unr
      pu(2,m     ,l) = uni
      pu(1,m+nlat,l) = usr
      pu(2,m+nlat,l) = usi
      pv(1,m     ,l) = vnr
      pv(2,m     ,l) = vni
      pv(1,m+nlat,l) = vsr
      pv(2,m+nlat,l) = vsi
   enddo ! m
enddo ! l
pz(1,2) = zsave ! Restore pz to absolute vorticity
end subroutine dv2uv


! =================
! SUBROUTINE MKTEND
! =================

subroutine mktend(d,t,z,tn,fu,fv,ke,ut,vt)
use legsym
implicit none

complex, intent(in) :: tn(nlon,nhpp)
complex, intent(in) :: fu(nlon,nhpp)
complex, intent(in) :: fv(nlon,nhpp)
complex, intent(in) :: ke(nlon,nhpp)
complex, intent(in) :: ut(nlon,nhpp)
complex, intent(in) :: vt(nlon,nhpp)

complex, intent(out) :: d(nesp/2)
complex, intent(out) :: t(nesp/2)
complex, intent(out) :: z(nesp/2)

integer :: l ! Loop index for latitude
integer :: m ! Loop index for zonal wavenumber m
integer :: w ! Loop index for spectral mode
integer :: e ! End index for w

complex :: fus,fua,fvs,fva,kes,kea,tns,tna,uts,uta,vts,vta

d(:) = (0.0,0.0) ! divergence
t(:) = (0.0,0.0) ! temperature
z(:) = (0.0,0.0) ! vorticity

do l = 1 , nhpp  ! process pairs of Nort-South latitudes
   w = 1
   do m = 1 , ntp1
      kes = ke(m,l) + ke(m+nlat,l) ; kea = ke(m,l) - ke(m+nlat,l)
      fvs = fv(m,l) + fv(m+nlat,l) ; fva = fv(m,l) - fv(m+nlat,l)
      fus = fu(m,l) + fu(m+nlat,l) ; fua = fu(m,l) - fu(m+nlat,l)
      uts = ut(m,l) + ut(m+nlat,l) ; uta = ut(m,l) - ut(m+nlat,l)
      vts = vt(m,l) + vt(m+nlat,l) ; vta = vt(m,l) - vt(m+nlat,l)
      tns = tn(m,l) + tn(m+nlat,l) ; tna = tn(m,l) - tn(m+nlat,l)
      e = w + ntp1 - m    ! vector of symmetric modes
      d(w:e:2) = d(w:e:2) + qq(w:e:2,l) * kes - qe(w:e:2,l) * fva + qx(w:e:2,l) * fus
      t(w:e:2) = t(w:e:2) + qe(w:e:2,l) * vta + qc(w:e:2,l) * tns - qx(w:e:2,l) * uts
      z(w:e:2) = z(w:e:2) + qe(w:e:2,l) * fua + qx(w:e:2,l) * fvs
      w = w + 1           ! vector of antisymmetric modes
      d(w:e:2) = d(w:e:2) + qq(w:e:2,l) * kea - qe(w:e:2,l) * fvs + qx(w:e:2,l) * fua
      t(w:e:2) = t(w:e:2) + qe(w:e:2,l) * vts + qc(w:e:2,l) * tna - qx(w:e:2,l) * uta
      z(w:e:2) = z(w:e:2) + qe(w:e:2,l) * fus + qx(w:e:2,l) * fva
      w = e + 1
   enddo ! m
enddo ! l
end subroutine mktend


! ================
! SUBROUTINE QTEND
! ================

subroutine qtend_old(q,qn,uq,vq)
use legsym
implicit none

complex, intent(in) :: qn(nlon,nhpp)
complex, intent(in) :: uq(nlon,nhpp)
complex, intent(in) :: vq(nlon,nhpp)

complex, intent(out) :: q(nesp/2)

integer :: l ! Loop index for latitude
integer :: m ! Loop index for zonal wavenumber m
integer :: w ! Loop index for spectral mode
integer :: e ! End index for w

complex :: qns,qna,uqs,uqa,vqs,vqa

q(:) = (0.0,0.0) ! humidity

do l = 1 , nhpp  ! process pairs of Nort-Souqh latitudes
   w = 1 
   do m = 1 , ntp1
      uqs = uq(m,l) + uq(m+nlat,l) ; uqa = uq(m,l) - uq(m+nlat,l)
      vqs = vq(m,l) + vq(m+nlat,l) ; vqa = vq(m,l) - vq(m+nlat,l)
      qns = qn(m,l) + qn(m+nlat,l) ; qna = qn(m,l) - qn(m+nlat,l)
      e = w + ntp1 - m    ! vector of symmetric modes
      q(w:e:2) = q(w:e:2) + qe(w:e:2,l) * vqa + qc(w:e:2,l) * qns - qx(w:e:2,l) * uqs 
      w = w + 1           ! vector of antisymmetric modes
      q(w:e:2) = q(w:e:2) + qe(w:e:2,l) * vqs + qc(w:e:2,l) * qna - qx(w:e:2,l) * uqa 
      w = e + 1 
   enddo ! m 
enddo ! l 
return
end subroutine qtend_old


! ================
! SUBROUTINE UV2DV
! ================

subroutine uv2dv(pu,pv,pd,pz)
use legsym
implicit none

complex, intent(in) :: pu(nlon,nhpp)
complex, intent(in) :: pv(nlon,nhpp)

complex, intent(out) :: pd(nesp/2)
complex, intent(out) :: pz(nesp/2)

integer :: l ! Loop index for latitude
integer :: m ! Loop index for zonal wavenumber m
integer :: w ! Loop index for spectral mode
integer :: e ! End index for w

complex :: zus,zua,zvs,zva

pd(:) = (0.0,0.0) ! divergence
pz(:) = (0.0,0.0) ! vorticity

do l = 1 , nhpp   ! process pairs of Nort-Souqh latitudes
   w = 1 
   do m = 1 , NTP1
      zus = pu(m,l) + pu(m+nlat,l) ; zua = pu(m,l) - pu(m+nlat,l) 
      zvs = pv(m,l) + pv(m+nlat,l) ; zva = pv(m,l) - pv(m+nlat,l) 
      e = w + ntp1 - m    ! vector of symmetric modes
      pz(w:e:2) = pz(w:e:2) + qe(w:e:2,l) * zua + qx(w:e:2,l) * zvs 
      pd(w:e:2) = pd(w:e:2) - qe(w:e:2,l) * zva + qx(w:e:2,l) * zus 
      w = w + 1           ! vector of antisymmetric modes
      pz(w:e:2) = pz(w:e:2) + qe(w:e:2,l) * zus + qx(w:e:2,l) * zva 
      pd(w:e:2) = pd(w:e:2) - qe(w:e:2,l) * zvs + qx(w:e:2,l) * zua 
      w = e + 1 
   enddo ! m 
enddo ! l 
end subroutine uv2dv


! ==================
! SUBROUTINE REG2ALT
! ==================

subroutine reg2alt(pr,klev)
use legsym
implicit none

integer, intent(in) :: klev
real, dimension(nlon,nlat,klev), intent(inout) :: pr(nlon,nlat,klev)
! local
real, dimension(nlon,nlat,klev) :: pa
integer :: jlat

do jlat = 1 , nlat / 2
  pa(:,2*jlat-1,:) = pr(:,jlat       ,:)
  pa(:,2*jlat  ,:) = pr(:,nlat-jlat+1,:)
enddo
pr = pa
end subroutine reg2alt


! ==================
! SUBROUTINE ALT2REG
! ==================

subroutine alt2reg(pa,klev)
use legsym
implicit none

integer, intent(in) :: klev
real, dimension(nlon,nlat,klev), intent(inout) :: pa
! local
real, dimension(nlon,nlat,klev) :: pr
integer :: jlat

do jlat = 1 , nlat / 2
  pr(:,jlat       ,:) = pa(:,2*jlat-1,:)
  pr(:,nlat-jlat+1,:) = pa(:,2*jlat  ,:)
enddo
pa = pr
end subroutine alt2reg


! ================
! SUBROUTINE ALTCS
! XW: Alternative Cross-Section
! ================

subroutine altcs(pcs)
use legsym
implicit none
real, dimension(nlat,nlev), intent(inout) :: pcs
real, dimension(nlat,nlev) :: pal
integer :: jlat

do jlat = 1 , nlat / 2
  pal(jlat       ,:) = pcs(2*jlat-1,:)
  pal(nlat-jlat+1,:) = pcs(2*jlat  ,:)
enddo
pcs = pal
end subroutine altcs


! =================
! SUBROUTINE ALTLAT 
! Alternative Latitude
! =================

subroutine altlat(pr,klat)
implicit none
integer, intent(in) :: klat
real, dimension(klat), intent(inout) :: pr  ! regular     grid
! local
real, dimension(klat) :: pa  ! alternating grid
integer :: jlat

do jlat = 1 , klat / 2
  pa(2*jlat-1) = pr(jlat       )
  pa(2*jlat  ) = pr(klat-jlat+1)
enddo
pr = pa
end subroutine altlat



! =======================================================
! SUBROUTINE INIGAU
!
! cleaned by Xinyu Wen on Apr/11/2017
! pure functions "ql" and "qld" are incorporated into the
! main pure function of "inigau", using "contains" block.
! And, gaussmod.f90 was pasted at the end of legsym.f90
!
! pure subroutine "inigau": 计算高斯纬度及权重的初始化
! =======================================================

pure &
subroutine inigau(klat,pz0,pzw)                             ! pz0 & pzw are (kind=8) reals !!!
   implicit none
   integer,                       intent(in)  :: klat       ! Number of Gaussian latitudes
   real(kind=8), dimension(klat), intent(out) :: pz0(klat)  ! Gaussian abscissas
   real(kind=8), dimension(klat), intent(out) :: pzw(klat)  ! Gaussian weights

   ! local
   integer      , parameter :: NITER = 50                   ! Maximum # of iterations
   real (kind=8), parameter :: PI    =  3.14159265358979_8  ! Pi
   real (kind=8), parameter :: ZEPS  =  1.0e-16             ! Convergence criterion
   integer                  :: jlat                         ! Latitudinal loop index
   integer                  :: jiter                        ! Iteration loop index
   real (kind=8)            :: z0,z1,z2,z3,z4,z5
   !real (kind=8)            :: ql,qld ! XW(2017/4/11): use "contains" to incorporate ql and qld inside

   ! Compute Gaussian abscissas & weights

   z0 = PI / (2*klat+1)
   z1 = 1.0_8 / (klat*klat*8)
   z4 = 2.0_8 / (klat*klat)

   do jlat = 1 , klat/2
      z2 = z0 * (2*jlat - 0.5_8)
      z2 = cos(z2 + z1 / tan(z2))
      do jiter = 1 , NITER
         z3 = ql(klat,z2) * qld(klat,z2)
         z2 = z2 - z3
         if (abs(z3) < ZEPS) exit ! converged
      enddo ! jiter
      z5 = ql(klat-1,z2) / sqrt(klat - 0.5_8)
      pz0(jlat) = z2
      pzw(jlat) = z4 * (1.0_8 - z2 * z2) / (z5 * z5)
      pz0(klat-jlat+1) = -z2
      pzw(klat-jlat+1) = pzw(jlat)
   enddo ! jlat

contains

! ===========
! FUNCTION QL
! ===========

   pure &
   function ql(k,p)
      implicit none
      integer      , intent(IN) :: k
      real (kind=8), intent(IN) :: p
      real (kind=8)             :: ql

      ! local
      real (kind=8) :: z0,z1,z2,z3,z4
      integer :: j

      z0 = acos(p)
      z1 = 1.0
      z2 = 0.0

      !XW(2017-4-8): avoid possibility of uninitialized use of z3 here,
      !I add "if (k.ge.0)" structure to avoid uninitializaion of z3"
      !if anything wrong with this, please remove the outter "if" structure.
      if (k.ge.0) then
         do j = k , 0 , -2
            z3 = z1 * cos(z0 * j)
            z2 = z2 + z3
            z4 = (k-j+1) * (k+j) * 0.5_8
            z1 = z1 * z4 / (z4 + (j-1))
         enddo ! j
         if (mod(k,2) == 0) z2 = z2 - 0.5_8 * z3
      end if

      z0 = sqrt(2.0_8)
      do j = 1 ,k
         z0 = z0 * sqrt(1.0_8 - 0.25_8 / (j*j))
      enddo ! j
      ql = z0 * z2
   end function ql

   ! ============
   ! FUNCTION QLD
   ! ============
   pure &
   function qld(k,p)
      implicit none
      integer      , intent(IN) :: k
      real (kind=8), intent(IN) :: p
      real (kind=8)             :: qld

      ! local
      real (kind=8) :: z
      !real (kind=8) :: ql

      z = p * ql(k,p) - sqrt((k + k + 1.0_8) / (k + k - 1.0_8)) * ql(k-1,p)
      qld = (p * p - 1.0_8) / (k * z)
   end function qld

end subroutine inigau
