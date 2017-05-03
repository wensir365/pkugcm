module dycore

!=========================================
!public ::  球谐变换
!           sp2fc_nlev
!           sp2fc_1lev
!           fc2sp_1lev
!           dv2uv_nlev
!           sp2fcdmu_1lev
!
!           傅立叶变换
!           gp2fc_acc
!              dfft2_acc, dfft3_acc, dfft4_acc, dfft8_acc
!           fc2gp_acc
!              ifft2_acc, ifft3_acc, ifft4_acc, ifft8_acc
!
!           非线性项的倾向计算
!           mktend_acc
!=========================================            

CONTAINS


! sp2fc_nlev与上面的sp2fc的唯一区别是增加垂向循环
! XW: for better performance in OpenACC
pure &
subroutine sp2fc_nlev(sp,fc,   qi, nlon,nlat,nlev,nhpp,ntp1,ncsp) ! Spectral to Fourier
implicit none
complex, dimension(ncsp,nlev),      intent(in ) :: sp   ! Coefficients of spherical harmonics
complex, dimension(nlon,nhpp,nlev), intent(out) :: fc   ! Fourier coefficients

real,    dimension(ncsp,nhpp),      intent(in ) :: qi
integer,                            intent(in ) :: nlon,nlat,nlev,nhpp,ntp1,ncsp

! local
integer :: l ! Loop index for latitude
integer :: m ! Loop index for zonal wavenumber m
integer :: w ! Index for spectral mode
integer :: e ! Index for last wavenumber
integer :: jlev
complex :: fs,fa

!$---acc kernels

fc = (0.0,0.0)

do jlev = 1, NLEV
   do l = 1 , nhpp
      w = 1  
      do m = 1 ,ntp1
         e = w + ntp1 - m
         fs = dot_product(qi(w  :e:2,l),sp(w  :e:2,jlev)) ! 对称模态
         fa = dot_product(qi(w+1:e:2,l),sp(w+1:e:2,jlev)) ! 非对称模态
         fc(m     ,l,jlev) = fs + fa
         fc(m+nlat,l,jlev) = fs - fa
         w = e + 1
      enddo ! m
   enddo ! l
end do

!$---acc end kernels

end subroutine sp2fc_nlev



pure &
subroutine sp2fc_1lev(sp,fc,   qi,nlon,nlat,nhpp,ntp1,ncsp) ! Spectral to Fourier
implicit none
complex, dimension(ncsp),      intent(in ) :: sp   ! Coefficients of spherical harmonics
complex, dimension(nlon,nhpp), intent(out) :: fc   ! Fourier coefficients

real,    dimension(ncsp,nhpp), intent(in ) :: qi
integer,                       intent(in ) :: nlon,nlat,nhpp,ntp1,ncsp

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
end subroutine sp2fc_1lev




pure &
subroutine fc2sp_1lev(fc,sp,   qc,nlon,nlat,nhpp,ntp1,ncsp)
implicit none
complex, dimension(nlon,nhpp), intent(in ) :: fc
complex, dimension(ncsp),      intent(out) :: sp

real,    dimension(ncsp,nhpp), intent(in ) :: qc
integer,                       intent(in ) :: nlon,nlat,nhpp,ntp1,ncsp

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
end subroutine fc2sp_1lev




! dv2uv_nlev与上面的dv2uv的唯一区别是增加垂向循环
! XW: for better performance in OpenACC
pure &
subroutine dv2uv_nlev(pd,pz,pu,pv,   qu,qv, nlon,nlat,nlev,nhpp,ntru,ntp1,ncsp,nesp, plavor)
implicit none

real, intent(in ) :: pd(2,nesp/2,nlev)    ! Spherical harmonics  of divergence
real, intent(in ) :: pz(2,nesp/2,nlev)    ! Spherical harmonics  of vorticity
real, intent(out) :: pu(2,nlon,nhpp,nlev) ! Fourier coefficients of u
real, intent(out) :: pv(2,nlon,nhpp,nlev) ! Fourier coefficients of v

real, intent(in)     :: qu(ncsp,nhpp),qv(ncsp,nhpp)
integer, intent(in)  :: nlon,nlat,nlev,nhpp,ntru,ntp1,ncsp,nesp
real, intent(in)     :: plavor

real :: pzz(2,nesp/2,nlev)                ! tmp pz
real :: zsave
real :: unr,uni,usr,usi,vnr,vni,vsr,vsi
real :: zdr,zdi,zzr,zzi

integer :: l      ! Loop index for latitude
integer :: m      ! Loop index for zonal wavenumber m
integer :: n      ! Loop index for total wavenumber n
integer :: w      ! Loop index for spectral mode
integer :: jlev   ! main loop

!$---acc kernels

pu = 0.0
pv = 0.0

pzz = pz
pzz(1,2,:) = pzz(1,2,:)-plavor

do jlev = 1, NLEV
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
         zdr = qu(w,l) * pd(1,w,jlev)  ! symmetric mode
         zdi = qu(w,l) * pd(2,w,jlev)
         zzr = qv(w,l) * pzz(1,w,jlev)
         zzi = qv(w,l) * pzz(2,w,jlev)
         unr = unr + zzr + zdi
         uni = uni + zzi - zdr
         usr = usr - zzr + zdi
         usi = usi - zzi - zdr
         zzr = qu(w,l) * pzz(1,w,jlev)
         zzi = qu(w,l) * pzz(2,w,jlev)
         zdr = qv(w,l) * pd(1,w,jlev)
         zdi = qv(w,l) * pd(2,w,jlev)
         vnr = vnr + zzi - zdr
         vni = vni - zzr - zdi
         vsr = vsr + zzi + zdr
         vsi = vsi - zzr + zdi
         w = w + 1
         zdr = qu(w,l) * pd(1,w,jlev)  ! antisymmetric mode
         zdi = qu(w,l) * pd(2,w,jlev)
         zzr = qv(w,l) * pzz(1,w,jlev)
         zzi = qv(w,l) * pzz(2,w,jlev)
         unr = unr + zzr + zdi
         uni = uni + zzi - zdr
         usr = usr + zzr - zdi
         usi = usi + zzi + zdr
         zzr = qu(w,l) * pzz(1,w,jlev)
         zzi = qu(w,l) * pzz(2,w,jlev)
         zdr = qv(w,l) * pd(1,w,jlev)
         zdi = qv(w,l) * pd(2,w,jlev)
         vnr = vnr + zzi - zdr
         vni = vni - zzr - zdi
         vsr = vsr - zzi - zdr
         vsi = vsi + zzr - zdi
         w = w + 1
      enddo
      if (n == ntp1) then         ! additional symmetric mode
         zdr = qu(w,l) * pd(1,w,jlev)  ! if (ntp1-m) is even
         zdi = qu(w,l) * pd(2,w,jlev)
         zzr = qv(w,l) * pzz(1,w,jlev)
         zzi = qv(w,l) * pzz(2,w,jlev)
         unr = unr + zzr + zdi
         uni = uni + zzi - zdr
         usr = usr - zzr + zdi
         usi = usi - zzi - zdr
         zzr = qu(w,l) * pzz(1,w,jlev)
         zzi = qu(w,l) * pzz(2,w,jlev)
         zdr = qv(w,l) * pd(1,w,jlev)
         zdi = qv(w,l) * pd(2,w,jlev)
         vnr = vnr + zzi - zdr
         vni = vni - zzr - zdi
         vsr = vsr + zzi + zdr
         vsi = vsi - zzr + zdi
         w = w + 1
      endif
      pu(1,m     ,l,jlev) = unr
      pu(2,m     ,l,jlev) = uni
      pu(1,m+nlat,l,jlev) = usr
      pu(2,m+nlat,l,jlev) = usi
      pv(1,m     ,l,jlev) = vnr
      pv(2,m     ,l,jlev) = vni
      pv(1,m+nlat,l,jlev) = vsr
      pv(2,m+nlat,l,jlev) = vsi
   enddo ! m
enddo ! l
end do

!$---acc end kernels

end subroutine dv2uv_nlev



pure &
subroutine sp2fcdmu_1lev(sp,fc,   qj,nlon,nlat,nhpp,ntp1,ncsp) ! Spectral to Fourier d/dmu
implicit none
complex, dimension(ncsp),      intent(in)  :: sp   ! Coefficients of spherical harmonics
complex, dimension(nlon,nhpp), intent(out) :: fc   ! Fourier coefficients

real,    dimension(ncsp,nhpp), intent(in ) :: qj
integer,                       intent(in ) :: nlon,nlat,nhpp,ntp1,ncsp

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
    fs = dot_product(qj(w  :e:2,l),sp(w  :e:2))    ! qj from legsym
    fa = dot_product(qj(w+1:e:2,l),sp(w+1:e:2))
    fc(m     ,l) = fa + fs
    fc(m+nlat,l) = fa - fs
    w = e + 1
  enddo ! m
enddo ! l
end subroutine sp2fcdmu_1lev




!=======================================================================
!           这两个关键过程调用安全 无任何依赖
! public:   gp2fc( x(N,lot), N, lot, trigs),  lot=NLAT*NLEV
!           fc2gp( x(N,lot), N, lot, trigs)
!           这两个subroutine都只需要use pumamod中的FFT基函数trigs(NLON)
!           而trigs已经在prolog中使用fftini初始化过了
!           在gridpoint中 频繁地调用了这两个子程序
! private:  fftini   基函数初始化 设置pumamod中的trigs 依赖use pumamod
!      pure dfft2, dfft3, dfft4, dfft8
!      pure ifft2, ifft3, ifft4, ifft8
! Comments by X. Wen, Peking Univ, Apr-22-2017
!=======================================================================


!     ================
!     SUBROUTINE GP2FC
!     ================

! XW (2017/4/8): can parallel along with lot!!! 
! lot:多少根纬线? (NLPP*NLEV); n:每根纬线上有多少点? (NLON)

pure  subroutine gp2fc_acc(a,n,lot,base)
      implicit none
      integer, intent(in) :: n, lot
      real, dimension(n,lot), intent(inout) :: a
      real, dimension(n    ), intent(in   ) :: base
      integer :: la,l

      interface
      pure  subroutine dfft8_acc(a,c,n,lot)
               implicit none
               integer, intent(in) :: n, lot
               real, dimension(n,lot), intent(in)  :: a
               real, dimension(n,lot), intent(out) :: c
            end subroutine dfft8_acc
      pure  subroutine dfft4_acc(a,trigs,n,lot,la)
               implicit none
               integer, intent(in) :: n, lot
               integer, intent(inout) :: la
               real, dimension(n), intent(in) :: trigs
               real, dimension(n,lot), intent(inout) :: a
            end subroutine dfft4_acc
      pure  subroutine dfft3_acc(a,trigs,n)
               implicit none
               integer, intent(in) :: n
               real, dimension(n), intent(in) :: trigs
               real, dimension(n), intent(inout) :: a
            end subroutine dfft3_acc
      pure  subroutine dfft2_acc(a,trigs,n)
               implicit none
               integer, intent(in) :: n
               real, dimension(n), intent(in) :: trigs
               real, dimension(n), intent(inout) :: a
            end subroutine dfft2_acc
      end interface

      call dfft8_acc(a,a,n,lot)
      la = n / 8
      do while (la >= 4)
         call dfft4_acc(a,base,n,lot,la)
      enddo

      if (la == 3) then
         do l = 1 , lot
            call dfft3_acc(a(1,l),base,n)
         enddo
      endif

      if (la == 2) then
         do l = 1 , lot
            call dfft2_acc(a(1,l),base,n)
         enddo
      endif

      end subroutine gp2fc_acc



!     ================
!     SUBROUTINE FC2GP
!     ================

pure  subroutine fc2gp_acc(a,n,lot,base)
      implicit none
      integer, intent(in) :: n, lot
      real, dimension(n,lot), intent(inout) :: a
      real, dimension(n    ), intent(in   ) :: base
      integer :: nf,la

      interface
      pure  subroutine ifft8_acc(a,c,n,lot)
               implicit none
               integer, intent(in) :: n, lot
               real, dimension(n,lot), intent(in)  :: a
               real, dimension(n,lot), intent(out) :: c
            end subroutine ifft8_acc
      pure  subroutine ifft4_acc(c,trigs,n,lot,la)
               implicit none
               integer, intent(in) :: n, lot
               integer, intent(inout) :: la
               real, dimension(n), intent(in) :: trigs
               real, dimension(n,lot), intent(inout) :: c
            end subroutine ifft4_acc
      pure  subroutine ifft3_acc(a,trigs,n,lot,la)
               implicit none
               integer, intent(in) :: n, lot
               integer, intent(inout) :: la
               real, dimension(n), intent(in) :: trigs
               real, dimension(n,lot), intent(inout) :: a
            end subroutine ifft3_acc
      pure  subroutine ifft2_acc(a,trigs,n,lot,la)
               implicit none
               integer, intent(in) :: n, lot
               integer, intent(inout) :: la
               real, dimension(n), intent(in) :: trigs
               real, dimension(n,lot), intent(inout) :: a
            end subroutine ifft2_acc
      end interface

      nf = n/8
      do while (nf >= 4)
         nf = nf/4
      enddo
      la = 1
      if (nf == 2) call ifft2_acc(a,base,n,lot,la)
      if (nf == 3) call ifft3_acc(a,base,n,lot,la)
      do while (la < n/8)
         call ifft4_acc(a,base,n,lot,la)
      enddo
      call ifft8_acc(a,a,n,lot)
      end subroutine fc2gp_acc



!     ================
!     SUBROUTINE DFFT2
!     ================

pure  subroutine dfft2_acc(a,trigs,n)
!$acc routine worker
      implicit none
      integer, intent(in) :: n
      real, dimension(n), intent(in) :: trigs
      real, dimension(n), intent(inout) :: a
      
      ! local
      real, dimension(n) :: c
      integer :: ja,jb
      real :: c1,s1,a1p3,a3m1
      integer :: i

      c(1) = a(1) + a(2)
      c(2) = 0.0

      ja = 3
      jb = n - 1

      do i=3,n-5,4
         c1 = trigs(ja  )
         s1 = trigs(ja+1)
         a1p3 = c1 * a(i+1) + s1 * a(i+3)
         a3m1 = c1 * a(i+3) - s1 * a(i+1)
         c(ja  ) = a(i) + a1p3
         c(jb  ) = a(i) - a1p3
         c(ja+1) = a3m1 + a(i+2)
         c(jb+1) = a3m1 - a(i+2)
         ja = ja + 2
         jb = jb - 2
      enddo

      c(ja  ) =  a(n-1)
      c(ja+1) = -a(n  )

      a = c
      end subroutine dfft2_acc


!     ================
!     SUBROUTINE DFFT3
!     ================

pure  subroutine dfft3_acc(a,trigs,n)
!$acc routine worker
      implicit none
      integer, intent(in) :: n
      real, dimension(n), intent(in) :: trigs
      real, dimension(n), intent(inout) :: a

      ! local
      real(kind=8), parameter :: SIN60 = 0.866025403784438D0
      real, dimension(n) :: c(n)
      integer :: ja,jb,jc
      real :: c1,s1, c2,s2, a1,b1, a2,b2, a3,b3
      integer :: i

      ja = 1              !  1
      jb = 2 * (n/3)  + 1 ! 65
      jc = jb             ! 65

      c(ja  ) = a(1) + a(2) + a(3)
      c(ja+1) = 0.0
      c(jb  ) = a(1) - 0.5 * (a(2) + a(3))
      c(jb+1) =      SIN60 * (a(3) - a(2))

      ja = 3         !  3, 5, 7, ... ,31
      jb = jb + 2    ! 67,69,71, ... ,95
      jc = jc - 2    ! 63,61,59, ... ,35

      do i = 4 , n-8 , 6 ! 88
         c1 = trigs(ja  )
         s1 = trigs(ja+1)
         c2 = trigs(ja+ja-1)
         s2 = trigs(ja+ja  )
         a1 = (c1*a(i+1)+s1*a(i+4))+(c2*a(i+2)+s2*a(i+5))
         b1 = (c1*a(i+4)-s1*a(i+1))+(c2*a(i+5)-s2*a(i+2))
         a2 = a(i  ) - 0.5 * a1
         b2 = a(i+3) - 0.5 * b1
         a3 = SIN60*((c1*a(i+1)+s1*a(i+4))-(c2*a(i+2)+s2*a(i+5)))
         b3 = SIN60*((c1*a(i+4)-s1*a(i+1))-(c2*a(i+5)-s2*a(i+2)))
         c(ja  ) = a(i  ) + a1
         c(ja+1) = a(i+3) + b1
         c(jb  ) = a2 + b3
         c(jb+1) = b2 - a3
         c(jc  ) = a2 - b3
         c(jc+1) =-b2 - a3
         ja = ja + 2
         jb = jb + 2
         jc = jc - 2
      enddo

      if (ja <= jc) then ! ja=33  jc=33
         c(ja  ) = a(n-2) + 0.5 * (a(n-1) - a(n)) ! 33
         c(ja+1) =       -SIN60 * (a(n-1) + a(n)) ! 34
      endif

      a = c
      end subroutine dfft3_acc


!     ================
!     SUBROUTINE DFFT4
!     ================

pure  subroutine dfft4_acc(a,trigs,n,lot,la)
      implicit none
      integer, intent(in) :: n, lot
      integer, intent(inout) :: la
      real, dimension(n), intent(in) :: trigs
      real, dimension(n,lot), intent(inout) :: a

      ! local
      real,parameter :: SIN45=sqrt(0.5)
      real, dimension(n,lot) :: c
      integer :: i1,i2,i3,i4,i5,i6,i7, ibase
      integer :: j1,j2,j3,j4,j5,j6,j7, j0,jink
      integer :: kb,kc,kd
      real :: c1,s1, c2,s2, c3,s3
      real :: a0p2, a1p3, a1m3, a1p5, a2p6, a3p7, a5m1, a6m2, a7m3
      real :: a0,a1,a2,a3, b0,b1,b2,b3
      integer :: i,j,k,l

      la = la / 4

      i1 = la
      i2 = la + i1
      i3 = la + i2
      i4 = la + i3
      i5 = la + i4
      i6 = la + i5
      i7 = la + i6

      j1 = n/2 - la
      j2 = n - la
      j3 = j1
      j5 = j1 + la

      do i=1,la
         do l=1,lot
         a0p2 = a(i   ,l) + a(i2+i,l)
         a1p3 = a(i1+i,l) + a(i3+i,l)
         c(   i,l) = a0p2 + a1p3
         c(j2+i,l) = a0p2 - a1p3
         c(j1+i,l) = a(   i,l) - a(i2+i,l)
         c(j5+i,l) = a(i3+i,l) - a(i1+i,l)
         enddo
      enddo

      jink = 2 * la
      j0 = la
      j1 = j1 + jink
      j2 = j2 - jink
      j3 = j3 - jink
      j4 = j0 + la
      j5 = j1 + la
      j6 = j2 + la
      j7 = j3 + la

      ibase=4*la

      do k=la,(n-4)/8,la
         kb=k+k
         kc=kb+kb
         kd=kc+kb
         c1=trigs(kb+1)
         s1=trigs(kb+2)
         c2=trigs(kc+1)
         s2=trigs(kc+2)
         c3=trigs(kd+1)
         s3=trigs(kd+2)

         i=ibase+1
         do j=1,la
            do l=1,lot
            a1p5 = c1 * a(i1+i,l) + s1 * a(i5+i,l)
            a2p6 = c2 * a(i2+i,l) + s2 * a(i6+i,l)
            a3p7 = c3 * a(i3+i,l) + s3 * a(i7+i,l)
            a5m1 = c1 * a(i5+i,l) - s1 * a(i1+i,l)
            a6m2 = c2 * a(i6+i,l) - s2 * a(i2+i,l)
            a7m3 = c3 * a(i7+i,l) - s3 * a(i3+i,l)
            a0 = a(i,l) + a2p6
            a2 = a(i,l) - a2p6
            a1 = a1p5 + a3p7
            a3 = a3p7 - a1p5
            b0 = a(i4+i,l) + a6m2
            b2 = a(i4+i,l) - a6m2
            b1 = a5m1 + a7m3
            b3 = a5m1 - a7m3
            c(j0+j,l) = a0+a1
            c(j2+j,l) = a0-a1
            c(j4+j,l) = b0+b1
            c(j6+j,l) = b1-b0
            c(j1+j,l) = a2+b3
            c(j3+j,l) = a2-b3
            c(j5+j,l) = a3+b2
            c(j7+j,l) = a3-b2
            enddo
            i=i+1
         enddo

         ibase=ibase+8*la
         j0 = j0 + jink
         j1 = j1 + jink
         j2 = j2 - jink
         j3 = j3 - jink
         j4 = j0 + la
         j5 = j1 + la
         j6 = j2 + la
         j7 = j3 + la
      end do

      if (j1 <= j2) then
         i=ibase+1
         do j=1,la
            do l=1,lot
            a1p3 = sin45 * (a(i1+i,l) + a(i3+i,l))
            a1m3 = sin45 * (a(i1+i,l) - a(i3+i,l))
            c(j0+j,l) =  a(   i,l) + a1m3
            c(j1+j,l) =  a(   i,l) - a1m3
            c(j4+j,l) = -a(i2+i,l) - a1p3
            c(j5+j,l) =  a(i2+i,l) - a1p3
            enddo
            i=i+1
         enddo
      endif

      if (la == 1) then
         do l=1,lot
            a(1,l) = c(1,l)
            a(2,l) = 0.0
            a(3:n,l) = c(2:n-1,l)
         enddo
      else
         a = c
      endif
      end subroutine dfft4_acc


!     ================
!     SUBROUTINE DFFT8
!     ================

pure  subroutine dfft8_acc(a,c,n,lot)
!$acc routine worker
      implicit none
      integer, intent(in) :: n, lot
      real, dimension(n,lot), intent(in)  :: a
      real, dimension(n,lot), intent(out) :: c

      ! local
      integer :: la, i0,i1,i2,i3,i4,i5,i6,i7
      real :: z, zsin45
      real, dimension(lot) :: a0p4,a1p5,a2p6,a3p7, a5m1,a7m3,a0m4,a6m2
      real, dimension(lot) :: a0p4p2p6, a1p5p3p7, a7m3p5m1, a7m3m5m1
      integer :: i
      character(len=200) :: tmp
      
      la = n / 8
      z  = 1.0 / n
      zsin45 = z * sqrt(0.5)

      do i=1, la
         i0 = i
         i1 = i0 + la
         i2 = i1 + la
         i3 = i2 + la
         i4 = i3 + la
         i5 = i4 + la
         i6 = i5 + la
         i7 = i6 + la

         a0p4 =  a(i0,:) + a(i4,:)
         a1p5 =  a(i1,:) + a(i5,:)
         a2p6 =  a(i2,:) + a(i6,:)
         a3p7 =  a(i3,:) + a(i7,:)
         a5m1 =  a(i5,:) - a(i1,:)
         a7m3 =  a(i7,:) - a(i3,:)
         a0m4 = (a(i0,:) - a(i4,:)) * z
         a6m2 = (a(i6,:) - a(i2,:)) * z

         a0p4p2p6 = a0p4 + a2p6
         a1p5p3p7 = a1p5 + a3p7
         a7m3p5m1 = (a7m3 + a5m1) * zsin45
         a7m3m5m1 = (a7m3 - a5m1) * zsin45

         c(i0,:) = z * (a0p4p2p6 + a1p5p3p7)
         c(i7,:) = z * (a0p4p2p6 - a1p5p3p7)
         c(i3,:) = z * (a0p4 - a2p6)
         c(i4,:) = z * (a3p7 - a1p5)
         c(i1,:) = a0m4 + a7m3m5m1
         c(i5,:) = a0m4 - a7m3m5m1
         c(i2,:) = a7m3p5m1 + a6m2
         c(i6,:) = a7m3p5m1 - a6m2
         
      enddo
      end subroutine dfft8_acc


!     ================
!     SUBROUTINE IFFT2
!     ================

pure  subroutine ifft2_acc(a,trigs,n,lot,la)
      implicit none
      integer, intent(in) :: n, lot
      integer, intent(inout) :: la
      real, dimension(n), intent(in) :: trigs
      real, dimension(n,lot), intent(inout) :: a

      ! local
      real, dimension(n,lot) :: c
      real :: c1,s1, amb,apb
      integer :: ia,ib
      integer :: j,l

      c(1,:) = 0.5 * a(1,:)
      c(2,:) = c(1,:)

      ia    =   3
      ib    = n-1

      do j = 3 , n-5 , 4
         c1 = trigs(ia  )
         s1 = trigs(ia+1)
         do l=1,lot
            amb = a(ia  ,l) - a(ib  ,l)
            apb = a(ia+1,l) + a(ib+1,l)
            c(j  ,l) = a(ia  ,l) + a(ib  ,l)
            c(j+2,l) = a(ia+1,l) - a(ib+1,l)
            c(j+1,l) = c1 * amb - s1 * apb
            c(j+3,l) = s1 * amb + c1 * apb
         enddo
         ia = ia + 2
         ib = ib - 2
      enddo
      c(n-1,:) =  a(ia  ,:)
      c(n  ,:) = -a(ia+1,:)

      la = 2
      a  = c
      end subroutine ifft2_acc


!     ================
!     SUBROUTINE IFFT3
!     ================

pure  subroutine ifft3_acc(a,trigs,n,lot,la)
      implicit none
      integer, intent(in) :: n, lot
      integer, intent(inout) :: la
      real, dimension(n), intent(in) :: trigs
      real, dimension(n,lot), intent(inout) :: a

      ! local
      real, parameter :: SIN60 = 0.866025403784438D0
      real, dimension(n,lot) :: c
      integer :: ia,ib,ic
      real :: c1,s1, c2,s2
      real :: hbpc,hbmc, sbpc,sbmc
      integer :: j,l

      ib = 2 * (n/3) + 1

      c(1,:) = 0.5 * a(1,:) + a(ib,:)
      c(2,:) = 0.5 * a(1,:) - 0.5 * a(ib,:) - SIN60 * a(ib+1,:)
      c(3,:) = 0.5 * a(1,:) - 0.5 * a(ib,:) + SIN60 * a(ib+1,:)

      ia = 3
      ic = ib - 2
      ib = ib + 2

      do j = 4 , n-8 , 6
         c1 = trigs(ia  )
         s1 = trigs(ia+1)
         c2 = trigs(ia+ia-1)
         s2 = trigs(ia+ia  )

         do l = 1 , lot
            hbpc = a(ia  ,l) - 0.5 * (a(ib  ,l) + a(ic  ,l))
            hbmc = a(ia+1,l) - 0.5 * (a(ib+1,l) - a(ic+1,l))
            sbmc = SIN60 * (a(ib  ,l) - a(ic  ,l))
            sbpc = SIN60 * (a(ib+1,l) + a(ic+1,l))

            c(j  ,l) = a(ia  ,l) + a(ib  ,l) + a(ic  ,l)
            c(j+3,l) = a(ia+1,l) + a(ib+1,l) - a(ic+1,l)
            c(j+1,l) = c1 * (hbpc-sbpc) - s1 * (hbmc+sbmc)
            c(j+4,l) = s1 * (hbpc-sbpc) + c1 * (hbmc+sbmc)
            c(j+2,l) = c2 * (hbpc+sbpc) - s2 * (hbmc-sbmc)
            c(j+5,l) = s2 * (hbpc+sbpc) + c2 * (hbmc-sbmc)
         enddo
         ia = ia + 2
         ib = ib + 2
         ic = ic - 2
      enddo

      c(n-2,:) = a(ia,:)
      c(n-1,:) =   0.5 * a(ia,:) - SIN60 * a(ia+1,:)
      c(n  ,:) = - 0.5 * a(ia,:) - SIN60 * a(ia+1,:)

      la = 3
      a  = c
      end subroutine ifft3_acc


!     ================
!     SUBROUTINE IFFT4
!     ================

pure  subroutine ifft4_acc(c,trigs,n,lot,la)
      implicit none
      integer, intent(in) :: n, lot
      integer, intent(inout) :: la
      real, dimension(n), intent(in) :: trigs
      real, dimension(n,lot), intent(inout) :: c

      ! local
      real, parameter :: SIN45=sqrt(0.5)
      real, dimension(n,lot) :: a
      integer :: m, kstop, kb, kc, kd
      integer :: i0,i1,i2,i3,i4,i5,i6,i7, iink
      integer ::    j1,j2,j3,j4,j5,j6,j7, jbase
      real :: c1,s1, c2,s2, c3,s3
      real :: a0p2,a0m2, a1p3,a1m3, a4p6,a4m6, a5p7,a5m7, a0p2m1p3,a4m6m5m7
      integer :: i,j,k,l

      if (la == 1) then
         a(1,:) = 0.5 * c(1,:)
         a(n,:) = 0.0
         a(2:n-1,:) = c(3:n,:)
      else
         a = c
      endif

      m=n/4
      kstop=(n-4)/8

      i1 = n/2 - la
      i2 = n   - la
      i5 = i1  + la

      j1 = la
      j2 = la+j1
      j3 = la+j2
      j4 = la+j3
      j5 = la+j4
      j6 = la+j5
      j7 = la+j6

      do i=1,la
      do l=1,lot
         c(   i,l) = a(i,l) + a(i2+i,l) + a(i1+i,l)
         c(j1+i,l) = a(i,l) - a(i2+i,l) - a(i5+i,l)
         c(j2+i,l) = a(i,l) + a(i2+i,l) - a(i1+i,l)
         c(j3+i,l) = a(i,l) - a(i2+i,l) + a(i5+i,l)
      enddo
      enddo

      iink  = 2 * la
      jbase = 4 * la + 1
      i0    = la
      i1    = i0 + n/2
      i2    = n - 3 * la
      i3    = i2 - n/2
      i4    = i0 + la
      i5    = i1 + la
      i6    = i2 + la
      i7    = i3 + la

      do k=la,kstop,la
         kb=k+k
         kc=kb+kb
         kd=kc+kb
         c1=trigs(kb+1)
         s1=trigs(kb+2)
         c2=trigs(kc+1)
         s2=trigs(kc+2)
         c3=trigs(kd+1)
         s3=trigs(kd+2)
         do i = 1 , la
            j = jbase
            do l=1,lot
               a0p2 = a(i0+i,l) + a(i2+i,l)
               a0m2 = a(i0+i,l) - a(i2+i,l)
               a1p3 = a(i1+i,l) + a(i3+i,l)
               a1m3 = a(i1+i,l) - a(i3+i,l)
               a4p6 = a(i4+i,l) + a(i6+i,l)
               a4m6 = a(i4+i,l) - a(i6+i,l)
               a5p7 = a(i5+i,l) + a(i7+i,l)
               a5m7 = a(i5+i,l) - a(i7+i,l)
 
               a0p2m1p3 = a0p2 - a1p3
               a4m6m5m7 = a4m6 - a5m7

               c(   j,l) = a0p2 + a1p3
               c(j4+j,l) = a4m6 + a5m7
               c(j2+j,l) = c2 * a0p2m1p3 - s2 * a4m6m5m7
               c(j6+j,l) = s2 * a0p2m1p3 + c2 * a4m6m5m7
               c(j1+j,l) = c1*(a0m2-a5p7)-s1*(a4p6+a1m3)
               c(j5+j,l) = s1*(a0m2-a5p7)+c1*(a4p6+a1m3)
               c(j3+j,l) = c3*(a0m2+a5p7)-s3*(a4p6-a1m3)
               c(j7+j,l) = s3*(a0m2+a5p7)+c3*(a4p6-a1m3)
            enddo
            jbase=jbase+1
         enddo
         i0 = i0 + iink
         i1 = i1 + iink
         i2 = i2 - iink
         i3 = i3 - iink
         i4 = i4 + iink
         i5 = i5 + iink
         i6 = i6 - iink
         i7 = i7 - iink
         jbase=jbase+7*la
      end do

      if (i1 <= i2) then
         do i=1,la
            j=jbase
            do l=1,lot
               c(   j,l)=a(i0+i,l)+a(i1+i,l)
               c(j1+j,l)=sin45*((a(i0+i,l)-a(i1+i,l))-(a(la+i0+i,l)+a(la+i1+i,l)))
               c(j2+j,l)=a(la+i1+i,l)-a(la+i0+i,l)
               c(j3+j,l)=-sin45*((a(i0+i,l)-a(i1+i,l))+(a(la+i0+i,l)+a(la+i1+i,l)))
            enddo
            jbase=jbase+1
         enddo
      endif

      la = la * 4
      end subroutine ifft4_acc


!     ================
!     SUBROUTINE IFFT8
!     ================

pure  subroutine ifft8_acc(a,c,n,lot)
!$acc routine worker
      implicit none
      integer, intent(in) :: n, lot
      real, dimension(n,lot), intent(in)  :: a
      real, dimension(n,lot), intent(out) :: c

      ! local
      real, parameter :: SQRT2 = 1.414213562373095D0
      integer :: la, i0,i1,i2,i3,i4,i5,i6,i7
      real, dimension(lot) :: a0p7,a0m7, a1p5,a1m5, a2p6,a2m6
      real, dimension(lot) :: a0p7p3,a0p7m3, a0m7p4,a0m7m4, a1m5p2p6,a1m5m2p6
      character(len=200) :: tmp
      integer :: i

      la = n / 8

      do i=1, la
         i0 = i
         i1 = i0 + la
         i2 = i1 + la
         i3 = i2 + la
         i4 = i3 + la
         i5 = i4 + la
         i6 = i5 + la
         i7 = i6 + la

         a0p7 = a(i0,:) + a(i7,:)
         a0m7 = a(i0,:) - a(i7,:)
         a1p5 = a(i1,:) + a(i5,:)
         a1m5 = a(i1,:) - a(i5,:)
         a2p6 = a(i2,:) + a(i6,:)
         a2m6 = a(i2,:) - a(i6,:)

         a0p7p3   = a0p7 + a(i3,:)
         a0p7m3   = a0p7 - a(i3,:)
         a0m7p4   = 2.0 * (a0m7 + a(i4,:))
         a0m7m4   = 2.0 * (a0m7 - a(i4,:))
         a1m5p2p6 = SQRT2 * (a1m5 + a2p6)
         a1m5m2p6 = SQRT2 * (a1m5 - a2p6)

         c(i0,:)  = 2.0 * (a0p7p3 + a1p5)
         c(i2,:)  = 2.0 * (a0p7m3 - a2m6)
         c(i4,:)  = 2.0 * (a0p7p3 - a1p5)
         c(i6,:)  = 2.0 * (a0p7m3 + a2m6)

         c(i1,:)  = a0m7m4 + a1m5m2p6
         c(i3,:)  = a0m7p4 - a1m5p2p6
         c(i5,:)  = a0m7m4 - a1m5m2p6
         c(i7,:)  = a0m7p4 + a1m5p2p6
      end do

      end subroutine ifft8_acc




! =================
! SUBROUTINE MKTEND
! =================
! XW(2017/5/1): 把单层改为对垂向循环

pure &
subroutine mktend_acc(d,t,z,tn,fu,fv,ke,ut,vt,   qq,qe,qc,qx, nlon,nlat,nlev,nhpp,ntp1,nesp,ncsp)
!$---acc routine worker
implicit none

complex, intent(in) :: tn(nlon,nhpp,nlev)
complex, intent(in) :: fu(nlon,nhpp,nlev)
complex, intent(in) :: fv(nlon,nhpp,nlev)
complex, intent(in) :: ke(nlon,nhpp,nlev)
complex, intent(in) :: ut(nlon,nhpp,nlev)
complex, intent(in) :: vt(nlon,nhpp,nlev)

complex, intent(out) :: d(nesp/2,nlev)
complex, intent(out) :: t(nesp/2,nlev)
complex, intent(out) :: z(nesp/2,nlev)

real,    intent(in)  :: qq(ncsp,nhpp),qe(ncsp,nhpp),qc(ncsp,nhpp)
complex, intent(in)  :: qx(ncsp,nhpp)
integer, intent(in)  :: nlon,nlat,nlev,nhpp,ntp1,nesp,ncsp

integer :: l      ! Loop index for latitude
integer :: m      ! Loop index for zonal wavenumber m
integer :: w      ! Loop index for spectral mode
integer :: e      ! End index for w
integer :: jlev   ! Main loop along with Z/level

complex :: fus,fua,fvs,fva,kes,kea,tns,tna,uts,uta,vts,vta


d = (0.0,0.0) ! divergence
t = (0.0,0.0) ! temperature
z = (0.0,0.0) ! vorticity

do jlev = 1, nlev
   do l = 1 , nhpp  ! process pairs of Nort-South latitudes
      w = 1
      do m = 1 , ntp1
         kes = ke(m,l,jlev) + ke(m+nlat,l,jlev) ; kea = ke(m,l,jlev) - ke(m+nlat,l,jlev)
         fvs = fv(m,l,jlev) + fv(m+nlat,l,jlev) ; fva = fv(m,l,jlev) - fv(m+nlat,l,jlev)
         fus = fu(m,l,jlev) + fu(m+nlat,l,jlev) ; fua = fu(m,l,jlev) - fu(m+nlat,l,jlev)
         uts = ut(m,l,jlev) + ut(m+nlat,l,jlev) ; uta = ut(m,l,jlev) - ut(m+nlat,l,jlev)
         vts = vt(m,l,jlev) + vt(m+nlat,l,jlev) ; vta = vt(m,l,jlev) - vt(m+nlat,l,jlev)
         tns = tn(m,l,jlev) + tn(m+nlat,l,jlev) ; tna = tn(m,l,jlev) - tn(m+nlat,l,jlev)
         e = w + ntp1 - m    ! vector of symmetric modes
         d(w:e:2,jlev) = d(w:e:2,jlev) + qq(w:e:2,l) * kes - qe(w:e:2,l) * fva + qx(w:e:2,l) * fus
         t(w:e:2,jlev) = t(w:e:2,jlev) + qe(w:e:2,l) * vta + qc(w:e:2,l) * tns - qx(w:e:2,l) * uts
         z(w:e:2,jlev) = z(w:e:2,jlev) + qe(w:e:2,l) * fua + qx(w:e:2,l) * fvs
         w = w + 1           ! vector of antisymmetric modes
         d(w:e:2,jlev) = d(w:e:2,jlev) + qq(w:e:2,l) * kea - qe(w:e:2,l) * fvs + qx(w:e:2,l) * fua
         t(w:e:2,jlev) = t(w:e:2,jlev) + qe(w:e:2,l) * vts + qc(w:e:2,l) * tna - qx(w:e:2,l) * uta
         z(w:e:2,jlev) = z(w:e:2,jlev) + qe(w:e:2,l) * fus + qx(w:e:2,l) * fva
         w = e + 1
      enddo ! m
   enddo ! l
end do

end subroutine mktend_acc

end module dycore
