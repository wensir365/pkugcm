!====================================
! PKUGCM MODULE "io"
! -----------------------------------
! Open, collect, and write outputs
! Xinyu Wen, Peking Univ, Apr/13/2017
!====================================

module io

   implicit none

   integer, parameter :: Nvar    = 7

   logical, dimension(Nvar)           :: won = (/.True.,.True.,.True.,.False.,.False.,.True.,.True./)
   integer, dimension(Nvar)           :: lev = (/ 1   , 10  , 10  , 10  , 10  , 10  ,  10 /)
   integer, dimension(Nvar)           :: fid = (/ 901 , 902 , 903 , 904 , 905 , 906 , 907 /)
   character(len=3),  dimension(Nvar) :: vid = (/"ps ","u  ","v  ","div","vor","t  ","z  "/)
   character(len=30), dimension(Nvar) :: desc= (/"Surface Pressure (mb)         ", &
                                                 "U Wind (m/s)                  ", &
                                                 "V Wind (m/s)                  ", &
                                                 "Div (m/s^2)                   ", &
                                                 "Vor (m/s^2)                   ", &
                                                 "Temperature (K)               ", &
                                                 "Geopotential Height (gpm)     " /)
   integer(kind=8) :: writecount=0, yearid=-1
   logical :: needctl = .FALSE.
end module io


!=====================
! io_write_output
!=====================

subroutine io_write_output
   use io
   use pumamod, only: sp,sd,sz,st,gu,gv,                 &
                      psurf,t0,ga,gascon,ww_scale,cv,ct, &
                      csq,sigma,sigmh,                   &
                      NLON,NLAT,NLEV,trigs,TAC
   implicit none

   ! local
   real, dimension(NLON,NLAT,NLEV)  :: x3d, y3d
   real, dimension(NLON,NLAT,0:NLEV):: x3dh,y3dh
   real, dimension(NLON,NLAT)       :: x2d, y2d
   real, dimension(0:NLEV)          :: sigmh2
   real, dimension(NLAT)            :: scal
   real                             :: offs
   integer :: jlon,jlat,jlev

   !--- creat new? or append to current?
   if (mod(writecount,360)==0) then     ! creat new
      if (writecount.gt.0) call io_flush_output
      yearid   = yearid + 1
      call io_open_output
   end if

   !--- Surface Pressure ---
   call sp2fc(sp,x2d)
   call fc2gp(x2d,NLON,NLAT,trigs)
   call alt2reg(x2d,1)
   x2d = exp(x2d)*psurf    ! psurf=1011mb
   x2d = x2d / 100.0       ! Pa to hPa/mb
   
   if (won(1)) write(fid(1)) x2d

   !---  U and V ---
   if (won(2).and.won(3)) then
   x3d   = reshape(gu(:,:),(/NLON,NLAT,NLEV/))     ! gu(NLONxNLAT,NLEV)
   y3d   = reshape(gv(:,:),(/NLON,NLAT,NLEV/))     ! gv(NLONxNLAT,NLEV)
   scal  = cv/sqrt(csq(:))
   offs  = 0.0

   do jlev = 1, NLEV
      do jlon = 1, NLON
         x3d(jlon,:,jlev) = x3d(jlon,:,jlev) * scal(:) + offs
         y3d(jlon,:,jlev) = y3d(jlon,:,jlev) * scal(:) + offs
      end do
   end do

   call alt2reg(x3d(:,:,:),NLEV)
   call alt2reg(y3d(:,:,:),NLEV)

   write(fid(2)) x3d
   write(fid(3)) y3d
   end if

   !--- Div and Vor ---
   if (won(4).and.won(5)) then
   do jlev = 1, NLEV
      ! div
      call sp2fc(sd(:,jlev),x3d(:,:,jlev))
      call fc2gp(x3d(:,:,jlev),NLON,NLAT,trigs)
      ! vor
      call sp2fc(sz(:,jlev),y3d(:,:,jlev))
      call fc2gp(y3d(:,:,jlev),NLON,NLAT,trigs)
   end do

   offs = 0.0
   x3d = x3d * ww_scale + offs
   y3d = y3d * ww_scale + offs
         
   call alt2reg(x3d,NLEV)
   call alt2reg(y3d,NLEV)

   write(fid(4)) x3d
   write(fid(5)) y3d
   end if

   !--- Temperature ---
   if (won(6)) then
   do jlev = 1, NLEV
      call sp2fc(st(:,jlev),x3d(:,:,jlev))
      call fc2gp(x3d(:,:,jlev),NLON,NLAT,trigs)
      x3d(:,:,jlev) = x3d(:,:,jlev) * ct + t0(jlev)*ct
   end do

   call alt2reg(x3d,NLEV)

   write(fid(6)) x3d
   end if

   !--- Geopotential at half ---
   if (won(6).and.won(7)) then
   x3d = x3d * gascon      ! depends on above temperature
   sigmh2(1:NLEV) = sigmh
   sigmh2(0)      = 0.0

   y3dh(:,:,NLEV) = 0.0    ! 这是地表位势 以后需要求精
   do jlev = NLEV-1, 1, -1
      y3dh (:,:,jlev)= y3dh(:,:,jlev+1)+x3d(:,:,jlev+1)*log(sigmh2(jlev+1)/sigmh2(jlev))
   end do

   ! 最顶层如果直接算就无穷远了 所以最顶层只算一半高度
   y3dh(:,:,0) = y3dh(:,:,1)+x3d(:,:,1)*log(sigmh2(1)/sigma(1))

   y3dh  = y3dh / ga  ! geopotential (m2/s2) to height (gpm)

   write(fid(7)) y3dh(:,:,0:NLEV-1)    ! 不写表面高度场
   end if

   !--- count ---
   writecount = writecount + 1

end subroutine


!=====================
! io_open_output
!=====================

subroutine io_open_output
   use io
   implicit none

   ! local
   integer :: i,f
   character(len=50) :: fn
   character(len=6)  :: yearstr

   write(yearstr,"(i4.4)") yearid
   !write(yearstr,"(i6.6)") yearid

   do i = 1, Nvar
      f  = fid(i)
      fn = "output_"//trim(vid(i))//"_"//trim(yearstr)//".bin"
      if (won(i)) then
         open(unit=f, file=trim(fn), access="stream", form="unformatted", status="replace")
         print *, "open new output: "//trim(fn)
      end if
   end do
end subroutine

!=====================
! io_flush_output
!=====================

subroutine io_flush_output
   use io
   implicit none
   integer :: i
   do i = 1, Nvar
      if (won(i)) close(fid(i))
   end do
end subroutine


!=====================
! io_close_output
!=====================

subroutine io_close_output
   use io
   use pumamod, only: NLON,NLAT,NLEV
   implicit none

   ! local
   integer :: i,f
   character(len=50) :: fn
   real :: dd, halfdd
   logical :: ifopen=.TRUE.

   dd     = 360.0/NLON
   halfdd = dd/2.0

   do i = 1, Nvar
      f = fid(i)
      ! close binary data file if it is still open
      inquire(f,OPENED=ifopen)
      if (ifopen) then
         print *, "closing ... ", f, vid(i), desc(i)
         close(f)
      end if

      if (needctl) then
      ! creat GrADS header file
      fn = "output_"//vid(i)
      open(unit=f, file=trim(fn)//".ctl", access="sequential", form="formatted", status="replace")
         write(f,*) "DSET  ^"//trim(fn)//"_%y4.bin"
         write(f,*) "TITLE "//trim(vid(i))
         write(f,*) "OPTIONS yrev zrev template"
         write(f,*) "UNDEF -999"
         write(f,*) "XDEF  ",NLON," LINEAR ",halfdd," ",dd
         write(f,*) "YDEF  ",NLAT," LINEAR ",(-90+halfdd)," ",dd
         if (vid(i)=="ps  ") then
            write(f,*) "ZDEF  1 LEVELS 1000"
         else
            !write(f,*) "ZDEF ",NLEV," LEVELS 1000 900 800 700 600 500 400 300 200 100"
            write(f,*) "ZDEF ",NLEV," LEVELS 950 850 750 650 550 450 350 250 150 50"
         end if
         write(f,*) "TDEF ",writecount," LINEAR 1jan0 1461mn"  ! 纯属为了迁就GrADS的弱智时间维定义 365.25day/360=1461min
         write(f,*) "VARS   1" 
         write(f,*) trim(vid(i))//"   ", lev(i), "   99   "//trim(desc(i))
         write(f,*) "ENDVARS"
      close(fid(i))
      end if

   end do
end subroutine

