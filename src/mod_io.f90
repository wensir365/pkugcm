module io
   implicit none

   integer, parameter :: Nvar = 6

   integer, dimension(Nvar)          :: lev = (/ 1    , 10   , 10   , 10   , 10   , 10   /)
   character(len=4), dimension(Nvar) :: vid = (/"ps  ","u   ","v   ","div ","vor ","t   "/)
   integer, dimension(Nvar)          :: fid = (/ 901  , 902  , 903  , 904  , 905  , 906  /)
                                      ! fid = 901 : Surface Pressure (mb)
                                      ! fid = 902 : U Wind (m/s)
                                      ! fid = 903 : V Wind (m/s)
                                      ! fid = 904 : Div (m/s^2)
                                      ! fid = 905 : Vor (m/s^2)
                                      ! fid = 906 : Temperature (K)

end module io

subroutine io_output
   use io
   use pumamod, only: NLON,NLAT,NLEV, gp,gu,gv,gd,gz,gt
   implicit none

   ! local
   real, dimension(NLON*NLAT)      :: g2d
   real, dimension(NLON*NLAT,NLEV) :: g3d

   g2d = gp
   call alt2reg(g2d,1)
   write(fid(1)) g2d       ! gp

   g3d = gu
   call alt2reg(g3d,NLEV)
   write(fid(2)) g3d       ! gu

   g3d = gv
   call alt2reg(g3d,NLEV)
   write(fid(3)) g3d       ! gv

   g3d = gd
   call alt2reg(g3d,NLEV)
   write(fid(4)) g3d       ! gd
   
   g3d = gz
   call alt2reg(g3d,NLEV)
   write(fid(5)) g3d       ! gz
   
   g3d = gt
   call alt2reg(g3d,NLEV)
   write(fid(6)) g3d       ! gt

end subroutine io_output

subroutine io_open_output
   use io
   implicit none
   ! local
   integer :: i,f
   character(len=50) :: fn

   do i = 1, Nvar
      f  = fid(i)
      fn = "output_"//trim(vid(i))//".bin"
      open(unit=f, file=trim(fn), access="stream", form="unformatted", status="replace")
   end do
end subroutine io_open_output

subroutine io_close_output
   use io
   implicit none
   ! local
   integer :: i,f
   character(len=50) :: fn

   do i = 1, Nvar
      f = fid(i)
      ! close binary data file
      close(f)
      ! creat GrADS header file
      fn = "output_"//vid(i)
      open(unit=f, file=trim(fn)//".ctl", access="sequential", form="formatted", status="replace")
         write(f,*) "DSET  ^"//trim(fn)//".bin"
         write(f,*) "TITLE "//trim(vid(i))
         write(f,*) "OPTIONS yrev"
         write(f,*) "UNDEF -999"
         write(f,*) "XDEF  64 LINEAR    2.8125  5.625"
         write(f,*) "YDEF  32 LINEAR  -87.1875  5.625"
         if (vid(i)=="ps  ") then
            write(f,*) "ZDEF  1 LEVELS 1000"
         else
            write(f,*) "ZDEF 10 LEVELS 1000 900 800 700 600 500 400 300 200 100"
         end if
         write(f,*) "TDEF 360 LINEAR 1jan0000 1dy"
         write(f,*) "VARS   1" 
         write(f,*) trim(vid(i))//"   ", lev(i), "  99  "//trim(vid(i))
         write(f,*) "ENDVARS"
      close(fid(i))
   end do
end subroutine io_close_output

