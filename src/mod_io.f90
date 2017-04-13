!====================================
! PKUGCM MODULE "io"
! -----------------------------------
! Creat, collect, and write outputs
! Xinyu Wen, Peking Univ, Apr/13/2017
!====================================

module io

   implicit none

   integer, parameter :: Nvar = 6

   integer, dimension(Nvar)           :: lev = (/ 1   , 10  , 10  , 10  , 10  , 10  /)
   integer, dimension(Nvar)           :: fid = (/ 901 , 902 , 903 , 904 , 905 , 906 /)
   character(len=3),  dimension(Nvar) :: vid = (/"ps ","u  ","v  ","div","vor","t  "/)
   character(len=30), dimension(Nvar) :: desc= (/"Surface Pressure (mb)         ", &
                                                 "U Wind (m/s)                  ", &
                                                 "V Wind (m/s)                  ", &
                                                 "Div (m/s^2)                   ", &
                                                 "Vor (m/s^2)                   ", &
                                                 "Temperature (K)               " /)
end module io


!=====================
! io_open_output
!=====================

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

   do i = 1, Nvar
      f = fid(i)
      ! close binary data file
      close(f)
      ! creat GrADS header file
      fn = "output_"//vid(i)
      open(unit=f, file=trim(fn)//".ctl", access="sequential", form="formatted", status="replace")
         write(f,*) "DSET  ^"//trim(fn)//".bin"
         write(f,*) "TITLE "//trim(vid(i))
         write(f,*) "OPTIONS yrev zrev"
         write(f,*) "UNDEF -999"
         write(f,*) "XDEF  ",NLON," LINEAR    2.8125  5.625"
         write(f,*) "YDEF  ",NLAT," LINEAR  -87.1875  5.625"
         if (vid(i)=="ps  ") then
            write(f,*) "ZDEF  1 LEVELS 1000"
         else
            write(f,*) "ZDEF ",NLEV," LEVELS 1000 900 800 700 600 500 400 300 200 100"
         end if
         write(f,*) "TDEF 360 LINEAR 1jan0000 1dy"
         write(f,*) "VARS   1" 
         write(f,*) trim(vid(i))//"   ", lev(i), "   99   "//trim(desc(i))
         write(f,*) "ENDVARS"
      close(fid(i))
   end do
end subroutine io_close_output

