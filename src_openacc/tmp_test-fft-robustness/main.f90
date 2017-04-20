program main

   implicit none

   integer, parameter :: N=64

   real, dimension(N) :: x,y
   integer :: i

   do i = 1, N
      x(i) = i
   end do
   y  = x

   !print "(16f8.2)", x
   !print *, "---"
   call gp2fc(x,n,1)
   call fc2gp(x,n,1)
   !print "(16f8.2)", x

   print "(8f8.2)", (x-y)

end program
