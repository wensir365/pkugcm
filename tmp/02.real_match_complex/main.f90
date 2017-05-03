program main

   implicit none
   real, allocatable, dimension(:) :: x, y
   complex :: z

   allocate(x(4))
   allocate(y(4))

   x = (/ 1, 2, 3, 4 /)
   call add(x,y)

   print *, "x=", x
   print *, "y=", y

   z = cmplx(1,2)
   print *, "z=", z
   print *, "real(z)=", real(z)
   print *, "realpart(z)=", real(z)
   print *, "imagpart(z)=", aimag(z)

end program

subroutine add(a,b)
   implicit none
   complex, dimension(2), intent(in) :: a
   complex, dimension(2), intent(out) :: b
   b = a+1
end subroutine add

