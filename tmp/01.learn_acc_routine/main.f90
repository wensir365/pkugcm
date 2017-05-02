module try
contains
subroutine add(x,y,z,N)
!$acc routine worker
   implicit none
   integer, intent(in) :: N
   real, dimension(N,N), intent(in) :: x,y
   real, dimension(N,N), intent(out) :: z
   real, dimension(N,N) :: c
   integer :: i,j
   do i = 1, N
      do j = 1, N
         !c(i,j) = sin(0.5*i)*cos(0.5*j)
         z(i,j) = x(i,j)+y(i,j)!*c(i,j)
      end do
   end do
end subroutine
end module


!==============================================================================
! Correct outputs should like this:
!==============================================================================
! max/min of a & b:    1.000000       -1.000000       0.9999962     -0.9999962    
!
! max/min of a+b:    1.928196       -1.928195    
! max/min of a-b:    1.928196       -1.928195    
! max/min of a*b:   0.9284945      -0.9284945    
!==============================================================================
! What I learn is:
! - subroutine只有放到contains里才有用 external subroutine无效
! - 单独call会给出错误结果 应放在至少一个循环的外套里 I know its stupid, but work!
!==============================================================================

program main
   use try
   implicit none
   integer, parameter :: M = 9999
   real, dimension(M,M) :: a, b, c_add, c_minus, c_times, c_nest
   integer :: i,j

   call setup(a,b,M)
   print *, "max/min of a & b:", maxval(a), minval(a), maxval(b), minval(b)
   print *, ""

!$acc kernels

   do i = 1, 1
      call add(a,b,c_add,M)
      call times(a,b,c_times,M)     ! 放在loop外面就会导致计算结果全是零 why?!
   end do
   call minus(a,b,c_minus,M)

!$acc end kernels

   call nest_times(a,b,c_nest,M)

   print *, "max/min of a+b:", maxval(c_add), minval(c_add)
   print *, "max/min of a-b:", maxval(c_minus), minval(c_minus)
   print *, "max/min of a*b:", maxval(c_times), minval(c_times)
   print *, "max/min of a*b (nested):", maxval(c_times), minval(c_times)

contains

subroutine setup(x,y,N)
!$acc routine worker
   implicit none
   integer, intent(in) :: N
   real, dimension(N,N), intent(inout) :: x,y
   integer :: i,j
   do i = 1, N
      do j = 1, N
         x(i,j) = sin(0.5*i)*cos(2.0*j)
         y(i,j) = sin(2.0*i)*cos(0.5*j)
      end do
   end do
end subroutine

function times_line(x,y,N)
   implicit none
   integer, intent(in) :: N
   real, dimension(N), intent(in) :: x,y
   real, dimension(N), intent(out) :: times_line
   integer :: j
!$acc kernels
   do j = 1, N
      times_line(j) = x(j)*y(j)
   end do
!$acc end kernels
end function times_line

subroutine nest_times(x,y,z,N)
   implicit none
   integer, intent(in) :: N
   real, dimension(N,N), intent(in) :: x,y
   real, dimension(N,N), intent(out) :: z
   real, dimension(N) :: xx,yy,zz
   integer :: i
   do i = 1, N
      z(i,:) = times_line(x(i,:),y(i,:),N)
   end do
end subroutine nest_times

subroutine minus(x,y,z,N)
   implicit none
   integer, intent(in) :: N
   real, dimension(N,N), intent(in) :: x,y
   real, dimension(N,N), intent(out) :: z
   integer :: i,j
!$acc kernels
   do i = 1, N
      do j = 1, N
         z(i,j) = x(i,j)-y(i,j)
      end do
   end do
!$acc end kernels
end subroutine

subroutine times(x,y,z,N)
!$acc routine worker
   implicit none
   integer, intent(in) :: N
   real, dimension(N,N), intent(in) :: x,y
   real, dimension(N,N), intent(out) :: z
   integer :: i,j
   do i = 1, N
      do j = 1, N
         z(i,j) = x(i,j)*y(i,j)
     end do
   end do
end subroutine

end program

