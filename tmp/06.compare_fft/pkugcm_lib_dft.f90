!module pkugcm_lib_dft
!
!	implicit none
!
!	real, parameter 	:: Pi = 3.1415926
!	complex, parameter	:: jj = (0,1)
!
!contains

subroutine dft1d_cplx(x,N,coef)
	implicit none
	integer, intent(in)			:: N
	complex, dimension(1:N), intent(in)	:: x
	complex, dimension(0:N-1), intent(out)	:: coef

	integer :: i,k
	real, parameter 	:: Pi = 3.1415926
	complex, parameter	:: jj = (0,1)

	do k = 0, N-1
		coef(k) = (0,0)
		do i = 1, N
			coef(k)	= coef(k) + x(i)*exp(-2.0*jj*Pi*k*i/N)
		end do
      coef(k) = coef(k)/real(N)
	end do
end subroutine

subroutine idft1d_cplx(coef,N,x)
	implicit none
	integer, intent(in)			:: N
	complex, dimension(0:N-1), intent(in)	:: coef
	complex, dimension(1:N), intent(out)	:: x

	integer :: i, k
	real, parameter 	:: Pi = 3.1415926
	complex, parameter	:: jj = (0,1)

	do i = 1, N
		x(i) = (0,0)
		do k = 0, N-1
			x(i) = x(i) + coef(k)*exp(2.0*jj*Pi*k*i/N)
		end do
	end do
end subroutine
	
!end module
