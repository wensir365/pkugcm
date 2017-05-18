program main
   implicit none

   integer, parameter :: N = 32
   real, dimension(N) :: xin, xc, xout, trig
   real, dimension(N) :: yin, yc, yout
   complex, dimension(N/2,N/2) :: trig_dft, trig_idft
   character(len=10) :: fmtstr = "(32f6.2)"

   ! init
   xin   = (/ 1, 2, 3, 4, 5, 6, 7, 8,   9, 10, 11, 12, 13, 14, 15, 16, &
              1, 2, 3, 4, 5, 6, 7, 8,   9, 10, 11, 12, 13, 14, 15, 16 /)
   yin   = xin


   ! gcm fft
   print *, "gcm FFT:"

   call fftini(N,trig)
   print fmtstr, xin

   call gp2fc(xin,N,1,trig)
   print fmtstr, xin
   call fc2gp(xin,N,1,trig)
   print fmtstr, xin

   call gp2fc(xin,N,1,trig)
   print fmtstr, xin
   call fc2gp(xin,N,1,trig)
   print fmtstr, xin

   ! pku dft
   print *, "------"
   print *, "pku DFT:"
   call dft1d_init(N/2,trig_dft,trig_idft)
   print fmtstr, yin

   call dft1d_cplx(yin,N/2,yc,trig_dft)
   print fmtstr, yc
   call idft1d_cplx(yc,N/2,yout,trig_idft)
   print fmtstr, yout

   call dft1d_cplx(yout,N/2,yc,trig_dft)
   print fmtstr, yc
   call idft1d_cplx(yc,N/2,yin,trig_idft)
   print fmtstr, yin

end program



pure &
subroutine dft1d_init(N,trig_dft,trig_idft)
   implicit none
   integer, intent(in) :: N
   complex, dimension(0:N-1,1:N), intent(out) :: trig_dft, trig_idft
   
   integer :: k,i
	real, parameter 	 :: Pi = 3.1415926
	complex, parameter :: jj = (0,1)

   do k = 0, N-1
      do i = 1, N
			trig_dft(k,i)  = exp(-2.0*jj*Pi*k*i/N)
			trig_idft(k,i) = exp( 2.0*jj*Pi*k*i/N)
      end do
   end do
end subroutine

pure &
subroutine dft1d_cplx(x,N,coef,trig)
	implicit none
	integer, intent(in)			:: N
	complex, dimension(1:N), intent(in)	:: x
	complex, dimension(0:N-1), intent(out)	:: coef
   complex, dimension(0:N-1,1:N), intent(in) :: trig

	integer :: i,k
	real, parameter 	:: Pi = 3.1415926
	complex, parameter	:: jj = (0,1)

	do k = 0, N-1
      !coef(k)  = dot_product(x,trig(k,:)) / real(N)
		coef(k) = (0,0)
		do i = 1, N
			!coef(k)	= coef(k) + x(i)*exp(-2.0*jj*Pi*k*i/N)
			coef(k)	= coef(k) + x(i)*trig(k,i)
		end do
      coef(k) = coef(k)/N
	end do
end subroutine

pure &
subroutine idft1d_cplx(coef,N,x,trig)
	implicit none
	integer, intent(in)			:: N
	complex, dimension(0:N-1), intent(in)	:: coef
	complex, dimension(1:N), intent(out)	:: x
   complex, dimension(0:N-1,1:N), intent(in) :: trig

	integer :: i, k
	real, parameter 	:: Pi = 3.1415926
	complex, parameter	:: jj = (0,1)

	do i = 1, N
      !x(i) = dot_product(coef,trig(:,i))
		x(i) = (0,0)
		do k = 0, N-1
			!x(i) = x(i) + coef(k)*exp(2.0*jj*Pi*k*i/N)
			x(i) = x(i) + coef(k)*trig(k,i)
		end do
	end do
end subroutine
	
