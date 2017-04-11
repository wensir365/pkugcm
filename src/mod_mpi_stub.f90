!     =======================================
!     mpimod_dummy.f90
!     ----------------
!     This module replaces <mpimod.f90> for 
!     single CPU runs
!
!     The module is shared by PUMA and PlaSim
!     =======================================

      subroutine mrdimensions
      return
      end

      subroutine mrdiff(p,d,n,l)
      integer, intent(in) :: n, l
      real, dimension(n), intent(in) :: p
      real, dimension(n), intent(in) :: d
      return
      end

      subroutine mrsum(k) ! sum up 1 integer
      integer, intent(in) :: k
      return
      end

      subroutine mrbci(k) ! broadcast 1 integer
      integer, intent(in) :: k
      return
      end

      subroutine mpbci(k) ! broadcast 1 integer
      integer, intent(in) :: k
      return
      end

      subroutine mpbcin(k,n) ! broadcast n integer
      integer, intent(in) :: n
      integer, dimension(n), intent(in) :: k
      return
      end

      subroutine mpbcr(p) ! broadcast 1 real
      integer, intent(in) :: p
      return
      end

      subroutine mpbcrn(p,n) ! broadcast n real
      integer, intent(in) :: n
      real, dimension(n), intent(in) :: p
      return
      end

      subroutine mpbcl(k) ! broadcast 1 logical
      logical, intent(in) :: k
      return
      end

      subroutine mpscin(k,n) ! scatter n integer
      integer, intent(in) :: n
      integer, dimension(n), intent(in) :: k
      return
      end

      subroutine mpscrn(p,n) ! scatter n real
      integer, intent(in) :: n
      real, dimension(n), intent(in) :: p
      return
      end

      subroutine mpscdn(p,n) ! scatter n double precision
      integer, intent(in) :: n
      real, dimension(n), intent(in) :: p
      return
      end

      subroutine mpscsp(pf,pp,klev) ! scatter spectral fields
      use pumamod
      real pf(NESP,klev)
      real pp(NSPP,klev)
      pp(1:NSPP,1:klev) = pf(1:NSPP,1:klev)
      return
      end

      subroutine mpscgp(pf,pp,klev) ! scatter gridpoint fields
      use pumamod
      real pf(NLON*NLAT,klev)
      real pp(NHOR,klev)
      pp(1:NHOR,1:klev) = pf(1:NHOR,1:klev)
      return
      end

      subroutine mpgasp(pf,pp,klev) ! gather spectral fields
      use pumamod
      real pf(NESP,klev)
      real pp(NSPP,klev)
      pf(1:NSPP,1:klev) = pp(1:NSPP,1:klev)
      return
      end

      subroutine mpgagp(pf,pp,klev) ! gather gridpoint fields
      use pumamod
      real pf(NHOR,klev)
      real pp(NHOR,klev)
      pf = pp
      return
      end

      subroutine mpgacs(pcs) ! gather cross sections
      return
      end

      subroutine mpgallsp(pf,pp,klev) ! gather spectral to all
      use pumamod
      real pf(NESP,klev)
      real pp(NSPP,klev)
      pf(1:NSPP,1:klev) = pp(1:NSPP,1:klev)
      return
      end

      subroutine mpsum(psp,klev) ! sum spectral fields
      return
      end

      subroutine mpsumsc(psf,psp,klev) ! sum & scatter spectral
      use pumamod
      real psf(NESP,klev)
      real psp(NSPP,klev)
      psp(1:NSPP,1:klev) = psf(1:NSPP,1:klev)
      return
      end

      subroutine mpsumr(pr,kdim) ! sum kdim reals
      return
      end subroutine mpsumr

      subroutine mpsumbcr(pr,kdim) ! sum & broadcast kdim reals
      return
      end

      subroutine mpstart ! initialization
      use pumamod
      npro = 1
      return
      end

      subroutine mpstop
      return
      end

      subroutine mpreadsp(ktape,p,kdim,klev)
      real p(kdim,klev)
      read (ktape) p
      return
      end

      subroutine mpreadgp(ktape,p,kdim,klev)
      real p(kdim,klev)
      read (ktape) p
      return
      end

      subroutine mpwritesp(ktape,p,kdim,klev)
      real p(kdim,klev)
      write (ktape) p
      return
      end

      subroutine mpwritegp(ktape,p,kdim,klev)
      real p(kdim,klev)
      write (ktape) p
      return
      end

      subroutine mpwritegph(ktape,p,kdim,klev,ihead)
      real :: p(kdim,klev)
      integer :: ihead(8)
      write (ktape) ihead
      write (ktape) p

      return
      end


      subroutine mpi_info(nprocess,pid)    ! get nproc and pid
      integer nprocess, pid
      nprocess = 1
      pid = 0
      return
      end subroutine mpi_info


      subroutine mpgetsp(yn,p,kdim,klev)
      character (len=*) :: yn
      real :: p(kdim,klev)
      call get_restart_array(yn,p,kdim,kdim,klev)
      return
      end subroutine mpgetsp


      subroutine mpgetgp(yn,p,kdim,klev)
      character (len=*) :: yn
      real :: p(kdim,klev)
      call get_restart_array(yn,p,kdim,kdim,klev)
      return
      end subroutine mpgetgp


      subroutine mpputsp(yn,p,kdim,klev)
      character (len=*) :: yn
      real :: p(kdim,klev)
      call put_restart_array(yn,p,kdim,kdim,klev)
      return
      end subroutine mpputsp


      subroutine mpputgp(yn,p,kdim,klev)
      character (len=*) :: yn
      real :: p(kdim,klev)
      call put_restart_array(yn,p,kdim,kdim,klev)
      return
      end subroutine mpputgp


      subroutine mpmaxval(p,kdim,klev,pmax)
      real :: p(kdim,klev)
      pmax = maxval(p(:,:))
      return
      end subroutine mpmaxval


      subroutine mpsumval(p,kdim,klev,psum)
      real :: p(kdim,klev)
      psum = sum(p(:,:))
      return
      end subroutine mpsumval

