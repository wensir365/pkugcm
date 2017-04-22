module pumamod

!*********************************************!
! Peking University General Circulation Model !
! Xinyu Wen                                   !
! April, 2017                                 !
!*********************************************!
! is based on                                 !
! Portable University Model of the Atmosphere !
! Version: 17.0   16-Feb-2011                 !
! Klaus Fraedrich                             !
! Frank Lunkeit  - Edilbert Kirk              !
! Frank Sielmann - Torben Kunz                !
! Hartmut Borth                               !
!*********************************************!
! Meteorologisches Institut                   !
! KlimaCampus - Universitaet Hamburg          !
!*********************************************!
! http://www.mi.uni-hamburg.de/puma           !
!*********************************************!

! ****************************************************************
! * Don't touch the following parameter definitions !            *
! ****************************************************************
integer, parameter :: PUMA   = 0              ! Model ID
integer, parameter :: PLASIM = 1              ! Model ID
integer, parameter :: PKUGCM = 2              ! Peking Univ GCM

integer        :: model      = PKUGCM
character(80)  :: modelver   = "prototype (Apr/15/2017)"


! +++++++++++++++++++++++
! Xinyu added global var
! +++++++++++++++++++++++
!integer (kind=8) :: Nbinary=0      ! number of new outputs in GrADS format
!integer          :: outformat=1    ! output format (0=default; 1=default+GrADS)
                                    ! see another key switch: "noutput"
logical           :: mloop = .FALSE.! check if in main loop (T), or in startup (F)
                                    ! derived from FFT module, to facilitate OpenACC for FFT
real, allocatable :: trigs(:)       ! XW: trigs(NLON) constant fourier base functions to be set in fftini
                                    ! derived from legsym module, to facilitate OpenACC for SH
real   , allocatable :: qi(:,:)     ! P(m,n) = Associated Legendre Polynomials  used in sp2fc
real   , allocatable :: qj(:,:)     ! Q(m,n) = Used for d/d(mu)                 used in sp2fcdmu
real   , allocatable :: qc(:,:)     ! P(m,n) * gwd                              used in fc2sp
real   , allocatable :: qe(:,:)     ! Q(mn,) * gwd / cos2                       used in mktend
real   , allocatable :: qq(:,:)     ! P(m,n) * gwd / cos2 * n * (n+1) / 2       used in mktend
real   , allocatable :: qu(:,:)     ! P(m,n) / (n*(n+1)) * m                    used in dv2uv
real   , allocatable :: qv(:,:)     ! Q(m,n) / (n*(n+1))                        used in dv2uv
complex, allocatable :: qx(:,:)     ! P(m,n) * gwd / cos2 * m                   used in mktend & uv2dv


!**************************************************************!
! The number of processes for processing on parallel machines  !
! NLAT/2 must be dividable by <npro>. npro can be set by the   !
! option -n <npro> when calling the puma executable            !
! This option is only available if the code is compiled with   !
! an mpi compiler.                                             !
!**************************************************************!
integer :: npro = 1

!**************************************************************!
! The horizontal resolution of PUMA is set by defining the     !
! number of latitudes <nlev> with the 1st. command line        !
! parameter and the number of levels with the 2nd. command     !
! parameter. A typical call for T42 is:                        !
! puma.x 64 10                                                 !
! which sets nlat=64 and nlev=10                               !
!**************************************************************!
integer :: nlat = 32

!example values:  32,  48,  64, 128,  192,  256,  512,  1024
!truncation:     T21, T31, T42, T85, T127, T170, T341,  T682

integer :: nlev = 10 

!*****************************************************!
! Grid related paramters, which are computed from the !
! command line arguments <nlat> and <nlev>            !
! Preset values are for T21 (nlat=32) and nlev=10     !
! ****************************************************!

integer :: nlem =    9 ! Levels - 1
integer :: nlep =   11 ! Levels + 1
integer :: nlsq =  100 ! Levels squared

integer :: nlon =   64 ! Longitudes = 2 * latitudes
integer :: nlah =   16 ! Half of latitudes
integer :: ntru =   21 ! (nlon-1) / 3
integer :: ntp1 =   22 ! ntru + 1
integer :: nzom =   44 ! Number of zonal modes
integer :: nrsp =  506 ! (ntru+1) * (ntru+2)
integer :: ncsp =  253 ! nrsp / 2
integer :: nspp =  506 ! nodes per process
integer :: nesp =  506 ! number of extended modes

integer :: nlpp =   32 ! Latitudes per process
integer :: nhpp =   16 ! Half latitudes per process
integer :: nhor = 2048 ! Horizontal part (nlon x nlat)
integer :: nugp = 2048 ! Horizontal total
integer :: npgp = 1024 ! Horizontal total packed words
                       ! 南极和北极合并为一个复数 i.e. packed word

integer :: nud  =    6 ! I/O unit for diagnostic output (*_diag)

!***********!
! filenames !
!***********!
character (256) :: puma_namelist       = "pkugcm_namelist"
character (256) :: puma_output         = "pkugcm_output"
character (256) :: puma_diag           = "pkugcm_diag"
character (256) :: puma_restart        = "pkugcm_restart"
character (256) :: puma_status         = "pkugcm_status"
character (256) :: efficiency_dat      = "efficiency.dat"
character (256) :: puma_sp_init        = "pkugcm_sp_init"

! *****************************************************************
! * For multiruns the instance number is appended to the filename *
! * e.g.: puma_namelist_1 puma_diag_1 etc. for instance # 1       *
! *****************************************************************

! XW (2017-04-07): change parameters to standard F90 format below
integer, parameter :: MAXLEV = 100            ! Maximum level dimension
integer, parameter :: NROOT  =   0            ! Master node

real, parameter :: PI = 3.141592653589793D0   ! Pi
real, parameter :: TWOPI = PI + PI            ! 2 Pi

real, parameter :: AKAP_EARTH   = 0.286       ! Kappa Earth (Poisson constant R/Cp)
real, parameter :: ALR_EARTH    = 0.0065      ! Lapse rate Earth
real, parameter :: GA_EARTH     = 9.80665     ! Gravity Earth (mean on NN)
real, parameter :: GASCON_EARTH = 287.0       ! Gas constant for dry air on Earth
real, parameter :: PSURF_EARTH  = 101100.0    ! Mean Surface pressure [Pa] on Earth
                                              ! Trenberth 1981, J. Geoph. Res., Vol.86, 5238-5246
real, parameter :: PLARAD_EARTH = 6371220.0   ! Earth radius

real, parameter :: PNU          = 0.02        ! Time filter
real, parameter :: PNU21        = 1.0-2.0*PNU ! Time filter 2

! *****************************************************************
! * EZ: Factor to multiply the spherical harmonic Y_(1,0) to get  *
! * the non-dimensional planetary vorticity 2 sin(phi). In PUMA   *
! * Y_(1,0) = sqrt(3/2)*sin(phi) (normalization factor 1/sqrt(2)).*
! * The time scale must be given by Tscale = 1/Omega              *   
! *****************************************************************

parameter(EZ     = 1.632993161855452D0) ! ez = 1 / sqrt(3/8)


! **************************************************************
! * Planetary parameters & Scales                              *
! * -----------------------------                              *
! * The Puma model is formulated in non-dimensional form with  * 
! * the planetary radius as length scale and the reciprocal of * 
! * the planetary rotation rate as time scale. The temperature * 
! * scale is given by the geopotential scale divided by the    * 
! * gas constant.                                              *  
! * For the time scale the length of the siderial day is used  *
! * as basic unit                                              *
! * The parameters are initialized for Earth settings. They    *
! * may be modified by the namelist file <puma_namelist>       *
! *                                                            *
! * The scales are derived internal quantities                 *
! **************************************************************

real :: omega        = TWOPI / 86164.0 ! Default scaling 
real :: sid_day      = 86164.0       ! Length of sideral day [sec] on Earth
real :: sol_day      = 86400.0       ! Length of solar   day [sec] on Earth
real :: plarad       = PLARAD_EARTH  ! Planetary radius [m] on Earth
real :: gascon       = GASCON_EARTH  ! Dry air gas consant [J/K kg] on Earth 
real :: akap         = AKAP_EARTH    ! Kappa [] on Earth
real :: alr          = ALR_EARTH     ! average lapse rate [K/km] on Earth
real :: ga           = GA_EARTH      ! Gravity [m/sec*sec] on Earth
real :: psurf        = PSURF_EARTH   ! Mean surface pressure for EARTH [Pa] 

real :: wwt          = 13713.4258    ! time scale 有些地方替换掉了 ww_time
real :: ww_time      = 0.0           ! time scale [sec] (day / 2 Pi)
real :: ww_scale     = 0.0           ! reciprocal of time scale [1/sec]
real :: cv           = 0.0           ! velocity scale [m/sec] on Earth
real :: ct           = 0.0           ! temperature scale [K] on Earth  

! **************************
! * Global Integer Scalars *
! **************************

logical :: lrestart =  .false. ! Existing "puma_restart" sets to .true.
logical :: lselect  =  .false. ! true: disable some zonal waves
logical :: lspecsel =  .false. ! true: disable some spectral modes

integer :: kick     =  1 ! kick > 0 initializes eddy generation
integer :: nafter   =  0 ! write data interval 0: controlled by nwpd
integer :: nwpd     =  1 ! number of writes per day
integer :: ncoeff   =  0 ! number of modes to print
integer :: ndel     =  6 ! ndel
integer :: ndiag    = 12 ! write diagnostics interval
integer :: newsr    =  0 ! 1: recalculate sr1 and sr2 after restart
!integer :: ngui     =  0 ! activate Graphical User Interface	!XW(Mar/25/2017) 0=off
integer :: nkits    =  3 ! number of initial timesteps
integer :: nlevt    =  9 ! tropospheric levels (set_vertical_grid)
integer :: noutput  =  2 ! global switch for output: 0=off; 1=on(puma format); 2=GrADS(pku format)
integer :: nwspini  =  1 ! write sp_init after initialization
integer :: nrun     =  0 ! if (nstop == 0) nstop = nstep + nrun
integer :: nstep1   =  0 ! start step (for cpu statistics)
integer :: nstep    = -1 ! current timestep step 0: 01-Jan-0001  00:00
integer :: nstop    =  0 ! finishing timestep
integer :: ntspd    =  0 ! one day = ntspd timesteps
integer :: mpstep   =  0 ! minutes per step 0 = automatic
integer :: ncu      =  0 ! check unit (debug output)
integer :: nwrioro  =  1 ! controls output of orography
integer :: nextout  =  0 ! 1: extended output (entropy production)
integer :: nruido   =  0 ! 1: global constant, temporal noise
!                          2: spatio-temporal noise
!                          3: spatio-temporal equator symmetric
integer :: nseedlen =  0 ! length of random seed (set by lib call)
integer :: nmonths  =  0 ! Simulation time (1 month =  30 days)
integer :: nyears   =  1 ! simulation time (1 year  = 360 days)
integer :: nsponge  =  0 ! 1: Create sponge layer
integer :: nhelsua  =  0 ! 1: Set up Held & Suarez T_R field
!                             instead of original PUMA T_R field
!                          2: Set up Held & Suarez T_R field
!                             instead of original PUMA T_R field
!                             AND use latitudinally varying
!                             heating timescale in PUMA (H&Z(94)),
!                             irrelevant for PumaPreProcessor (ppp)
!                          3: Use latitudinally varying
!                             heating timescale in PUMA (H&Z(94)),
!                             irrelevant for PumaPreProcessor (ppp)
integer  :: ndiagp  = 0 ! 0/1 switch for grid point diabatic heating 
integer  :: nconv   = 0 ! 0/1 switch for convecive heating
integer  :: nvg     = 0 ! type of vertical grid
                        ! 0 = linear
                        ! 1 = Scinocca & Haynes
                        ! 2 = Polvani & Kushner
integer  :: nenergy = 0 ! energy diagnostics (on/off 1/0)
integer  :: nentropy= 0 ! entropy diagnostics (on/off 1/0)
integer  :: ndheat  = 0 ! energy recycling (on/off 1/0)

integer  :: nradcv = 0  ! use two restoration fields


! ***********************
! * Global Real Scalars *
! ***********************

real :: alpha  =     1.0         ! Williams filter factor
real :: alrs  =      0.0         ! stratospheric lapse rate [K/m]
real :: delt                     ! normalized timestep
real :: delt2                    ! 2 * delt
real :: dtep   =    60.0         ! delta T equator <-> pole  [K]
real :: dtns   =   -70.0         ! delta T   north <-> south [K]
real :: dtrop  = 12000.0         ! Tropopause height [m]
real :: dttrp  =     2.0         ! Tropopause smoothing [K]
real :: dtzz   =    10.0         ! delta(Theta)/H additional lapserate in
                                 ! Held & Suarez T_R field
real :: orofac =    1.0          ! factor to scale the orograpy
real :: plavor =    EZ           ! planetary vorticity
real :: psmean = PSURF_EARTH     ! Mean of Ps on Earth
real :: rotspd =     1.0         ! rotation speed 1.0 = normal Earth rotation
real :: sigmax =  6.0e-7         ! sigma for top half level
real :: spstep =  0.0            ! seconds per step 0 = automatic
real :: diffts = 21600.0         ! diffusion time scale [sec]
real :: tac    =   360.0         ! length of annual cycle [days] (0 = no cycle)
real :: pac    =     0.0         ! phase of the annual cycle [days]
real :: tgr    =   288.0         ! Ground Temperature in mean profile [K]
real :: dvdiff =     0.0         ! vertical diffusion coefficient [m2/s]
!                                ! dvdiff =0. means no vertical diffusion
real :: disp   =     0.0         ! noise dispersion
real :: tauta  =    40.0         ! heating timescale far from surface
real :: tauts  =     4.0         ! heating timescale close to surface
real :: pspon  = 50.             ! apply sponge layer where p < pspon
!                                ! pressure [Pa]
real :: dcsponge = 0.5 / 86400.0 ! damping coefficient for sponge layer [1/sec]

! **************************
! * Global Spectral Arrays *
! **************************

real, allocatable :: sd(:,:)     ! Spectral Divergence
!real, allocatable :: sdd(:,:)    ! Difference between instances         ! 干掉之 
real, allocatable :: st(:,:)     ! Spectral Temperature
!real, allocatable :: std(:,:)    ! Difference between instances         ! 干掉之
real, allocatable :: st1(:,:)    ! Spectral Temperature at t-1 (for NEXTOUT == 1)
real, allocatable :: st2(:,:)    ! Spectral Temperature at t-2 (for NEXTOUT == 1)
real, allocatable :: sz(:,:)     ! Spectral Vorticity
!real, allocatable :: szd(:,:)    ! Difference between instances         ! 干掉之
real, allocatable :: sp(:)       ! Spectral Pressure (ln Ps)
!real, allocatable :: spd(:)      ! Difference between instances         ! 干掉之
real, allocatable :: sq(:,:)     ! For compatibility with PlaSim
real, allocatable :: sp1(:)      ! Spectral Pressure at t-1 (for NEXTOUT == 1)
real, allocatable :: sp2(:)      ! Spectral Pressure at t-2 (for NEXTOUT == 1)
real, allocatable :: so(:)       ! Spectral Orography
real, allocatable :: sr1(:,:)    ! Spectral Restoration Temperature
real, allocatable :: sr2(:,:)    ! Spectral Restoration Temperature

real, allocatable :: sdp(:,:)    ! Spectral Divergence  Partial
real, allocatable :: stp(:,:)    ! Spectral Temperature Partial
real, allocatable :: szp(:,:)    ! Spectral Vorticity   Partial
real, allocatable :: spp(:)      ! Spectral Pressure    Partial
real, allocatable :: sop(:)      ! Spectral Orography   Partial
real, allocatable :: srp1(:,:)   ! Spectral Restoration Partial
real, allocatable :: srp2(:,:)   ! Spectral Restoration Partial

real, allocatable :: sdt(:,:)    ! Spectral Divergence  Tendency
real, allocatable :: stt(:,:)    ! Spectral Temperature Tendency
real, allocatable :: szt(:,:)    ! Spectral Vorticity   Tendency
real, allocatable :: spt(:)      ! Spectral Pressure    Tendency

real, allocatable :: sdm(:,:)    ! Spectral Divergence  Minus
real, allocatable :: stm(:,:)    ! Spectral Temperature Minus
real, allocatable :: szm(:,:)    ! Spectral Vorticity   Minus
real, allocatable :: spm(:)      ! Spectral Pressure    Minus

real, allocatable :: sak(:)      ! Hyper diffusion
real, allocatable :: srcn(:)     ! 1.0 / (n * (n+1))
real, allocatable :: span(:)     ! Pressure for diagnostics
real, allocatable :: spnorm(:)   ! Factors for output normalization

integer, allocatable :: nindex(:)  ! Holds wavenumber
integer, allocatable :: nscatsp(:) ! Used for reduce_scatter op
integer, allocatable :: nselzw(:)  ! Enable/disable selected zonal waves
integer, allocatable :: nselsp(:)  ! Enable/disable slected spectral modes

! ***************************
! * Global Gridpoint Arrays *
! ***************************

real, allocatable :: gd(:,:)     ! Divergence
real, allocatable :: gt(:,:)     ! Temperature
real, allocatable :: gz(:,:)     ! Vorticity
real, allocatable :: gu(:,:)     ! u * cos(phi)
real, allocatable :: gv(:,:)     ! v * cos(phi)
real, allocatable :: gp(:)       ! Ln(Ps)
real, allocatable :: gq(:,:)     ! For compatibilty with PlaSim
real, allocatable :: gfu(:,:)    ! Term Fu in Primitive Equations
real, allocatable :: gfv(:,:)    ! Term Fv in Primitive Equations
real, allocatable :: gut(:,:)    ! Term u * T
real, allocatable :: gvt(:,:)    ! Term v * T
real, allocatable :: gke(:,:)    ! Kinetic energy u * u + v * v
real, allocatable :: gpj(:)      ! d(Ln(Ps)) / d(mu)
real, allocatable :: rcsq(:)     ! 1 / cos2(phi)
real, allocatable :: ruido(:,:,:)! noise (nlon,nlat,nlev)
real, allocatable :: ruidop(:,:) ! noise partial (nhor,nlev)
real, allocatable :: gtdamp(:,:) ! 3D reciprocal damping times [1/sec] 
                                 ! for relaxation in grid point space 
                                 ! for radiative restoration temperature 
                                 ! (e.g. for Held&Suarez)
real, allocatable :: gr1(:,:)    ! constant radiative restoration time scale
real, allocatable :: gr2(:,:)    ! variable radiative restoration time scale
real, allocatable :: gtdampc(:,:)! the same as gtdamp, but for convective  
                                 ! restoration temperature
real, allocatable :: gr1c(:,:)   ! constant convective restoration time scale
real, allocatable :: gr2c(:,:)   ! variable convective restoration time scale

! *********************
! * Diagnostic Arrays *
! *********************

integer, allocatable :: ndil(:)  ! Set diagnostics level

real, allocatable :: csu(:,:)    ! Cross section u [m/s]
real, allocatable :: csv(:,:)    ! Cross section v [m/s]
real, allocatable :: cst(:,:)    ! Cross section T [Celsius]

real,allocatable :: denergy(:,:) ! energy diagnostics
real,allocatable :: dentropy(:,:)! entropy diagnostics

! *******************
! * Latitude Arrays *
! *******************

character (3),allocatable :: chlat(:) ! label for latitudes
real (kind=8),allocatable :: sid(:)   ! sin(phi)
real (kind=8),allocatable :: cid(:)   ! cos(phi)   !XW(2017/4/15)
real (kind=8),allocatable :: gwd(:)   ! Gaussian weight (phi)
real, allocatable :: csq(:)           ! cos2(phi)
real, allocatable :: rcs(:)           ! 1/cos(phi)

! ****************
! * Level Arrays *
! ****************

real, allocatable :: t0(:)            ! reference temperature
real, allocatable :: t0d(:)           ! vertical t0 gradient
real              :: taur(MAXLEV)     ! tau R [sec]
real              :: tauf(MAXLEV)     ! tau F [sec]
real, allocatable :: damp(:)          ! 1.0 / (2 Pi * taur)
real, allocatable :: fric(:)          ! 1.0 / (2 Pi * tauf)

real, allocatable :: bm1(:,:,:)
real, allocatable :: dsigma(:)
real, allocatable :: rdsig(:)         ! 每一层的sigma厚度的一半
real, allocatable :: sigma(:)         ! full level sigma
real, allocatable :: sigmh(:)         ! half level sigma
real, allocatable :: tkp(:)
real, allocatable :: c(:,:)           ! XW: C matrix, used in calcgp
real, allocatable :: xlphi(:,:)       ! matrix Lphi (g)
real, allocatable :: xlt(:,:)         ! matrix LT (tau)

! ******************
! * Parallel Stuff *
! ******************

!integer :: myworld = 0                   ! MPI variable
!integer :: mpinfo  = 0                   ! MPI variable
integer :: mypid   = 0                   ! My Process Id
!character(80), allocatable :: ympname(:) ! Processor name

real(kind=8)    :: tmstart = 0.0         ! CPU time at start
real(kind=8)    :: tmstop  = 0.0         ! CPU time at stop

! **********************
! * Multirun variables *
! **********************

!integer :: mrworld =  0   ! MPI communication
!integer :: mrinfo  =  0   ! MPI info
!integer :: mrpid   = -1   ! MPI instance id
!integer :: mrnum   =  0   ! MPI number of instances
!integer :: mintru  =  0   ! Lowest resolution of all instances
!integer :: mrdim   =  0   ! Exchange dimension  (min. NRSP)
!integer :: nsync   =  0   ! Synchronization on or off    ! 特殊的同步要求 _namelist required 早晚干掉它
!integer, allocatable :: mrtru(:) ! Truncations of members

!real    :: syncstr  =  0.0 ! Coupling strength (0 .. 1)  ! 特殊的同步要求 _namelist required 早晚干掉它
!real    :: syncsecs =  0.0 ! Coupling time [sec]


! ***************
! * Random seed *
! ***************

integer              :: seed(8) = 0 ! settable in namelist
integer, allocatable :: meed(:)     ! machine dependent seed
real                 :: ganext = 0.0! y part of gaussian noise

end module pumamod


! ***************** !
! * MODULE PPPMOD * !
! * 删除所有PPP部分: prepmod (module); ppp_def_ini; ppp_def_real; ppp_read_i; ppp_read_r
! ***************** !


! *********************
! * PROGRAM PUMA_MAIN *
! *********************

program puma_main
use pumamod

! ***********
! * History *
! ***********

! 1972 - W. Bourke:
!        An efficient one-level primitive equation spectral model
!        Mon. Weath. Rev., 100, pp. 683-689

! 1975 - B.J. Hoskins and A.J. Simmons: 
!        A multi-layer spectral model and the semi-implicit method
!        Qart. J. R. Met. Soc., 101, pp. 637-655

! 1993 - I.N. James and J.P. Dodd:
!        A Simplified Global Circulation Model
!        Users' Manual, Dept. of Meteorology, University of Reading

! 1998 - Klaus Fraedrich, Edilbert Kirk, Frank Lunkeit
!        Portable University Model of the Atmosphere
!        DKRZ Technical Report No. 16

! 2009 - PUMA Version 16.0
!        http://www.mi.uni-hamburg.de/puma

! ******************
! * Recent Changes *
! ******************

! 10-Jun-2002 - Puma Workshop - Documentation of subroutine SPECTRAL
! 04-Jul-2002 - Frank Lunkeit - Annual cycle
! 08-Jul-2002 - Edilbert Kirk - Factor for rotation speed
! 25-Sep-2002 - Puma Workshop - Documentation of subroutine CALCGP
! 11-Nov-2002 - Edilbert Kirk - Add Orography to output file
! 26-Feb-2003 - Edilbert Kirk - Read preprocessed initial file
! 07-Sep-2004 - Edilbert Kirk - Graphical User Interface
! 23-Aug-2006 - Torben Kunz   - Held & Suarez forcing
! 23-Aug-2006 - Torben Kunz   - new spacing schemes of sigma levels
! 23-Aug-2006 - Edilbert Kirk - individual selection of zonal waves
! 23-Aug-2006 - Edilbert Kirk - optimized Legendre trasnformation module
! 19-Feb-2007 - Edilbert Kirk - new flexible restart I/O
! 15-Sep-2009 - Edilbert Kirk - static arrays replaced by allocatable
! 15-Sep-2009 - Frank Lunkeit - diagnostics for entropy production
! 27-Sep-2010 - Edilbert Kirk - cleaned up ruido routines
! 15-Mar-2017 - Xinyu Wen     - Cleanup & reorganization for PKU OpenMP version

! Let's Rock Here !
! "pass" implies 经过完整检查和简化，并用implicit none最后封装
! 直接移除开始阶段的mpstart & setfilenames 和最后阶段的mpstop


open(nud,file=puma_diag)   ! pass (simple)  ! opendiag简化而来
call read_resolution       ! pass (simple)  ! 从命令行读取2个参数: NLAT, NLEV
call resolution            ! pass (simple)  ! 设置SP和GP的分辨率参数
call allocate_arrays       ! pass (simple)  ! 分配SP和GP所有动态数组的内存空间
call prolog                ! pass (complex) ! 初始化 ......
call master                ! pass (complex) ! Key: gridpoint & spectral
call epilog                ! pass (simple)  ! 终止化 关闭_output; 写restart; 显示时间信息

print *, "STOP Normally!!!"
end program puma_main


! ==========================
! SUBROUTINE READ_RESOLUTION
! ==========================
subroutine read_resolution
   use pumamod
   implicit none

   character (80) :: ylat
   character (80) :: ylev
   call get_command_argument(1,ylat)
   call get_command_argument(2,ylev)
   read(ylat,*) nlat
   read(ylev,*) nlev
end subroutine


! =====================
! SUBROUTINE RESOLUTION
! =====================
subroutine resolution
   use pumamod
   implicit none

   ! 垂直方向
   nlem = nlev - 1
   nlep = nlev + 1
   nlsq = nlev * nlev

   ! 水平方向
   nlon = nlat + nlat                  ! Longitudes
   nlah = nlat / 2
   nlpp = nlat
   nhpp = nlat / 2                     ! N-S数据对
   nhor = nlon * nlat
   nugp = nlon * nlat
   npgp = nlon * nlat / 2

   ! 谱空间
   ntru = (nlon - 1) / 3               ! 截断波数，比如：21 for T21, 42 for T42
   ntp1 = ntru + 1
   nzom = ntp1 * 2
   nrsp = (ntru + 1) * (ntru + 2)      ! 所有SP谱空间系数，Pn重复一遍有ntru+1个点 (22x23=506 for T21)
   ncsp = nrsp / 2                     !     ^----------- 的右一半，左一半与之共轭
   !nspp = (nrsp + npro - 1) / npro    ! 分配到每个processor上的SP谱空间mode
   !nesp = nspp * npro                 ! 扩展的所有SP谱空间系数
   !nesp = nesp + 3 - mod(nesp-1,4)    ! 506+3-mod(506-1,4) = 506+3-1 = 508
                                       ! 在前面此值为506, why? I guess it should be "nrsp" too!
   nspp = nrsp
   nesp = nrsp
end subroutine resolution


! ******************************
! * SUBROUTINE ALLOCATE_ARRAYS *
! ******************************
subroutine allocate_arrays
use pumamod
implicit none

! SP全模态 (nesp)
allocate(st(nesp,nlev))   ; st(:,:)  = 0.0 ! Spectral Temperature
allocate(sd(nesp,nlev))   ; sd(:,:)  = 0.0 ! Spectral Divergence
allocate(sz(nesp,nlev))   ; sz(:,:)  = 0.0 ! Spectral Vorticity
allocate(sp(nesp))        ; sp(:)    = 0.0 ! Spectral Pressure (ln Ps)
allocate(so(nesp))        ; so(:)    = 0.0 ! Spectral Orography
allocate(sr1(nesp,nlev))  ; sr1(:,:) = 0.0 ! Spectral Restoration Temperature
allocate(sr2(nesp,nlev))  ; sr2(:,:) = 0.0 ! Spectral Restoration Temperature

! 每个CPU分配的SP部分模态 (nspp)
allocate(stp(nspp,nlev))  ; stp(:,:) = 0.0 ! Spectral Temperature Partial
allocate(sdp(nspp,nlev))  ; sdp(:,:) = 0.0 ! Spectral Divergence  Partial
allocate(szp(nspp,nlev))  ; szp(:,:) = 0.0 ! Spectral Vorticity   Partial
allocate(spp(nspp))       ; spp(:)   = 0.0 ! Spectral Pressure    Partial
allocate(sop(nspp))       ; sop(:)   = 0.0 ! Spectral Orography   Partial
allocate(srp1(nspp,nlev)) ; srp1(:,:)= 0.0 ! Spectral Restoration Partial
allocate(srp2(nspp,nlev)) ; srp2(:,:)= 0.0 ! Spectral Restoration Partial

allocate(stt(nspp,nlev))  ; stt(:,:) = 0.0 ! Spectral Temperature Tendency
allocate(sdt(nspp,nlev))  ; sdt(:,:) = 0.0 ! Spectral Divergence  Tendency
allocate(szt(nspp,nlev))  ; szt(:,:) = 0.0 ! Spectral Vorticity   Tendency
allocate(spt(nspp))       ; spt(:)   = 0.0 ! Spectral Pressure    Tendency

allocate(stm(nspp,nlev))  ; stm(:,:) = 0.0 ! Spectral Temperature Minus
allocate(sdm(nspp,nlev))  ; sdm(:,:) = 0.0 ! Spectral Divergence  Minus
allocate(szm(nspp,nlev))  ; szm(:,:) = 0.0 ! Spectral Vorticity   Minus
allocate(spm(nspp))       ; spm(:)   = 0.0 ! Spectral Pressure    Minus

allocate(sak(nesp))       ; sak(:)   = 0.0 ! Hyper diffusion
allocate(srcn(nesp))      ; srcn(:)  = 0.0 ! 1.0 / (n * (n+1))
allocate(span(nesp))      ; span(:)  = 0.0 ! Pressure for diagnostics
allocate(spnorm(nesp))    ; spnorm(:)= 0.0 ! Factors for output normalization

allocate(nindex(nesp))    ; nindex(:)  = ntru ! Holds wavenumber
allocate(nscatsp(npro))   ; nscatsp(:) = nspp ! Used for reduce_scatter op    ! WHAT???
allocate(nselzw(0:ntru))  ; nselzw(:)  =    1 ! Enable selected zonal waves
allocate(nselsp(ncsp))    ; nselsp(:)  =    1 ! Enable slected spectral modes

! 格点空间 (nhor=NLONxNLAT)
allocate(gt(nhor,nlev))   ; gt(:,:)  = 0.0 ! Temperature
allocate(gd(nhor,nlev))   ; gd(:,:)  = 0.0 ! Divergence
allocate(gz(nhor,nlev))   ; gz(:,:)  = 0.0 ! Vorticity
allocate(gu(nhor,nlev))   ; gu(:,:)  = 0.0 ! u * cos(phi)
allocate(gv(nhor,nlev))   ; gv(:,:)  = 0.0 ! v * sin(phi)
allocate(gp(nhor))        ; gp(:)    = 0.0 ! Ln(Ps)

allocate(gfu(nhor,nlev))  ; gfu(:,:) = 0.0 ! Term Fu in Primitive Equations
allocate(gfv(nhor,nlev))  ; gfv(:,:) = 0.0 ! Term Fv in Primitive Equations
allocate(gut(nhor,nlev))  ; gut(:,:) = 0.0 ! Term u * T
allocate(gvt(nhor,nlev))  ; gvt(:,:) = 0.0 ! Term v * T
allocate(gke(nhor,nlev))  ; gke(:,:) = 0.0 ! Kinetic energy u * u + v * v
allocate(gpj(nhor))       ; gpj(:)   = 0.0 ! d(Ln(Ps)) / d(mu)

allocate(rcsq(nhor))      ; rcsq(:)  = 0.0 ! 1 / cos2(phi)

! 切片用
allocate(ndil(nlev))      ; ndil(:)  = 0
allocate(csu(nlat,nlev))  ; csu(:,:) = 0.0
allocate(csv(nlat,nlev))  ; csv(:,:) = 0.0
allocate(cst(nlat,nlev))  ; cst(:,:) = 0.0

! 与南北纬度有关的参数
allocate(chlat(nlat))     ; chlat(:) = '   '
allocate(sid(nlat))       ; sid(:)   = 0.0      ! sin(phi)
allocate(gwd(nlat))       ; gwd(:)   = 0.0      ! Gaussian weight (phi)
allocate(csq(nlat))       ; csq(:)   = 0.0      ! cos2(phi)
allocate(rcs(nlat))       ; rcs(:)   = 0.0      ! 1/cos(phi)

! 与垂直层结有关的参数
allocate(t0(nlev))        ; t0(:)     = 250.0   ! reference temperature
allocate(t0d(nlev))       ; t0d(:)    =   0.0   ! vertical t0 gradient
allocate(damp(nlev))      ; damp(:)   =   0.0   ! 1.0 / (2 Pi * taur)
allocate(fric(nlev))      ; fric(:)   =   0.0   ! 1.0 / (2 Pi * tauf)
allocate(dsigma(nlev))    ; dsigma(:) =   0.0
allocate(rdsig(nlev))     ; rdsig(:)  =   0.0
allocate(sigma(nlev))     ; sigma(:)  =   0.0
allocate(sigmh(nlev))     ; sigmh(:)  =   0.0
allocate(tkp(nlev))       ; tkp(:)    =   0.0
allocate(c(nlev,nlev))    ; c(:,:)    =   0.0
allocate(xlphi(nlev,nlev)); xlphi(:,:) = 0.0    ! matrix Lphi (g)
allocate(xlt(nlev,nlev))  ; xlt(:,:)   = 0.0    ! matrix LT (tau)
allocate(bm1(nlev,nlev,0:NTRU)) ; bm1(:,:,:)  = 0.0

!if (mrnum == 2) then
!   allocate(std(nesp,nlev))   ; std(:,:)  = 0.0   ! instance之间的T差别
!   allocate(sdd(nesp,nlev))   ; sdd(:,:)  = 0.0   ! instance之间的div差别
!   allocate(szd(nesp,nlev))   ; szd(:,:)  = 0.0   ! instance之间的vor差别
!   allocate(spd(nesp     ))   ; spd(:  )  = 0.0   ! instance之间的Ps差别
!endif
end subroutine allocate_arrays


! =================
! SUBROUTINE PROLOG
! =================
! XW(2017/4/16)
! 以下有"?"的部分
! 没经过understanding & implicit none检查
!------------------
! Calling tree:
! - restart_ini         : 读入restart file的变量名列表
! - inigau              : 计算Gauss纬度权重gwd及sin纬度sid
!   ql & qld            : Legendre and Associated Legendre functions
! - inilat              : 设置csq, rcs, chlat
! - fftini              : 我把mod_fft中的FFT基函数初始化放到这里来 只需执行一次
! - legpri              : 只打印一个表格显示lat, csq, gwd
! - readnl              : 读入_namelist并显示出来，修正基本的尺度因子
! - ppp_interface       : 读入ppp产生的输出文本并与_namelist对比看是否一致 (ALL DELETED)
!   ppp_define_ini
!   ppp_define_real
!   ppp_read_i
!   ppp_read_r
! - initpm (complex)    : 准备用于快速计算用的各种垂向数组
!   select_zonal_waves    如果T21的22个纬向波数(m)都打开(=1)，那么就关闭滤波；否则打开，有波数损失
!   select_spectral_modes 如果T21的所有经向波数(m,n)都打开，那么就关闭滤波；否则打开，有模态损失
!   set_vertical_grid     设置sigma levels (3 methods selectable)
!   sponge              : 计算模式最顶层使用Rayleigh摩擦的系数fric(NLEV)
! - initsi              : 半隐时间方案初始化
! - altlat              : csq(reg2d) to csq(alt2d)
! - initrandom          : 设随机数种子
! - initruido           : 分配ruido,ruidop数组
!   mpscgp(nuido,nuidop): 分发nuido
! - legini              : SP初始化, 计算Legendre Polynomials和其它SP空间参数
! - if(lrestart)        : 最终来到最重要的初始化分支！
!   = 1                 : --- Restart
!     read_atmos_restart: 如果使用restart，读入restartfile变量，随机数种子，并设置水球理想
!     noise             : 随机数种子，并设置水球理想
!     setzt             : 设置水球理想T回复场（sr1=年平均；sr2=季节循环）
!   = 0                 : --- Init from zero
!   initfd              : GP空间初始化(key)
! ?   read_surf         : 读入surf变量（地表高度等）
! ?   read_vargp        : 读入其它变量
!     setzt             : 设置水球理想T回复场（sr1=年平均；sr2=季节循环）
!     printprofile      : 打印垂直profile表格
! - printseed           : 打印随机种子(three sources: 1 namelist; 2 clock initialized; 3 restart file)
! - ntomin              : 时间步数 to 年月日时分
! - mastercpu_time      : 系统时间，现在只向月初对齐，以后可改为向2000年对齐
! - io_open_output      : 新输出格式(GrADS)初始化

subroutine prolog
use pumamod
implicit none

real     :: zsig(nlon*nlat)   ! 实际只用前几个数保存sigma值，只借用2D壳而已，便于Service记录
integer  :: istep,iyear,imonth,iday,ihour,imin  ! 临时保存时间信息,用于向_output写入第一个2D记录zsig（sigma信息）
integer  :: jlat              ! loop var

!if (mypid == NROOT) then
   call mastercpu_time(tmstart)   ! XW(2017/4/9): replace cpu_time(tmstart) to accurate OpenMP time cost
   write(nud, '("********************************************************")')
   write(nud, '("* Peking University General Circulation Model (PKUGCM) *")')
   write(nud, '("* based on PUMA from University of Hamburg             *")')
   write(nud, '("********************************************************")')
   write(nud, '("  NTRU =",i4,"  NLEV =",i4,"  NLON = ",i4,"   NLAT =",i4 )') NTRU,NLEV,NLON,NLAT
   write(nud, '("********************************************************")')

   call restart_ini(lrestart,puma_restart)
   !call fftini(NLON) XW(2017/4/11): I added this line incorrectly. Affect MPI.
   call inigau(NLAT,sid,gwd)
   call inilat
   call fftini(NLON)
   call legpri
   call readnl
   !call ppp_interface  完全不必要
   call initpm

print *, "LSELECT =", lselect       ! check if zonal-filter work after calling initpm
print *, "LSPECSEL = ", lspecsel    ! check if sp-truncator work after calling initpm

   call initsi
   call altlat(csq,NLAT) ! csq -> alternating grid
   !XW(Mar/25/2017) to remove GUI: if (ngui > 0) call guistart
   if (nrun == 0 .and. nstop  > 0) nrun = nstop-nstep
   if (nrun == 0) nrun = ntspd * (nyears * 360 + nmonths * 30)
   call initrandom     ! set random seed
!endif ! (mypid == NROOT)

call initruido      ! allocate ruido arrays

! scatter gridpoint arrays
! 虽然是scatter，但这里并不能写成ruidop=ruido
! 因为ruido(nlon,nlat,nlev), ruidop(nhor,nlev)两者维度rank不同
! 可见通过mpscgp只是达到了方便快捷地transform数组的作用
if (nruido > 0) ruidop = reshape(ruido,(/nhor,nlev/)) !call mpscgp(ruido,ruidop,NLEV)


! ***********************
! * broadcast & scatter *
! * 全部删除            *
!   broadcast integer
!   broadcast logical
!   broadcast real
!   broadcast integer arrays
!   broadcast real arrays
!   scatter integer arrays
! ***********************

do jlat = 1 , NLPP
   rcsq(1+(jlat-1)*NLON:jlat*NLON) = 1.0 / csq(jlat)
enddo

call legini(nlat,nlpp,nesp,nlev,plavor,sid,gwd)

if (lrestart) then
   call read_atmos_restart
   !if (mypid == NROOT) then
      if (kick > 10) call noise(kick-10)
      if (newsr > 0) call setzt
   !endif
else
   call initfd
endif

!if (mypid == NROOT) then
   call printseed ! either namelist, clock initialized or from restart file
!endif

!     broadcast spectral arrays  删除4行
!     scatter spectral arrays    
!     以下7行不能删除，否则会导致全零场的积分僵化，scatter本质是复制操作，可简化为
sdp   = sd  ! call mpscsp(sd,sdp,NLEV)
stp   = st  ! call mpscsp(st,stp,NLEV)
szp   = sz  ! call mpscsp(sz,szp,NLEV)
srp1  = sr1 ! call mpscsp(sr1,srp1,NLEV)
srp2  = sr2 ! call mpscsp(sr2,srp2,NLEV)
spp   = sp  ! call mpscsp(sp,spp,1)
sop   = so  ! call mpscsp(so,sop,1)

!
!     initialize energy and entropy diagnostics
!
if(nenergy > 0) then
 allocate(denergy(NHOR,9))
 denergy(:,:)=0.0
endif
if(nentropy > 0) then
 allocate(dentropy(NHOR,9))
 dentropy(:,:)=0.0
endif
if(ndheat > 1) then
 open(9,file=efficiency_dat,form='formatted')
endif
!
!     write first service record containing sigma coordinates
!
!if (mypid == NROOT) then
   if (noutput == 1) then
      istep = nstep
      if (istep > 0) istep = istep + nafter ! next write after restart
      open(40,file=puma_output,form='unformatted')
      call ntomin(istep,imin,ihour,iday,imonth,iyear)
      zsig(1:nlev) = sigmh(:)    ! zsig(nlon,nlat)只用了前10个数来保存level的sigma值，其它无用
      zsig(nlev+1:) = 0.0
      write(40) 333,0,iyear*10000+imonth*100+iday,0,nlon,nlat,nlev,ntru
      write(40) zsig
   endif

   if (noutput == 2) then
      call io_open_output
   end if

!endif
end subroutine prolog


!===========================!
! SUBROUTINE MASTER         !
! Calling tree              !
! - makebm    P: 计算B矩阵  !
! ? minvers
! ?   ludcmp
! ?   lubksb
! - gridpoint P: GP空间计算 !
! - spectral   : SP空间计算 !
! - outsp      : Output SP  !
! - outgp      : Output GP  !
! - diag
! ?   prisp
! ?      wrspam
! ?   xsect
! ?      wrzs
! ?   energy
! ?      rmssp
! ?      powerspec
! ?      powerprint
! - checkunit
!===========================!
subroutine master
use pumamod

! ***************************
! * short initial timesteps *
! ***************************


ikits = nkits
do jkits = 1 , ikits
   delt  = spstep * ww_scale / (2**nkits)
   delt2 = delt + delt
   call gridpoint
   call makebm
   call spectral
   nkits = nkits - 1
enddo

delt  = spstep * ww_scale
delt2 = delt + delt
call makebm

nstep1 = nstep ! remember 1.st timestep

! Xinyu added, for explicitly tell compiler
! the leapfrog block in subroutine spectral can be paralled
mloop = .TRUE.

do jstep = 1 , nrun

   nstep = nstep + 1

!  ************************************************************
!  * calculation of non-linear quantities in grid point space *
!  ************************************************************

   call gridpoint
   !if (mypid == NROOT) then
      !if (mod(nstep,nafter)==0 .and. noutput==1) call outsp
      !if (mod(nstep,ndiag )==0) call diag
      !if (ncu > 0) call checkunit
      if (mod(nstep,ndiag )==0)  print "(i10,2f10.4)", nstep, maxval(gu), minval(gu)
   !endif

!  ******************************
!  * adiabatic part of timestep *
!  ******************************

   call spectral
   !if (mod(nstep,nafter)==0 .and. noutput==1) call outgp

   !***********************
   ! XW: pku output format
   !***********************
   if (mod(nstep,nafter)==0 .and. noutput==2) call io_write_output
enddo

end subroutine master


!     =================
!     SUBROUTINE EPILOG
!     =================

      subroutine epilog
      use pumamod
      implicit none

      real    :: tmrun, zspy, zypd

      close(40)   ! close default "_output" file

      ! XW(2017/4/12): close new output files
      if (noutput==2) call io_close_output 

!     write restart file

      !if (mypid == NROOT) then
         call restart_prepare(puma_status)
         sp(1) = psmean ! save psmean
         call put_restart_integer('nstep'   ,nstep   )
         call put_restart_integer('nlat'    ,NLAT    )
         call put_restart_integer('nlon'    ,NLON    )
         call put_restart_integer('nlev'    ,NLEV    )
         call put_restart_integer('nrsp'    ,NRSP    )

!        Save current random number generator seed

         call random_seed(get=meed)
         call put_restart_array('seed',meed,nseedlen,nseedlen,1)
         call put_restart_array('ganext',ganext,1,1,1)

         call put_restart_array('sz' ,sz ,NRSP,NESP,NLEV)
         call put_restart_array('sd' ,sd ,NRSP,NESP,NLEV)
         call put_restart_array('st' ,st ,NRSP,NESP,NLEV)
         call put_restart_array('sr1',sr1,NRSP,NESP,NLEV)
         call put_restart_array('sr2',sr2,NRSP,NESP,NLEV)
         call put_restart_array('sp' ,sp ,NRSP,NESP,   1)
         call put_restart_array('so' ,so ,NRSP,NESP,   1)
         if (nruido > 0) then
            call put_restart_array('ruido',ruido,nugp,nugp,nlev)
         endif
      !endif

      ! 以下写sp和gp数组的函数，貌似是mpi_stub里的，实际都可以直接跨过
      ! 今后把mpi的这些函数都去掉，直接调用它们调用的mod_restart里的put_restart_array

!     write sp arrays

      call put_restart_array('szm',szm,NSPP,NSPP,NLEV)            !mpputsp('szm',szm,NSPP,NLEV)
      call put_restart_array('sdm',sdm,NSPP,NSPP,NLEV)            !mpputsp('sdm',sdm,NSPP,NLEV)
      call put_restart_array('stm',stm,NSPP,NSPP,NLEV)            !mpputsp('stm',stm,NSPP,NLEV)
      call put_restart_array('spm',spm,NSPP,NSPP,   1)            !mpputsp('spm',spm,NSPP,   1)

!     write gridpoint arrays

      if (allocated(gr1)) then
         call put_restart_array('gr1',gr1,nhor,nhor,nlev)         !mpputgp('gr1',gr1,nhor,nlev)
      endif
      if (allocated(gr2)) then
         call put_restart_array('gr2',gr2,nhor,nhor,nlev)         !mpputgp('gr2',gr2,nhor,nlev)
      endif
      if (allocated(gtdamp)) then
         call put_restart_array('gtdamp',gtdamp,nhor,nhor,nlev)   !mpputgp('gtdamp',gtdamp,nhor,nlev)
      endif

      if (allocated(gr1c)) then
         call put_restart_array('gr1c',gr1c,nhor,nhor,nlev)       !mpputgp('gr1c',gr1c,nhor,nlev)
      endif
      if (allocated(gr2c)) then
         call put_restart_array('gr2c',gr2c,nhor,nhor,nlev)       !mpputgp('gr2c',gr2c,nhor,nlev)
      endif
      if (allocated(gtdampc)) then
         call put_restart_array('gtdampc',gtdampc,nhor,nhor,nlev) !mpputgp('gtdampc',gtdampc,nhor,nlev)
      endif

      !if (mypid == NROOT) then 
!        Get resource stats from function resources in file pumax.c

         !ires = 1
         !XW/Mar-23-2017: this line need to call a function in pumax.c
         !ires = nresources(zut,zst,imem,ipr,ipf,isw,idr,idw)

         call mastercpu_time(tmstop)   ! XW(2017/4/9): replace cpu_time(tmstop)
         tmrun = tmstop - tmstart

         if (nstep > nstep1) then 
            zspy = tmrun * 360.0 * real(ntspd) / (nstep - nstep1) ! sec / siy
            zypd = (24.0 * 3600.0 / zspy)                         ! siy / day
            write(nud, '(/,"****************************************")')
            write(nud, '("* Total CPU time      : ", f10.3," sec *")') tmrun
            write(nud, '("****************************************")')
            if (zspy < 600.0) then
               write(nud,'("* Seconds per sim year: ",i6,9x,"*")') nint(zspy)
            else if (zspy < 900000.0) then
               write(nud,'("* Minutes per sim year  ",i6,9x,"*")') nint(zspy/60.0)
            else
               write(nud,'("* Days per sim year:    ",i6,5x,"*")') nint(zspy/sol_day)
            endif
            write(nud,'("* Sim years per day   :",i7,9x,"* <-- Running Speed Index")') nint(zypd)
            write(nud,'("****************************************")')
         endif
      !endif
      end subroutine epilog


!     =============================
!     SUBROUTINE READ_ATMOS_RESTART
!     =============================

      subroutine read_atmos_restart
      use pumamod
      implicit none

      integer :: k = 0
      integer :: ktmp

!     read scalars and full spectral arrays

      !if (mypid == NROOT) then
         call get_restart_integer('nstep',nstep)
         call get_restart_array('seed',meed,nseedlen,nseedlen,1)
         call get_restart_array('ganext',ganext,1,1,1)
         call get_restart_array('sz' ,sz ,NRSP,NESP,NLEV)
         call get_restart_array('sd' ,sd ,NRSP,NESP,NLEV)
         call get_restart_array('st' ,st ,NRSP,NESP,NLEV)
         call get_restart_array('sr1',sr1,NRSP,NESP,NLEV)
         call get_restart_array('sr2',sr2,NRSP,NESP,NLEV)
         call get_restart_array('sp' ,sp ,NRSP,NESP,   1)
         call get_restart_array('so' ,so ,NRSP,NESP,   1)
         if (nruido > 0) then
            call get_restart_array('ruido',ruido,nugp,nugp,nlev)
         endif
         psmean = sp(1)
         sp(1)  = 0.0
         call random_seed(put=meed)
      !endif

      !call mpbci(nstep)     ! broadcast current timestep
      !call mpbcr(psmean)    ! broadcast mean surface pressure

!     read and scatter spectral arrays

      call get_restart_array('szm',szm,NSPP,NSPP,NLEV)   !mpgetsp('szm',szm,NSPP,NLEV)
      call get_restart_array('sdm',sdm,NSPP,NSPP,NLEV)   !mpgetsp('sdm',sdm,NSPP,NLEV)
      call get_restart_array('stm',stm,NSPP,NSPP,NLEV)   !mpgetsp('stm',stm,NSPP,NLEV)
      call get_restart_array('spm',spm,NSPP,NSPP,   1)   !mpgetsp('spm',spm,NSPP,   1)

!     allocate, read and scatter gridpoint arrays

      call varseek('gr1',ktmp)
      !call mpbci(ktmp)
      if (ktmp > 0) then 
         allocate(gr1(nhor,nlev))
         call get_restart_array('gr1',gr1,nhor,nhor,nlev)         !mpgetgp('gr1',gr1,nhor,nlev)
      endif

      call varseek('gr2',ktmp)
      !call mpbci(ktmp)
      if (ktmp > 0) then 
         allocate(gr2(nhor,nlev))
         call get_restart_array('gr2',gr2,nhor,nhor,nlev)         !mpgetgp('gr2',gr2,nhor,nlev)
      endif

      call varseek('gtdamp',ktmp)
      !call mpbci(ktmp)
      if (ktmp > 0) then 
         allocate(gtdamp(nhor,nlev))
         call get_restart_array('gtdamp',gtdamp,nhor,nhor,nlev)   !mpgetgp('gtdamp',gtdamp,nhor,nlev)
      endif

      call varseek('gr1c',ktmp)
      !call mpbci(ktmp)
      if (ktmp > 0) then 
         allocate(gr1c(nhor,nlev))
         call get_restart_array('gr1c',gr1c,nhor,nhor,nlev)       !mpgetgp('gr1c',gr1c,nhor,nlev)
      endif

      call varseek('gr2c',ktmp)
      !call mpbci(ktmp)
      if (ktmp > 0) then 
         allocate(gr2c(nhor,nlev))
         call get_restart_array('gr2c',gr2c,nhor,nhor,nlev)       !mpgetgp('gr2c',gr2c,nhor,nlev)
      endif

      call varseek('gtdampc',ktmp)
      !call mpbci(ktmp)
      if (ktmp > 0) then 
         allocate(gtdampc(nhor,nlev))
         call get_restart_array('gtdampc',gtdampc,nhor,nhor,nlev) !mpgetgp('gtdampc',gtdampc,nhor,nlev)
      endif

      end subroutine read_atmos_restart


!     =================
!     SUBROUTINE INITFD
!     =================

      subroutine initfd
      use pumamod
      implicit none

      ! local
      integer  :: iread1, iread2, iread3, iread4
      integer  :: iread121, iread122, iread123, iread124, iread125, iread126

      if (nkits < 1) nkits = 1

!     Look for start data and read them if there

      call read_surf(129,so,    1,iread1)
      call read_surf(134,sp,    1,iread2)
      call read_surf(121,sr1,NLEV,iread3)
      call read_surf(122,sr2,NLEV,iread4)
      call read_vargp(123,NLEV,iread123)
      !if (mypid == NROOT .and. iread123 == 0) then
      if (iread123 == 0) then
         if (nhelsua > 1) then
            write(nud,*) "*** ERROR no *_surf_0123.sra file for Held&Suarez"
            stop
         endif
      endif
   
      if (ndiagp > 0) then  
         call read_vargp(121,NLEV,iread121)
         call read_vargp(122,NLEV,iread122)
         if (.not. allocated(gtdamp)) then
            call read_vargp(123,NLEV,iread123)
         endif
         !if (mypid == NROOT) then
            if (iread121==0 .or. iread122==0 .or. iread123==0) then
               write(nud,*) "*** ERROR not all fields (121,122,123) for grid point heating found"
               stop
            endif
         !endif
      endif

      if (nconv > 0) then
         call read_vargp(124,NLEV,iread124)
         call read_vargp(125,NLEV,iread125)
         call read_vargp(126,NLEV,iread126)
         !if (mypid == NROOT) then
            if (iread124==0 .or. iread125==0 .or. iread126==0) then
               write(nud,*) "*** ERROR not all fields (124,125,126) for convective heating found"
               stop
            endif
         !endif
      endif

      !if (mypid == NROOT) then
         if (iread1==0 .or. iread2==0 .or. iread3==0 .or. iread4==0) then
            call setzt ! setup for aqua-planet
         else
            psmean = psurf * exp(spnorm(1) * sp(1)) 
            sp(1)  = 0.0
            so(:) = so(:) / (cv * cv) ! descale from [m2/s2]
            sr1(:,:) = sr1(:,:) / ct  ! descale from [K]
            sr2(:,:) = sr2(:,:) / ct  ! descale from [K]
            sr1(1,:) = sr1(1,:) - t0(:) * sqrt(2.0) ! subtract profile
            write(nud,'(a,f8.2,a)') ' Mean of Ps = ',0.01*psmean, '[hPa]'
         endif
      !endif

!     Add initial noise if wanted


      !if (mypid == NROOT) then
         call printprofile
         if (kick > 10) then
            call noise(kick-10)
         else
            call noise(kick)
         endif
      !endif ! (mypid == NROOT)

      !call mpscsp(sp,spm,1)     ! replace mpscsp with
      spm = sp

      !if (mypid == NROOT) then
          st(1,:) = sr1(1,:)
         stm(1,:) = sr1(1,:)
          sz(3,:) = plavor
         szm(3,:) = plavor
      !endif
      end subroutine initfd


!     =================
!     SUBROUTINE READNL
!     =================

      subroutine readnl
      use pumamod
      implicit none

!     This workaround is necessaray, because allocatable arrays are
!     not allowed in namelists for FORTRAN versions < F2003

      integer, parameter :: MAXSELZW = 85
      integer, parameter :: MAXSELSP = ((MAXSELZW+1) * (MAXSELZW+2)) / 2
      integer :: nselect(0:MAXSELZW) = 1      ! NSELECT can be used up tp T42
      integer :: nspecsel(MAXSELSP)  = 1      ! Default setting: all modes active
      integer :: ndl(MAXLEV)         = 0      ! Diagnostics off
      real    :: sigmah(MAXLEV)      = 0.0    ! Half level sigma
      real    :: t0k(MAXLEV)         = 250.0  ! Reference temperature

      integer :: itru, icsp, ilev
      integer :: ios

      ! pkugcm_namelist was derived and extended from puma_namelist as below
      namelist /pkugcm_nl/ &
        akap    , alpha   , alr     , alrs    , diffts  , disp    &
      , dtep    , dtns    , dtrop   , dttrp   , dtzz    , dvdiff  &
      , ga      , gascon  &
      , kick    , mpstep  , nafter  , ncoeff  , nconv   , ncu     &
      , ndel    , ndheat  , ndiag   , ndiagp  , ndl     , nenergy &
      , nentropy, newsr   , nextout , nhelsua , nkits   &
      , nlevt   , nmonths , noutput , nradcv  , nruido  , nrun    &
      , nselect , nspecsel, nsponge , nstep   , nstop   &
      , ntspd   , nvg     , nwpd    , nwspini , nyears  &
      , orofac  , pac     , plarad  , pspon   , psurf   &
      , rotspd  , seed    , sid_day , sigmah  , sigmax  , dcsponge&
      , spstep  , t0k     , tauf    , taur    &
      , tac     , tauta   , tauts   , tgr     , ww_time

      open(13,file=puma_namelist,iostat=ios)
      if (ios == 0) then
         read (13,pkugcm_nl)
         close(13)
      endif

!--- modify basic scales according to namelist 

      if (ww_time < 1.0) ww_time = sid_day / TWOPI ! time scale
      ww_scale = 1.0 / ww_time      ! reciprocal of time scale 1/Omega
      cv    = plarad*ww_scale       ! velocity scale (velocity at the equator)
      ct    = cv*cv/gascon    ! temperature scale from hydrostatic equation 
      if (mpstep == 0 .and. spstep == 0.0) then ! automatic timestep
         mpstep = (60 * 32) / nlat              ! 60 min for T21
         spstep = mpstep * 60.0
      endif
      if (spstep == 0.0) spstep = 60.0 * mpstep
      if (ntspd == 0) ntspd = int(sol_day / spstep) !XW(2017-4-7): add int()
      
      nafter = ntspd                             ! daily output
      if (nwpd > 0 .and. nwpd <= ntspd) then
         nafter = ntspd / nwpd
      endif
      if (ndiag  < 1) ndiag  = ntspd * 10       ! every 10th. day

      !if (syncsecs > 0.0) syncstr = spstep / syncsecs
      !if (syncstr  > 1.0) syncstr = 1.0

      write(nud,pkugcm_nl)

      itru = ntru
      if (itru > MAXSELZW) itru = MAXSELZW
      icsp = ncsp
      if (icsp > MAXSELSP) icsp = MAXSELSP
      ilev = nlev
      if (ilev > MAXLEV)   ilev = MAXLEV

      nselzw(0:itru) = nselect(0:itru)  ! Copy values to allocated array
      nselsp(1:icsp) = nspecsel(1:icsp) 
      ndil(1:ilev)   = ndl(1:ilev)
      sigmh(1:ilev)  = sigmah(1:ilev)
      t0(1:ilev)     = t0k(1:ilev)

      end subroutine readnl

      
!     =============================
!     SUBROUTINE SELECT_ZONAL_WAVES
!     =============================

      subroutine select_zonal_waves
      use pumamod
      implicit none

      if (sum(nselzw(:)) /= NTP1) then ! some wavenumbers disabled
         lselect = .true.
      endif
      return
      end subroutine select_zonal_waves


!     ================================
!     SUBROUTINE SELECT_SPECTRAL_MODES
!     ================================

      subroutine select_spectral_modes
      use pumamod
      implicit none

      if (sum(nselsp(:)) /= NCSP) then ! some modes disabled
         lspecsel = .true.
      endif
      return
      end subroutine select_spectral_modes


!     =====================
!     * SET VERTICAL GRID *
!     =====================

      subroutine set_vertical_grid
      use pumamod
      implicit none

      ! local
      real :: zsigtran, zsigmin
      integer :: inl
      integer :: jlev

      if (sigmh(NLEV) /= 0.0) return ! Already read in from namelist puma

      ! nvg==1: Scinocca & Haynes sigma levels
      !
      if (nvg == 1) then

         if (nlevt >= NLEV) then      ! Security check for 'nlevt'
            write(nud,*) '*** ERROR *** nlevt >= NLEV'
            write(nud,*) 'Number of levels (NLEV): ',NLEV
            write(nud,*) 'Number of tropospheric levels (nlevt): ',nlevt
         endif   

!     troposphere: linear spacing in sigma
!     stratosphere: linear spacing in log(sigma)
!     after (see their Appendix):
!     Scinocca, J. F. and P. H. Haynes (1998): Dynamical forcing of
!        stratospheric planetary waves by tropospheric baroclinic eddies.
!        J. Atmos. Sci., 55 (14), 2361-2392

!     Here, zsigtran is set to sigma at dtrop (tropopause height for
!     construction of restoration temperature field). If tgr=288.15K,
!     ALR=0.0065K/km and dtrop=11.km, then zsigtran=0.223 (=0.1 in
!     Scinocca and Haynes (1998)).
!     A smoothing of the transition between linear and logarithmic
!     spacing, as noted in Scinocca and Haynes (1998), is not yet
!     implemented.

         zsigtran = (1. - alr * dtrop / tgr)**(ga/(gascon*alr))
         zsigmin = 1. - (1. - zsigtran) / real(nlevt)

         do jlev=1,NLEV
            if (jlev == 1) then
               sigmh(jlev) = SIGMAX
            elseif (jlev > 1 .and. jlev < NLEV - nlevt) then
               sigmh(jlev) = exp((log(SIGMAX) - log(zsigtran))         &
     &             / real(NLEV - nlevt - 1) * real(NLEV - nlevt - jlev) &
     &             + log(zsigtran))
            elseif (jlev >= NLEV - nlevt .and. jlev < NLEV - 1) then
               sigmh(jlev) = (zsigtran - zsigmin) / real(nlevt - 1)    &
     &                        * real(NLEV - 1 - jlev) + zsigmin
            elseif (jlev == NLEV - 1) then
               sigmh(jlev) = zsigmin
            elseif (jlev == NLEV) then
               sigmh(jlev) = 1.
            endif
         enddo
      end if
     
      ! nvg==2: Polvani & Kushner sigma levels
      !
      if (nvg == 2) then
         inl = int(real(NLEV)/(1.0 - sigmax**(1.0/5.0)))
         do jlev=1,NLEV
            sigmh(jlev) = (real(jlev + inl - NLEV) / real(inl))**5
         enddo
      end if

      ! Default (nvg == 0) : equidistant sigma levels
      !
      if (nvg==0) then
         do jlev = 1 , NLEV
            sigmh(jlev) = real(jlev) / real(NLEV)
         enddo
      end if
      end subroutine set_vertical_grid


!     =================
!     SUBROUTINE INITPM
!     =================

      subroutine initpm
      use pumamod
      implicit none

      real (kind=8) :: radea,zakk,zzakk
      real :: zsigb           ! sigma_b for Held & Suarez frictional
!                               and heating timescales
      ! local
      real     :: zrsq2       ! 临时变量=1/sqrt(2), 后续会不断变号
      real     :: zsq         ! 临时变量 只保存某整数jn的近似平方=jn*(jn+1)
      integer  :: jdelh, jlev
      integer  :: ji,jm,jn,jr,jw

      radea  = plarad         ! planet radius in high precision
      plavor = EZ * rotspd * omega * ww_time ! planetary vorticity

!     *************************************************************
!     * carries out all initialisation of model prior to running. *
!     * major sections identified with comments.                  *
!     * this s/r sets the model parameters and all resolution     *
!     * dependent quantities.                                     *
!     *************************************************************

      if (lrestart) nkits=0

!     ****************************************************
!     * Check for enabling / disabling zonal wavenumbers *
!     ****************************************************

      call select_zonal_waves
      if (npro == 1) call select_spectral_modes

!     *********************
!     * set vertical grid *
!     *********************

      call set_vertical_grid

      dsigma(1     ) = sigmh(1)
      dsigma(2:NLEV) = sigmh(2:NLEV) - sigmh(1:NLEM)

      rdsig(:) = 0.5 / dsigma(:)

      sigma(1     ) = 0.5 * sigmh(1)
      sigma(2:NLEV) = 0.5 * (sigmh(1:NLEV-1) + sigmh(2:NLEV))

!     Initialize profile of tau R if not set in namelist

      if (taur(NLEV) == 0.0) then
         do jlev = 1 , NLEV
            taur(jlev) = sid_day * 50.0 * atan(1.0 - sigma(jlev))
            if (taur(jlev) > 30.0 * sid_day) taur(jlev) = 30.0 * sid_day
         enddo
      endif

!     Initialize profile of tau F if not set in namelist

      if (tauf(NLEV) == 0.0) then
         do jlev = 1 , NLEV
            if (sigma(jlev) > 0.8) then
               tauf(jlev) = exp(10.0 * (1.0 - sigma(jlev))) / 2.718 * sid_day
            endif
         enddo
      endif

!     Compute 1.0 / (2 Pi * tau) for efficient use in calculations
!     A day is 2 Pi in non dimensional units using omega as scaling

      where (taur(1:NLEV) > 0.0)
!        damp(1:NLEV) = ww_time / taur(1:NLEV)
         damp(1:NLEV) = wwt     / taur(1:NLEV)
      endwhere

      where (tauf(1:NLEV) > 0.0)
!         fric(1:NLEV) = ww_time / tauf(1:NLEV)
          fric(1:NLEV) = wwt     / tauf(1:NLEV)
      endwhere

      if (nsponge == 1) call sponge


!     annual cycle period and phase in timesteps

      if (tac > 0.0) tac = TWOPI / (ntspd * tac)
      pac = pac * ntspd

!     compute internal diffusion parameter

      jdelh = ndel/2
      if (diffts > 0.0) then
!        zakk = ww_scale*(radea**ndel)/(TWOPI*diffts/sol_day &
         zakk = 1.0/wwt*(radea**ndel)/(TWOPI*diffts/sol_day &
              * ((NTRU*(NTRU+1.))**jdelh))
      else
         zakk = 0.0
      endif
      zzakk = zakk / (1.0/wwt*(radea**ndel))

!     set coefficients which depend on wavenumber

      zrsq2 = 1.0 / sqrt(2.0)

      jr =-1
      jw = 0
      do jm=0,NTRU
         do jn=jm,NTRU
            jr=jr+2
            ji=jr+1
            jw=jw+1
            nindex(jr)=jn
            nindex(ji)=jn
            spnorm(jr)=zrsq2
            spnorm(ji)=zrsq2
            zsq = jn * (jn+1)
            if (jn > 0) then
               srcn(jr) = 1.0 / zsq
               srcn(ji) = srcn(jr)
            endif
            sak(jr) = -zzakk * zsq**jdelh
            sak(ji) = sak(jr)
         enddo
         zrsq2=-zrsq2
      enddo

! finally make temperatures dimensionless

      dtns  = dtns    / ct
      dtep  = dtep    / ct
!     dttrp = dttrp   / ct
      t0(:) = t0(:) / ct

!     print out

      write(nud,8120)
      write(nud,8000)
      write(nud,8010) NLEV
      write(nud,8020) NTRU
      write(nud,8030) NLAT
      write(nud,8040) NLON
      if (zakk == 0.0) then
         write(nud,8060)
      else
         write(nud,8070) ndel
         write(nud,8080)
         write(nud,8090) zakk,ndel
         write(nud,8100) diffts
      endif
      write(nud,8110) PNU
      write(nud,8000)
      write(nud,8120)
      return

 8000 format('*****************************************************')
 8010 format('* NLEV = ',i6,'   Number of levels                  *')
 8020 format('* NTRU = ',i6,'   Triangular truncation             *')
 8030 format('* NLAT = ',i6,'   Number of latitudes               *')
 8040 format('* NLON = ',i6,'   Number of longitues               *')
 8060 format('*                 No lateral dissipation            *')
 8070 format('* ndel = ',i6,'   Lateral dissipation               *')
 8080 format('* on vorticity, divergence and temperature          *')
 8090 format('* with diffusion coefficient = ',e13.4,' m**',i1,'/s *')
 8100 format('* e-folding time for smallest scale is ',f7.0,' sec  *')
 8110 format('* Robert time filter with parameter PNU =',f8.3,'   *')
 8120 format(/)
      end subroutine initpm


!     =================
!     SUBROUTINE MAKEBM
!     =================

      subroutine makebm
      use pumamod
      implicit none

      ! local
      real :: zdeltsq, zaq
      integer :: jn, jlev, jlev1, jlev2

      zdeltsq = delt * delt

      do jlev1 = 1 , NLEV
         do jlev2 = 1 , NLEV
            zaq = zdeltsq * (t0(jlev1) * dsigma(jlev2)&
     &          + dot_product(xlphi(:,jlev1),xlt(jlev2,:)))
            bm1(jlev2,jlev1,1:NTRU) = zaq
         enddo
      enddo

      do jn=1,NTRU
         do jlev = 1 , NLEV
            bm1(jlev,jlev,jn) = bm1(jlev,jlev,jn) + 1.0 / (jn*(jn+1))
         enddo
         call minvers(bm1(1,1,jn),NLEV)
      enddo
      return
      end subroutine makebm


!     =================
!     SUBROUTINE INITSI
!     =================

      subroutine initsi
      use pumamod
      implicit none

!     **********************************************
!     * Initialisation of the Semi Implicit scheme *
!     **********************************************

      ! local
      real, dimension(NLEV)      :: zalp,   zh
      real, dimension(NLEV,NLEV) :: ztautk, ztaudt
      real :: zfctr
      integer :: i,j,ilev,jlev

      tkp(:) = akap * t0(:)
      t0d(1:NLEM) = t0(2:NLEV) - t0(1:NLEM)

      zalp(2:NLEV) = log(sigmh(2:NLEV)) - log(sigmh(1:NLEM))

      xlphi(:,:) = 0.0
      xlphi(1,1) = 1.0
      do jlev = 2 , NLEV
         xlphi(jlev,jlev) = 1.0 - zalp(jlev)*sigmh(jlev-1)/dsigma(jlev)
         xlphi(jlev,1:jlev-1) = zalp(jlev)
      enddo

      do jlev = 1 , NLEV
         c(jlev,:) = xlphi(:,jlev) * (dsigma(jlev) / dsigma(:))
      enddo

!     ***********************   tkp(i) = t0(i) * AKAP
!     * matrix xlt - part 1 *
!     ***********************

      do jlev = 1 , NLEV
         ztautk(:,jlev) = tkp(jlev) * c(:,jlev)
      enddo

!     *********************   dsigma(i) = sigmh(i) - sigmh(i-1)
!     * matrix xlt part 2 *   rdsig (i) = 0.5 / dsigma(i)
!     *********************

      ztaudt(1,1)      = 0.5 * t0d(1) * (sigmh(1) - 1.0)
      ztaudt(2:NLEV,1) = 0.5 * t0d(1) *  dsigma(2:NLEV)

      do j= 2 , NLEV
         do i = 1 , j-1
            ztaudt(i,j) =  dsigma(i) * rdsig(j) &
            * (t0d(j-1) * (sigmh(j-1)-1.0) + t0d(j) * (sigmh(j)-1.0))
         enddo
            ztaudt(j,j) =  0.5                  &
            * (t0d(j-1) *  sigmh(j-1)      + t0d(j) * (sigmh(j)-1.0))
         do i = j+1 , NLEV
            ztaudt(i,j) =  dsigma(i) * rdsig(j) &
            * (t0d(j-1) *  sigmh(j-1)      + t0d(j) *  sigmh(j)     )
         enddo
      enddo

      xlt(:,:) = ztautk(:,:) + ztaudt(:,:)

!     xlt finished

      zfctr=0.001*cv*cv/ga
      do jlev=1,NLEV
         zh(jlev) = dot_product(xlphi(:,jlev),t0(:)) * zfctr
      enddo

!     **********************************
!     * write out vertical information *
!     **********************************

      ilev = min(NLEV,5)
      write(nud,9001)
      write(nud,9002)
      write(nud,9003)
      write(nud,9002)
      do jlev=1,NLEV
        write(nud,9004) jlev,sigma(jlev),t0(jlev)*ct,zh(jlev)
      enddo
      write(nud,9002)
      write(nud,9001)

!     matrix c

      write(nud,9012)
      write(nud,9013) 'c',(jlev,jlev=1,ilev)
      write(nud,9012)
      do jlev=1,NLEV
        write(nud,9014) jlev,(c(i,jlev),i=1,ilev)
      enddo
      write(nud,9012)
      write(nud,9001)

!     matrix xlphi

      write(nud,9012)
      write(nud,9013) 'xlphi',(jlev,jlev=1,ilev)
      write(nud,9012)
      do jlev=1,NLEV
        write(nud,9014) jlev,(xlphi(i,jlev),i=1,ilev)
      enddo
      write(nud,9012)
      write(nud,9001)
      return
 9001 format(/)
 9002 format(33('*'))
 9003 format('* Lv *    Sigma Basic-T  Height *')
 9004 format('*',i3,' * ',3f8.3,' *')
 9012 format(69('*'))
 9013 format('* Lv * ',a5,i7,4i12,' *')
 9014 format('*',i3,' * ',5f12.8,' *')
      end subroutine initsi


!     =====================
!     SUBROUTINE INITRANDOM
!     =====================

      subroutine initrandom
      use pumamod
      implicit none

      integer :: i, clock

!     Set random number generator seed

      call random_seed(size=nseedlen)
      allocate(meed(nseedlen))

!     Take seed from namelist parameter 'SEED' ?

      if (seed(1) /= 0) then
         meed(:) = 0
         i = nseedlen
         if (i > 8) i = 8
         meed(1:i) = seed(1:i)
      else
         call system_clock(count=clock)
         meed(:) = clock + 37 * (/(i,i=1,nseedlen)/)
      endif
      call random_seed(put=meed)
      return
      end subroutine initrandom


!     ====================
!     SUBROUTINE PRINTSEED
!     ====================

      subroutine printseed
      use pumamod
      implicit none

      integer :: i

      write (nud,9020)
      write (nud,9010)
      do i = 1 , nseedlen
         write (nud,9000) i,meed(i)
      enddo
      write (nud,9010)
      write (nud,9020)
      return
 9000 format('* seed(',i1,') = ',i10,' *')
 9010 format('************************')
 9020 format(/)
      end subroutine printseed


!     ====================
!     SUBROUTINE INITRUIDO
!     ====================

      subroutine initruido
      use pumamod
      implicit none
      if (nruido > 0) then
         allocate(ruido(nlon,nlat,nlev))
         allocate(ruidop(nhor,nlev))
         ruido = 77
         ruidop = 88
      endif
      return
      end subroutine initruido


!     ====================
!     SUBROUTINE STEPRUIDO
!     ====================

      subroutine stepruido
      use pumamod
      implicit none

      ! local
      real :: zr
      real :: gasdev
      integer :: jlon,jlat,jlev

      !if (mypid == NROOT) then
         if (nruido == 1) then         ! ruido = 全场均一随机数
            zr = disp*gasdev()
            ruido(:,:,:) = zr
         elseif (nruido == 2) then     ! ruido = 全场不均一随机数
            do jlev=1,NLEV
               do jlat=1,NLAT
                  do jlon=1,NLON
                     ruido(jlon,jlat,jlev) = disp*gasdev()
                  enddo
               enddo
            enddo
         elseif (nruido == 3) then     ! ruido = 全场南北半球对称的不均一随机数
            do jlev=1,NLEV
               do jlat=1,NLAT,2
                  do jlon=1,NLON
                     ruido(jlon,jlat  ,jlev) = disp*gasdev()
                     ruido(jlon,jlat+1,jlev) = ruido(jlon,jlat,jlev)
                  enddo
               enddo
            enddo
         endif
      !endif ! (mypid == NROOT)

      !call mpscgp(ruido,ruidop,NLEV)
      ruidop = reshape(ruido,(/nhor,nlev/))
      return
      end subroutine stepruido


!     ==================
!     SUBROUTINE MINVERS
!     ==================

      subroutine minvers(a,n)
      implicit none
      integer, intent(in) :: n
      real, dimension(n,n), intent(inout) :: a

      ! local
      real, dimension(n,n) :: b
      integer, dimension(n)   :: indx
      integer :: j

      b = 0.0
      do j = 1 , n
         b(j,j) = 1.0
      enddo
      call ludcmp(a,n,indx)
      do j = 1 , n
         call lubksb(a,n,indx,b(1,j))
      enddo
      a = b
      return
      end subroutine minvers


!     =================
!     SUBROUTINE LUBKSB
!     =================

      subroutine lubksb(a,n,indx,b)
      dimension a(n,n),b(n),indx(n)
      k = 0
      do i = 1 , n
         l    = indx(i)
         sum  = b(l)
         b(l) = b(i)
         if (k > 0) then
            do j = k , i-1
               sum = sum - a(i,j) * b(j)
            enddo
         else if (sum /= 0.0) then
            k = i
         endif
         b(i) = sum
      enddo

      do i = n , 1 , -1
         sum = b(i)
         do j = i+1 , n
            sum = sum - a(i,j) * b(j)
         enddo
         b(i) = sum / a(i,i)
      enddo
      return
      end subroutine lubksb


!     =================
!     SUBROUTINE LUDCMP
!     =================

      subroutine ludcmp(a,n,indx)
      dimension a(n,n),indx(n),vv(n)

      d = 1.0
      vv = 1.0 / maxval(abs(a),2)

      do 19 j = 1 , n
         do i = 2 , j-1
            a(i,j) = a(i,j) - dot_product(a(i,1:i-1),a(1:i-1,j))
         enddo
         aamax = 0.0
         do i = j , n
            if (j > 1) &
     &      a(i,j) = a(i,j) - dot_product(a(i,1:j-1),a(1:j-1,j))
            dum = vv(i) * abs(a(i,j))
            if (dum .ge. aamax) then
               imax = i
               aamax = dum
            endif
         enddo
         if (j .ne. imax) then
            do 17 k = 1 , n
               dum = a(imax,k)
               a(imax,k) = a(j,k)
               a(j,k) = dum
   17       continue
            d = -d
            vv(imax) = vv(j)
         endif
         indx(j) = imax
         if (a(j,j) == 0.0) a(j,j) = tiny(a(j,j))
         if (j < n) a(j+1:n,j) = a(j+1:n,j) / a(j,j)
   19 continue
      return
      end subroutine ludcmp


!     =============================
!     SUBROUTINE FILTER_ZONAL_WAVES
!     =============================

      subroutine filter_zonal_waves(pfc)
      use pumamod
      implicit none

      real, dimension(2,NLON/2,NLPP) :: pfc
      integer :: jlat

      do jlat = 1 , NLPP
         pfc(1,1:NTP1,jlat) = pfc(1,1:NTP1,jlat) * nselzw(:)
         pfc(2,1:NTP1,jlat) = pfc(2,1:NTP1,jlat) * nselzw(:)
      enddo

      return
      end subroutine filter_zonal_waves
      

!     ================================
!     SUBROUTINE FILTER_SPECTRAL_MODES
!     ================================

      subroutine filter_spectral_modes
      use pumamod
      implicit none

      integer :: j,k,m,n

      j =  0
      k = -1
      do m = 0 , NTRU
         do n = m , NTRU
            k = k + 2
            j = j + 1
            if (nselsp(j) == 0) then
               spp(k:k+1  ) = 0.0
               sdp(k:k+1,:) = 0.0
               stp(k:k+1,:) = 0.0
               spt(k:k+1  ) = 0.0
               sdt(k:k+1,:) = 0.0
               stt(k:k+1,:) = 0.0
               spm(k:k+1  ) = 0.0
               sdm(k:k+1,:) = 0.0
               stm(k:k+1,:) = 0.0
              srp1(k:k+1,:) = 0.0
              srp2(k:k+1,:) = 0.0
               if (n < NTRU) then
                  szp(k+2:k+3,:) = 0.0
                  szt(k+2:k+3,:) = 0.0
                  szm(k+2:k+3,:) = 0.0
               endif
            endif
         enddo
      enddo

      return
      end subroutine filter_spectral_modes
      

!     ================
!     SUBROUTINE NOISE
!     ================

      subroutine noise(kickval)
      use pumamod
      implicit none

!     kickval = -1 : read ln(ps) from puma_sp_init
!     kickval =  0 : model runs zonally symmetric with no eddies
!     kickval =  1 : add white noise to ln(Ps) asymmetric hemispheres
!     kickval =  2 : add white noise to ln(Ps) symmetric to the equator
!     kickval =  3 : force mode(1,2) of ln(Ps) allowing reproducable runs
!     kickval =  4 : add white noise to symmetric zonal wavenumbers 7 of ln(Ps)

      integer :: kickval
      integer :: jsp, jsp1, jn, jm
      integer :: jr, ji, ins
      real    :: zr, zi, zscale, zrand
      integer :: iostat

      zscale = 0.000001         ! amplitude of noise
      zr     = 0.001            ! kickval=3 value for mode(1,2) real
      zi     = 0.0005           ! kickval=3 value for mode(1,2) imag

      select case (kickval)
      case (-1)
         open(71, file=puma_sp_init,form='unformatted',iostat=iostat)
         if (iostat /= 0) then
            write(nud,*) ' *** kick=-1: needs file <',trim(puma_sp_init),'> ***'
            stop
         endif
         read(71,iostat=iostat) sp(:)
         if (iostat /= 0) then
            write(nud,*) ' *** error reading file <',trim(puma_sp_init),'> ***'
            stop
         endif
         close(71)
         write(nud,*) 'initial ln(ps) field read from <',trim(puma_sp_init),'>'
         return
      case (0)                  ! do nothing
      case (1)
         jsp1=2*NTP1+1
         do jsp=jsp1,NRSP
            call random_number(zrand)
            !if (mrpid > 0) zrand = zrand + mrpid * 0.01
            sp(jsp)=sp(jsp)+zscale*(zrand-0.5)
         enddo
         write(nud,*) 'white noise added'
      case (2)
         jr=2*NTP1-1
         do jm=1,NTRU
            do jn=jm,NTRU
               jr=jr+2
               ji=jr+1
               if (mod(jn+jm,2) == 0) then
                  call random_number(zrand)
                  !if (mrpid > 0) zrand = zrand + mrpid * 0.01
                  sp(jr)=sp(jr)+zscale*(zrand-0.5)
                  sp(ji)=sp(ji)+zscale*(zrand-0.5)
               endif
            enddo
         enddo
         write(nud,*) 'symmetric white noise added'
      case (3)
         sp(2*NTP1+3) = sp(2*NTP1+3) + zr
         sp(2*NTP1+4) = sp(2*NTP1+4) + zi
         write(nud,*) 'mode(1,2) of ln(Ps) set to (',sp(2*NTP1+3),',',sp(2*NTP1+4),')'
      case (4)
         jr=2*NTP1-1
         do jm=1,NTRU
            do jn=jm,NTRU
               jr=jr+2
               ji=jr+1
               if (mod(jn+jm,2) == 0 .and. jm == 7) then
                  call random_number(zrand)
                  sp(jr)=sp(jr)+zscale*(zrand-0.5)
                  sp(ji)=sp(ji)+zscale*(zrand-0.5)
               endif
            enddo
         enddo
         write(nud,*) 'symmetric zonal wavenumbers 7 of ln(Ps) perturbed',   &
     &        'with white noise.'
      case default
         write(nud,*) 'Value ',kickval  ,' for kickval not implemented.'
         stop
      end select

      if (nwspini == 1) then
         open(71, file=puma_sp_init, form='unformatted')
         write(71) sp(:)
         close(71)
      endif

      return
      end subroutine noise


!     ================
!     SUBROUTINE SETZT
!     ================
      subroutine setzt
      use pumamod
      implicit none

!     *************************************************************
!     * Set up the restoration temperature fields sr1 and sr2     *
!     * for aqua planet conditions.                               *
!     * The temperature at sigma = 1 is <tgr>, entered in kelvin. *
!     * The lapse rate of ALR K/m is assumed under the tropopause *
!     * and zero above. The tropopause is defined by <dtrop>.     *
!     * The smoothing ot the tropopause depends on <dttrp>.       *
!     ************************************************************* 

      real, dimension(NLEV) :: ztrs   ! Mean profile
      real, dimension(NLEV) :: zfac
      real :: zsqrt2, zsqrt6, zsqrt04
      real :: zsigprev
      real :: ztp, ztpm, ztpp, ztprev, ztps, zttrop, zzp, zzpp, zzprev
      integer :: jlev

      sr1(:,:) = 0.0 ! NESP,NLEV
      sr2(:,:) = 0.0 ! NESP,NLEV

!     Temperatures in [K]

      zsigprev = 1.0  ! sigma value
      ztprev   = tgr  ! Temperature [K]
      zzprev   = 0.0  ! Height      [m]

      do jlev = NLEV , 1 , -1   ! from bottom to top of atmosphere
        zzp=zzprev+(gascon*ztprev/ga)*log(zsigprev/sigma(jlev))
        ztp=tgr-dtrop*alr ! temperature at tropopause
        ztp=ztp+sqrt((.5*alr*(zzp-dtrop))**2+dttrp**2)
        ztp=ztp-.5*alr*(zzp-dtrop)
        ztpm=.5*(ztprev+ztp)
        zzpp=zzprev+(gascon*ztpm/ga)*log(zsigprev/sigma(jlev))
        ztpp=tgr-dtrop*alr
        ztpp=ztpp+sqrt((.5*alr*(zzpp-dtrop))**2+dttrp**2)
        ztpp=ztpp-.5*alr*(zzpp-dtrop)
        ztrs(jlev)=ztpp
        zzprev=zzprev+(.5*(ztpp+ztprev)*gascon/ga)*log(zsigprev/sigma(jlev))
        ztprev=ztpp
        zsigprev=sigma(jlev)
      enddo

      do jlev=1,NLEV
         ztrs(jlev)=ztrs(jlev)/ct
      enddo

!******************************************************************
! loop to set array zfac - this controls temperature gradients as a
! function of sigma in tres. it is a sine wave from one at
! sigma = 1 to zero at stps (sigma at the tropopause) .
!******************************************************************
! first find sigma at dtrop
!
      zttrop=tgr-dtrop*alr
      ztps=(zttrop/tgr)**(ga/(alr*gascon))
!
! now the latitudinal variation in tres is set up ( this being in terms
! of a deviation from t0 which is usually constant with height)
!
      zsqrt2  = sqrt(2.0)
      zsqrt04 = sqrt(0.4)
      zsqrt6  = sqrt(6.0)
      do 2100 jlev=1,NLEV
        zfac(jlev)=sin(0.5*PI*(sigma(jlev)-ztps)/(1.-ztps))
        if (zfac(jlev).lt.0.0) zfac(jlev)=0.0
        sr1(1,jlev)=zsqrt2*(ztrs(jlev)-t0(jlev))
        sr2(3,jlev)=(1./zsqrt6)*dtns*zfac(jlev)
        sr1(5,jlev)=-2./3.*zsqrt04*dtep*zfac(jlev)
 2100 continue
      write(nud,*) '**************************************************'
      write(nud,*) '* Restoration Temperature set up for aqua planet *'
      write(nud,*) '**************************************************'
      return
      end subroutine setzt


!     =======================
!     SUBROUTINE PRINTPROFILE
!     =======================

      subroutine printprofile
      use pumamod
      implicit none

      integer :: jlev
      real :: zt

!     **********************************
!     * write out vertical information *
!     **********************************

      write(nud,9001)
      write(nud,9002)
      write(nud,9003)
      write(nud,9002)

      do jlev=1,NLEV
         zt = (sr1(1,jlev)/sqrt(2.0) + t0(jlev)) * ct
         if (tauf(jlev) > 0.1) then
            write(nud,9004) jlev,sigma(jlev),zt,taur(jlev),tauf(jlev)
         else
            write(nud,9005) jlev,sigma(jlev),zt,taur(jlev)
         endif
      enddo

      write(nud,9002)
      write(nud,9001)
      return
 9001 format(/)
 9002 format(46('*'))
 9003 format('* Lv *    Sigma Restor-T tauR [s]  tauF [s]  *')
 9004 format('*',i3,' * ',f8.3,f9.3,2e10.3,' *')
 9005 format('*',i3,' * ',f8.3,f9.3,e10.3,'         - *')
      end subroutine printprofile


!     ====================
!     SUBROUTINE READ_SURF
!     ====================

      subroutine read_surf(kcode,psp,klev,kread)
      use pumamod

      logical :: lexist
      integer :: kread
      integer :: ihead(8)
      character(len=256) :: yfilename
      real :: psp(NESP,klev)
      real :: zgp(NUGP,klev)
      real :: zpp(NHOR,klev)

      kread = 0
      !if (mypid == NROOT) then
         if (NLAT < 1000) then
         write(yfilename,'("N",I3.3,"_surf_",I4.4,".sra")') NLAT,kcode
         else
         write(yfilename,'("N",I4.4,"_surf_",I4.4,".sra")') NLAT,kcode
         endif
         inquire(file=yfilename,exist=lexist)
      !endif
      !call mpbcl(lexist)
      if (.not. lexist) return

      !if (mypid == NROOT) then
         open(65,file=yfilename,form='formatted')
         write(nud,*) 'Reading file <',trim(yfilename),'>'
         do jlev = 1 , klev
            read (65,*) ihead(:)
            read (65,*) zgp(:,jlev)
         enddo
         close(65)
         if (kcode == 134) then
            write(nud,*) "Converting Ps to LnPs"
            zscale   = log(100.0) - log(psurf) ! Input [hPa] / PSURF [Pa]
            zgp(:,:) = log(zgp(:,:)) + zscale
         endif
         call reg2alt(zgp,klev)
      !endif ! (mypid == NROOT)

      !call mpscgp(zgp,zpp,klev)
      zpp = zgp

      call gp2fc(zpp,NLON,NLPP*klev,trigs)
      do jlev = 1 , klev
         call fc2sp(zpp(1,jlev),psp(1,jlev))
      enddo
      !call mpsum(psp,klev)
      kread = 1
      return
      end subroutine read_surf


!     =====================
!     SUBROUTINE READ_VARGP
!     =====================

      subroutine read_vargp(kcode,klev,kread)
      use pumamod
    
      logical :: lexist
      integer :: ihead(8)
      character(len=256) :: yfilename
      real :: zgp(NUGP,klev)

      kread = 0
      !if (mypid == NROOT) then
         if (NLAT < 1000) then
         write(yfilename,'("N",I3.3,"_surf_",I4.4,".sra")') NLAT,kcode
         else
         write(yfilename,'("N",I4.4,"_surf_",I4.4,".sra")') NLAT,kcode
         endif
         inquire(file=yfilename,exist=lexist)
      !endif
      !call mpbcl(lexist)
      if (.not. lexist) then
         !if (mypid == NROOT) then
            write(nud,*) 'File <',trim(yfilename),'> not found'
         !endif
         return
      endif

      !if (mypid == NROOT) then
         open(65,file=yfilename,form='formatted')
         write(nud,*) 'Reading file <',trim(yfilename),'>'
         do jlev = 1 , klev
            read (65,*) ihead(:)
            read (65,*) zgp(:,jlev)
         enddo
         close(65)
         call reg2alt(zgp,klev)
      !endif ! (mypid == NROOT)

      select case(kcode)
         case(121)
            !--- non-dimensionalize and shift const radiative rest. temp.
            !if (mypid == NROOT) then
               zgp(:,:) = zgp(:,:)/ct
               do jhor = 1,nugp
                  zgp(jhor,:) = zgp(jhor,:) - t0(:)
               enddo
            !endif
            allocate(gr1(nhor,klev))
            !if (mypid == NROOT) then
               write(nud,*) 'Field gr1 allocated'
            !endif
            !call mpscgp(zgp,gr1,klev)
            gr1 = zgp

         case(122)
            !--- non-dimensionalize variable. radiative rest. temp.
            !if (mypid == NROOT) then
               zgp(:,:) = zgp(:,:)/ct
            !endif
            allocate(gr2(nhor,klev))
            !if (mypid == NROOT) then
               write(nud,*) 'Field gr2 allocated'
            !endif
            !call mpscgp(zgp,gr2,klev)
            gr2 = zgp

         case(123)
            !--- non-dimensionalize radiative relaxation time scale
            !if (mypid == NROOT) then
               zgp(:,:) = zgp(:,:)/ww_scale
            !endif
            allocate(gtdamp(nhor,klev))
            !if (mypid == NROOT) then
               write(nud,*) 'Field gtdamp allocated'
            !endif
            !call mpscgp(zgp,gtdamp,klev)
            gtdamp = zgp

         case(124)
            !--- non-dimensionalize and shift const. convective rest. temp.
            !if (mypid == NROOT) then
               zgp(:,:) = zgp(:,:)/ct
               do jhor = 1,nugp
                  zgp(jhor,:) = zgp(jhor,:) - t0(:)
               enddo
            !endif
            allocate(gr1c(nhor,klev))
            !if (mypid == NROOT) then
               write(nud,*) 'Field gr1c allocated'
            !endif
            !call mpscgp(zgp,gr1c,klev)
            gr1c = zgp

         case(125)
            !--- non-dimensionalize variable. convective rest. temp.
            !if (mypid == NROOT) then
               zgp(:,:) = zgp(:,:)/ct
            !endif
            allocate(gr2c(nhor,klev))
            !if (mypid == NROOT) then
               write(nud,*) 'Field gr2c allocated'
            !endif
            !call mpscgp(zgp,gr2c,klev)
            gr2c = zgp

         case(126)
            !--- non-dimensionalize convective relaxation time scale
            !if (mypid == NROOT) then
               zgp(:,:) = zgp(:,:)/ww_scale
            !endif
            allocate(gtdampc(nhor,klev))
            !if (mypid == NROOT) then
               write(nud,*) 'Field gtdampc allocated'
            !endif
            !call mpscgp(zgp,gtdampc,klev)
            gtdampc = zgp

      end select
      kread = 1
      return
      end subroutine read_vargp


!     ===============
!     SUBROUTINE DIAG
!     ===============

      subroutine diag
      use pumamod
      if (noutput > 0 .and. mod(nstep,ndiag) == 0) then
         if (ncoeff > 0) call prisp
         call xsect
      endif
      call energy
      return
      end subroutine diag


!     ================
!     SUBROUTINE PRISP
!     ================

      subroutine prisp
      use pumamod

      character(30) :: title

      scale = 100.0
      title = 'Vorticity [10-2]'
      do 100 jlev=1,NLEV
         if (ndil(jlev).ne.0) call wrspam(sz(1,jlev),jlev,title,scale)
  100 continue

      title = 'Divergence [10-2]'
      do 200 jlev=1,NLEV
         if (ndil(jlev).ne.0) call wrspam(sd(1,jlev),jlev,title,scale)
  200 continue

      scale = 1000.0
      title = 'Temperature [10-3]'
      do 300 jlev=1,NLEV
         if (ndil(jlev).ne.0) call wrspam(st(1,jlev),jlev,title,scale)
  300 continue

      title = 'Pressure [10-3]'
      call wrspam(sp,0,title,scale)

      return
      end subroutine prisp


!     ====================
!     SUBROUTINE POWERSPEC
!     ====================

      subroutine powerspec(pf,pspec)
      use pumamod
      real :: pf(2,NCSP)
      real :: pspec(NTP1)

      do j = 1 , NTP1
         pspec(j) = 0.5 * (pf(1,j) * pf(1,j) + pf(2,j) * pf(2,j))
      enddo

      j = NTP1 + 1
      do m = 2 , NTP1
         do l = m , NTP1
            pspec(l) = pspec(l) + pf(1,j) * pf(1,j) + pf(2,j) * pf(2,j)
            j = j + 1
         enddo
      enddo
      return
      end subroutine powerspec


!     =====================
!     SUBROUTINE POWERPRINT
!     =====================

      subroutine powerprint(text,pspec)
      use pumamod
      character(3) :: text
      real :: pspec(NTP1)

      zmax = maxval(pspec(:))
      if (zmax <= 1.0e-20) return
      zsca = 10 ** (4 - int(log10(zmax)))
      write(nud,1000) text,(int(pspec(j)*zsca),j=2,13)
      return
 1000 format('* Power(',a3,') ',i8,11i5,' *')
      end subroutine powerprint


!     ==============
!     FUNCTION RMSSP
!     ==============

      function rmssp(pf)
      use pumamod
      real pf(NESP,NLEV)

      zsum = 0.0
      do jlev = 1 , NLEV
         zsum = zsum + dsigma(jlev)&
     &        * (dot_product(pf(1:NZOM,jlev),pf(1:NZOM,jlev)) * 0.5&
     &        +  dot_product(pf(NZOM+1:NRSP,jlev),pf(NZOM+1:NRSP,jlev)))
      enddo
      rmssp = zsum
      return
      end function rmssp


!     =================
!     SUBROUTINE ENERGY
!     =================

      subroutine energy
      use pumamod

      parameter (idim=5) ! Number of scalars for GUI timeseries

!     calculates various global diagnostic quantities
!     remove planetary vorticity so sz contains relative vorticity

      real :: spec(NTP1)
      real :: ziso(idim) !XW(2017-4-7): remove "kind=4" for ziso

      sz(3,:) = sz(3,:) - plavor

!    ***********************************************
!     calculate means - zpsitot rms vorticity
!                       zchitot rms divergence
!                       ztmptot rms temperature
!                       ztotp  ie+pe potential energy
!                       zamsp mean surface pressure
!     ***********************************************

      zsqrt2 = sqrt(2.0)
      zamsp  = 1.0 + span(1) / zsqrt2
      zst    = dot_product(dsigma(:),st(1,:)) / zsqrt2
      ztout1 = dot_product(dsigma(:),t0(:))

      ztout2 = 0.0
      zst2b  = 0.0
      ztoti  = 0.0
      do jlev = 1 , NLEV
         ztout2 = ztout2 + dsigma(jlev) * t0(jlev) * t0(jlev)
         zst2b  = zst2b  + dsigma(jlev) * t0(jlev) * st(1,jlev)
         ztoti  = ztoti + dsigma(jlev)&
     &          * (dot_product(span(1:NZOM),st(1:NZOM,jlev)) * 0.5&
     &          +  dot_product(span(NZOM+1:NRSP),st(NZOM+1:NRSP,jlev)))
      enddo

      ztotp = dot_product(span(1:NZOM),so(1:NZOM)) * 0.5&
     &      + dot_product(span(NZOM+1:NRSP),so(NZOM+1:NRSP))&
     &      + so(1)/zsqrt2 + (zamsp*ztout1+ztoti+zst) / akap

      zpsitot = sqrt(rmssp(sz))
      zchitot = sqrt(rmssp(sd))
      ztmptot = sqrt(rmssp(st)+ztout2+zst2b*zsqrt2)

      ziso(1) = ct * (spnorm(1) * st(1,NLEV) + t0(NLEV)) - 273.16 ! T(NLEV) [C]
      ziso(2) = ww_scale * zchitot * 1.0e6
      ziso(3) = ztmptot
      ziso(4) = ztotp
      ziso(5) = sz(3,2)
      !XW(Mar/25/2017) to remove GUI: call guiput("SCALAR" // char(0) ,ziso,idim,1,1)

!     restore sz to absolute vorticity

      sz(3,:) = sz(3,:) + plavor

      if (mod(nstep,ndiag) /= 0) return ! was called for GUI only
      write(nud,9001)
      write(nud,9002) nstep,zpsitot,zchitot,ztmptot,ztotp,zamsp
      write(nud,9002)
      write(nud,9011) (j,j=1,12)
      write(nud,9012)
      call powerspec(span,spec)
      call powerprint('Pre',spec)
      call powerspec(sz(1,NLEV),spec)
      call powerprint('Vor',spec)
      call powerspec(sd(1,NLEV),spec)
      call powerprint('Div',spec)
      call powerspec(st(1,NLEV),spec)
      call powerprint('Tem',spec)
      return
 9001 format(/,'     nstep     rms z       rms d       rms t       &
     & pe+ie       msp')
 9002 format(i10,4x,4g12.5,g15.8)
!9009 format('*',75(' '),' *')
!9010 format('* Power(',a,') ',7e9.2,' *')
 9011 format('* Wavenumber ',i8,11i5,' *')
 9012 format('',78('*'))
      end subroutine energy


!     =================
!     SUBROUTINE NTOMIN
!     =================

      subroutine ntomin(kstep,imin,ihou,iday,imon,iyea)
      use pumamod
      implicit none
      integer, intent(inout) :: kstep, imin, ihou, iday, imon, iyea
      integer :: istep

      istep = kstep                          ! day [0-29] month [0-11]
      if (istep .lt. 0) istep = 0            ! min [0-59] hour  [0-23]
      imin = mod(istep,ntspd) * 1440 / ntspd ! minutes of current day
      ihou = imin / 60                       ! hours   of current day
      imin = imin - ihou * 60                ! minutes of current hour
      iday = istep / ntspd                   ! days    in this run
      imon = iday / 30                       ! months  in this run
      iday = iday - imon * 30                ! days    of current month
      iyea = imon / 12                       ! years   in this run
      imon = imon - iyea * 12                ! month   of current year
      iday = iday + 1
      imon = imon + 1
      iyea = iyea + 1
      return
      end subroutine ntomin


!     =================
!     SUBROUTINE NTODAT
!     =================

      subroutine ntodat(istep,datch)
      character(18) :: datch
      character(3) :: mona(12)
      data mona /'Jan','Feb','Mar','Apr','May','Jun',&
     &           'Jul','Aug','Sep','Oct','Nov','Dec'/
      call ntomin(istep,imin,ihou,iday,imon,iyea)
      write(datch,20030) iday,mona(imon),iyea,ihou,imin
20030 format(i2,'-',a3,'-',i4.4,2x,i2,':',i2.2)
      end subroutine ntodat


!     =================
!     SUBROUTINE WRSPAM
!     =================

      subroutine wrspam(ps,klev,title,scale)
      use pumamod
!
      dimension ps(NRSP)
      character(30) :: title
      character(18) :: datch

!     cab(i)=real(scale*sqrt(ps(i+i-1)*ps(i+i-1)+ps(i+i)*ps(i+i)))

      call ntodat(nstep,datch)
      write(nud,'(1x)')
      write(nud,20000)
      write(nud,20030) datch,title,klev
      write(nud,20000)
      write(nud,20020) (i,i=0,9)
      write(nud,20000)
      write(nud,20100) (cab(i),i=1,10)
      write(nud,20200) (cab(i),i=NTRU+2,NTRU+10)
      write(nud,20300) (cab(i),i=2*NTRU+2,2*NTRU+9)
      write(nud,20400) (cab(i),i=3*NTRU+1,3*NTRU+7)
      write(nud,20000)
      write(nud,'(1x)')

20000 format(78('*'))
20020 format('* n * ',10i7,' *')
20030 format('*   * ',a18,2x,a30,'  Level ',i2,11x,'*')
20100 format('* 0 *',f8.2,9f7.2,' *')
20200 format('* 1 *',8x,9f7.2,' *')
20300 format('* 2 *',15x,8f7.2,' *')
20400 format('* 3 *',22x,7f7.2,' *')
      contains
      function cab(i)
         cab = scale * sqrt(ps(i+i-1)*ps(i+i-1)+ps(i+i)*ps(i+i))
      end function cab
      end subroutine wrspam


!     ===============
!     SUBROUTINE WRZS
!     ===============

      subroutine wrzs(zs,title,scale)
      use pumamod
!
      dimension zs(NLAT,NLEV)
      character(30) :: title
      character(18) :: datch

      ip = NLAT / 16
      ia = ip/2
      ib = ia + 7 * ip
      id = NLAT + 1 - ia
      ic = id - 7 * ip

      call ntodat(nstep,datch)
      write(nud,'(1x)')
      write(nud,20000)
      write(nud,20030) datch,title
      write(nud,20000)
      write(nud,20020) (chlat(i),i=ia,ib,ip),(chlat(j),j=ic,id,ip)
      write(nud,20000)
      do 200 jlev = 1 , NLEV
         write(nud,20100) jlev,((int(zs(i,jlev)*scale)),i=ia,ib,ip),&
     &                       ((int(zs(j,jlev)*scale)),j=ic,id,ip),jlev
  200 continue
      write(nud,20000)
      write(nud,'(1x)')

20000 format(78('*'))
20020 format('* Lv * ',16(1x,a3),' * Lv *')
20030 format('*    * ',a18,2x,a30,20x,'*')
20100 format('* ',i2,' * ',16i4,' * ',i2,' *')
      end subroutine wrzs


!     ================
!     SUBROUTINE XSECT
!     ================

      subroutine xsect
      use pumamod
      character(30) :: title

      scale = 10.0
      title = 'Zonal Wind [0.1 m/s]'
      call wrzs(csu,title,scale)
      title = 'Meridional Wind [0.1 m/s]'
      call wrzs(csv,title,scale)
      scale = 1.0
      title = 'Temperature [C]'
      call wrzs(cst,title,scale)
      return
      end subroutine xsect


!     ==================
!     SUBROUTINE WRITESP
!     ==================

      subroutine writesp(kunit,pf,kcode,klev,pscale,poff)
      use pumamod
      real    :: pf(NRSP)
      real    :: zf(NRSP)
      integer :: ihead(8)

      call ntomin(nstep,nmin,nhour,nday,nmonth,nyear)

      ihead(1) = kcode
      ihead(2) = klev
      ihead(3) = nday + 100 * nmonth + 10000 * nyear
      ihead(4) = nmin + 100 * nhour
      ihead(5) = NRSP
      ihead(6) = 1
      ihead(7) = 1
      ihead(8) = 0

!     normalize ECHAM compatible and scale to physical dimensions

      zf(:) = pf(:) * spnorm(1:NRSP) * pscale
      zf(1) = zf(1) + poff ! Add offset if necessary
      write(kunit) ihead
      write(kunit) zf

      return
      end subroutine writesp


!     ================
!     SUBROUTINE OUTSP
!     ================

      subroutine outsp
      use pumamod
      real zsr(NESP)

      if (nwrioro == 1) then
         call writesp(40,so,129,0,cv*cv,0.0)
         nwrioro = 0
      endif

      if (nextout == 1) then
         call writesp(40,sp2,40,0,1.0,log(psmean))
         call writesp(40,sp1,41,0,1.0,log(psmean))
         do jlev = 1,NLEV
            call writesp(40,st2(1,jlev),42,jlev,ct,t0(jlev)*ct)
         enddo
         do jlev = 1,NLEV
            call writesp(40,st1(1,jlev),43,jlev,ct,t0(jlev)*ct)
         enddo
      endif

!     ************
!     * pressure *

      call writesp(40,sp,152,0,1.0,log(psmean))

!     ***************
!     * temperature *

      do jlev = 1 , NLEV
         call writesp(40,st(1,jlev),130,jlev,ct,t0(jlev)*ct)
      enddo

!     ********************
!     * res. temperature *

      zampl = cos((real(nstep)-pac)*tac)
      do jlev = 1 , NLEV
         zsr(:)=sr1(:,jlev)+sr2(:,jlev)*zampl
         call writesp(40,zsr,154,jlev,ct,t0(jlev)*ct)
      enddo

!     **************
!     * divergence *

      do jlev = 1 , NLEV
         call writesp(40,sd(1,jlev),155,jlev,ww_scale,0.0)
      enddo

!     *************
!     * vorticity *

      do jlev = 1 , NLEV
         zsave = sz(3,jlev)
         sz(3,jlev) = sz(3,jlev) - plavor
         call writesp(40,sz(1,jlev),138,jlev,ww_scale,0.0)
         sz(3,jlev) = zsave
      enddo

      return
      end subroutine outsp


!     ==================
!     SUBROUTINE WRITEGP
!     ==================

      subroutine writegp(kunit,pf,kcode,klev)
      use pumamod
      real :: pf(NHOR)
      real :: zf(NUGP)
      integer :: ihead(8)

      !call mpgagp(zf,pf,1)
      zf = pf

      !if (mypid == NROOT) then 
         call alt2reg(zf,1)
         call ntomin(nstep,nmin,nhour,nday,nmonth,nyear)
   
         ihead(1) = kcode
         ihead(2) = klev 
         ihead(3) = nday + 100 * nmonth + 10000 * nyear
         ihead(4) = nmin + 100 * nhour
         ihead(5) = NLON 
         ihead(6) = NLAT 
         ihead(7) = 1
         ihead(8) = 0

         write(kunit) ihead
         write(kunit) zf
      !endif

      return
      end subroutine writegp


!     ================
!     SUBROUTINE OUTGP
!     ================

      subroutine outgp
      use pumamod
      real zhelp(NHOR)
      real, dimension(NUGP) :: x2d
!     
!     energy diagnostics
!   
      if(nenergy > 0) then
       do je=1,9
        jcode=300+je
        zhelp(:)=denergy(:,je)
        call writegp(40,zhelp,jcode,0)
       enddo
      endif
      if(nentropy > 0) then
       do je=1,9
        jcode=310+je
        zhelp(:)=dentropy(:,je)
        call writegp(40,zhelp,jcode,0)
       enddo
      endif

      end subroutine outgp


!     ====================
!     SUBROUTINE CHECKUNIT
!     ====================

      subroutine checkunit
      use pumamod
      implicit none

      write(ncu,1000) nstep,'sp(  1  )',sp(1),sp(1)*spnorm(1)+log(psmean)
      write(ncu,1000) nstep,'st(  1,1)',st(1,1),st(1,1)*spnorm(1)*ct+t0(1)*ct
      write(ncu,1000) nstep,'sd(  1,1)',sd(1,1),sd(1,1)*spnorm(1)*ww_scale
      write(ncu,1000) nstep,'sz(  1,1)',sz(1,1),sz(1,1)*spnorm(1)*ww_scale

      write(ncu,1000) nstep,'st(  1,NLEV)',st(1,NLEV),st(1,NLEV)*spnorm(1)*ct+t0(5)*ct
      write(ncu,1000) nstep,'sd(  1,NLEV)',sd(1,NLEV),sd(1,NLEV)*spnorm(1)*ww_scale
      write(ncu,1000) nstep,'sz(  1,NLEV)',sz(1,NLEV),sz(1,NLEV)*spnorm(1)*ww_scale

      if (100 < NRSP) then
      write(ncu,1000) nstep,'sp(100  )',sp(100),sp(100)*spnorm(100)
      write(ncu,1000) nstep,'st(100,NLEV)',st(100,NLEV),st(100,NLEV)*spnorm(100)*ct
      write(ncu,1000) nstep,'sd(100,NLEV)',sd(100,NLEV),sd(100,NLEV)*spnorm(100)*ww_scale
      write(ncu,1000) nstep,'sz(100,NLEV)',sz(100,NLEV),sz(100,NLEV)*spnorm(100)*ww_scale
      endif

      return
 1000 format(i5,1x,a,1x,2f14.7)
      end subroutine checkunit


!     =====================
!     * SUBROUTINE LEGPRI *
!     =====================

      subroutine legpri
      use pumamod
      implicit none

      ! local
      integer :: jlat
      real :: zalat

      write(nud,231)
      write(nud,232)
      write(nud,233)
      write(nud,232)
      do jlat = 1 , NLAT
         zalat = asin(sid(jlat))*180.0/PI
         write(nud,234) jlat,zalat,csq(jlat),gwd(jlat)
      end do 
      write(nud,232)
      write(nud,231)
      return
  231 format(/)
  232 format(37('*'))
  233 format('*  No *   Lat *       csq    weight *')
  234 format('*',i4,' *',f6.1,' *',2f10.4,' *')
      end subroutine legpri


!     =================
!     SUBROUTINE INILAT
!     =================

      subroutine inilat
      use pumamod
      implicit none

      real (kind=8) :: zcsq
      integer :: jlat, ideg

      do jlat = 1 , NLAT
         zcsq       = 1.0 - sid(jlat) * sid(jlat)
         csq(jlat)  = zcsq
         rcs(jlat)  = 1.0 / sqrt(zcsq)
      enddo
      do jlat = 1 , NLAT/2
         ideg = nint(180.0/PI * asin(sid(jlat)))
         write(chlat(jlat),'(i2,a1)') ideg,'N'
         write(chlat(NLAT+1-jlat),'(i2,a1)') ideg,'S'
      enddo
      return
      end subroutine inilat


!     ====================
!     SUBROUTINE GRIDPOINT
!     Key callings:
!     - calcgp              : 超级密集计算               纯数组计算 没有调用任何其它函数或子程序
!     - mktend (mod_legsym) : Compute Nonlinear terms    纯数组计算 没有调用任何其它函数或子程序 但"use legsym"
!
!     Others: <--- All neglected for OpenACC testing
!     - filter_zonal_waves  : 纬向滤波
!     - stepruido           : 运行中增加随机性
!       gasdev                stepruido调用了这个Gauss随机数发生器
!     - altcs               : 诊断cs切片时转换alt grid
!
!     MPI: mpsumsc; mpgagp; mpsum; mpgacs <--- All removed
!
!     本质上此子程序的输入: sd,st,sz,sp
!           中间计算并修改: gd,gt,gz,gp,gpj,gu,gv,gut,gvt,gke,gfu,gfv
!      最终得到4变量的倾向: spf,sdf,szf,stf (local)
!             再复制给全局: spt,sdt,szt,stt (pumamod)
!     ====================

      subroutine gridpoint
      use pumamod, only: sd,st,sz,sp, sdt,stt,szt,spt,   &  ! 前4个是总体的最根本input; 后四个是最根本output
                         gd,gt,gz,gp,gpj,gu,gv,          &  ! 这些变量都会被update
                         gut,gvt,gke,gfu,gfv,            &  ! 这些变量都会被update
                         NLON,NLAT,NLEV,NLPP,NHOR,NRSP,NESP,trigs
      implicit none

      ! 在calcgp中重要的传递变量
      real, dimension(NLON,NLPP,NLEV) :: gtn    ! tempurature
      real, dimension(NHOR)           :: gvpp   ! 动量的垂直积分 最终决定Ps的倾向spf
      real, dimension(NLON,NLPP)      :: gpmt   ! Ps

      ! 此子程序最终计算出的4变量的倾向 最终要赋值给sdt,stt,szt,spt
      real, dimension(NESP,NLEV)      :: sdf    ! tendency of div
      real, dimension(NESP,NLEV)      :: stf    ! tendency of temperature
      real, dimension(NESP,NLEV)      :: szf    ! tendency of vor
      real, dimension(NESP)           :: spf    ! tendency of surface pressure

      ! 用于最后diag时用的
      real, dimension(NLON,NLAT)      :: zgp
      real, dimension(NHOR)           :: zgpp
      real, dimension(NLAT,NLEV)      :: zcs  !XW(2017-4-7): remove "(kind=4)" for zcs and zsp
      real, dimension(NRSP)           :: zsp

      real :: sec
      integer :: jlon, jlat, jlev

      !$OMP parallel

      !$OMP sections
      !$OMP section
      do jlev = 1 , NLEV
         call sp2fc(sd(:,jlev),gd(:,jlev))
      end do
      !$OMP section
      do jlev = 1 , NLEV
         call sp2fc(st(:,jlev),gt(:,jlev))
      end do
      !$OMP section
      do jlev = 1 , NLEV
         call sp2fc(sz(:,jlev),gz(:,jlev))
      enddo

      !$OMP section
      call sp2fc(sp,gp)                                          ! lnPs
      call sp2fcdmu(sp,gpj)                                      ! d(lnPs) / d(mu)
      do jlev = 1 , NLEV
         call dv2uv(sd(:,jlev),sz(:,jlev),gu(:,jlev),gv(:,jlev)) ! div,vor->ucos(phi),vcos(phi)
      enddo
      !$OMP end sections

      !---$OMP single

      ! 纬向FFT滤波
      !if (lselect) then
      !   call filter_zonal_waves(gp)
      !   call filter_zonal_waves(gpj)
      !   do jlev = 1 , NLEV
      !      call filter_zonal_waves(gu(:,jlev))
      !      call filter_zonal_waves(gv(:,jlev))
      !      call filter_zonal_waves(gd(:,jlev))
      !      call filter_zonal_waves(gt(:,jlev))
      !      call filter_zonal_waves(gz(:,jlev))
      !   enddo
      !endif

      ! 可暂时关闭diag信息
      ! 取cross-section: gu,gv,gt每lev每lat第一个经度
      !if (mod(nstep,ndiag) == 0) then
      !  do jlev = 1 , NLEV
      !    do jlat = 1 , NLPP
      !      sec = cv / sqrt(csq(jlat))
      !      csu(jlat,jlev) = gu(1+(jlat-1)*NLON,jlev) * sec
      !      csv(jlat,jlev) = gv(1+(jlat-1)*NLON,jlev) * sec
      !      cst(jlat,jlev) =(gt(1+(jlat-1)*NLON,jlev) + t0(jlev))*ct-273.16
      !    enddo
      !  enddo
      !endif

      !---$OMP end single

      !$OMP do collapse(2)
      do jlat = 1 , NLPP
         do jlon = 1 , NLON-1 , 2
           gpmt(jlon  ,jlat) = -gp(jlon+1+(jlat-1)*NLON) * ((jlon-1)/2)
           gpmt(jlon+1,jlat) =  gp(jlon  +(jlat-1)*NLON) * ((jlon-1)/2)
         end do
      end do
      !$OMP end do

      !$OMP sections
      !$OMP section
      call fc2gp(gu ,NLON,NLPP*NLEV,trigs)
      !$OMP section
      call fc2gp(gv ,NLON,NLPP*NLEV,trigs)
      !$OMP section
      call fc2gp(gt ,NLON,NLPP*NLEV,trigs)
      !$OMP section
      call fc2gp(gd ,NLON,NLPP*NLEV,trigs)
      !$OMP section
      call fc2gp(gz ,NLON,NLPP*NLEV,trigs)
      !$OMP section
      call fc2gp(gpj,NLON,NLPP,trigs)
      !$OMP section
      call fc2gp(gpmt,NLON,NLPP,trigs)
      !$OMP end sections

      !$OMP single
      call calcgp(gtn,gpmt,gvpp)

      gut(:,:) = gu(:,:) * gt(:,:)
      gvt(:,:) = gv(:,:) * gt(:,:)
      gke(:,:) = gu(:,:) * gu(:,:) + gv(:,:) * gv(:,:)
      !$OMP end single

      !$OMP sections
      !$OMP section
      call gp2fc(gtn ,NLON,NLPP*NLEV,trigs)
      !$OMP section
      call gp2fc(gut ,NLON,NLPP*NLEV,trigs)
      !$OMP section
      call gp2fc(gvt ,NLON,NLPP*NLEV,trigs)
      !$OMP section
      call gp2fc(gfv ,NLON,NLPP*NLEV,trigs)
      !$OMP section
      call gp2fc(gfu ,NLON,NLPP*NLEV,trigs)
      !$OMP section
      call gp2fc(gke ,NLON,NLPP*NLEV,trigs)
      !$OMP section
      call gp2fc(gvpp,NLON,NLPP     ,trigs)
      !$OMP end sections

      !$OMP single
      call fc2sp(gvpp,spf)

      ! 纬向滤波
      !if (lselect) then
      !   call filter_zonal_waves(gvpp)
      !   do jlev = 1 , NLEV
      !      call filter_zonal_waves(gtn(:,:,jlev))
      !      call filter_zonal_waves(gut(:,jlev))
      !      call filter_zonal_waves(gvt(:,jlev))
      !      call filter_zonal_waves(gfv(:,jlev))
      !      call filter_zonal_waves(gfu(:,jlev))
      !      call filter_zonal_waves(gke(:,jlev))
      !   enddo
      !endif
      !$OMP end single

      !$OMP do
      do jlev = 1 , NLEV
         call mktend(sdf(:,jlev),stf(:,jlev),szf(:,jlev),   &  ! these are THE 3 intent(out) variables
                     gtn(:,:,jlev),gfu(:,jlev),gfv(:,jlev), &  ! there are 6 intent(in) variables 
                     gke(:,jlev),gut(:,jlev),gvt(:,jlev))
      enddo
      !$OMP end do
      !$OMP end parallel

      ! 实际上未加 暂时注释掉
      !if (nruido > 0) call stepruido      ! 在运行中是否还要加随机扰动? 这里计算ruido数组

      ! 这是关键--要传递给spectral的4变量倾向
      spt   = spf    !call mpsumsc(spf,spt,1)
      sdt   = sdf    !call mpsumsc(sdf,sdt,NLEV)
      szt   = szf    !call mpsumsc(szf,szt,NLEV)
      stt   = stf    !call mpsumsc(stf,stt,NLEV)

! 为OpenACC测试暂时关闭 打开第一列叹号就是原来的优化带注释
!
!      if (mod(nstep,ndiag) == 0) then
!         call fc2gp(gp,NLON,NLPP,trigs)
!         zgpp(:) = exp(gp)                ! LnPs -> Ps
!         !call mpgagp(zgp,zgpp,1)         ! zgp = Ps (full grid)
!         zgp = reshape(zgpp,(/NLON,NLAT/))! zgp = Ps (full grid)
!      !XW(Mar/25/2017) to remove GUI:
!         !if (ngui > 0) then
!         !   call guips(zgp,psmean)        
!         !   call guigv("GU" // char(0),gu)
!         !   call guigv("GV" // char(0),gv)
!         !   call guigt(gt)
!         !endif
!         zgpp(:) =  zgpp(:) - 1.0         ! Mean(LnPs) = 0  <-> Mean(Ps) = 1
!         call gp2fc(zgpp,NLON,NLPP,trigs)
!         call fc2sp(zgpp,span)
!
!         !call mpsum(span,1)              ! span = Ps spectral
!         !call mpgacs(csu)
!         !call mpgacs(csv)
!         !call mpgacs(cst)
!         !if (mypid == NROOT) then
!            call altcs(csu)
!            call altcs(csv)
!            call altcs(cst)
!      !XW(Mar/25/2017) to remove GUI:
!            !if (ngui > 0) then
!            !   zcs(:,:) = csu(:,:)
!            !   call guiput("CSU"  // char(0) ,zcs ,NLAT,NLEV,1)
!            !   zcs(:,:) = csv(:,:)
!            !   call guiput("CSV"  // char(0) ,zcs ,NLAT,NLEV,1)
!            !   zcs(:,:) = cst(:,:)
!            !   call guiput("CST"  // char(0) ,zcs ,NLAT,NLEV,1)
!            !   zsp(:) = span(1:NRSP)
!            !   call guiput("SPAN" // char(0) ,zsp ,NCSP,-NTP1,1)
!            !endif
!         !endif
!      endif
      return
      end subroutine gridpoint


!     =================
!     SUBROUTINE CALCGP
!     高密度计算模块
!     深挖scalable computing
!     纯数组计算 没有调用任何其它函数和子程序
!     intent(in   ) : gpm,gpj, gu,gv,gd,gz,gt
!     intent(inout) : gfu,gfv
!     intent(  out) : gtn,gvp
!     =================
      subroutine calcgp(gtn,gpm,gvp)            ! gpm进 gtn/gvp出
      use pumamod, only: gfu, gfv,           &  ! 注意这两个变量也要被修改
                         gu,gv,gd,gz,gt,gpj, &  ! Others just intent(in)
                         t0d,tkp,rdsig,rcsq,c,dsigma,sigmh,akap,ruidop, NLEV,NLEM,NHOR,NRUIDO
      implicit none

!     Comments by Torben Kunz and Guido Schroeder

!     Compute nonlinear tendencies in grid point space.
!     Hoskins and Simmons 1975 (Q.J.R.Meteorol.Soc.,101,637-655) (HS75)

!     For terms calculated in this routine, see HS75, eqs. (8)-(10) and
!     appendix I:
!     - script Fu, Fv as contributions to script D: gl. arrays gfu, gfv
!     - script T: returned as gtn
!     - script P: returned as gvp


!     parameters (in)
!     ---------------

!     gpm --              d(ln(ps)) / d(lambda)

!     parameters (out)
!     ---------------

!     gtn -- temperature tendency
!     gvp -- vertical integral of (u,v) * grad(ln(ps))

!     global arrays variable in time
!     ------------------------------

!     gfu, gfv -- terms Fu, Fv in primitive equations,
!                 see HS75 (eqs. (1), (2))
!     gu, gv   -- components u, v of horizontal velocity vector
!     gd       -- divergence D
!     gz       -- absolute vorticity
!     gt       -- temperature deviation T'

!     global arrays constant in time
!     ------------------------------

!     t0d   -- reference temperature difference between two adjacent
!              full levels
!     tkp   -- reference temperature times kappa (global parameter AKAP)
!     rdsig -- 1 / (2 * dsigma)
!     rcsq  -- 1 / (1 - mu^2) 

!     notations used in subsequent comments
!     -------------------------------------

!     aINTb(A)dsigma :<=> the integral of A over the interval [a,b]
!                         with respect to sigma

      real, intent(out) :: gtn(NHOR,NLEV)
      real, intent(in ) :: gpm(NHOR)
      real, intent(out) :: gvp(NHOR)

      ! local
      real zsdotp(NHOR,NLEM),zsumd(NHOR),zsumvp(NHOR),zsumvpm(NHOR)
      real ztpta(NHOR),ztptb(NHOR)
      real zvgpg(NHOR,NLEV)
      real gtd(NHOR,NLEM)
      real gud(NHOR,NLEM)
      real gvd(NHOR,NLEM)

      integer :: jlev,jlej

!$acc kernels

!     1.
!     1.1 zvgpg: (u,v) * grad(ln(ps))

      do jlev = 1 , NLEV
         zvgpg(:,jlev) = rcsq  * (gu(:,jlev)*gpm(:)+gv(:,jlev)*gpj(:))
      enddo

!     1.2 Calculate vertical integral of A = D + (u,v) * grad(ln(ps)),
!         separated into divergence and ln(ps) advection.
!         zsumd  : 0INT1(D)dsigma
!         gvp    : 0INT1[(u,v) * grad ln(ps)]dsigma
!         zsdotp : 0INTsigma(A)dsigma

      zsumd = dsigma(1) * gd(:,1)
      gvp   = dsigma(1) * zvgpg(:,1)
      zsdotp(:,1) = zsumd + gvp

      do jlev = 2 , NLEM
         zsumd = zsumd + dsigma(jlev) * gd(:,jlev)
         gvp   = gvp   + dsigma(jlev) * zvgpg(:,jlev)
         zsdotp(:,jlev) = zsumd + gvp
      enddo

      zsumd = zsumd + dsigma(NLEV) * gd(:,NLEV)
      gvp   = gvp   + dsigma(NLEV) * zvgpg(:,NLEV)

!     2. Calculate vertical velocity and vertical advection terms
!        on half levels.

      do jlev = 1 , NLEM
         zsdotp(:,jlev) = (sigmh(jlev) * (zsumd+gvp) - zsdotp(:,jlev))
      enddo

      gtd(:,:) = zsdotp(:,:) * (gt(:,2:NLEV) - gt(:,1:NLEM))
      gud(:,:) = zsdotp(:,:) * (gu(:,2:NLEV) - gu(:,1:NLEM))
      gvd(:,:) = zsdotp(:,:) * (gv(:,2:NLEV) - gv(:,1:NLEM))

!     3. Calculate nonlinear contributions to temperature tendency and
!        nonlinear terms Fu, Fv as used in vorticity and
!        divergence equation.

!     3.1 top level:

!     3.1.1 zsumvp: 0INTsigma[(u,v) * grad(ln(ps))]dsigma

      zsumvp = zvgpg(:,1) * dsigma(1)

!     3.1.2 Calculation of gtn, gfv and gfu as for inner levels (3.2),
!           but somewhat simplified:
!           a) For the top level the following equation holds in the
!              discretized form: (1/sigma)*0INTsigma(A)dsigma == A
!              (HS75, second equation following eq. (7)). Therefore,
!              (3.2.3) simplifies to -kappa*T' * D and (3.2.4) vanishes.
!           b) Vertical advection terms (gtd, gud, gvd (see section 2)
!              and vertical T0 advection (3.2.6)) vanish at upper
!              boundary (sigma == 0).

      gtn(:,1) = (1.0-akap) * gt(:,1) * gd(:,1) - rdsig(1) * (gtd(:,1) &
               + t0d(1) * (sigmh(1)*gvp-zsumvp))

      gfv(:,1) = - gu(:,1)*gz(:,1) - gpj(:)*gt(:,1) - rdsig(1)*gvd(:,1)
      gfu(:,1) =   gv(:,1)*gz(:,1) - gpm(:)*gt(:,1) - rdsig(1)*gud(:,1)

!     3.2 inner levels:

      do jlev = 2 , NLEM

!        3.2.1 ztpta: (1/sigma)*0INTsigma(A-D)dsigma
!              ztptb: (1/sigma)*0INTsigma(A)dsigma
!              Matrix c contains factors for discretized integration, see
!              HS75 (second equation following eq. (7)).

         ztpta = c(1,jlev) *  zvgpg(:,1)
         ztptb = c(1,jlev) * (zvgpg(:,1) + gd(:,1))

         do jlej = 2 , jlev
            ztpta = ztpta + c(jlej,jlev) *  zvgpg(:,jlej)
            ztptb = ztptb + c(jlej,jlev) * (zvgpg(:,jlej) + gd(:,jlej))
         enddo

         zsumvpm = zsumvp
         zsumvp  = zsumvp + zvgpg(:,jlev) * dsigma(jlev)

!        3.2.2 D * T' 

         gtn(:,jlev) = gt(:,jlev) * gd(:,jlev)

!        3.2.3 kappa*T' *
!                    [(u,v)*grad(ln(ps)) - (1/sigma)*0INTsigma(A)dsigma]

         gtn(:,jlev) = gtn(:,jlev)                                      &
     &       + akap * gt(:,jlev) * (zvgpg(:,jlev) - ztptb)

!        3.2.4 kappa*T0 *
!                  [(u,v)*grad(ln(ps)) - (1/sigma)*0INTsigma(A-D)dsigma]

         gtn(:,jlev) = gtn(:,jlev)                                      &
     &       + tkp(jlev) * (zvgpg(:,jlev) - ztpta)

!        3.2.5 Calculate vertical T' advection on full levels by
!              averaging two half level advection terms (gtd, calculated
!              in section 2).

!        and

!        3.2.6 Calculate vertical T0 advection on full levels by
!              averaging two half level advection terms.

         gtn(:,jlev) = gtn(:,jlev)                                      &
     &       - rdsig(jlev) * (gtd(:,jlev) + gtd(:,jlev-1)               &
     &         +(sigmh(jlev)   * gvp - zsumvp)  * t0d(jlev)            &
     &         +(sigmh(jlev-1) * gvp - zsumvpm) * t0d(jlev-1))

!        3.2.7 terms Fv, Fu, see HS75 (equations following eq. (5));
!              vertical advection terms interpolated to full levels by
!              averaging two half level advection terms.

         gfv(:,jlev) = - gu(:,jlev)*gz(:,jlev) - gpj(:)*gt(:,jlev)      &
     &                 - rdsig(jlev)*(gvd(:,jlev) + gvd(:,jlev-1))

         gfu(:,jlev) =   gv(:,jlev)*gz(:,jlev) - gpm(:)*gt(:,jlev)      &
     &                 - rdsig(jlev)*(gud(:,jlev) + gud(:,jlev-1))
      enddo

!     3.3 bottom level

!     3.3.1 ztpta, ztptb: see 3.2.1

      ztpta = c(1,NLEV) *  zvgpg(:,1)
      ztptb = c(1,NLEV) * (zvgpg(:,1) + gd(:,1))

      do jlej = 2 , NLEV
         ztpta = ztpta + c(jlej,NLEV) *  zvgpg(:,jlej)
         ztptb = ztptb + c(jlej,NLEV) * (zvgpg(:,jlej) + gd(:,jlej))
      enddo

!     3.3.2 Calculation of gtn, gfv and gfu as for inner levels (3.2),
!           but somewhat simplified:
!           Vertical advection terms (gtd, gud, gvd (see section 2) and 
!           vertical T0 advection (3.2.6)) vanish at
!           lower boundary (sigma == 1).

      gtn(:,NLEV) = gt(:,NLEV) * gd(:,NLEV)                            &
     &            + akap*gt(:,NLEV)*(zvgpg(:,NLEV)-ztptb)              &
     &            + tkp(NLEV)*(zvgpg(:,NLEV)-ztpta)                    &
     &            - rdsig(NLEV) * (gtd(:,NLEM)                         &
     &            + t0d(NLEM)*(sigmh(NLEM)*gvp-zsumvp))

      gfv(:,NLEV) = -gu(:,NLEV) * gz(:,NLEV) - gpj(:) * gt(:,NLEV)      &
     &              - rdsig(NLEV) * gvd(:,NLEM)
      gfu(:,NLEV) =  gv(:,NLEV) * gz(:,NLEV) - gpm(:) * gt(:,NLEV)      &
     &              - rdsig(NLEV) * gud(:,NLEM)

!     3.3.3 Add gaussian noise to T (controlled by nruido)

      if (nruido > 0) gtn(:,:) = gtn(:,:) + ruidop(:,:)

!$acc end kernels

      end subroutine calcgp


!     ===================
!     SUBROUTINE SPECTRAL
!     Calling tree:
!     - filter_spectral_modes : SP滤波
!     - mkenerdiag            : 诊断Energy
!     - mkentrodiag           : 诊断Entropy
!     - heatgp                : GP空间中 运行中时间维的加热 理想的Newtonian Cooling
!     - diagp                 : GP空间中 运行中时间维的加热 来自文件的读入数据
!     - vdiff                 : 计算垂直方向的耗散 (div, vor, T)
!     - mkdheat               : 对耗散动能的能量回收，加热周围环境
!     MPI: mpgallsp, mpsumbcr, mrdiff
!
!     Key Info
!     Input  : intent(inout)  : 倾向spt,sdt,szt,stt; 前一时刻spm,sdm,szm,stm; 当前时刻spp,sdp,szp,stp
!              intent(in   )  : 其它pumamod引入的变量都只是引用 不做任何修改
!     Output : intent(out  )  : 新时刻sp,sd,sz,st <--- 就是最新的当前时刻spp,sdp,szp,stp
!     ===================

      subroutine spectral
      use pumamod, only: sp,sd,sz,st, spt,sdt,szt,stt, spm,sdm,szm,stm, spp,sdp,szp,stp, &
                         sak,sop,srp1,srp2,nindex,srcn,damp,fric,t0,xlphi,xlt,bm1,dsigma,&
                         nstep,pac,tac,nhelsua,plavor,pnu,alpha,delt,delt2,mloop, NSPP,NLEV
      implicit none

!*    Add adiabatic and diabatic tendencies - perform leapfrog

!     The adiabatic tendencies are added using the semi implicit scheme
!     Hoskins and Simmons 1975 (Q.J.R.Meteorol.Soc.,101,637-655) (HS75)
!     To compare the code directly with HS75 the following notes might
!     be helpful (in addition to the comments below):

!     Name rule for global arrays <abc>:
!     a : representation (s=spectral, g=grid, z=local)
!     b : variable (p=ln(ps), d=divergence, z=vorticity, t=temperature)
!     c : modifier (m=previous timestep, p=present timestep, t=tendency)

!     global arrays variable in time
!     ------------------------------

!     spt - pressure    tendency HS75 (10)
!     sdt - divergence  tendency HS75 ( 8)
!     szt - vorticity   tendency
!     stt - temperature tendency HS75 ( 9)

!     spm - pressure    at previous timestep
!     sdm - divergence  at previous timestep
!     szm - vorticity   at previous timestep
!     stm - temperature at previous timestep

!     spp - pressure    at present timestep
!     sdp - divergence  at present timestep
!     szp - vorticity   at present timestep
!     stp - temperature at present timestep

!     global arrays constant in time
!     ------------------------------

!     sak(NSPP)      -     = hyper diffusion
!     sop(NSPP)      - g*  = orography as geopotential
!     srp1(NSPP,NLEV) - Tr = radiative equilibrium temperature (annual mean)
!     srp2(NSPP,NLEV) - Tr = radiative equilibrium temperature (annual cycle)
!     nindex(NSPP)   - n   = total wavenumber n for spectral modes
!     srcn(NSPP)   - 1/Cn  = 1.0 / (n * (n+1))
!     damp(NLEV)  1/tau R  = time constant for newtonian cooling
!     fric(NLEV)  1/tau F  = time constant for Rayleigh friction

      real zpm(NSPP)      ! new spm
      real zdm(NSPP,NLEV) ! new sdm
      real zzm(NSPP,NLEV) ! new szm
      real ztm(NSPP,NLEV) ! new stm

      real zwp(NSPP)      ! timefilter delta pm
      real zwd(NSPP,NLEV) ! timefilter delta sd
      real zwz(NSPP,NLEV) ! timefilter delta sz
      real zwt(NSPP,NLEV) ! timefilter delta st

      real zsrp(NSPP)     ! restoring temperature (mean + annual cycle)
                          ! 此变量在计算牛顿冷却时用的临时变量

      real zgt(NSPP,NLEV) ! work array 临时工作数组

      ! For OpenACC: 以下变量似乎都没有使用
      real,allocatable :: zstte(:,:,:) ! temp. tendencies for energy diag.
      real,allocatable :: zszte(:,:,:) ! vort. tendencies for energy recycling
      real,allocatable :: zsdte(:,:,:) ! div. tendencies for energy recycling
      real,allocatable :: zdps(:)      ! surf pressure for energy diag

      real,allocatable :: zsp(:)       ! surf pressure spectral
      real,allocatable :: zspf(:)      ! surf pressure spectral
      real,allocatable :: zspt(:)      ! surf pressure tendency 

      real,allocatable :: zst(:,:)     ! temperature for entropy diagnostics
      real,allocatable :: zstt(:,:)    ! tem. tendencies for entropy diag.
      real,allocatable :: ztgp(:,:)    ! 
      real,allocatable :: zdtgp(:,:)   ! 
      real,allocatable :: zsum1(:)
      real,allocatable :: zgw(:)

      ! Xinyu added
      real     :: z0, za, zb, zm, zq
      real     :: zampl, zcp, zsum3, zt, ztp, zztm
      integer  :: jhor, jsp, jn, jlon, jlat, jlev

!$acc kernels

!     0. Special code for experiments with mode filtering

      ! 暂时关闭 for OpenACC testing
      !if (lspecsel) call filter_spectral_modes

!     1. Initialize local arrays

      zpm(:)   = spp(:)
      zdm(:,:) = sdp(:,:)
      zzm(:,:) = szp(:,:)
      ztm(:,:) = stp(:,:)
!
!     allocate diagnostic arrays if needed
!

      ! 暂时关闭 for OpenACC testing
      ! checked, they are ALL 0 :-)
      !if(nenergy > 0 .or. nentropy > 0 .or. ndheat > 0) then
      ! allocate(zstte(NSPP,NLEV,3))
      !endif
      !if(ndheat > 0) then
      ! allocate(zszte(NSPP,NLEV,2))
      ! allocate(zsdte(NSPP,NLEV,2))
      !endif
!
!     allocate and compute surface pressure if needed
!
      ! 暂时关闭 for OpenACC testing
      ! checked, they are ALL 0 :-)
      !if(nenergy > 0 .or. nentropy > 0 .or. ndheat > 0) then
      ! allocate(zspt(NSPP))
      ! allocate(zsp(NSPP))
      !endif

!     2. Calculate divergence on timelevel t (sdt) HS75 (17)
!        which will replace the divergence tendency sdt
!        (semi implicit scheme)

!        The vertical scheme has being changed to the ECMWF scheme
!        (see e.g. Simmons and Burridge 1981, Mon.Wea.Rev.,109,758-766).
!        in this scheme, matrix xlphi (g) differs from that in HS75.

!         z0 : reference temperature To
!         zq : 1.0 / Cn
!         zt : xlphi * script T - To * script P
!         zm : xlphi * T        + To * ln(Ps)(t-dt)

!         (note that phi is needed in HS75 (17) and, therefore,
!         the surface geopotential phi* [sop] is added

      do jlev=1,NLEV
         z0 = t0(jlev)
         do jsp=1,NSPP
            zq = srcn(jsp)                   ! 1.0 / (n * (n + 1))
            zt = dot_product(xlphi(:,jlev),stt(jsp,:)) - z0 * spt(jsp)
            zm = dot_product(xlphi(:,jlev),stm(jsp,:)) + z0 * spm(jsp)
            za = sdt(jsp,jlev) * zq
            zb = sdm(jsp,jlev) * zq
            zgt(jsp,jlev) = zb + delt * (za + zm + sop(jsp) + zt * delt)
         enddo
      enddo

!     bm1 is the invers of matrix (1/cn I+B dt**2) (lhs HS75 (17))

      do jlev = 1 , NLEV
         do jsp = 1 , NSPP
            jn = nindex(jsp)            ! total wavenumber n
            sdt(jsp,jlev) = dot_product(zgt(jsp,:),bm1(:,jlev,jn))
         enddo
      enddo

!     3. Calculate surface pressure tendency -ln(ps) HS75 (15)

      do jlev = 1 , NLEV
         spt(:) = spt(:) + dsigma(jlev) * sdt(:,jlev)
      enddo

!     4. Calculate temperature tendency   HS75 (14)

      do jlev = 1 , NLEV
       do jsp = 1 , NSPP
        stt(jsp,jlev)=stt(jsp,jlev)-dot_product(xlt(:,jlev),sdt(jsp,:))
       enddo
      enddo

!     5. Add tendencies

      spp(:)   = spm(:) - delt2 * spt(:)     ! ln(ps) :  spt = -ln(ps) tendency
      sdp(:,:) =   2.0 * sdt(:,:) - sdm(:,:) ! div    :  sdt = sdm + delt * tend.
      szp(:,:) = delt2 * szt(:,:) + szm(:,:) ! vor
      stp(:,:) = delt2 * stt(:,:) + stm(:,:) ! t

      ! 暂时关闭 for OpenACC testing
      ! checked, they are ALL 0 :-)
      !if(nenergy > 0) then
      ! zspt(:)=-spt(:)
      ! call mkenerdiag(stm,stt,spm,zspt,denergy(:,1))
      !endif
      !if(nentropy > 0) then
      ! call mkentrodiag(stm,stt,spm,dentropy(:,1))
      !endif

!     6. Calculate newtonian cooling, friction and biharmonic diffusion
!        (srp - stp) * damp = (Tr' -T') / tau R = newtonian cooling
!        srp1 = annual mean  component
!        srp2 = annual cycle component
!        sak  = diffusion
!        fric = friction
!        zampl = annual cycle

      zampl = cos((real(nstep)-pac)*tac)

      if (nhelsua == 0 .or. nhelsua == 1) then
         do jlev=1,NLEV
            zsrp(:)=srp1(:,jlev)+srp2(:,jlev)*zampl
            sdt(:,jlev) =  sdp(:,jlev) * (sak(1:NSPP) - fric(jlev))
            szt(:,jlev) =  szp(:,jlev) * (sak(1:NSPP) - fric(jlev))
            stt(:,jlev) = (zsrp(:) - stp(:,jlev)) * damp(jlev) + stp(:,jlev) * sak(1:NSPP)
      ! 暂时关闭 for OpenACC testing
      ! checked, nenergy = ndheat = ndiagp = 0
      !      if(nenergy > 0) then
      !       zstte(:,jlev,2)=(zsrp(:)-stp(:,jlev))*damp(jlev)
      !       zstte(:,jlev,3)=stp(:,jlev)*sak(1:NSPP)
      !      endif
      !      if(ndheat > 0) then
      !       zsdte(:,jlev,1) =  -sdp(:,jlev) * fric(jlev)
      !       zszte(:,jlev,1) =  -szp(:,jlev) * fric(jlev)
      !       zsdte(:,jlev,2) =  sdp(:,jlev) * sak(1:NSPP)
      !       zszte(:,jlev,2) =  szp(:,jlev) * sak(1:NSPP)
      !      endif
         enddo
      !elseif (nhelsua == 2 .or. nhelsua == 3 .or. ndiagp > 0) then     ! checked: ndiagp=0
      !   if (ndiagp == 0) then
      !      call heatgp(zampl)  ! stt(:,:) = Newtonian cooling
      !   else
      !      call diagp(zampl)  ! stt(:,:) = Newtonian cooling
      !   endif
      !   if(nenergy > 0) then
      !    zstte(:,:,2)=stt(:,:)
      !   endif
      !   do jlev=1,NLEV
      !      sdt(:,jlev) = sdp(:,jlev) * (sak(1:NSPP) - fric(jlev))
      !      szt(:,jlev) = szp(:,jlev) * (sak(1:NSPP) - fric(jlev))
      !      stt(:,jlev) = stt(:,jlev) + stp(:,jlev) * sak(1:NSPP)
      !      if(nenergy > 0) then
      !       zstte(:,jlev,3)=stp(:,jlev)*sak(1:NSPP)
      !      endif
      !      if(ndheat > 0) then
      !       zsdte(:,jlev,1) =  -sdp(:,jlev) * fric(jlev)
      !       zszte(:,jlev,1) =  -szp(:,jlev) * fric(jlev)
      !       zsdte(:,jlev,2) =  sdp(:,jlev) * sak(1:NSPP)
      !       zszte(:,jlev,2) =  szp(:,jlev) * sak(1:NSPP)
      !      endif
      !   enddo
      endif

!        Conserve ln(ps) by forcing mode(0,0) to zero
!        Correct vorticity by canceling the friction and diffusion
!        applied to planetary vorticity
!        Only root node processes the first NSPP modes

      ! 暂时关闭 for OpenACC testing
      ! checked, they are ALL 0 :-)
      !if(nenergy > 0) then
      ! zspt(:)=0.
      ! call mkenerdiag(stp,zstte(:,:,2),spp,zspt,denergy(:,2))
      ! call mkenerdiag(stp,zstte(:,:,3),spp,zspt,denergy(:,3))
      !endif
      !if(nentropy > 0) then
      ! call mkentrodiag(stp,zstte(:,:,2),spp,dentropy(:,2))
      ! call mkentrodiag(stp,zstte(:,:,3),spp,dentropy(:,3))
      !endif
      !if(nenergy > 0 .or. nentropy > 0 .or. ndheat > 0) then
      ! zsp(:)=spp(:)
      ! zstte(:,:,1)=stt(:,:)
      !endif

      !if (mypid == NROOT) then
         spp(1) = 0.0
         spp(2) = 0.0
         szt(3,:) = szt(3,:) + plavor * (fric(:) - sak(3))
      ! 暂时关闭 for OpenACC testing
      ! checked, ndheat=0 :-)
      !   if(ndheat > 0) then
      !    zszte(3,:,1) = zszte(3,:,1) + plavor * fric(:)
      !    zszte(3,:,2) = zszte(3,:,2) - plavor * sak(3)
      !   endif
      !endif
!
!     6b) call for vertical diffusion
!

      ! 暂时关闭 for OpenACC testing
      ! checked, dvdiff=0.0 :-)
      ! 但是我认为这个if语句的条件很成问题 fix it later
      !if(dvdiff > 0.0) call vdiff(stp,szp,sdp,stt,szt,sdt)

!
!     recycle kin energy dissipation
! 

      ! 暂时关闭 for OpenACC testing
      ! checked, ndheat=nenergy=nentropy=0
      !if(ndheat > 0) then
      ! call mkdheat(zszte(:,:,1),zszte(:,:,2),     &
      !              zsdte(:,:,1),zsdte(:,:,2),zsp)
      !endif


      !if(nenergy > 0 .or. nentropy > 0) then
      ! zstte(:,:,1)=stt(:,:)-zstte(:,:,1)
      !endif
      !if(nenergy > 0) then
      ! call mkenerdiag(stp,zstte(:,:,1),zsp,zspt,denergy(:,4))
      !endif
      !if(nentropy > 0) then
      ! zstte(:,:,1)=stt(:,:)-zstte(:,:,1)
      ! call mkentrodiag(stp,zstte(:,:,1),zsp,dentropy(:,4))
      !endif
      !if(nenergy > 0 .or. nentropy > 0) then
      ! zstte(:,:,1)=0.
      ! zspt(:)=(spp(:)-zsp(:))/delt2
      !endif
      !if(nenergy > 0) then
      ! call mkenerdiag(stp,zstte(:,:,1),zsp,zspt,denergy(:,8))
      !endif
      !if(nentropy > 0) then
      ! call mkentrodiag(stp,zstte(:,:,1),zsp,dentropy(:,8))
      !endif

!
!     diagnostics of efficiency
!

   ! 暂时关闭 for OpenACC testing
   ! checked, ndheat=0
   !   if(ndheat > 1) then
   !    zcp=gascon/akap
   !    allocate(zst(NESP,NLEV))
   !    allocate(zstt(NESP,NLEV))
   !    allocate(zspf(NESP))
   !    allocate(ztgp(NHOR,NLEV))
   !    allocate(zdtgp(NHOR,NLEV))
   !    allocate(zdps(NHOR))
   !    allocate(zsum1(4))
   !    allocate(zgw(NHOR))
   !    jhor=0
   !    do jlat=1,NHPP
   !     do jlon=1,NLON*2
   !      jhor=jhor+1
   !      zgw(jhor)=gwd(jlat)
   !     enddo
   !    enddo
   !    zst  = stp    !call mpgallsp(zst,stp,NLEV)
   !    zstt = stt    !call mpgallsp(zstt,stt,NLEV)
   !    zspf = zsp    !call mpgallsp(zspf,zsp,1)
   !    do jlev = 1 , NLEV
   !       call sp2fc(zst(1,jlev),ztgp(1,jlev))
   !       call sp2fc(zstt(1,jlev),zdtgp(1,jlev))
   !    enddo
   !    call sp2fc(zspf,zdps)
   !    call fc2gp(ztgp,NLON,NLPP*NLEV,trigs)
   !    call fc2gp(zdtgp,NLON,NLPP*NLEV,trigs)
   !    call fc2gp(zdps,NLON,NLPP,trigs)
   !    zdps(:)=psurf*exp(zdps(:))
   !    zsum1(:)=0.
   !    do jlev=1,NLEV
   !     ztgp(:,jlev)=ct*(ztgp(:,jlev)+t0(jlev))
   !     zdtgp(:,jlev)=ct*ww_scale*zdtgp(:,jlev)
   !     zsum1(1)=zsum1(1)+SUM(zdtgp(:,jlev)*zgw(:)                      &
   !  &                       *zcp*zdps(:)/ga*dsigma(jlev)               &
   !  &                       ,mask=(zdtgp(:,jlev) >= 0.))
   !     zsum1(2)=zsum1(2)+SUM(zdtgp(:,jlev)*zgw(:)                      &
   !  &                       *zcp*zdps(:)/ga*dsigma(jlev)               &
   !  &                       ,mask=(zdtgp(:,jlev) < 0.))
   !     zsum1(3)=zsum1(3)+SUM(zdtgp(:,jlev)/ztgp(:,jlev)*zgw(:)         &
   !  &                       *zcp*zdps(:)/ga*dsigma(jlev)               &
   !  &                       ,mask=(zdtgp(:,jlev) >= 0.))
   !     zsum1(4)=zsum1(4)+SUM(zdtgp(:,jlev)/ztgp(:,jlev)*zgw(:)         &
   !  &                       *zcp*zdps(:)/ga*dsigma(jlev)               &
   !  &                       ,mask=(zdtgp(:,jlev) < 0.))
   !    enddo
   !    zsum3=SUM(zgw(:))
   !    !call mpsumbcr(zsum1,4)
   !    !call mpsumbcr(zsum3,1)
   !    zsum1(:)=zsum1(:)/zsum3
   !    !if(mypid == NROOT) then
   !      ztp=zsum1(1)/zsum1(3)
   !      zztm=zsum1(2)/zsum1(4)
   !      write(9,*) zsum1(:),zsum1(1)/zsum1(3),zsum1(2)/zsum1(4),(ztp-zztm)/ztp
   !    !endif
   !    deallocate(zst)
   !    deallocate(zstt)
   !    deallocate(zspf)
   !    deallocate(ztgp)
   !    deallocate(zdps)
   !    deallocate(zdtgp)
   !    deallocate(zsum1)
   !    deallocate(zgw)
   !   endif

!     7. Add newtonian cooling, friction and diffusion tendencies

      sdp(:,:) = sdp(:,:) + delt2 * sdt(:,:)
      szp(:,:) = szp(:,:) + delt2 * szt(:,:)
      stp(:,:) = stp(:,:) + delt2 * stt(:,:)

!     11. Coupling for synchronization runs

      !if (mrnum == 2 .and. nsync > 0) then
      !   call mrdiff(stp,std,NESP,NLEV)         ! 我认为完全可以去掉 !
      !   call mrdiff(sdp,sdd,NESP,NLEV)
      !   call mrdiff(szp,szd,NESP,NLEV)
      !   call mrdiff(spp,spd,NESP,   1)
      !   stp(:,:) = stp(:,:) + syncstr * std(:,:)
      !   sdp(:,:) = sdp(:,:) + syncstr * sdd(:,:)
      !   szp(:,:) = szp(:,:) + syncstr * szd(:,:)
      !   spp(:  ) = spp(:  ) + syncstr * spd(:  )
      !endif

!     8. Apply Robert Asselin time filter (not for short initial timesteps)
!        d(t) = pnu * f(t-1) + pnu * f(t+1) - 2 * pnu * f(t)

   ! nkits 是初始化小起步时的标志 程序一开始nkits=3 进行3个小起步后nkits=0并进入mainloop
   ! mainloop中nkits==0 所以以下if块在mainloop中执行 在mainloop之前不执行 故必须保留
   !   if (nkits == 0) then 
      if (mloop) then 
         zwp(:)   = pnu * (spm(:)   + spp(:)   - 2.0 * zpm(:)  )
         zwd(:,:) = pnu * (sdm(:,:) + sdp(:,:) - 2.0 * zdm(:,:))
         zwz(:,:) = pnu * (szm(:,:) + szp(:,:) - 2.0 * zzm(:,:))
         zwt(:,:) = pnu * (stm(:,:) + stp(:,:) - 2.0 * ztm(:,:))

!        Add Robert-Asselin-Williams filter value to f(t)

         spm(:)   = zpm(:)   + alpha * zwp(:)
         sdm(:,:) = zdm(:,:) + alpha * zwd(:,:)
         szm(:,:) = zzm(:,:) + alpha * zwz(:,:)
         stm(:,:) = ztm(:,:) + alpha * zwt(:,:)

!        Add filter value to f(t+1)
         
         spp(:)   = spp(:)   - (1.0 - alpha) * zwp(:)
         sdp(:,:) = sdp(:,:) - (1.0 - alpha) * zwd(:,:)
         szp(:,:) = szp(:,:) - (1.0 - alpha) * zwz(:,:)
         stp(:,:) = stp(:,:) - (1.0 - alpha) * zwt(:,:)
      endif

   ! 暂时关闭 for OpenACC testing
   ! checked, nenergy=nentropy=0
   !   if (nenergy > 0 .or. nentropy > 0) then
   !    zstte(:,:,1)=(stm(:,:)-ztm(:,:))/delt2
   !    zspt(:)=(spm(:)-zpm(:))/delt2
   !   endif
   !   if(nenergy > 0) then
   !    call mkenerdiag(ztm,zstte(:,:,1),zpm,zspt,denergy(:,9))
   !   endif
   !   if (nentropy > 0) then
   !    call mkentrodiag(ztm,zstte(:,:,1),zpm,dentropy(:,9))
   !   endif

!     9. Save spectral arrays for extended output

   ! 暂时关闭 for OpenACC testing
   ! checked, nextout=0, 用于entropy诊断输出的开关 默认off
   !   if (nextout == 1) then
   !      if (mod(nstep,nafter) == nafter - 2) then
   !         if (.not. allocated(st2)) allocate(st2(nesp,nlev))
   !         st2(:,:) = st(:,:)
   !         if (.not. allocated(sp2)) allocate(sp2(nesp))
   !         sp2(:) = sp(:)
   !      endif
   !      if (mod(nstep,nafter) == nafter - 1) then
   !         if (.not. allocated(st1)) allocate(st1(nesp,nlev))
   !         st1(:,:) = st(:,:)
   !         if (.not. allocated(sp1)) allocate(sp1(nesp))
   !         sp1(:) = sp(:)
   !      endif
   !   endif

!     10. Gather spectral modes from all processes

      sp = spp    !call mpgallsp(sp,spp,   1)
      sd = sdp    !call mpgallsp(sd,sdp,NLEV)
      sz = szp    !call mpgallsp(sz,szp,NLEV)
      st = stp    !call mpgallsp(st,stp,NLEV)

   ! 暂时关闭 for OpenACC testing
   ! checked, nenergy=nentropy=ndheat=0
   !   if(nenergy > 0 .or. nentropy > 0) then
   !    deallocate(zstte)
   !   endif
   !   if(ndheat > 0) then
   !    deallocate(zszte)
   !    deallocate(zsdte)
   !   endif
   !   if(nenergy > 0 .or. nentropy > 0 .or. ndheat > 0) then
   !    deallocate(zsp)
   !    deallocate(zspt)
   !   endif

!$acc end kernels

      end subroutine spectral


!     ================
!     SUBROUTINE DIAGP
!     ================

      subroutine diagp(zampl)
      use pumamod
      implicit none

      real :: zstf(NESP,NLEV)
      real :: zgr12(NHOR,NLEV)
      real :: zgtt(NHOR,NLEV)
      real :: gr12(NHOR,NLEV)
      real :: gr12c(NHOR,NLEV)
     

      real :: gdtmp(NHOR)

      real :: zampl
      integer :: jlev

      !--- transform temperature and divergence to grid point space
      st = stp    !call mpgallsp(st,stp,NLEV)
      if (nconv > 0) then
         sd = sdp !call mpgallsp(sd,sdp,NLEV)
      endif
      do jlev=1,NLEV
         call sp2fc(st(1,jlev)   ,gt(1,jlev)   )
         if (nconv > 0) then
            call sp2fc(sd(1,jlev)   ,gd(1,jlev)   )
         endif
      enddo
      call fc2gp(gt   ,NLON,NLPP*NLEV,trigs)
      if (nconv > 0) then
         call fc2gp(gd   ,NLON,NLPP*NLEV,trigs)
      endif


      !--- radiative temperature tendencies 
      gr12(:,:) = gr1(:,:) + gr2(:,:)*zampl
      zgtt(:,:) = (gr12(:,:) - gt(:,:))*gtdamp(:,:)

      !--- add convective temperature tendencies
      if (nconv > 0) then
         gdtmp(:) = gd(:,nlev)
         do jlev = 1,nlev
            where (gdtmp < 0.0)
               gr12c(:,jlev) = gr1c(:,jlev) + gr2c(:,jlev)*zampl
               zgtt(:,jlev)  = zgtt(:,jlev) + (gr12c(:,jlev) - gt(:,jlev))*gtdampc(:,jlev)
            endwhere
         enddo
      endif

      !--- transform temperature tendencies to spectral space
      call gp2fc(zgtt ,NLON,NLPP*NLEV,trigs)
      do jlev=1,NLEV
         call fc2sp(zgtt(1,jlev),zstf(1,jlev))
      enddo
      stt   = zstf   !call mpsumsc(zstf,stt,NLEV)

      return
      end subroutine diagp


!     =================
!     SUBROUTINE HEATGP
!     =================

      subroutine heatgp(zampl)
      use pumamod
      implicit none

      real :: zsr12(NESP,NLEV)
      real :: zsrp12(NSPP,NLEV)
      real :: zstf(NESP,NLEV)
      real :: zgr12(NHOR,NLEV)
      real :: zgtt(NHOR,NLEV)

      real :: zampl
      integer :: jlev

      zsrp12(:,:)=srp1(:,:)+srp2(:,:)*zampl
      zsr12 = zsrp12    !call mpgallsp(zsr12,zsrp12,NLEV)
      st    = stp       !call mpgallsp(st,stp,NLEV)
      do jlev=1,NLEV
         call sp2fc(zsr12(1,jlev),zgr12(1,jlev))
         call sp2fc(st(1,jlev)   ,gt(1,jlev)   )
      enddo
      call fc2gp(zgr12,NLON,NLPP*NLEV,trigs)
      call fc2gp(gt   ,NLON,NLPP*NLEV,trigs)

!     Newtonian cooling

      zgtt(:,:) = (zgr12(:,:) - gt(:,:)) * gtdamp(:,:)

      call gp2fc(zgtt ,NLON,NLPP*NLEV,trigs)
      do jlev=1,NLEV
         call fc2sp(zgtt(1,jlev),zstf(1,jlev))
      enddo
      stt   = zstf   !call mpsumsc(zstf,stt,NLEV)

      return
      end subroutine heatgp


!     ================
!     SUBROUTINE VDIFF
!     ================

      subroutine vdiff(pt,pz,pd,ptt,pzt,pdt)
      use pumamod
      implicit none
!
      real, parameter :: ztref=250.0

      real pt(NSPP,NLEV),pz(NSPP,NLEV),pd(NSPP,NLEV)
      real ptt(NSPP,NLEV),pzt(NSPP,NLEV),pdt(NSPP,NLEV)
      real ztn(NSPP,NLEV),zzn(NSPP,NLEV),zdn(NSPP,NLEV)
      real zebs(NLEM)
      real zskap(NLEV),zskaph(NLEV)
      real zkdiff(NLEM)

      real :: zdelt, zkonst1, zkonst2
      integer :: jlp,jlev,jlep,jlem
!
      zdelt=delt2/ww_scale
      zkonst1=ga*zdelt/gascon
      zkonst2=zkonst1*ga/gascon
!
      zskap(:)=sigma(:)**akap
      zskaph(:)=sigmh(:)**akap
!
!     1) modified diffusion coefficents
!
      do jlev=1,NLEM
       jlp=jlev+1
       zkdiff(jlev)=zkonst2*sigmh(jlev)*sigmh(jlev)/(ztref*ztref)     &
     &             *dvdiff/(sigma(jlp)-sigma(jlev))
      enddo
!
!     2. semi implicit scheme
!
!     2a momentum
!
!     top layer elimination
!
      zebs(1)=zkdiff(1)/(dsigma(1)+zkdiff(1))
      zdn(:,1)=dsigma(1)*pd(:,1)/(dsigma(1)+zkdiff(1))
      zzn(:,1)=dsigma(1)*pz(:,1)/(dsigma(1)+zkdiff(1))
!
!     middle layer elimination
!
      do jlev=2,NLEM
       jlem=jlev-1
       zebs(jlev)=zkdiff(jlev)/(dsigma(jlev)+zkdiff(jlev)               &
     &           +zkdiff(jlem)*(1.-zebs(jlem)))
       zdn(:,jlev)=(pd(:,jlev)*dsigma(jlev)+zkdiff(jlem)*zdn(:,jlem))   &
     &            /(dsigma(jlev)+zkdiff(jlev)                           &
     &            +zkdiff(jlem)*(1.-zebs(jlem)))
       zzn(:,jlev)=(pz(:,jlev)*dsigma(jlev)+zkdiff(jlem)*zzn(:,jlem))   &
     &            /(dsigma(jlev)+zkdiff(jlev)                           &
     &            +zkdiff(jlem)*(1.-zebs(jlem)))
      enddo
!
!     bottom layer elimination
!
      zdn(:,NLEV)=(pd(:,NLEV)*dsigma(NLEV)+zkdiff(NLEM)*zdn(:,NLEM))    &
     &           /(dsigma(NLEV)+zkdiff(NLEM)*(1.-zebs(NLEM)))
      zzn(:,NLEV)=(pz(:,NLEV)*dsigma(NLEV)+zkdiff(NLEM)*zzn(:,NLEM))    &
     &           /(dsigma(NLEV)+zkdiff(NLEM)*(1.-zebs(NLEM)))
!
!     back-substitution
!
      do jlev=NLEM,1,-1
       jlep=jlev+1
       zdn(:,jlev)=zdn(:,jlev)+zebs(jlev)*zdn(:,jlep)
       zzn(:,jlev)=zzn(:,jlev)+zebs(jlev)*zzn(:,jlep)
      enddo
!
!     tendencies
!
      pdt(:,1:NLEV)=pdt(:,1:NLEV)+(zdn(:,1:NLEV)-pd(:,1:NLEV))/delt2
      pzt(:,1:NLEV)=pzt(:,1:NLEV)+(zzn(:,1:NLEV)-pz(:,1:NLEV))/delt2
!
!     2c potential temperature
!
      do jlev=1,NLEM
       zkdiff(jlev)=zkdiff(jlev)*zskaph(jlev)
      enddo
!
!     semi implicit scheme
!
!     top layer elimination
!
      zebs(1)=zkdiff(1)/(dsigma(1)+zkdiff(1)/zskap(1))
      ztn(:,1)=dsigma(1)*pt(:,1)/(dsigma(1)+zkdiff(1)/zskap(1))
!
!     middle layer elimination
!
      do jlev=2,NLEM
       jlem=jlev-1
       zebs(jlev)=zkdiff(jlev)/(dsigma(jlev)+(zkdiff(jlev)              &
     &           +zkdiff(jlem)*(1.-zebs(jlem)/zskap(jlem)))/zskap(jlev))
       ztn(:,jlev)=(pt(:,jlev)*dsigma(jlev)                             &
     &             +zkdiff(jlem)/zskap(jlem)*ztn(:,jlem))               &
     &            /(dsigma(jlev)+(zkdiff(jlev)                          &
     &             +zkdiff(jlem)*(1.-zebs(jlem)/zskap(jlem)))           &
     &             /zskap(jlev))
      enddo
!
!     bottom layer elimination
!
      ztn(:,NLEV)=(pt(:,NLEV)*dsigma(NLEV)                              &
     &            +zkdiff(NLEM)*ztn(:,NLEM)/zskap(NLEM))                &
     &           /(dsigma(NLEV)+zkdiff(NLEM)/zskap(NLEV)                &
     &                         *(1.-zebs(NLEM)/zskap(NLEM)))
!
!     back-substitution
!
      do jlev=NLEM,1,-1
       jlep=jlev+1
       ztn(:,jlev)=ztn(:,jlev)+zebs(jlev)*ztn(:,jlep)/zskap(jlep)
      enddo
!
!     tendencies
!
      ptt(:,1:NLEV)=ptt(:,1:NLEV)+(ztn(:,1:NLEV)-pt(:,1:NLEV))/delt2
!
      return
      end subroutine vdiff


!     =================
!     SUBROUTINE GASDEV
!     =================
!     Gaussian noise generator with zero mean and unit variance.

      function gasdev() 
      use pumamod
      implicit none
      real :: gasdev
      real :: fr, vx, vy, ra

      if (ganext == 0.0) then
         ra = 2.0
         do while (ra >= 1.0 .or. ra < 1.0e-20)
            call random_number(vx)
            call random_number(vy)
            vx = 2.0 * vx - 1.0
            vy = 2.0 * vy - 1.0
            ra = vx * vx + vy * vy
         enddo
         fr = sqrt(-2.0 * log(ra) / ra)
         gasdev = vx * fr
         ganext = vy * fr
      else
         gasdev = ganext
         ganext = 0.0
      endif
      end function gasdev


!     =================
!     SUBROUTINE SPONGE
!     =================

      subroutine sponge
      use pumamod
      implicit none

      real :: zp
      integer :: jlev

!     This introduces a simple sponge layer to the highest model levels
!     by applying Rayleigh friction there, according to
!     Polvani & Kushner (2002, GRL), see their appendix.

      write(nud,*)
      write(nud,9991)
      write(nud,9997)
      write(nud,9991)
      write(nud,9996)
      write(nud,9991)
      do jlev=1,NLEV
         zp = sigma(jlev)*psurf
         if (zp < pspon) then
            fric(jlev) = (dcsponge * sol_day &
                       * ((pspon - zp) / pspon)**2) / TWOPI
         endif

!        some output
         if (zp > pspon) then
            if (fric(jlev) == 0) then
               write(nud,9992) jlev
            else
               write(nud,9993) jlev, fric(jlev)*TWOPI
            endif
         else
            if (fric(jlev) == 0) then
               write(nud,9994) jlev
            else
               write(nud,9995) jlev, fric(jlev)*TWOPI
            endif
         endif
      enddo
      write(nud,9991)
      write(nud,*)
      return
 9991 format(33('*'))
 9992 format('*',i4,' * ',7('-'),' *               *')
 9993 format('*',i4,' * ',f7.4,' *               *')
 9994 format('*',i4,' * ',7('-'),' *',' SPONGE        *')
 9995 format('*',i4,' * ',f7.4,' *',' SPONGE        *')
 9996 format('*  Lv * [1/day] *               *')
 9997 format('* Rayleigh damping coefficients *')
      end subroutine sponge


!     =====================
!     SUBROUTINE MKENERDIAG
!     =====================

      subroutine mkenerdiag(pst,pstt,psp,pspt,penergy)
      use pumamod
      implicit none
!
      real :: pst(NSPP,NLEV),pstt(NSPP,NLEV)
      real :: psp(NSPP),pspt(NSPP)
      real :: penergy(NHOR)
!
      real :: zsttf(NESP,NLEV),zstf(NESP,NLEV)
      real :: zsptf(NESP),zspf(NESP)
      real :: zgtt(NHOR,NLEV),zgt(NHOR,NLEV)
      real :: zgps(NHOR),zgpst(NHOR)
      real :: ztm(NHOR)

      real :: zcp, zdelt
      integer :: jlev
!
      zcp=gascon/akap
      zdelt=delt2/ww_scale
!
      zsttf = pstt   !call mpgallsp(zsttf,pstt,NLEV)
      zstf  = pst    !call mpgallsp(zstf,pst,NLEV)
      zsptf = pspt   !call mpgallsp(zsptf,pspt,1)
      zspf  = psp    !call mpgallsp(zspf,psp,1)

      do jlev=1,NLEV
       call sp2fc(zsttf(:,jlev),zgtt(:,jlev))
       call sp2fc(zstf(:,jlev),zgt(:,jlev))
      enddo
      call sp2fc(zsptf,zgpst)
      call sp2fc(zspf,zgps)
      call fc2gp(zgtt,NLON,NLPP*NLEV,trigs)
      call fc2gp(zgt,NLON,NLPP*NLEV,trigs)
      call fc2gp(zgps,NLON,NLPP,trigs)
      call fc2gp(zgpst,NLON,NLPP,trigs)
      zgpst(:)=psurf*(exp(zgps(:)+delt2*zgpst(:))-exp(zgps(:)))/zdelt
      zgps(:)=psurf*exp(zgps(:))+zdelt*zgpst(:)
      zgtt(:,:)=ct*ww_scale*zgtt(:,:)
      do jlev=1,NLEV
       zgt(:,jlev)=ct*(zgt(:,jlev)+t0(jlev))
      enddo
!
      ztm(:)=0.
      penergy(:)=0.
      do jlev=1,NLEV
       ztm(:)=ztm(:)+zgt(:,jlev)*dsigma(jlev)
       penergy(:)=penergy(:)+zgtt(:,jlev)*dsigma(jlev)
      enddo
      penergy(:)=ztm(:)*zcp*zgpst(:)/ga+penergy(:)*zcp*zgps(:)/ga
!
      end subroutine mkenerdiag


!     ======================
!     SUBROUTINE MKENTRODIAG
!     ======================

      subroutine mkentrodiag(pst,pstt,psp,pentropy)
      use pumamod
      implicit none
!
      real :: pst(NSPP,NLEV),pstt(NSPP,NLEV)
      real :: psp(NSPP)
      real :: pentropy(NHOR)
!
      real :: zsttf(NESP,NLEV),zstf(NESP,NLEV)
      real :: zspf(NESP)
      real :: zgtt(NHOR,NLEV),zgt(NHOR,NLEV)
      real :: zgps(NHOR)

      real :: zcp
      integer :: jlev
!
      zcp=gascon/akap
!
      zsttf = pstt   !call mpgallsp(zsttf,pstt,NLEV)
      zstf  = pst    !call mpgallsp(zstf,pst,NLEV)
      zspf  = psp    !call mpgallsp(zspf,psp,1)

      do jlev=1,NLEV
       call sp2fc(zsttf(:,jlev),zgtt(:,jlev))
       call sp2fc(zstf(:,jlev),zgt(:,jlev))
      enddo
      call sp2fc(zspf,zgps)
      call fc2gp(zgtt,NLON,NLPP*NLEV,trigs)
      call fc2gp(zgt,NLON,NLPP*NLEV,trigs)
      call fc2gp(zgps,NLON,NLPP,trigs)
      zgps(:)=psurf*exp(zgps(:))
      zgtt(:,:)=ct*ww_scale*zgtt(:,:)
      do jlev=1,NLEV
       zgt(:,jlev)=ct*(zgt(:,jlev)+t0(jlev))
      enddo
!
      pentropy(:)=0.
      do jlev=1,NLEV
       pentropy(:)=pentropy(:)+zgtt(:,jlev)*dsigma(jlev)/zgt(:,jlev)
      enddo
      pentropy(:)=pentropy(:)*zcp*zgps(:)/ga                               
!
      return
      end subroutine mkentrodiag


!     ==================
!     SUBROUTINE MKDHEAT
!     ==================

      subroutine mkdheat(zszt1,zszt2,zsdt1,zsdt2,zsp)
      use pumamod
      implicit none
!
!     'recycle' kin. energy loss by heating the environment
!
!     zszt1/zsdt1 : vorticity/divergence tendency due to friction
!     zszt2/zsdt2 : vorticity/divergence tendency fue to diffusion  
!     zp          : surface pressure
!
      real zszt1(NSPP,NLEV),zszt2(NSPP,NLEV)
      real zsdt1(NSPP,NLEV),zsdt2(NSPP,NLEV)
      real zsp(NSPP)
      real zp(NHOR)
!
      real zsd(NESP,NLEV),zsz(NESP,NLEV)
      real zspf(NESP),zspt(NSPP)
      real zsdp(NSPP,NLEV),zszp(NSPP,NLEV)
      real zu(NHOR,NLEV),zun(NHOR,NLEV),zdu1(NHOR,NLEV),zdu2(NHOR,NLEV)
      real zv(NHOR,NLEV),zvn(NHOR,NLEV),zdv1(NHOR,NLEV),zdv2(NHOR,NLEV)
      real zdtdt1(NHOR,NLEV),zdtdt2(NHOR,NLEV),zdtdt3(NHOR,NLEV)
!
      real zdtdt(NHOR,NLEV),zdekin(NHOR,NLEV)
!
      real zsde(NSPP,NLEV),zsdef(NESP,NLEV)
      real zstt(NSPP,NLEV),zstf(NESP,NLEV)
      real zstt1(NSPP,NLEV),zstf1(NESP,NLEV),zstt3(NSPP,NLEV)
      real zstt2(NSPP,NLEV),zstf2(NESP,NLEV),zstf3(NESP,NLEV)

      real :: zcp, zdelt
      integer :: jlev
!
!     some constants
!
      zdelt=delt2/ww_scale     ! timestep in s
      zcp=gascon/akap    ! heat capacity
!
!     'recycle' friction
!
!     a) gather the 'partial' field of z and d, and make u and v 
!        at old time level
!
      zsdp(:,:)=sdp(:,:)
      zszp(:,:)=szp(:,:)
      zsd   = zsdp   !call mpgallsp(zsd,zsdp,NLEV)
      zsz   = zszp   !call mpgallsp(zsz,zszp,NLEV)
      do jlev = 1 , NLEV
         call dv2uv(zsd(1,jlev),zsz(1,jlev),zu(1,jlev),zv(1,jlev))
      enddo
      call fc2gp(zu,NLON,NLPP*NLEV,trigs)
      call fc2gp(zv,NLON,NLPP*NLEV,trigs)
!
!     b) add fricton tendencies and create new u and v
!
      zsdp(:,:)=sdp(:,:)+zsdt1(:,:)*delt2
      zszp(:,:)=szp(:,:)+zszt1(:,:)*delt2
      zsd   = zsdp   !call mpgallsp(zsd,zsdp,NLEV)
      zsz   = zszp   !call mpgallsp(zsz,zszp,NLEV)
      do jlev = 1 , NLEV
         call dv2uv(zsd(1,jlev),zsz(1,jlev),zun(1,jlev),zvn(1,jlev))
      enddo
      call fc2gp(zun,NLON,NLPP*NLEV,trigs)
      call fc2gp(zvn,NLON,NLPP*NLEV,trigs)
!
!     c) compute temperature tendency
!
      do jlev=1,NLEV
       zu(:,jlev)=cv*zu(:,jlev)*SQRT(rcsq(:))
       zv(:,jlev)=cv*zv(:,jlev)*SQRT(rcsq(:))
       zun(:,jlev)=cv*zun(:,jlev)*SQRT(rcsq(:))
       zvn(:,jlev)=cv*zvn(:,jlev)*SQRT(rcsq(:))
       zdu1(:,jlev)=zun(:,jlev)-zu(:,jlev)
       zdv1(:,jlev)=zvn(:,jlev)-zv(:,jlev)
       zdtdt1(:,jlev)=-(zun(:,jlev)*zun(:,jlev)                         &
     &                 -zu(:,jlev)*zu(:,jlev)                           &
     &                 +zvn(:,jlev)*zvn(:,jlev)                         &
     &                 -zv(:,jlev)*zv(:,jlev))*0.5/zdelt/zcp     
      enddo

!
!     'recycle' momentum diffusion  
!    
!     a) add tendencies and create new u and v and get surface pressure
!
!
      zsdp(:,:)=sdp(:,:)+zsdt2(:,:)*delt2
      zszp(:,:)=szp(:,:)+zszt2(:,:)*delt2
      zsd   = zsdp   !call mpgallsp(zsd,zsdp,NLEV)
      zsz   = zszp   !call mpgallsp(zsz,zszp,NLEV)
      zspf  = zsp    !call mpgallsp(zspf,zsp,1)
      do jlev = 1 , NLEV
         call dv2uv(zsd(1,jlev),zsz(1,jlev),zun(1,jlev),zvn(1,jlev))
      enddo
      call fc2gp(zun,NLON,NLPP*NLEV,trigs)
      call fc2gp(zvn,NLON,NLPP*NLEV,trigs)
      call sp2fc(zspf,zp)
      call fc2gp(zp,NLON,NLPP,trigs)
      zp(:)=psurf*exp(zp(:))
!
!     b) compute loss of kinetic energy
!        (note: only the global average change of kin. e. is 'lost'
!         the other changes are just diffusion)
!
      do jlev = 1 , NLEV
       zun(:,jlev)=cv*zun(:,jlev)*SQRT(rcsq(:))
       zvn(:,jlev)=cv*zvn(:,jlev)*SQRT(rcsq(:))
       zdu2(:,jlev)=zun(:,jlev)-zu(:,jlev)
       zdv2(:,jlev)=zvn(:,jlev)-zv(:,jlev)
       zdekin(:,jlev)=(zun(:,jlev)*zun(:,jlev)                          &
     &                -zu(:,jlev)*zu(:,jlev)                            &
     &                +zvn(:,jlev)*zvn(:,jlev)                          &
     &                -zv(:,jlev)*zv(:,jlev))*0.5/zdelt                 &
     &               *zp(:)/ga*dsigma(jlev)   
      enddo
!
!     c) get the global average and transform it back
!
      call gp2fc(zdekin,NLON,NLPP*NLEV,trigs)
      do jlev=1,NLEV
       call fc2sp(zdekin(:,jlev),zsdef(:,jlev))
      enddo
      zsde  = zsdef  !call mpsumsc(zsdef,zsde,NLEV)
      zsdef = zsde   !call mpgallsp(zsdef,zsde,NLEV)
      zsdef(2:NESP,:)=0.
      do jlev = 1 , NLEV
         call sp2fc(zsdef(1,jlev),zdekin(1,jlev))
      enddo
      call fc2gp(zdekin,NLON,NLPP*NLEV,trigs)
!
!     d) compute temperature tendency
!
      do jlev=1,NLEV
       zdtdt2(:,jlev)=-zdekin(:,jlev)*ga/zp(:)/dsigma(jlev)/zcp
       zdtdt3(:,jlev)=-(zdu1(:,jlev)*zdu2(:,jlev)                       &
     &                 +zdv1(:,jlev)*zdv2(:,jlev))/zdelt/zcp
      enddo
!
      zdtdt1(:,:)=zdtdt1(:,:)/ct/ww_scale
      zdtdt2(:,:)=zdtdt2(:,:)/ct/ww_scale
      zdtdt3(:,:)=zdtdt3(:,:)/ct/ww_scale
!
      call gp2fc(zdtdt1,NLON,NLPP*NLEV,trigs)
      call gp2fc(zdtdt2,NLON,NLPP*NLEV,trigs)
      call gp2fc(zdtdt3,NLON,NLPP*NLEV,trigs)
      do jlev=1,NLEV
       call fc2sp(zdtdt1(:,jlev),zstf1(:,jlev))
       call fc2sp(zdtdt2(:,jlev),zstf2(:,jlev))
       call fc2sp(zdtdt3(:,jlev),zstf3(:,jlev))
      enddo
      zstt1 = zstf1  !call mpsumsc(zstf1,zstt1,NLEV)
      zstt2 = zstf2  !call mpsumsc(zstf2,zstt2,NLEV)
      zstt3 = zstf3  !call mpsumsc(zstf3,zstt3,NLEV)
!
!     add the temprature tendencies
!
      stt(:,:)=stt(:,:)+zstt1(:,:)+zstt2(:,:)+zstt3(:,:)
!
!     energy diagnostics
!
      if(nenergy > 0) then
       zspt(:)=0.
       call mkenerdiag(stp,zstt1,zsp,zspt,denergy(:,5))
       call mkenerdiag(stp,zstt2,zsp,zspt,denergy(:,6))
       call mkenerdiag(stp,zstt3,zsp,zspt,denergy(:,7))
      endif
      if(nentropy > 0) then
       call mkentrodiag(stp,zstt1,zsp,dentropy(:,5))
       call mkentrodiag(stp,zstt2,zsp,dentropy(:,6))
       call mkentrodiag(stp,zstt3,zsp,dentropy(:,7))
      endif
!
      return
      end subroutine mkdheat


!     =================
!     SUBROUTINE MKEKIN
!     XW deleted:
!      subroutine mkekin(zszp,zsdp,zp,zekin)
!      subroutine mkekin2(zszp,zsdp,zspp,zekin)
!     =================

!     =================
!     SUBROUTINE MKEPOT
!     XW deleted:
!     subroutine mkepot(zstp,zp,zepot)
!     subroutine mkepot2(zstp,zspp,zepot)
!     =================


!========================
! Compute mastercpu time
! XW(2017/4/9)
!========================
subroutine mastercpu_time(x)
   implicit none
   real(kind=8), intent(out) :: x
   ! local
   integer, dimension(8) :: v
   call date_and_time(values=v)
   ! SO FAR, only count from the begining of THIS month
   ! 目前只以当月1日为起点进行计时，以后要改进为绝对零点，考虑闰年问题即可
   x = 86400.0*v(3)+3600.0*v(5)+60.0*v(6)+v(7)+v(8)/1000.0
end subroutine mastercpu_time

