# Change Log

## 2017-Mar-27: Add MPI function

Add mpimod.f90 and make minor modification on makefile, so that I can deliver MPI runs using 4 CPUs. Seems the excutable with MPI still can be run with 1 single CPU, like this: "./puma_mpi.x 32 10". Normally, the MPI run can be done by "mpiexec -np 4 puma_mpi.x 32 10".

The running speed, under the unif of "simulate model years per day",  was tested on my Intel i5 CPU with 4 processors, with gfortran as the compilor and -O3 as the optimal argument.

| `gfortran -O3` |   T21 |  T42 |
|:---------------|------:|-----:|
| 1 CPU (No MPI) |  4373 |  299 |
| 4 CPU (MPI)    | 13918 |  782 |

Note that the 1st number here (4373) is 1755 tested on pkuclimate.club, an 5-year-old Intel i3 Xeon CPU.

## 2017-Mar-25: Add Post-Processing (pp)

Add directory "pp" for the post-processing excutable "burn7.x" and 2 associated namelist, one for 3 surface variables, another for 10 multi-level variables. To compile burn7.x you should apt install libnetcdf-cxx-legacy-dev first.

## 2017-Mar-25: Remove GUI-related code
Remove guimod_stub.f90, and comment out all the GUI-related code in puma.f90 with a description "XW(Mar/25/2017) to remove GUI:". Please note I did not delete any lines, just comment out. The number of source code were shinked from 7 to 6.

## 2017-Mar-23: Make PUMA running

Successfully make PUMA running on 1CPU, with fake MPI and GUI interface (stub). 
I changed "nresources" line in puma.f90, which might be related to retrieving system time used in pumax.c.
They are marked with "XW".

## 2017-Mar-17: Make PLASIM17 running with MOST

### PUMA 4core 1year T21 run

	wensir@himalaya ~/model/puma17 $ ./most.x 
	mpif90 -c -O3 -ffpe-trap=invalid,zero,overflow -ffpe-summary=none -finit-real=snan  puma.f90
	mpif90 -c -O3 -ffpe-trap=invalid,zero,overflow -ffpe-summary=none -finit-real=snan  mpimod.f90
	mpif90 -c -O3 -ffpe-trap=invalid,zero,overflow -ffpe-summary=none -finit-real=snan  fftmod.f90
	mpif90 -c -O3 -ffpe-trap=invalid,zero,overflow -ffpe-summary=none -finit-real=snan  guimod.f90
	mpicc -c -O3 pumax.c
	mpif90 -c -O3 -ffpe-trap=invalid,zero,overflow -ffpe-summary=none -finit-real=snan  legini.f90
	as -o legfast32.o legfast32.s
	mpif90 -c -O3 -ffpe-trap=invalid,zero,overflow -ffpe-summary=none -finit-real=snan  restartmod.f90
	mpif90 -c -O3 -ffpe-trap=invalid,zero,overflow -ffpe-summary=none -finit-real=snan  gaussmod.f90
	mpif90 -o puma.x -O3 -ffpe-trap=invalid,zero,overflow -ffpe-summary=none -finit-real=snan 
               mpimod.o fftmod.o guimod.o pumax.o legini.o legfast32.o puma.o restartmod.o gaussmod.o 
               -L/usr/lib/X11 -lX11

	=== Success: Launched process most_puma_run ===


### PUMA 1core 1year T21 run

	wensir@himalaya ~/model/puma17 $ ./most.x 
	gfortran -c -O3 -ffpe-trap=invalid,zero,overflow -ffpe-summary=none -finit-real=snan  puma.f90
	gfortran -c -O3 -ffpe-trap=invalid,zero,overflow -ffpe-summary=none -finit-real=snan  mpimod_stub.f90
	gfortran -c -O3 -ffpe-trap=invalid,zero,overflow -ffpe-summary=none -finit-real=snan  fftmod.f90
	gfortran -c -O3 -ffpe-trap=invalid,zero,overflow -ffpe-summary=none -finit-real=snan  guimod.f90
	gcc   -c -O3 -I /usr/lib/X11/include pumax.c
	gfortran -c -O3 -ffpe-trap=invalid,zero,overflow -ffpe-summary=none -finit-real=snan  legini.f90
	as -o legfast32.o legfast32.s
	gfortran -c -O3 -ffpe-trap=invalid,zero,overflow -ffpe-summary=none -finit-real=snan  restartmod.f90
	gfortran -c -O3 -ffpe-trap=invalid,zero,overflow -ffpe-summary=none -finit-real=snan  gaussmod.f90
	gfortran -o puma.x -O3 -ffpe-trap=invalid,zero,overflow -ffpe-summary=none -finit-real=snan 
                 mpimod_stub.o fftmod.o guimod.o pumax.o legini.o legfast32.o puma.o restartmod.o gaussmod.o 
                 -L/usr/lib/X11 -lX11

	=== Success: Launched process most_puma_run ===


### Configure and Compile Source Code

	wensir@himalaya ~/model/puma17 $ ./configure.sh 
	Found FORTRAN-90 compiler     at: /usr/bin/gfortran
	gfortran version  5.4
	Found Xlib (X11)              at: /usr/lib/X11
	Found C compiler              at: /usr/bin/gcc
	Found C++ compiler            at: /usr/bin/c++
	Found MPI FORTRAN-90 compiler at: /usr/bin/mpif90
	Found MPI C compiler          at: /usr/bin/mpicc
	Found MPI C++ compiler        at: /usr/bin/mpicxx
	Found GNU assembler           at: /usr/bin/as
	Fast Legendre Transformation    : ACTIVE (Linux)
	gfortran -o f90check.x f90check.f90
	gcc -o cc_check.x cc_check.c

	System info for <himalaya>
	Architecture: Linux himalaya 4.4.0-53-generic #74-Ubuntu SMP Fri Dec 2 15:59:10 UTC 2016 x86_64 GNU/Linux
	Endian format             : little endian
	FORTRAN control word size : 4 bytes
	FORTRAN integer size      : 4 bytes
	FORTRAN real    size      : 4 bytes
	C       int     size      : 4 bytes
	C       float   size      : 4 bytes
	C       long    size      : 8 bytes
	C  long long    size      : 8 bytes
	FORTRAN Compiler: gfortran
	C       Compiler: gcc
	gcc -o most.x most.c -I/usr/lib/X11/include -lm -L/usr/lib/X11 -lX11

	configuration complete - run <most.x>


### Follow README_UBUNTU:

- apt install libx11-dev
- apt install openmpi-bin openmpi-common openmpi-doc libopenmpi-dev
