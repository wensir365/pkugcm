# Makefile for PUMA without GUI
# Xinyu Wen, Peking Univ, Mar/23/2017

#=====================
# 1. Compiler & Flags
#=====================
# GNU gfortran
F90    = gfortran
FLAGS  = -O3 -fopenmp #-Wall #-pg 
#OTHERS = -ffree-line-length-none
#USELESS= -fdefault-real-8 -fbacktrace -ffpe-trap=underflow,denormal

# PGI pgfortran
#F90    = pgfortran
#FLAGS  = -fast

# Intel ifort
#F90    = ifort
#FLAGS  = -fast
#OTHERS = -parallel
#OTHERS = -qopenmp
#USELESS= -parallel

#=====================
# 2. if MPI?
#=====================
# No
#MPIMOD	= mod_mpi_stub

#====================
# 3. Compile
#====================
OBJ = pkugcm.o mod_fft.o mod_legsym.o mod_restart.o mod_io.o

%.o : %.f90
	${F90} $(FLAGS) $(OTHERS) -c $<

puma.x:	$(OBJ)
	${F90} $(FLAGS) -o pkugcm.x $(OBJ)

pkugcm.o:      pkugcm.f90
mod_fft.o:     mod_fft.f90
mod_legsym.o:  mod_legsym.f90
mod_restart.o: mod_restart.f90
mod_io.o:      mod_io.f90

#====================
# 4. Clean
#====================
clean:
	rm -f *.o *.mod

cleano:
	rm -f *.o

cleanmod:
	rm -f *.mod
 
cleanall:
	rm -f *.o *.mod
	rm -f puma.x

