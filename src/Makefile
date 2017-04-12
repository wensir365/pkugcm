# Makefile for PUMA without GUI
# Xinyu Wen, Peking Univ, Mar/23/2017

#=====================
# 1. Compiler & Flags
#=====================
# GNU gfortran
F90    = gfortran
#FLAGS  = -O3 -fdefault-real-8 -ffree-line-length-none #-Wall #-pg
FLAGS  = -O3 -ffree-line-length-none -fopenmp #-Wall #-pg 

# PGI pgfortran
#F90    = pgfortran
#FLAGS  = -r8 -O2 ‑Ktrap=fp,ovf,unf,inexact #-Kieee

#=====================
# 2. if MPI?
#=====================
# No
#MPIMOD	= mod_mpi_stub

# Yes
F90	= mpif90
FLAGS  = -O3
MPIMOD	= mod_mpi

#====================
# 3. Compile
#====================
OBJ=${MPIMOD}.o mod_fft.o puma.o mod_legsym.o mod_restart.o

%.o : %.f90
	${F90} $(FLAGS) -c $<

puma.x:	$(OBJ)
	${F90} $(FLAGS) -o puma.x $(OBJ)

puma.o:		puma.f90
mod_fft.o:	mod_fft.f90
mod_legsym.o:	mod_legsym.f90
mod_restart.o:	mod_restart.f90
${MPIMOD}.o:	${MPIMOD}.f90 puma.o

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

