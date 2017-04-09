# Makefile for PUMA without GUI
# Xinyu Wen, Peking Univ, Mar/23/2017

#=====================
# 1. Compiler & Flags
#=====================
# GNU gfortran
F90    = gfortran
#FLAGS  = -O3 -fdefault-real-8 -ffree-line-length-none #-Wall #-pg
FLAGS  = -O3 -ffree-line-length-none -fopenmp -pg #-Wall #-pg 

# PGI pgfortran
#F90    = pgfortran
#FLAGS  = -r8 -O2 ‑Ktrap=fp,ovf,unf,inexact #-Kieee

#=====================
# 2. if MPI?
#=====================
# No
MPIMOD	= mpimod_stub
FFTMOD	= fftmod

# Yes
#F90	= mpif90
#FLAGS  = -O2 -Wall
#MPIMOD	= mpimod
#FFTMOD	= fftmod

#====================
# 3. Compile
#====================
OBJ=${MPIMOD}.o ${FFTMOD}.o puma.o legsym.o restartmod.o gaussmod.o

%.o : %.f90
	${F90} $(FLAGS) -c $<

puma.x:	$(OBJ)
	${F90} $(FLAGS) -o puma.x $(OBJ)

puma.o:		puma.f90
${FFTMOD}.o:	${FFTMOD}.f90
${MPIMOD}.o:	${MPIMOD}.f90 puma.o
legsym.o:	legsym.f90
restartmod.o:	restartmod.f90
gaussmod.o:	gaussmod.f90

#====================
# 4. Clean
#====================
clean:
	rm -f *.o *.mod

