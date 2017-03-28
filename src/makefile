# Makefile for PUMA without MPI and GUI
# Xinyu Wen, Peking Univ, Mar/23/2017

F90=gfortran
MPIMOD=mpimod_stub
FFTMOD=fftmod

# Makefile for PUMA with MPI, without GUI
# Xinyu Wen, Peking Univ, Mar/27/2017

#F90=mpif90
#MPIMOD=mpimod
#FFTMOD=fftmod

OBJ=${MPIMOD}.o ${FFTMOD}.o puma.o legsym.o restartmod.o gaussmod.o

%.o : %.f90
	${F90} -c -O3 $<

puma.x:	$(OBJ)
	${F90} -o puma.x $(OBJ)

puma.o:		puma.f90
${FFTMOD}.o:	${FFTMOD}.f90
${MPIMOD}.o:	${MPIMOD}.f90 puma.o
legsym.o:	legsym.f90
restartmod.o:	restartmod.f90
gaussmod.o:	gaussmod.f90

clean:
	rm -f *.o *.mod
