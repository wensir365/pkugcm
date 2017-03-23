# Makefile for PUMA without GUI and MPI
# Xinyu Wen, Peking Univ, Mar/23/2017

F90=gfortran

%.o : %.f90
	${F90} -c -O3 $<

MPIMOD=mpimod_stub
GUIMOD=guimod_stub
FFTMOD=fftmod

OBJ=${MPIMOD}.o ${FFTMOD}.o ${GUIMOD}.o puma.o legsym.o restartmod.o gaussmod.o

puma.x:	$(OBJ)
	${F90} -o puma.x $(OBJ)
	rm -f *.o *.mod
	cp puma.x ../run

puma.o:		puma.f90
${FFTMOD}.o:	${FFTMOD}.f90
${GUIMOD}.o:	${GUIMOD}.f90
${MPIMOD}.o:	${MPIMOD}.f90 puma.o
legsym.o:	legsym.f90
restartmod.o:	restartmod.f90
gaussmod.o:	gaussmod.f90

clean:
	rm -f *.o *.mod
