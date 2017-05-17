# Xinyu Wen, Peking University, Apr/3/2017
# 1. make
# 2. get ready "ppp_namelist"
# 3. get ready "N032_surf_0129.sra"
# 4. run "./ppp.x" to produce another 4 sra files (134, 121, 122, 123)

F90	= gfortran  	# name of the FORTRAN compiler
OPTS	= -O3  		# FORTRAN options

%.o : %.f90
	$(F90) -c $(OPTS) $<

OBJ=ppp.o legsym.o fftmod.o gaussmod.o
ppp.x:	$(OBJ)
	$(F90) -o ppp.x $(OPTS) $(OBJ)

gaussmod.o:	./gaussmod.f90
legsym.o:	./legsym.f90
fftmod.o:	./fftmod.f90
ppp.o:		./ppp.f90

clean:
	rm -f *.o *.mod
