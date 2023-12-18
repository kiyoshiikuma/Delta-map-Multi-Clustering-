F90	= ifort
#FFLAGS	= -mkl -warn all -check all -fpe0 -CB -traceback -g -debug parallel
FFLAGS	= -mkl -xHOST -O3 -ipo -no-prec-div -qopt-matmul -qopenmp -heap-arrays
#FFLAGS	= -mkl -O0 -heap-arrays
OBJ	= Utils.o LBird.o main.o
INCLUDE	= -I/home/ichiki/usr/local/include -I/home/ichiki/usr/local/Healpix_3.31/include
HEALPIXL= -L/home/ichiki/usr/local/Healpix_3.31/lib -lhealpix
CFITSIO = -L/home/ichiki/usr/local/lib64 -lcfitsio
NLOPTLIB= -L/usr/local/lib64 -lnlopt
LIBS	= $(HEALPIXL) $(CFITSIO) $(NLOPTLIB) 

default: minimize
minimize: $(OBJ)
	$(F90) $(OBJ) $(FFLAGS) $(LIBS) $(INCLUDE) -o $@

%.o: %.f90
	$(F90) $(FFLAGS) $(LIBS) $(INCLUDE) -c $*.f90
clean:
	-rm -f *.o *.a *.d core *.mod *~
