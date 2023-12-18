F90	= ifort 
FFLAGS	= -mkl -qopenmp -heap_arrays -warn all
#FFLAGS	= -warn all -mkl -traceback -check all

OBJ	= SH_mod.o simulatemaps.o main.o 
INCLUDE	= -I/home/ichiki/usr/local/Healpix_3.31/include
HEALPIXL= -L/home/ichiki/usr/local/Healpix_3.31/lib -lhealpix
CFITSIO = -L/home/ichiki/usr/local/lib64 -lcfitsio
LAPACKL = 
LIBS	=  $(HEALPIXL) $(CFITSIO) $(LAPACKL)

default: kk2011

kk2011:	$(OBJ)
	$(F90) $(OBJ) $(FFLAGS) $(LIBS) $(INCLUDE) -o $@

%.o: %.f90
	$(F90) $(FFLAGS) $(LIBS) $(INCLUDE) -c $*.f90

clean:
	-rm -f *.o *.a *.d core *.mod *~
