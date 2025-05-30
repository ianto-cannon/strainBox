HDF5_INC=$(HOME)/drops/boxStrain/src/lib/hdf5_install/include
HDF5_LIB=$(HOME)/drops/boxStrain/src/lib/hdf5_install/lib
LD_LIBRARY_PATH=$HDF5_LIB:$LD_LIBRARY_PATH
FC	= mpif90
LD	= $(FC)
RM	= /bin/rm -f
OLEVEL	= -O3
DBG := #-g -fbounds-check -Wall -fcheck=all -std=f2008ts 
FOPTS	= -mcmodel=large -fconvert=big-endian  -ffixed-line-length-140 -fno-align-commons -cpp -I$(HDF5_INC)
LIB = -L$(HDF5_LIB)  -Wl,-rpath,$(HDF5_LIB) -llapack  -lhdf5_fortran -lhdf5 -lfftw3 -lfftw3f 
#FOPTS	= -align none -mcmodel medium -warn all
FFLAGS	= $(FOPTS) $(OLEVEL) $(DBG)

FLWOBJS = \
./strainBox.f90 

OBJS	= $(FLWOBJS)
EXEC    =  ./strainBox

$(EXEC):	$(OBJS)
	$(LD) $(FFLAGS) $(OBJS) -o $@ $(LIB)

clean:
	$(RM) $(EXEC)
	rm -f *.mod

.SUFFIXES: .o

.f90.o:
	$(FC)  -c $(FFLAGS) $<
