HDF5_INC=$(HOME)/drops/boxStrain/src/lib/hdf5_install/include
HDF5_LIB=$(HOME)/drops/boxStrain/src/lib/hdf5_install/lib
LD_LIBRARY_PATH=$HDF5_LIB:$LD_LIBRARY_PATH

FFTW_INC =/home/spack/opt/spack/linux-centos7-broadwell/gcc-12.2.0/fftw-3.3.10-aw4yurpymgyijtqutxl46sk52h2gnbo6/include/
FFTW_LIB =/home/spack/opt/spack/linux-centos7-broadwell/gcc-12.2.0/fftw-3.3.10-aw4yurpymgyijtqutxl46sk52h2gnbo6/lib/

FC	= mpif90
LD	= $(FC)
RM	= /bin/rm -f
OLEVEL	= -O3
DBG := #-g -fbounds-check -Wall -fcheck=all -std=f2008ts 
FOPTS	= -mcmodel=large -fconvert=big-endian  -ffixed-line-length-140 -fno-align-commons -cpp -I$(HDF5_INC) -I$(FFTW_INC)
LIB = -L$(HDF5_LIB) -Wl,-rpath,$(HDF5_LIB) -llapack  -lhdf5_fortran -lhdf5 -L$(FFTW_LIB) -lfftw3 -lfftw3f 
#FOPTS	= -align none -mcmodel medium -warn all
FFLAGS	= $(FOPTS) $(OLEVEL) $(DBG)

FLWOBJS = \
./strainBox.f90  ./fftw3.f03

OBJS	= $(FLWOBJS)
EXEC    =  ./strainBox

$(EXEC):	$(OBJS)
	$(LD) $(FFLAGS) $(OBJS) -o $@ $(LIB)

fftw3.o: $(FFTW_INC)/fftw3.f03
	$(FC) $(FFLAGS) -c $(FFTW_INC)/fftw3.f03

fftw3.mod: $(FFTW_INC)/fftw3.f03
	# Built automatically when compiling fftw3.o

clean:
	$(RM) $(EXEC)
	rm -f *.mod

.SUFFIXES: .o

.f90.o:
	$(FC)  -c $(FFLAGS) $<

%.o: %.f03
	$(FC) $(FFLAGS) -c $<
