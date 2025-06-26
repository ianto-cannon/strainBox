HDF5_INC = $(HOME)/drops/boxStrain/src/lib/hdf5_parallel/include
HDF5_LIB = $(HOME)/drops/boxStrain/src/lib/hdf5_parallel/lib

# Compiler and flags
FC = mpif90
DBG = #-O0 -Wall -fcheck=all -g #-fsanitize=address -fno-omit-frame-pointer
FFLAGS = -mcmodel=large -fconvert=big-endian -ffixed-line-length-140 -fno-align-commons -cpp \
         -I$(HDF5_INC) -O3 $(DBG)

# Libraries
LIBS = -L$(HDF5_LIB) -Wl,-rpath,$(HDF5_LIB) \
       -llapack -lhdf5_fortran -lhdf5 -lfftw3 -lfftw3f

# Sources and objects
OBJS = main.o
EXE = spectrum

# Rules
all: $(EXE)

$(EXE): $(OBJS)
	$(FC) $(FFLAGS) -o $@ $^ $(LIBS)

%.o: %.f90
	$(FC) $(FFLAGS) -c $<

clean:
	rm -f $(EXE) *.mod *.o

