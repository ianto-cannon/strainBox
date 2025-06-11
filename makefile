HDF5_INC = $(HOME)/drops/boxStrain/src/lib/hdf5_install/include
HDF5_LIB = $(HOME)/drops/boxStrain/src/lib/hdf5_install/lib

# Compiler and flags
FC = mpif90
FFLAGS = -mcmodel=large -fconvert=big-endian -ffixed-line-length-140 -fno-align-commons -cpp \
         -I$(HDF5_INC) -g \
         -O0 -Wall -fcheck=all -fsanitize=address -fno-omit-frame-pointer
         #-O3 
# Libraries
LIBS = -L$(HDF5_LIB) -Wl,-rpath,$(HDF5_LIB) \
       -llapack -lhdf5_fortran -lhdf5 -lfftw3 -lfftw3f

# Sources and executable
SRC = ./strainBox.f90
EXE = strainBox

# Rules
all: $(EXE)

$(EXE): $(SRC)
	$(FC) $(FFLAGS) $^ -o $@ $(LIBS)

clean:
	rm -f $(EXE) *.mod *.o
