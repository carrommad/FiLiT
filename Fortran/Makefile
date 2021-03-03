SHELL=/bin/bash
FC=mpifort
LINK=mpifort

FFLAGS=-Dplplot -DNETCDF -cpp -ffree-line-length-none -g -fbacktrace -fdump-core #-qopenmp
LDFLAGS=-ffree-line-length-none -g -fbacktrace -fdump-core #-qopenmp 

# PSFT-crmadrid
LIBS=-L/home/crmadrid/apps/netcdf/lib -lnetcdff
MOD=-I/home/crmadrid/apps/netcdf/include

# CICA
#LIBS=-L/home/software/libraries/netcdf-fortran-4.5.3/lib -lnetcdff
#MOD=-I/home/software/libraries/netcdf-fortran-4.5.3/include


PROGRAM=tracer.exe
OBJ=general_mod.o\
field_line_mod.o\
field_line_tracer.o


.SUFFIXES:
.SUFFIXES: .f90 .o

.f90.o:
	$(FC) -c $(MOD) $(FFLAGS) $<

$(PROGRAM): $(OBJ)
	$(LINK) $(OBJ) $(LDFLAGS) $(LIBS) -o $(PROGRAM)

clean:
	rm -f *.o *.mod *.log *~ core

