.SUFFIXES:
.SUFFIXES: .F90 .f90 .o .s .c

# Compilers and flags
CC=mpicc
FC=mpif90
FCOPT=-O0 -march=native
OFFLOAD=-fopenacc -foffload=nvptx-none="-march=sm_70"

# Directories
LIBDIR=../lib
INCDIR=../include
MODDIR=../modules
IMPLDIR=./impl  # Adding the impl directory

# Include and library paths
INCLUDES=-I$(LIBDIR) -I$(INCDIR) -I$(MODDIR) -I$(IMPLDIR)
LIBS=-L$(LIBDIR) -lpsb_util -lpsb_ext -lpsb_base -lopenblas -lmetis

# Source files
FOBJS= psb_i_oacc_vect_mod.o psb_d_oacc_vect_mod.o \
	psb_oacc_mod.o  psb_d_oacc_csr_mat_mod.o \
	impl/psb_d_oacc_csr_vect_mv.o 

# Library name
LIBNAME=libpsb_openacc.a

# Rules
all: $(LIBNAME)

$(LIBNAME): $(FOBJS)
	ar cr $(LIBNAME) $(FOBJS)
	/bin/cp -p $(LIBNAME) $(LIBDIR)

clean:
	/bin/rm -fr *.o $(LIBNAME) *.mod impl/*.o

.f90.o:
	$(FC) $(FCOPT) $(OFFLOAD) $(INCLUDES) -c $< -o $@

.c.o:
	$(CC) -c $< -o $@

.F90.o:
	$(FC) $(FCOPT) $(OFFLOAD) $(INCLUDES) -c $< -o $@

.F90.s:
	$(FC) $(FCOPT) $(INCLUDES) -c -S $< -o $@
