TOPDIR=../..
include $(TOPDIR)/Make.inc

# Directories
LIBDIR=$(TOPDIR)/lib/
PSBLIBDIR=$(TOPDIR)/lib/
PSBINCDIR=$(TOPDIR)/include
PSBMODDIR=$(TOPDIR)/modules
INCDIR=$(TOPDIR)/include
MODDIR=$(TOPDIR)/modules
EXEDIR=./runs

# Libraries
PSBLAS_LIB= -L$(LIBDIR) -L$(PSBLIBDIR) -lpsb_util -lpsb_ext -lpsb_base -lpsb_openacc -lopenblas -lmetis
LDLIBS=$(PSBGPULDLIBS)

# Includes
FINCLUDES=$(FMFLAG)$(MODDIR) $(FMFLAG)$(INCDIR) $(FMFLAG). $(FMFLAG)$(PSBMODDIR) $(FMFLAG)$(PSBINCDIR) $(LIBRSB_DEFINES)

# Compiler flags
FFLAGS=-O0 -march=native -fopenacc -foffload=nvptx-none="-march=sm_70"
CFLAGS=-O0 -march=native

# Source files
SRCS=test.F90 vectoacc.F90
CSRC=timers.c

# Object files
OBJS=$(SRCS:.F90=.o) $(CSRC:.c=.o)

# Default rule
all: dir

dir:
	@if test ! -d $(EXEDIR); then mkdir $(EXEDIR); fi

# Pattern rule for creating executables
%: %.o timers.o
	$(FC) $(FFLAGS) $^ -o $@ $(FINCLUDES) $(PSBLAS_LIB) $(LDLIBS)
	/bin/mv $@ $(EXEDIR)

# Compilation rules
%.o: %.F90
	$(FC) $(FFLAGS) $(FINCLUDES) -c $< -o $@

%.o: %.c
	$(CC) $(CFLAGS) $(FINCLUDES) -c $< -o $@

clean:
	/bin/rm -fr *.o *.mod $(EXEDIR)/* 

# Phony targets
.PHONY: all dir clean
