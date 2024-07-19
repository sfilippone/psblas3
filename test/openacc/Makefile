TOPDIR=../..
include $(TOPDIR)/Make.inc

LIBDIR=$(TOPDIR)/lib/
PSBLIBDIR=$(TOPDIR)/lib/
PSBINCDIR=$(TOPDIR)/include
PSBMODDIR=$(TOPDIR)/modules
INCDIR=$(TOPDIR)/include
MODDIR=$(TOPDIR)/modules
EXEDIR=./runs

PSBLAS_LIB= -L$(LIBDIR) -L$(PSBLIBDIR) -lpsb_openacc -lpsb_base -lpsb_ext -lpsb_util -lpsb_krylov -lpsb_prec -lopenblas -lmetis
LDLIBS=$(PSBGPULDLIBS)

FINCLUDES=$(FMFLAG)$(MODDIR) $(FMFLAG)$(INCDIR) $(FMFLAG). $(FMFLAG)$(PSBMODDIR) $(FMFLAG)$(PSBINCDIR) $(LIBRSB_DEFINES)

FFLAGS=-O0 -march=native -fopenacc -foffload=nvptx-none="-march=sm_70"
CFLAGS=-O0 -march=native

SRCS=vectoacc.F90 datavect.F90
CSRC=timers.c

OBJS=$(SRCS:.F90=.o) $(CSRC:.c=.o)

all: dir $(OBJS)
	$(FC) $(FFLAGS) $(OBJS) -o datavect $(FINCLUDES) $(PSBLAS_LIB) $(LDLIBS)
	/bin/mv datavect $(EXEDIR)

dir:
	@if test ! -d $(EXEDIR); then mkdir $(EXEDIR); fi

%: %.o timers.o
	$(FC) $(FFLAGS) $^ -o $@ $(FINCLUDES) $(PSBLAS_LIB) $(LDLIBS)
	/bin/mv $@ $(EXEDIR)

%.o: %.F90
	$(FC) $(FFLAGS) $(FINCLUDES) -c $< -o $@

%.o: %.c
	$(CC) $(CFLAGS) $(FINCLUDES) -c $< -o $@

clean:
	/bin/rm -fr *.o *.mod $(EXEDIR)/*

.PHONY: all dir clean
