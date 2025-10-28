TOPDIR=../..
include $(TOPDIR)/Make.inc

LIBDIR=$(TOPDIR)/lib/
PSBLIBDIR=$(TOPDIR)/lib/
PSBINCDIR=$(TOPDIR)/include
PSBMODDIR=$(TOPDIR)/modules
INCDIR=$(TOPDIR)/include
MODDIR=$(TOPDIR)/modules
EXEDIR=./runs

PSBLAS_LIB= -L$(LIBDIR) -L$(PSBLIBDIR) -lpsb_openacc -lpsb_ext -lpsb_util -lpsb_linsolve -lpsb_prec -lpsb_base 
LDLIBS=$(PSBGPULDLIBS)

FINCLUDES= $(FMFLAG). $(FMFLAG)$(MODDIR) $(FMFLAG)$(INCDIR) $(FMFLAG)$(PSBMODDIR) $(FMFLAG)$(PSBINCDIR) $(LIBRSB_DEFINES)

FFLAGS=-O3 -march=native $(FCOPENACC)  -DPSB_OPENACC -DPSB_HAVE_CUDA
CFLAGS=-O3 -march=native

VTC=vectoacc.o
DVT=datavect.o
CSRC=timers.c

OBJS=$(SRCS:.F90=.o) $(CSRC:.c=.o)

all: dir dpdegenmv

#$(OBJS)
#	$(FC) $(FFLAGS) $(OBJS) -o datavect $(FINCLUDES) $(PSBLAS_LIB) $(LDLIBS)
#	/bin/mv datavect $(EXEDIR)

dir:
	@if test ! -d $(EXEDIR); then mkdir $(EXEDIR); fi

%: %.o timers.o
	$(FC) $(FFLAGS) $^ -o $@ $(FINCLUDES) $(PSBLAS_LIB) $(LDLIBS)
	/bin/mv $@ $(EXEDIR)

%.o: %.F90
	$(FC) $(FFLAGS) $(FINCLUDES) -c $< -o $@

%.o: %.c
	$(CC) $(CFLAGS) $(FINCLUDES) -c $< -o $@


dpdegenmv.o:  psb_d_pde3d_mod.o
dpdegenmv: dpdegenmv.o psb_d_pde3d_mod.o
	$(FLINK)  $(FCOPENACC) dpdegenmv.o  psb_d_pde3d_mod.o -o dpdegenmv -fuse-ld=mold  $(PSBLAS_LIB) $(LDLIBS)
	/bin/mv dpdegenmv $(EXEDIR)

clean:
	/bin/rm -fr *.o *.mod $(EXEDIR)/*

.PHONY: all dir clean
