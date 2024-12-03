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

FINCLUDES=$(FMFLAG)$(MODDIR) $(FMFLAG)$(INCDIR) $(FMFLAG). $(FMFLAG)$(PSBMODDIR) $(FMFLAG)$(PSBINCDIR) $(LIBRSB_DEFINES)
FDEFINES=-DHAVE_OPENACC
FFLAGS=-O3 -march=native -fopenacc -foffload=nvptx-none="-march=sm_75" -DHAVE_OPENACC
CFLAGS=-O3 -march=native

VTC=vectoacc.o
DVT=datavect.o
CSRC=timers.c

OBJS=$(SRCS:.F90=.o) $(CSRC:.c=.o)
DPGOBJS=dpdegenmv.o #data_input.o

all: dir psb_d_oacc_pde3d dpdegenmv 
dpdegenmv: dir 


#$(OBJS)
#	$(FC) $(FFLAGS) $(OBJS) -o datavect $(FINCLUDES) $(PSBLAS_LIB) $(LDLIBS)
#	/bin/mv datavect $(EXEDIR)

dir:
	@if test ! -d $(EXEDIR); then mkdir $(EXEDIR); fi
#dpdegenmv.o spdegenmv.o: data_input.o
dpdegenmv: $(DPGOBJS) 
	$(FLINK) $(LOPT) -fopenacc $(DPGOBJS) -o dpdegenmv $(FINCLUDES) $(PSBLAS_LIB) $(LDLIBS) 
	/bin/mv dpdegenmv $(EXEDIR) 

%: %.o timers.o
	$(FC) $(FFLAGS) $^ -o $@ $(FINCLUDES) $(PSBLAS_LIB) $(LDLIBS)
	/bin/mv $@ $(EXEDIR)

%.o: %.F90
	$(FC) $(FFLAGS) $(FINCLUDES) -c $< -o $@

%.o: %.c
	$(CC) $(CFLAGS) $(FINCLUDES) -c $< -o $@

psb_d_oacc_pde3d:
	mpifort -fallow-argument-mismatch -frecursive -g -O3 -frecursive  -I../../modules/ -I. -DOPENACC -DHAVE_LAPACK -DHAVE_FLUSH_STMT -DLPK8 -DIPK4 -DMPI_MOD   -c psb_d_oacc_pde3d.F90 -o psb_d_oacc_pde3d.o
	$(FLINK) -fopenacc -DOPENACC psb_d_oacc_pde3d.o -o psb_d_oacc_pde3d $(PSBLAS_LIB) $(LDLIBS)
	/bin/mv psb_d_oacc_pde3d $(EXEDIR)

clean:
	/bin/rm -fr *.o *.mod $(EXEDIR)/*

.PHONY: all dir clean
