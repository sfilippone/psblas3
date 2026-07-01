INSTALLDIR=../../..
INCDIR=$(INSTALLDIR)/include/
MODDIR=$(INSTALLDIR)/modules/
include $(INCDIR)/Make.inc.psblas
#
# Libraries used
#
LIBDIR=$(INSTALLDIR)/lib/
PSBLAS_LIB= -L$(LIBDIR) -lpsb_util -lpsb_linsolve -lpsb_prec -lpsb_cuda -lpsb_ext -lpsb_base
LDLIBS=$(PSBGPULDLIBS)

FINCLUDES=$(FMFLAG)$(MODDIR) $(FMFLAG).

TOBJS=psb_spmv_test.o

EXEDIR=./runs

all: runsd psb_spmv_kernel

runsd:
	(if test ! -d runs ; then mkdir runs; fi)

psb_spmv_test.o: psb_spmv_test.f90
	$(FC) $(FCOPT) -cpp $(FINCLUDES) $(FDEFINES) $(FCUDEFINES) -c $< -o $@

psb_spmv_kernel: $(TOBJS)
	$(FLINK) $(LOPT) $(TOBJS) -o psb_spmv_kernel $(PSBLAS_LIB) $(LDLIBS)
	/bin/mv psb_spmv_kernel $(EXEDIR)

clean:
	/bin/rm -f $(TOBJS) $(TOBJS_API) *$(.mod) $(EXEDIR)/spmv_overlap $(EXEDIR)/psb_spmv_kernel

lib:
	(cd ../../; make library)
verycleanlib:
	(cd ../../; make veryclean)