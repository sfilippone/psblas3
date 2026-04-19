INSTALLDIR=../../..
INCDIR=$(INSTALLDIR)/include/
MODDIR=$(INSTALLDIR)/modules/
include $(INCDIR)/Make.inc.psblas
#
# Libraries used
#
LIBDIR=$(INSTALLDIR)/lib/
PSBLAS_LIB= -L$(LIBDIR) -lpsb_util -lpsb_linsolve -lpsb_prec -lpsb_base
LDLIBS=$(PSBLDLIBS)

FINCLUDES=$(FMFLAG)$(MODDIR) $(FMFLAG).

TOBJS=psb_spmv_overlap_test.o spmv_overlap.o

EXEDIR=./runs

all: runsd spmv_overlap

runsd:
	(if test ! -d runs ; then mkdir runs; fi)

spmv_overlap: $(TOBJS)
	$(FLINK) $(LOPT) $(TOBJS) -o spmv_overlap $(PSBLAS_LIB) $(LDLIBS)
	/bin/mv spmv_overlap $(EXEDIR)
	/bin/cp -f $(EXEDIR)/spmv_overlap $(EXEDIR)/psb_spmv_overlap_test

clean:
	/bin/rm -f $(TOBJS) $(TOBJS_API) *$(.mod) $(EXEDIR)/spmv_overlap $(EXEDIR)/psb_spmv_overlap_test

lib:
	(cd ../../; make library)
verycleanlib:
	(cd ../../; make veryclean)