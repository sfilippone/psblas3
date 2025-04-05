INSTALLDIR=../..
INCDIR=$(INSTALLDIR)/include/
MODDIR=$(INSTALLDIR)/modules/
include $(INCDIR)/Make.inc.psblas
#
# Libraries used
#
LIBDIR=$(INSTALLDIR)/lib/
PSBLAS_LIB= -L$(LIBDIR) -lpsb_util -lpsb_prec -lpsb_base
LDLIBS=$(PSBLDLIBS)

FINCLUDES=$(FMFLAG)$(MODDIR) $(FMFLAG).

OBJS=psb_spmm_test.o   

EXEDIR=./runs

all: runsd psb_spmm_test

runsd:
	(if test ! -d runs ; then mkdir runs; fi)

psb_spmm_test: $(OBJS)
	$(FLINK) $(LOPT) $(OBJS) -o psb_spmm_test $(PSBLAS_LIB) $(LDLIBS) 
	/bin/mv psb_spmm_test $(EXEDIR) 

clean: 
	/bin/rm -f $(OBJS)\
        *$(.mod) $(EXEDIR)/psb_spmm_test 

lib:
	(cd ../../; make library)
verycleanlib:
	(cd ../../; make veryclean)

