INSTALLDIR=../..
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

TOBJS=psb_comm_test.o

EXEDIR=./runs

all: runsd psb_comm_test 

runsd:
	(if test ! -d runs ; then mkdir runs; fi)

psb_comm_test: $(TOBJS)
	$(FLINK) $(LOPT) $(TOBJS) -o psb_comm_test $(PSBLAS_LIB) $(LDLIBS)
	/bin/mv psb_comm_test $(EXEDIR)

clean:
	/bin/rm -f $(TOBJS) $(TOBJS_API) *$(.mod) $(EXEDIR)/psb_comm_test

lib:
	(cd ../../; make library)
verycleanlib:
	(cd ../../; make veryclean)