#
# Libraries used
#
INSTALLDIR=../..
INCDIR=$(INSTALLDIR)/include/
MODDIR=$(INSTALLDIR)/modules/
include $(INCDIR)/Make.inc.psblas
LIBDIR=$(INSTALLDIR)/lib/
PSBLAS_LIB= -L$(LIBDIR) -lpsb_util -lpsb_krylov -lpsb_prec -lpsb_base
LDLIBS=$(PSBLDLIBS)

FINCLUDES=$(FMFLAG)$(MODDIR) $(FMFLAG).


DTOBJS=d_file_spmv.o
STOBJS=s_file_spmv.o
DPGOBJS=pdgenspmv.o
DVECOBJS=vecoperation.o
EXEDIR=./runs

all: runsd pdgenspmv
#d_file_spmv s_file_spmv pdgenspmv vecoperation

runsd:
	(if test ! -d runs ; then mkdir runs; fi)

d_file_spmv: $(DTOBJS)
	$(FLINK) $(LOPT) $(DTOBJS) -o d_file_spmv $(PSBLAS_LIB) $(LDLIBS)
	/bin/mv d_file_spmv $(EXEDIR)

pdgenspmv: $(DPGOBJS)
	$(FLINK) $(LOPT) $(DPGOBJS) -o pdgenspmv $(PSBLAS_LIB) $(LDLIBS)
	/bin/mv pdgenspmv $(EXEDIR)


s_file_spmv: $(STOBJS)
	$(FLINK) $(LOPT) $(STOBJS) -o s_file_spmv $(PSBLAS_LIB) $(LDLIBS)
	/bin/mv s_file_spmv $(EXEDIR)

vecoperation: $(DVECOBJS)
	$(FLINK) $(LOPT) $(DVECOBJS) -o vecoperation $(PSBLAS_LIB) $(LDLIBS)
	/bin/mv vecoperation $(EXEDIR)

clean:
	/bin/rm -f $(DBOBJSS) $(DPGOBJS) $(DTOBJS) $(STOBJS) $(DVECOBJS)

lib:
	(cd ../../; make library)
verycleanlib:
	(cd ../../; make veryclean)
