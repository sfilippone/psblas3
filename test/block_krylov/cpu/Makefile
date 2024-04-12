INSTALLDIR=../..
INCDIR=$(INSTALLDIR)/include/
MODDIR=$(INSTALLDIR)/modules/
include $(INCDIR)/Make.inc.psblas
#
# Libraries used
#
LIBDIR=$(INSTALLDIR)/lib/
PSBLAS_LIB= -L$(LIBDIR) -lpsb_util -lpsb_krylov -lpsb_prec -lpsb_base
LDLIBS=$(PSBLDLIBS)

FINCLUDES=$(FMFLAG)$(MODDIR) $(FMFLAG). 

DBFOBJS=getp.o psb_dbf_sample.o

EXEDIR=./runs

all: runsd psb_dbf_sample

runsd:
	(if test ! -d runs ; then mkdir runs; fi)

psb_dbf_sample.o: getp.o

psb_dbf_sample: $(DBFOBJS)
	$(FLINK) $(LOPT) $(DBFOBJS) -o psb_dbf_sample $(PSBLAS_LIB) $(LDLIBS)
	/bin/mv psb_dbf_sample $(EXEDIR)

.f90.o:
	$(MPFC) $(FCOPT) $(FINCLUDES) $(FDEFINES) -c $<

clean: 
	/bin/rm -f $(DBFOBJS)\
        *$(.mod) $(EXEDIR)/psb_*f_sample 

lib:
	(cd ../../; make library)
verycleanlib:
	(cd ../../; make veryclean)
