INSTALLDIR=../../..
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

DBFOBJS=psb_dbf_sample.o
DPDEGEN=psb_dpde_gen.o

EXEDIR=./runs

all: runsd psb_dbf_sample psb_dpde_gen

runsd:
	(if test ! -d runs ; then mkdir runs; fi)

psb_dbf_sample: runsd $(DBFOBJS)
	$(FLINK) $(LOPT) $(DBFOBJS) -o psb_dbf_sample $(PSBLAS_LIB) $(LDLIBS)
	/bin/mv psb_dbf_sample $(EXEDIR)

psb_dpde_gen: runsd $(DPDEGEN)
	$(FLINK) $(LOPT) $(DPDEGEN) -o psb_dpde_gen $(PSBLAS_LIB) $(LDLIBS)
	/bin/mv psb_dpde_gen $(EXEDIR)

clean: 
	/bin/rm -f $(DBFOBJS) $(DPDEGEN)\
        *$(.mod) $(EXEDIR)/psb_*

lib:
	(cd ../../; make library)
verycleanlib:
	(cd ../../; make veryclean)
