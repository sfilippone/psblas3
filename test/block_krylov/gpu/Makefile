TOPDIR=../../..
include $(TOPDIR)/Make.inc
#
# Libraries used
#
LIBDIR=$(TOPDIR)/lib/
PSBLIBDIR=$(TOPDIR)/lib/
OPTDIR=$(LIBDIR)
PSBINCDIR=$(TOPDIR)/include
PSBMODDIR=$(TOPDIR)/modules
PSBLAS_LIB= -L$(LIBDIR) -L$(PSBLIBDIR) $(LCUDA) -lpsb_ext -lpsb_util -lpsb_krylov -lpsb_prec -lpsb_base
INCDIR=$(TOPDIR)/include
MODDIR=$(TOPDIR)/modules

LDLIBS=$(PSBLDLIBS)
LDLIBSG=$(PSBGPULDLIBS)

FINCLUDES=$(FMFLAG)$(MODDIR) $(FMFLAG)$(INCDIR) $(FMFLAG). $(FMFLAG)$(PSBMODDIR) $(FMFLAG)$(PSBINCDIR) $(LIBRSB_DEFINES)

DBFOBJS=getp.o psb_dbf_sample.o
DCOMP=psb_dbf_compare.o
DPDEGEN=psb_dpde_gen.o
EXEDIR=./runs

all: dir psb_dbf_sample psb_dbf_compare psb_dpde_gen

dir:
	(if test ! -d $(EXEDIR); then mkdir $(EXEDIR); fi)

psb_dbf_sample.o: getp.o

psb_dbf_sample: dir $(DBFOBJS)
	$(FLINK) $(LOPT) $(DBFOBJS) -fopenmp -o psb_dbf_sample  $(FINCLUDES) $(PSBLAS_LIB) $(LDLIBSG)
	/bin/mv psb_dbf_sample $(EXEDIR)

psb_dbf_compare: $(DCOMP)
	$(FLINK) $(LOPT) $(DCOMP) -fopenmp -o psb_dbf_compare  $(FINCLUDES) $(PSBLAS_LIB) $(LDLIBSG)
	/bin/mv psb_dbf_compare $(EXEDIR)

psb_dpde_gen: dir $(DPDEGEN)
	$(FLINK) $(LOPT) $(DPDEGEN) -fopenmp -o psb_dpde_gen  $(FINCLUDES) $(PSBLAS_LIB) $(LDLIBSG)
	/bin/mv psb_dpde_gen $(EXEDIR)

clean: 
	/bin/rm -f $(DBFOBJS) $(DCOMP) $(DPDEGEN) *$(.mod) $(EXEDIR)/psb_*

lib:
	(cd ../../; make library)
verycleanlib:
	(cd ../../; make veryclean)
