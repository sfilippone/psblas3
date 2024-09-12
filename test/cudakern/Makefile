TOPDIR=../..
include $(TOPDIR)/Make.inc
#
# Libraries used
#
LIBDIR=$(TOPDIR)/lib/
PSBLIBDIR=$(TOPDIR)/lib/
OPTDIR=$(LIBDIR)
PSBINCDIR=$(TOPDIR)/include
PSBMODDIR=$(TOPDIR)/modules
PSBLAS_LIB= -L$(LIBDIR) -L$(PSBLIBDIR) $(LCUDA) -lpsb_ext -lpsb_util -lpsb_base 
INCDIR=$(TOPDIR)/include
MODDIR=$(TOPDIR)/modules

LDLIBS=$(PSBGPULDLIBS) 

FINCLUDES=$(FMFLAG)$(MODDIR) $(FMFLAG)$(INCDIR) $(FMFLAG). $(FMFLAG)$(PSBMODDIR)  $(FMFLAG)$(PSBINCDIR) $(LIBRSB_DEFINES)

DPGOBJS=dpdegenmv.o data_input.o
SPGOBJS=spdegenmv.o data_input.o
EXEDIR=./runs

all: dir pgen 
pgen: dpdegenmv spdegenmv 
dpdegenmv spdegenmv: dir 
dir:
	(if test ! -d $(EXEDIR); then mkdir $(EXEDIR); fi)
dpdegenmv.o spdegenmv.o: data_input.o
dpdegenmv: $(DPGOBJS) 
	$(FLINK) $(LOPT) $(DPGOBJS) -fopenmp -o dpdegenmv $(FINCLUDES) $(PSBLAS_LIB) $(LDLIBS) 
	/bin/mv dpdegenmv $(EXEDIR) 
spdegenmv: $(SPGOBJS)
	$(FLINK) $(LOPT) $(SPGOBJS) -o spdegenmv $(PSBLAS_LIB) $(LDLIBS) 
	/bin/mv spdegenmv $(EXEDIR) 

clean: 
	/bin/rm -f $(DTOBJS) $(STOBJS) $(DPGOBJS) $(SPGOBJS) $(ZTOBJS) $(CTOBJS)  \
	$(EXEDIR)/dpdegenmv $(EXEDIR)/spdegenmv *mod

lib:
	(cd ../../; make library)
verycleanlib:
	(cd ../../; make veryclean)
