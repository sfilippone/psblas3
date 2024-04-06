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
PSBLAS_LIB= -L$(LIBDIR) -L$(PSBLIBDIR) $(LCUDA) -lpsb_ext -lpsb_util -lpsb_base 
INCDIR=$(TOPDIR)/include
MODDIR=$(TOPDIR)/modules

LDLIBS=$(PSBGPULDLIBS) 

FINCLUDES=$(FMFLAG)$(MODDIR) $(FMFLAG)$(INCDIR) $(FMFLAG). $(FMFLAG)$(PSBMODDIR)  $(FMFLAG)$(PSBINCDIR) $(LIBRSB_DEFINES)

DPMGOBJS=dpdegenmm.o
EXEDIR=./runs

all: dir pgen file
pgen: dpdegenmm
file: 
dpdegenmm: dir
dir:
	(if test ! -d $(EXEDIR); then mkdir $(EXEDIR); fi)

dpdegenmm: $(DPMGOBJS)
	$(FLINK) $(LOPT) $(DPMGOBJS) -fopenmp -o dpdegenmm $(FINCLUDES) $(PSBLAS_LIB) $(LDLIBS) 
	/bin/mv dpdegenmm $(EXEDIR) 

clean: 
	/bin/rm -f $(DPMGOBJS) $(EXEDIR)/dpdegenmm

lib:
	(cd ../../; make library)
verycleanlib:
	(cd ../../; make veryclean)
