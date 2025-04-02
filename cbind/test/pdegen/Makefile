TOP=../../..
include $(TOP)/Make.inc
LIBDIR=$(TOP)/lib
INCLUDEDIR=$(TOP)/include
MODDIR=$(TOP)/modules/
HERE=../..

FINCLUDES=$(FMFLAG). $(FMFLAG)$(HERE) $(FMFLAG)$(MODDIR)
CINCLUDES=-I. -I$(HERE) -I$(INCLUDEDIR)

PSBC_LIBS= -L$(LIBDIR) -lpsb_cbind
PSB_LIBS=-lpsb_util -lpsb_linsolve -lpsb_prec -lpsb_base -L$(LIBDIR)

#
# Compilers and such
#

EXEDIR=./runs

all: pdegen3dc 

pdegen3dc: pdegen3dc.o
	$(MPFC) pdegen3dc.o  -o pdegen3dc $(PSBC_LIBS) $(PSB_LIBS) $(PSBLDLIBS) -lm -lgfortran
	/bin/mv pdegen3dc $(EXEDIR)


.f90.o:
	$(MPFC) $(FCOPT) $(FINCLUDES) $(FDEFINES) -c $<
.c.o:
	$(MPCC) $(CCOPT) $(CINCLUDES) $(CDEFINES) -c $<


clean:
	/bin/rm -f  pdegen3dc.o  $(EXEDIR)/pdegen3dc
verycleanlib:
	(cd ../..; make veryclean)
lib:
	(cd ../../; make library)

tests: all
	cd runs ; ./pdegen3dc < pdegen3d.inp
