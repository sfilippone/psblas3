TOP=../../..
include $(TOP)/Make.inc
LIBDIR=$(TOP)/lib
INCLUDEDIR=$(TOP)/include
HERE=../..

FINCLUDES=$(FMFLAG). $(FMFLAG)$(HERE) $(FMFLAG)$(INCLUDEDIR) 
CINCLUDES=-I. -I$(HERE) -I$(INCLUDEDIR)

PSBC_LIBS= -L$(LIBDIR) -lpsb_cbind
PSB_LIBS=-lpsb_prec -lpsb_base -L$(LIBDIR)
#    -lpsb_krylov_cbind

#
# Compilers and such
#
CCOPT= -g
#FINCLUDES=$(FMFLAG)$(LIBDIR) $(FMFLAG).
CINCLUDES=-I$(LIBDIR) $(FIFLAG)$(INCLUDEDIR) $(FIFLAG)$(PSBLAS_INCDIR)

EXEDIR=./runs

all: ppdec

ppdec: ppdec.o
	$(MPFC) ppdec.o  -o ppdec $(PSBC_LIBS) $(PSB_LIBS) $(PSBLDLIBS) -lm -lgfortran 
# \
#	-lifcore -lifcoremt -lguide -limf -lirc -lintlc -lcxaguard -L/opt/intel/fc/10.0.023/lib/ -lm
	/bin/mv ppdec $(EXEDIR)

.f90.o:
	$(MPFC) $(FCOPT) $(FINCLUDES) $(FDEFINES) -c $<
.c.o:
	$(MPCC) $(CCOPT) $(CINCLUDES) $(CDEFINES) -c $<


clean: 
	/bin/rm -f  ppdec.o  $(EXEDIR)/ppdec
verycleanlib: 
	(cd ../..; make veryclean)
lib:
	(cd ../../; make library)

tests: all
	cd runs ; ./ppdec < ppde.inp


