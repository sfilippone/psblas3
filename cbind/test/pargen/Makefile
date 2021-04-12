TOP=../../..
include $(TOP)/Make.inc
LIBDIR=$(TOP)/lib
INCLUDEDIR=$(TOP)/include
MODDIR=$(TOP)/modules/
HERE=../..

FINCLUDES=$(FMFLAG). $(FMFLAG)$(HERE) $(FMFLAG)$(MODDIR)
CINCLUDES=-I. -I$(HERE) -I$(INCLUDEDIR)

PSBC_LIBS= -L$(LIBDIR) -lpsb_cbind
PSB_LIBS=-lpsb_util -lpsb_krylov -lpsb_prec -lpsb_base -L$(LIBDIR)

#
# Compilers and such
#

EXEDIR=./runs

all: ppdec build

ppdec: ppdec.o
	$(MPFC) ppdec.o  -o ppdec $(PSBC_LIBS) $(PSB_LIBS) $(PSBLDLIBS) -lm -lgfortran
	/bin/mv ppdec $(EXEDIR)

build: build.o
		$(MPFC) build.o  -o build $(PSBC_LIBS) $(PSB_LIBS) $(PSBLDLIBS) -lm -lgfortran
		/bin/mv build $(EXEDIR)
# \
#	-lifcore -lifcoremt -lguide -limf -lirc -lintlc -lcxaguard -L/opt/intel/fc/10.0.023/lib/ -lm

.f90.o:
	$(MPFC) $(FCOPT) $(FINCLUDES) $(FDEFINES) -c $<
.c.o:
	$(MPCC) $(CCOPT) $(CINCLUDES) $(CDEFINES) -c $<


clean:
	/bin/rm -f  ppdec.o build.o  $(EXEDIR)/ppdec
verycleanlib:
	(cd ../..; make veryclean)
lib:
	(cd ../../; make library)

tests: all
	cd runs ; ./ppdec < ppde.inp
