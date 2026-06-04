TOP=../../..
include $(TOP)/Make.inc
LIBDIR=$(TOP)/lib
INCLUDEDIR=$(TOP)/include
MODDIR=$(TOP)/modules/
HERE=../..

FINCLUDES=$(FMFLAG). $(FMFLAG)$(HERE) $(FMFLAG)$(MODDIR) $(CUDA_INCLUDES)
CINCLUDES=-I. -I$(HERE) -I$(INCLUDEDIR) $(CUDA_INCLUDES)

PSBC_LIBS= -L$(LIBDIR) -lpsb_cbind
# PSB_LIBS=-lpsb_util -lpsb_linsolve -lpsb_prec -lpsb_base -L$(LIBDIR) $(PSBGPULDLIBS)
PSB_LIBS=-lpsb_ext -lpsb_util -lpsb_linsolve -lpsb_prec -lpsb_base -L$(LIBDIR) $(PSBGPULDLIBS) 


#
# Compilers and such
#

EXEDIR=./runs

all: gputest pdedgen3dc

gputest: gputest.o
	$(FLINK) gputest.o  -o gputest $(PSBC_LIBS) \
	$(PSB_LIBS) \
	-lm 
	/bin/mv gputest $(EXEDIR)

pdedgen3dc: pdegen3dc.o
	$(FLINK) pdegen3dc.o  -o pdedgen3dc $(PSBC_LIBS) $(PSB_LIBS) -lstdc++ -lm 
	/bin/mv pdedgen3dc $(EXEDIR)

.f90.o:
	$(MPFC) $(FCOPT) $(FINCLUDES) $(FDEFINES) -c $<
.c.o:
	$(MPCC) $(CCOPT) $(CINCLUDES) $(CDEFINES) -c $<


clean:
	/bin/rm -f  gputest.o  $(EXEDIR)/gputest
	/bin/rm -f  pdedgen3dc.o  $(EXEDIR)/pdedgen3dc
verycleanlib:
	(cd ../..; make veryclean)
lib:
	(cd ../../; make library)

tests: all
	cd runs ; ./gputest
	cd runs ; ./pdedgen3dc < runs/pdegen3d.inp > pdegen3d.out