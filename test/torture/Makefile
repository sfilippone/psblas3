include ../../Make.inc
#
# Libraries used
#
LIBDIR=../../lib/
PSBLAS_LIB= -L$(LIBDIR) -lpsb_util -lpsb_krylov -lpsb_prec -lpsb_base 
LDLIBS=$(PSBLDLIBS)
#
# Compilers and such
#
CCOPT= -g
FINCLUDES=$(FMFLAG)$(LIBDIR) $(FMFLAG).

PSBTOBJS=psb_mvsv_tester.o psbtf.o
EXEDIR=./runs


all: psbtf

psbtf: $(PSBTOBJS)
	$(F90LINK)  $(PSBTOBJS)  -o psbtf $(PSBLAS_LIB) $(LDLIBS)
	/bin/mv psbtf $(EXEDIR)

psbtf.o: psb_mvsv_tester.o 
.f90.o:
	$(MPF90) $(F90COPT) $(FINCLUDES) $(FDEFINES) -c $<


clean: 
	/bin/rm -f $(PSBTOBJS) ppde.o spde.o  $(EXEDIR)/ppde 
verycleanlib: 
	(cd ../..; make veryclean)
lib:
	(cd ../../; make library)



