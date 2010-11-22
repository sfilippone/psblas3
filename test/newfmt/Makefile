include ../../Make.inc
#
# Libraries used
#
LIBDIR=../../lib/
OPTDIR=../../opt
PSBLAS_LIB= -L$(LIBDIR) -lpsb_util -lpsb_krylov -lpsb_prec -lpsb_base -L$(OPTDIR) -lopt
LDLIBS=$(PSBLDLIBS)
#
# Compilers and such
#
CCOPT= -g
FINCLUDES=$(FMFLAG)$(LIBDIR) $(FMFLAG). $(FMFLAG)$(OPTDIR)


EXEDIR=./runs

all: ppde spde

ppde: ppde.o 
	$(F90LINK)  ppde.o -o ppde $(PSBLAS_LIB) $(LDLIBS)
	/bin/mv ppde $(EXEDIR)


spde: spde.o 
	$(F90LINK) spde.o  -o spde $(PSBLAS_LIB) $(LDLIBS)
	/bin/mv spde $(EXEDIR)


clean: 
	/bin/rm -f ppde.o spde.o  $(EXEDIR)/ppde 
verycleanlib: 
	(cd ../..; make veryclean)
lib:
	(cd ../../; make library)



