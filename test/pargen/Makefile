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


EXEDIR=./runs


ppde: ppde.o  
	$(F90LINK) ppde.o  -o ppde $(PSBLAS_LIB) $(LDLIBS)
	/bin/mv ppde $(EXEDIR)

tpde: tpde.o  
	$(F90LINK) tpde.o -pg  -o tpde $(PSBLAS_LIB) $(LDLIBS)
	/bin/mv tpde $(EXEDIR)

.f90.o:
	$(MPF90) $(F90COPT) $(FINCLUDES) $(FDEFINES) -c $<


clean: 
	/bin/rm -f ppde.o   $(EXEDIR)/ppde 
verycleanlib: 
	(cd ../..; make veryclean)
lib:
	(cd ../../; make library)



