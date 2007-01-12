include ../../Make.inc
#
# Libraries used
#
LIBDIR=../../lib/
PSBLAS_LIB= -L$(LIBDIR) -lpsb_util -lpsb_krylov -lpsb_prec -lpsb_base 

#
# Compilers and such
#
CCOPT= -g
INCDIRS=-I$(LIBDIR)


EXEDIR=./runs


ppde90: ppde90.o  
	$(F90LINK) ppde90.o  -o ppde90 $(PSBLAS_LIB) $(LDLIBS)
	/bin/mv ppde90 $(EXEDIR)


.f90.o:
	$(MPF90) $(F90COPT) $(INCDIRS) -c $<


clean: 
	/bin/rm -f ppde90.o   $(EXEDIR)/ppde90 
verycleanlib: 
	(cd ../..; make veryclean)
lib:
	(cd ../../; make library)



