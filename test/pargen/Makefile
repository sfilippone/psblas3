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

all: ppde spde

ppde: ppde.o  psb_d_csc_impl.o psb_d_csc_mat_mod.o
	$(F90LINK) ppde.o  psb_d_csc_impl.o psb_d_csc_mat_mod.o -o ppde $(PSBLAS_LIB) $(LDLIBS)
	/bin/mv ppde $(EXEDIR)

psb_d_csc_impl.o ppde.o: psb_d_csc_mat_mod.o

spde: spde.o  
	$(F90LINK) spde.o  -o spde $(PSBLAS_LIB) $(LDLIBS)
	/bin/mv spde $(EXEDIR)

.f90.o:
	$(MPF90) $(F90COPT) $(FINCLUDES) $(FDEFINES) -c $<


clean: 
	/bin/rm -f ppde.o spde.o  $(EXEDIR)/ppde 
verycleanlib: 
	(cd ../..; make veryclean)
lib:
	(cd ../../; make library)



