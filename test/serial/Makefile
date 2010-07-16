include ../../Make.inc
#
# Libraries used
#
LIBDIR=../../lib/
PSBLAS_LIB= -L$(LIBDIR) -lpsb_base
#-lpsb_util -lpsb_krylov -lpsb_prec -lpsb_base 
LDLIBS=$(PSBLDLIBS)
#
# Compilers and such
#
CCOPT= -g
FINCLUDES=$(FMFLAG)$(LIBDIR) $(FMFLAG).


EXEDIR=./runs

all: d_coo_matgen d_matgen 
d_coo_matgen: d_coo_matgen.o
	$(F90LINK) $(LINKOPT) d_coo_matgen.o  -o d_coo_matgen  $(PSBLAS_LIB) $(LDLIBS)
	/bin/mv d_coo_matgen $(EXEDIR)	
d_matgen: d_matgen.o psb_d_cxx_mat_mod.o psb_d_cxx_impl.o
	$(F90LINK) $(LINKOPT) d_matgen.o psb_d_cxx_mat_mod.o psb_d_cxx_impl.o -o d_matgen  $(PSBLAS_LIB) $(LDLIBS)
	/bin/mv d_matgen $(EXEDIR)	
d_matgen.o: psb_d_cxx_mat_mod.o 

#ppde spde

ppde: ppde.o  
	$(F90LINK) ppde.o  -o ppde $(PSBLAS_LIB) $(LDLIBS)
	/bin/mv ppde $(EXEDIR)

spde: spde.o  
	$(F90LINK) spde.o  -o spde $(PSBLAS_LIB) $(LDLIBS)
	/bin/mv spde $(EXEDIR)

.f90.o:
	$(MPF90) $(F90COPT) $(FINCLUDES) $(FDEFINES) -c $<


clean: 
	/bin/rm -f d_coo_matgen.o d_matgen.o tpg.o ppde.o spde.o \
	 psb_d_cxx_mat_mod.o psb_d_cxx_impl.o *$(.mod) $(EXEDIR)/ppde 
verycleanlib: 
	(cd ../..; make veryclean)
lib:
	(cd ../../; make library)



