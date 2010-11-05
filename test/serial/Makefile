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
psb_d_cxx_impl.o d_matgen.o: psb_d_cxx_mat_mod.o 
psb_d_cyy_impl.o d_matgen.o: psb_d_cyy_mat_mod.o 
psb_d_rsb_impl.o d_matgen.o: psb_d_rsb_mat_mod.o  rsb_mod.o 
d_matgen: d_matgen.o psb_d_cxx_mat_mod.o psb_d_cxx_impl.o psb_d_cyy_mat_mod.o psb_d_cyy_impl.o  psb_d_czz_mat_mod.o psb_d_rsb_mat_mod.o
	$(F90LINK) $(LINKOPT) d_matgen.o psb_d_cxx_mat_mod.o psb_d_cxx_impl.o psb_d_cyy_mat_mod.o psb_d_cyy_impl.o  psb_d_czz_mat_mod.o  psb_d_rsb_mat_mod.o -o d_matgen  $(PSBLAS_LIB) $(LDLIBS)
	/bin/mv d_matgen $(EXEDIR)	

.f90.o:
	$(MPF90) $(F90COPT) $(FINCLUDES) $(FDEFINES) -c $<


clean: 
	/bin/rm -f d_coo_matgen.o d_matgen.o  \
	 psb_d_czz_mat_mod.o   \
	 psb_d_rsb_mat_mod.o   \
	 psb_d_cyy_mat_mod.o psb_d_cyy_impl.o  \
	 psb_d_cxx_mat_mod.o psb_d_cxx_impl.o *$(.mod) 
verycleanlib: 
	(cd ../..; make veryclean)
lib:
	(cd ../../; make library)



