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

# FIXME: martone will clean up this file from RSB_EXTRA
RSB_EXTRA=-lgomp

all: d_coo_matgen d_matgen 
d_coo_matgen: d_coo_matgen.o rsb_mod.o
	$(F90LINK) $(LINKOPT) d_coo_matgen.o  -o d_coo_matgen  $(PSBLAS_LIB) $(LDLIBS)
	/bin/mv d_coo_matgen $(EXEDIR)	
psb_d_cxx_impl.o d_matgen.o: psb_d_cxx_mat_mod.o 
psb_d_cyy_impl.o d_matgen.o: psb_d_cyy_mat_mod.o 
d_matgen.o: psb_d_czz_mat_mod.o 
psb_d_rsb_mat_mod.o:  rsb_mod.o 
psb_d_rsb_impl.o d_matgen.o: psb_d_rsb_mat_mod.o  rsb_mod.o 
d_matgen: d_matgen.o psb_d_cxx_mat_mod.o psb_d_cxx_impl.o psb_d_cyy_mat_mod.o psb_d_cyy_impl.o  psb_d_czz_mat_mod.o psb_d_rsb_mat_mod.o
	$(F90LINK) $(LINKOPT) d_matgen.o psb_d_cxx_mat_mod.o psb_d_cxx_impl.o psb_d_cyy_mat_mod.o psb_d_cyy_impl.o  psb_d_czz_mat_mod.o  psb_d_rsb_mat_mod.o -o d_matgen  $(PSBLAS_LIB) $(LDLIBS) $(RSB_EXTRA)
	/bin/mv d_matgen $(EXEDIR)	

check: all
	cd runs && echo 5 | ./d_matgen

clean: 
	/bin/rm -f d_coo_matgen.o d_matgen.o  \
	 psb_d_czz_mat_mod.o   \
	 psb_d_rsb_mat_mod.o   \
	 psb_d_cyy_mat_mod.o psb_d_cyy_impl.o  \
	 psb_d_cxx_mat_mod.o psb_d_cxx_impl.o *$(.mod)  \
	 rsb_mod.o
verycleanlib: 
	(cd ../..; make veryclean)
lib:
	(cd ../../; make library)



