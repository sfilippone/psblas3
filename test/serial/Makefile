#
# Libraries used
#
BASEDIR=../..
INCDIR=$(BASEDIR)/include/
include $(INCDIR)/Make.inc.psblas
LIBDIR=$(BASEDIR)/lib/
PSBLAS_LIB= -L$(LIBDIR) -lpsb_util -lpsb_krylov -lpsb_prec -lpsb_base
LDLIBS=$(PSBLDLIBS)

FINCLUDES=$(FMFLAG)$(INCDIR) $(FMFLAG).


EXEDIR=./runs


all: d_matgen 

psb_d_cxx_impl.o d_matgen.o: psb_d_cxx_mat_mod.o 

d_matgen: d_matgen.o psb_d_cxx_mat_mod.o psb_d_cxx_impl.o   
	$(F90LINK) $(LINKOPT) d_matgen.o psb_d_cxx_mat_mod.o psb_d_cxx_impl.o \
	 -o d_matgen  $(PSBLAS_LIB) $(LDLIBS) 
	/bin/cp -p $(CPUPDFLAG) d_matgen $(EXEDIR)	
#	/bin/mv d_matgen $(EXEDIR)	

check: all
	cd runs && echo 5 | ./d_matgen

clean: 
	/bin/rm -f  d_matgen.o  \
	 psb_d_cxx_mat_mod.o psb_d_cxx_impl.o *$(.mod)  
verycleanlib: 
	(cd ../..; make veryclean)
lib:
	(cd ../../; make library)



