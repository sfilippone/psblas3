include ../Make.inc
#
# Libraries used
#
LIBDIR=../lib/
PSBLAS_LIB= -L$(LIBDIR) -lpsb_base
#-lpsb_util -lpsb_krylov -lpsb_prec -lpsb_base 
LDLIBS=$(PSBLDLIBS)
#
# Compilers and such
#
CCOPT= -g
FINCLUDES=$(FMFLAG)$(LIBDIR) $(FMFLAG).


EXEDIR=./runs

all: psb_d_ell_impl.o psb_d_ell_mat_mod.o
psb_d_ell_impl.o: psb_d_ell_mat_mod.o

clean: 
	/bin/rm -f  psb_d_ell_impl.o psb_d_ell_mat_mod.o *$(.mod)

verycleanlib: 
	(cd ../..; make veryclean)
lib:
	(cd ../../; make library)



