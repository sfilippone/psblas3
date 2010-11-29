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

OBJS=psb_d_ell_impl.o psb_d_ell_mat_mod.o rsb_mod.o psb_d_rsb_mat_mod.o
all: libopt.a
libopt.a: $(OBJS)
	ar cur libopt.a $(OBJS)
psb_d_ell_impl.o: psb_d_ell_mat_mod.o
psb_d_rsb_mat_mod.o:  rsb_mod.o

clean: 
	/bin/rm -f  $(OBJS) *$(.mod) 

verycleanlib: 
	(cd ../..; make veryclean)
lib:
	(cd ../../; make library)



