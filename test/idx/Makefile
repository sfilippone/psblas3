INSTALLDIR=../..
INCDIR=$(INSTALLDIR)/include
MODDIR=$(INSTALLDIR)/modules/
include $(INCDIR)/Make.inc.psblas
#
# Libraries used
LIBDIR=$(INSTALLDIR)/lib
PSBLAS_LIB= -L$(LIBDIR) -lpsb_util -lpsb_krylov -lpsb_prec -lpsb_base 
LDLIBS=$(PSBLDLIBS)
#
# Compilers and such
#
CCOPT= -g
FINCLUDES=$(FMFLAG)$(MODDIR) $(FMFLAG).


EXEDIR=./

all: tryidxijk 

tryidxijk: tryidxijk.o 
	$(FLINK)  tryidxijk.o -o tryidxijk $(PSBLAS_LIB) $(LDLIBS)
	/bin/mv tryidxijk $(EXEDIR)

clean: 
	/bin/rm -f tryidxijk.o  *$(.mod) \
	$(EXEDIR)/tryidxijk $(EXEDIR)/psb_d_pde3d
verycleanlib: 
	(cd ../..; make veryclean)
lib:
	(cd ../../; make library)



