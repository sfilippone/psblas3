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


EXEDIR=./runs

all: exed psb_d_pde3d 

exed:
	(if test ! -d $(EXEDIR) ; then mkdir $(EXEDIR); fi)

psb_d_pde3d: psb_d_pde3d.o 
	$(FLINK)  psb_d_pde3d.o -o psb_d_pde3d $(PSBLAS_LIB) $(LDLIBS)
	/bin/mv psb_d_pde3d $(EXEDIR)



clean: 
	/bin/rm -f psb_d_pde3d.o  *$(.mod) $(EXEDIR)/psb_d_pde3d 
verycleanlib: 
	(cd ../..; make veryclean)
lib:
	(cd ../../; make library)



