BASEDIR=../..
INCDIR=$(BASEDIR)/include
MODDIR=$(BASEDIR)/modules/
include $(INCDIR)/Make.inc.psblas
#
# Libraries used
LIBDIR=$(BASEDIR)/lib
PSBLAS_LIB= -L$(LIBDIR) -lpsb_util -lpsb_krylov -lpsb_prec -lpsb_base 
LDLIBS=$(PSBLDLIBS)
#
# Compilers and such
#
CCOPT= -g
FINCLUDES=$(FMFLAG)$(MODDIR) $(FMFLAG).


EXEDIR=./runs

all: runsd hello pingpong tsum tsum1

runsd:
	(if test ! -d runs ; then mkdir runs; fi)

hello: hello.o 
	$(FLINK)  hello.o -o hello $(PSBLAS_LIB) $(LDLIBS)
	/bin/mv hello $(EXEDIR)

pingpong: pingpong.o 
	$(FLINK)  pingpong.o -o pingpong $(PSBLAS_LIB) $(LDLIBS)
	/bin/mv pingpong $(EXEDIR)
tsum: tsum.o 
	$(FLINK)  tsum.o -o tsum $(PSBLAS_LIB) $(LDLIBS)
	/bin/mv tsum $(EXEDIR)
tsum1: tsum1.o 
	$(FLINK)  tsum1.o -o tsum1 $(PSBLAS_LIB) $(LDLIBS)
	/bin/mv tsum1 $(EXEDIR)



clean: 
	/bin/rm -f hello.o pingpong.o tsum.o tsum1.o
	$(EXEDIR)/hello 
verycleanlib: 
	(cd ../..; make veryclean)
lib:
	(cd ../../; make library)



