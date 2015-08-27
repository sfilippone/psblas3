BASEDIR=../..
INCDIR=$(BASEDIR)/include
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
FINCLUDES=$(FMFLAG)$(INCDIR) $(FMFLAG).


EXEDIR=./runs

all: hello pingpong

hello: hello.o 
	$(F90LINK)  hello.o -o hello $(PSBLAS_LIB) $(LDLIBS)
	/bin/mv hello $(EXEDIR)

pingpong: pingpong.o 
	$(F90LINK)  pingpong.o -o pingpong $(PSBLAS_LIB) $(LDLIBS)
	/bin/mv pingpong $(EXEDIR)



clean: 
	/bin/rm -f hello.o pingpong.o
	$(EXEDIR)/hello 
verycleanlib: 
	(cd ../..; make veryclean)
lib:
	(cd ../../; make library)



