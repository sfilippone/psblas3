include ../../Make.inc
#
# Libraries used
#
LIBDIR=../../lib/
OPTDIR=../../opt
PSBLAS_LIB= -L$(LIBDIR) -lpsb_util -lpsb_krylov -lpsb_prec -lpsb_base -L$(OPTDIR) -lopt\
	-lspgpu -L/usr/local/spgpu/lib -lcudart -lcublas -L/usr/local/cuda/lib64
LDLIBS=$(PSBLDLIBS)
#
# Compilers and such
#
CCOPT= -g
FINCLUDES=$(FMFLAG)$(LIBDIR) $(FMFLAG). $(FMFLAG)$(OPTDIR)


EXEDIR=./runs

all: ppde 

ppde: ppde.o 
	$(F90LINK)  ppde.o -o ppde $(PSBLAS_LIB) $(LDLIBS) -lgomp
	/bin/mv ppde $(EXEDIR)

clean: 
	/bin/rm -f ppde.o spde.o  $(EXEDIR)/ppde 
verycleanlib: 
	(cd ../..; make veryclean)
lib:
	(cd ../../; make library)



