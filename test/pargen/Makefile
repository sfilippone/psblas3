include ../../Make.inc
#
# Libraries used
#
LIBDIR=../../lib/
PSBLAS_LIB= -L$(LIBDIR) -lpsb_util -lpsb_krylov -lpsb_prec -lpsb_base 
LDLIBS=$(PSBLDLIBS)
#
# Compilers and such
#
CCOPT= -g
FINCLUDES=$(FMFLAG)$(LIBDIR) $(FMFLAG).


EXEDIR=./runs

all: ppde3d spde3d ppde2d spde2d


ppde3d: ppde3d.o 
	$(F90LINK)  ppde3d.o -o ppde3d $(PSBLAS_LIB) $(LDLIBS)
	/bin/mv ppde3d $(EXEDIR)

ppde3d_ext: ppde3d_ext.o 
	$(F90LINK)  ppde3d_ext.o -o ppde3d_ext $(PSBLAS_LIB) $(LDLIBS)
	/bin/mv ppde3d_ext $(EXEDIR)


spde3d: spde3d.o 
	$(F90LINK) spde3d.o  -o spde3d $(PSBLAS_LIB) $(LDLIBS)
	/bin/mv spde3d $(EXEDIR)

ppde2d: ppde2d.o 
	$(F90LINK)  ppde2d.o -o ppde2d $(PSBLAS_LIB) $(LDLIBS)
	/bin/mv ppde2d $(EXEDIR)


spde2d: spde2d.o 
	$(F90LINK) spde2d.o  -o spde2d $(PSBLAS_LIB) $(LDLIBS)
	/bin/mv spde2d $(EXEDIR)


clean: 
	/bin/rm -f ppde3d.o spde3d.o  ppde2d.o spde2d.o ppde3d_ext.o  \
	$(EXEDIR)/ppde3d $(EXEDIR)/spde3d $(EXEDIR)/ppde2d $(EXEDIR)/spde2d 
verycleanlib: 
	(cd ../..; make veryclean)
lib:
	(cd ../../; make library)



