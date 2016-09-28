BASEDIR=../..
INCDIR=$(BASEDIR)/include
include $(INCDIR)/Make.inc.psblas
#
# Libraries used
PFUNIT = /opt/pfunit/pfunit-coarrays-last
FFLAGS += -I$(INCDIR) -I$(PFUNIT)/mod -ISuites
LIBDIR=$(BASEDIR)/lib
PSBLAS_LIB= -L/opencoarrays6.2 -L$(LIBDIR) -lcaf_mpi -lpsb_util -lpsb_base -llapack -lblas
LDLIBS=$(PSBLDLIBS)
#
# Compilers and such
#
CCOPT= -g
FINCLUDES=$(FMFLAG)$(INCDIR) $(FMFLAG).


all: test_psb_halo

%: %.x
	mpirun -np 4 ./$^	
%.x:%.o driver.o
	$(FCFLAGS)$(F90LINK) -g $(LOPT) -o $@ $^ -L$(PFUNIT)/lib $(PSBLAS_LIB) $(LDLIBS) -DUSE_CAF -lpfunit -lmpi
#Create .F90 file
%.F90: %.pf
	$(PFUNIT)/bin/pFUnitParser.py $< $@ -lmpi
#Create .o file
%.o: %.F90
	$(FC) -g -DUSE_PFUNIT -DUSE_CAF -c $(FFLAGS) $(FPPFLAGS) $^ $(PFUNIT)/include/driver.F90 -L$(PFUNIT)/lib $(PSBLAS_LIB) $(LDLIBS) -lpfunit -lmpi



clean: 
	/bin/rm -f *.F90 *.o *.mod
verycleanlib: 
	(cd ../..; make veryclean)
lib:
	(cd ../../; make library)

