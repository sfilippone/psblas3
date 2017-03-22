BASEDIR=../..
INCDIR=$(BASEDIR)/include
include $(INCDIR)/Make.inc.psblas
#
# Libraries used
PFUNIT = /opt/pfunit/pfunit-marzo2017
FFLAGS += -I$(INCDIR) -I$(PFUNIT)/mod -ISuites
LIBDIR=$(BASEDIR)/lib
PSBLAS_LIB= -L/opencoarrays6.2 -L$(LIBDIR) -lcaf_mpi -lpsb_util -lpsb_base -llapack -lblas
LDLIBS=$(PSBLDLIBS)
#
# Compilers and such
#
CCOPT= -g
FINCLUDES=$(FMFLAG)$(INCDIR) $(FMFLAG).


all: prova

prova.x: test_psb_dmatdist.o test_psb_ihalo.o test_psb_dhalo.o test_psb_shalo.o test_psb_chalo.o test_psb_zhalo.o test_psb_reduce_nrm2.o test_psb_sum.o test_psb_max.o test_psb_amx.o test_psb_min.o test_psb_amn.o test_psb_broadcast.o test_psb_caf.o driver.o
	$(FCFLAGS)$(F90LINK) -g $(LOPT) -o $@ $^ -L$(PFUNIT)/lib $(PSBLAS_LIB) $(LDLIBS) -DUSE_CAF -lpfunit -lmpi
%: %.x
	mpirun -np 8 ./$^	
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

