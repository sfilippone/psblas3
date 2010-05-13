include ../../Make.inc
#
# Libraries used
#
LIBDIR=../../lib/
PSBLAS_LIB= -L$(LIBDIR) -lpsb_base 
LDLIBS=$(PSBLDLIBS)
#
# Compilers and such
#
CCOPT= -g
FINCLUDES=$(FMFLAG)$(LIBDIR) $(FMFLAG).

PSBTOBJS=psbtf.o psb_mvsv_tester.o  \
	psb_s_mvsv_tester.o psb_d_mvsv_tester.o psb_c_mvsv_tester.o \
	psb_z_mvsv_tester.o 
EXEDIR=./runs


all: psbtf
psbtf.o: psb_mvsv_tester.o
psb_mvsv_tester.o: psb_s_mvsv_tester.o psb_d_mvsv_tester.o psb_c_mvsv_tester.o \
	psb_z_mvsv_tester.o

psbtf: $(PSBTOBJS)
	$(F90LINK)  $(PSBTOBJS)  -o psbtf $(PSBLAS_LIB) $(LDLIBS)
	/bin/mv psbtf $(EXEDIR)

psbtf.o: psb_mvsv_tester.o 


psbt1:  psbt1.o
	$(F90LINK) psbt1.o  -o psbt1 $(PSBLAS_LIB) $(LDLIBS)
	/bin/mv psbt1 $(EXEDIR)


.f90.o:
	$(MPF90) $(F90COPT) $(FINCLUDES) $(FDEFINES) -c $<


clean: 
	/bin/rm -f $(PSBTOBJS) ppde.o spde.o  $(EXEDIR)/ppde 
verycleanlib: 
	(cd ../..; make veryclean)
lib:
	(cd ../../; make library)



