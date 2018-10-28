INSTALLDIR=../..
INCDIR=$(INSTALLDIR)/include/
MODDIR=$(INSTALLDIR)/modules/
include $(INCDIR)/Make.inc.psblas
LIBDIR=$(INSTALLDIR)/lib/
PSBLAS_LIB= -L$(LIBDIR) -lpsb_util -lpsb_krylov -lpsb_prec -lpsb_base
LDLIBS=$(PSBLDLIBS)
CCOPT= -g
FINCLUDES=$(FMFLAG)$(MODDIR) $(FMFLAG).

PSBTOBJS=psbtf.o psb_mvsv_tester.o  \
	psb_s_mvsv_tester.o psb_d_mvsv_tester.o psb_c_mvsv_tester.o \
	psb_z_mvsv_tester.o 
EXEDIR=./runs


all: runsd psbtf

runsd:
	(if test ! -d runs ; then mkdir runs; fi)

psbtf.o: psb_mvsv_tester.o
psb_mvsv_tester.o: psb_s_mvsv_tester.o psb_d_mvsv_tester.o psb_c_mvsv_tester.o \
	psb_z_mvsv_tester.o

psbtf: $(PSBTOBJS)
	$(FLINK)  $(PSBTOBJS)  -o psbtf $(PSBLAS_LIB) $(LDLIBS)
	/bin/mv psbtf $(EXEDIR)

psbtf.o: psb_mvsv_tester.o 


.f90.o:
	$(MPFC) $(FCOPT) $(FINCLUDES) $(FDEFINES) -c $<


clean: 
	/bin/rm -f $(PSBTOBJS) ppde.o spde.o  $(EXEDIR)/ppde 
verycleanlib: 
	(cd ../..; make veryclean)
lib:
	(cd ../../; make library)



