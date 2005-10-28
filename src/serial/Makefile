include ../../Make.inc


FOBJS = psb_cest.o psb_dcoins.o psb_dcsdp.o psb_dcsmm.o psb_dcsmv.o \
	psb_dcsnmi.o psb_dcsprt.o psb_dcsrws.o psb_dcssm.o psb_dcssv.o \
	psb_dfixcoo.o psb_dipcoo2csr.o psb_dipcsr2coo.o psb_dneigh.o \
	psb_dnumbmm.o psb_drwextd.o psb_dspgtdiag.o psb_dspgtrow.o \
	psb_dspinfo.o psb_dspscal.o psb_dsymbmm.o psb_dtransp.o \
	string_impl.o psb_dipcoo2csc.o 


INCDIRS = -I ../../lib -I .
LIBDIR = ../../lib

lib: auxd cood csrd jadd f77d dpd lib1
	$(AR) $(LIBDIR)/$(LIBNAME) $(FOBJS)
	$(RANLIB) $(LIBDIR)/$(LIBNAME)


lib1:	$(FOBJS)


auxd:
	(cd aux; make lib)

cood:
	(cd coo; make lib)

csrd:
	(cd csr; make lib)

jadd:
	(cd jad; make lib)

dpd:
	(cd dp; make lib)

f77d:
	(cd f77; make lib)

clean:
	/bin/rm -f $(FOBJS)
	(cd aux; make clean)
	(cd coo; make clean)
	(cd csr; make clean)
	(cd jad; make clean)
	(cd dp; make  clean)
	(cd f77; make clean)

veryclean:
	/bin/rm -f $(FOBJS)
	(cd aux; make veryclean)
	(cd coo; make veryclean)
	(cd csr; make veryclean)
	(cd jad; make veryclean)
	(cd dp; make  veryclean)
	(cd f77; make veryclean)
