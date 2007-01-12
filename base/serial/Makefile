include ../../Make.inc


FOBJS = psb_cest.o psb_dcoins.o psb_dcsdp.o psb_dcsmm.o psb_dcsmv.o \
	psb_dcsnmi.o psb_dcsprt.o psb_dcsrws.o psb_dcssm.o psb_dcssv.o \
	psb_dfixcoo.o psb_dipcoo2csr.o psb_dipcsr2coo.o psb_dneigh.o \
	psb_dnumbmm.o psb_drwextd.o psb_dspgtdiag.o psb_dspgtblk.o \
	psb_dspscal.o psb_dsymbmm.o psb_dtransp.o \
	psb_dipcoo2csc.o psb_dspgetrow.o lsame.o psb_zspgetrow.o\
	psb_zcsmm.o psb_zcsmv.o  psb_zspgtdiag.o psb_zspgtblk.o\
	psb_zcsnmi.o psb_zcsrws.o psb_zcssm.o psb_zcssv.o psb_zcsdp.o\
	psb_zfixcoo.o psb_zipcoo2csr.o psb_zipcsr2coo.o psb_zipcoo2csc.o \
	psb_zcoins.o psb_zcsprt.o psb_zneigh.o psb_ztransp.o psb_ztransc.o\
	psb_zrwextd.o psb_zsymbmm.o psb_znumbmm.o psb_zspscal.o\
	psb_getifield.o psb_setifield.o psb_update_mod.o

LIBDIR = ..
MODDIR = ../modules
INCDIRS = -I $(LIBDIR) -I $(MODDIR) -I . 

lib: auxd cood csrd jadd f77d dpd lib1
	$(AR) $(LIBDIR)/$(LIBNAME) $(FOBJS)
	$(RANLIB) $(LIBDIR)/$(LIBNAME)


lib1:	$(FOBJS)

psb_dcoins.o psb_zcoins.o: psb_update_mod.o

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
	/bin/rm -f $(FOBJS) *$(.mod)
	(cd aux; make clean)
	(cd coo; make clean)
	(cd csr; make clean)
	(cd jad; make clean)
	(cd dp; make  clean)
	(cd f77; make clean)

veryclean: clean 
