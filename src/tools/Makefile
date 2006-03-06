include ../../Make.inc

FOBJS = psb_dallc.o psb_dasb.o psb_dcsrp.o psb_cdprt.o \
	psb_dfree.o psb_dgelp.o psb_dins.o \
	psb_cdall.o psb_cdalv.o psb_dscasb.o psb_dsccpy.o \
	psb_cddec.o psb_cdfree.o psb_dscins.o psb_cdovr.o \
	psb_dscren.o psb_cdrep.o psb_dspalloc.o psb_dspasb.o \
	psb_dspcnv.o psb_dspfree.o psb_dspins.o psb_dsprn.o \
	psb_glob_to_loc.o psb_ialloc.o psb_iasb.o \
	psb_ifree.o psb_iins.o psb_loc_to_glob.o

MPFOBJS = psb_cdovrbld.o psb_dspovr.o

INCDIRS = -I ../../lib -I .
LIBDIR = ../../lib

lib: mpfobjs $(FOBJS)
	$(AR) $(LIBDIR)/$(LIBNAME) $(MPFOBJS) $(FOBJS)
	$(RANLIB) $(LIBDIR)/$(LIBNAME)


mpfobjs: 
	(make $(MPFOBJS) F90="$(MPF90)" FC="$(MPF90)" FCOPT="$(F90COPT)")


clean:
	/bin/rm -f $(MPFOBJS) $(FOBJS)

veryclean: clean
