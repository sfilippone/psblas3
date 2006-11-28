include ../../Make.inc

FOBJS = psb_dallc.o psb_dasb.o psb_dcsrp.o psb_cdprt.o \
	psb_dfree.o psb_dgelp.o psb_dins.o \
	psb_cdall.o psb_cdalv.o psb_cd_inloc.o psb_cdcpy.o \
	psb_cddec.o psb_cdfree.o psb_cdins.o  \
	psb_cdren.o psb_cdrep.o psb_cdtransfer.o psb_get_overlap.o\
	psb_dspalloc.o psb_dspasb.o \
	psb_dspcnv.o psb_dspfree.o psb_dspins.o psb_dsprn.o \
	psb_glob_to_loc.o psb_ialloc.o psb_iasb.o \
	psb_ifree.o psb_iins.o psb_loc_to_glob.o\
	psb_zallc.o psb_zasb.o psb_zfree.o psb_zins.o \
	psb_zspalloc.o psb_zspasb.o psb_zspcnv.o psb_zspfree.o\
	psb_zspins.o psb_zsprn.o psb_zcdovr.o psb_zgelp.o

MPFOBJS = psb_dcdovrbld.o psb_dsphalo.o psb_zcdovrbld.o psb_zsphalo.o psb_cdasb.o \
 psb_dcdovr.o psb_zcdovr.o

INCDIRS = -I ../../lib -I .
LIBDIR = ../../lib

lib: mpfobjs $(FOBJS)
	$(AR) $(LIBDIR)/$(LIBNAME) $(MPFOBJS) $(FOBJS)
	$(RANLIB) $(LIBDIR)/$(LIBNAME)


mpfobjs: 
	(make $(MPFOBJS) F90="$(MPF90)" FC="$(MPF90)" FCOPT="$(F90COPT)")

psb_glob_to_loc.o: psb_glob_to_loc.f90
	$(F90) $(F90COPT) $(INCDIRS) -c $< $(F90INLINE)
clean:
	/bin/rm -f $(MPFOBJS) $(FOBJS)

veryclean: clean
