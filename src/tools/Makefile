include ../../Make.inc

FOBJS = psb_dallc.o psb_dasb.o psb_dcsrp.o psb_descprt.o \
	psb_dfree.o psb_dgelp.o psb_dins.o \
	psb_dscall.o psb_dscalv.o psb_dscasb.o psb_dsccpy.o \
	psb_dscdec.o psb_dscfree.o psb_dscins.o psb_dscov.o \
	psb_dscren.o psb_dscrep.o psb_dspalloc.o psb_dspasb.o \
	psb_dspcnv.o psb_dspfree.o psb_dspins.o psb_dsprn.o \
	psb_dspupdate.o psb_glob_to_loc.o psb_ialloc.o psb_iasb.o \
	psb_ifree.o psb_iins.o psb_loc_to_glob.o

MPFOBJS = psb_descasb.o psb_dcsrovr.o

INCDIRS = ../../lib .

lib: mpfobjs $(FOBJS)


mpfobjs: 
	(make $(MPFOBJS) F90="$(MPF90)" FC="$(MPF90)" FCOPT="$(F90COPT)")


clean:
	/bin/rm -f $(MPFOBJS) $(FOBJS)
