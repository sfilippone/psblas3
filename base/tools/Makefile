include ../../Make.inc

FOBJS = psb_sallc.o psb_sasb.o  \
	psb_sfree.o psb_sins.o \
	psb_dallc.o psb_dasb.o \
	psb_dfree.o psb_dins.o \
	psb_cdals.o psb_cdalv.o psb_cd_inloc.o psb_cdins.o psb_cdprt.o \
	psb_cdren.o psb_cdrep.o psb_get_overlap.o psb_cd_lstext.o\
	psb_dspalloc.o psb_dspasb.o \
	psb_dspfree.o psb_dspins.o psb_dsprn.o \
	psb_sspalloc.o psb_sspasb.o \
	psb_sspfree.o psb_sspins.o psb_ssprn.o\
	psb_glob_to_loc.o psb_ialloc.o psb_iasb.o \
	psb_ifree.o psb_iins.o psb_loc_to_glob.o\
	psb_zallc.o psb_zasb.o psb_zfree.o psb_zins.o \
	psb_zspalloc.o psb_zspasb.o psb_zspfree.o\
	psb_zspins.o psb_zsprn.o \
	psb_cspalloc.o psb_cspasb.o psb_cspfree.o\
	psb_callc.o psb_casb.o psb_cfree.o psb_cins.o \
	psb_cspins.o psb_csprn.o psb_map.o psb_inter_desc.o
#	\
#	psb_s_inter_desc.o psb_s_inter_desc_noidx.o \
#	psb_d_inter_desc.o psb_d_inter_desc_noidx.o \
#	psb_c_inter_desc.o psb_c_inter_desc_noidx.o \
#	psb_z_inter_desc.o psb_z_inter_desc_noidx.o 

MPFOBJS = psb_ssphalo.o psb_csphalo.o psb_dsphalo.o psb_zsphalo.o psb_icdasb.o \
	psb_dcdbldext.o psb_zcdbldext.o psb_scdbldext.o psb_ccdbldext.o

LIBDIR=..
MODDIR=../modules
FINCLUDES=$(FMFLAG)$(LIBDIR) $(FMFLAG)$(MODDIR) $(FMFLAG). $(FIFLAG)$(MODDIR)

lib: mpfobjs $(FOBJS)
	$(AR) $(LIBDIR)/$(LIBNAME) $(MPFOBJS) $(FOBJS)
	$(RANLIB) $(LIBDIR)/$(LIBNAME)


mpfobjs: 
	(make $(MPFOBJS) F90="$(MPF90)" FC="$(MPF90)" FCOPT="$(F90COPT)")

clean:
	/bin/rm -f $(MPFOBJS) $(FOBJS)

veryclean: clean
