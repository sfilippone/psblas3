include ../../Make.inc

FOBJS = psb_cdall.o psb_cdals.o psb_cdalv.o psb_cd_inloc.o psb_cdins.o psb_cdprt.o \
	psb_cdren.o psb_cdrep.o psb_get_overlap.o psb_cd_lstext.o psb_cd_remap.o\
	psb_cdcpy.o psb_cd_reinit.o psb_cd_switch_ovl_indxmap.o psb_cd_renum_block.o \
	psb_dspalloc.o psb_dspasb.o psb_d_remap.o \
	psb_dspfree.o psb_dspins.o psb_dsprn.o \
	psb_sspalloc.o psb_sspasb.o \
	psb_sspfree.o psb_sspins.o psb_ssprn.o\
	psb_glob_to_loc.o psb_loc_to_glob.o\
	psb_iallc.o psb_iasb.o psb_ifree.o psb_iins.o \
	psb_lallc.o psb_lasb.o psb_lfree.o psb_lins.o \
	psb_sallc.o psb_sasb.o psb_sfree.o psb_sins.o \
	psb_dallc.o psb_dasb.o psb_dfree.o psb_dins.o \
	psb_callc.o psb_casb.o psb_cfree.o psb_cins.o \
	psb_zallc.o psb_zasb.o psb_zfree.o psb_zins.o \
	psb_mallc_a.o psb_masb_a.o psb_mfree_a.o psb_mins_a.o \
	psb_eallc_a.o psb_easb_a.o psb_efree_a.o psb_eins_a.o \
	psb_sallc_a.o psb_sasb_a.o psb_sfree_a.o psb_sins_a.o \
	psb_dallc_a.o psb_dasb_a.o psb_dfree_a.o psb_dins_a.o \
	psb_callc_a.o psb_casb_a.o psb_cfree_a.o psb_cins_a.o \
	psb_zallc_a.o psb_zasb_a.o psb_zfree_a.o psb_zins_a.o \
	psb_zspalloc.o psb_zspasb.o psb_zspfree.o\
	psb_zspins.o psb_zsprn.o \
	psb_cspalloc.o psb_cspasb.o psb_cspfree.o\
	psb_cspins.o psb_csprn.o psb_cd_set_bld.o \
	psb_s_map.o psb_d_map.o psb_c_map.o psb_z_map.o \
	psb_s_par_csr_spspmm.o psb_d_par_csr_spspmm.o psb_c_par_csr_spspmm.o psb_z_par_csr_spspmm.o \
	psb_s_glob_transpose.o psb_d_glob_transpose.o psb_c_glob_transpose.o psb_z_glob_transpose.o \
	psb_cgetelem.o psb_dgetelem.o psb_sgetelem.o psb_zgetelem.o   
#	psb_lallc.o psb_lasb.o psb_lfree.o psb_lins.o \

MPFOBJS = psb_icdasb.o psb_ssphalo.o  psb_dsphalo.o psb_csphalo.o psb_zsphalo.o  \
	psb_dcdbldext.o psb_zcdbldext.o psb_scdbldext.o psb_ccdbldext.o

LIBDIR=..
INCDIR=..
MODDIR=../modules
FINCLUDES=$(FMFLAG). $(FMFLAG)$(MODDIR) $(FMFLAG)$(INCDIR) $(FIFLAG)$(MODDIR)

lib: mpfobjs $(FOBJS)
	$(AR) $(LIBDIR)/$(LIBNAME) $(MPFOBJS) $(FOBJS)
	$(RANLIB) $(LIBDIR)/$(LIBNAME)


mpfobjs:
	(make $(MPFOBJS) FC="$(MPFC)")

clean:
	/bin/rm -f $(MPFOBJS) $(FOBJS)

veryclean: clean
