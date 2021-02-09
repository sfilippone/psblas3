include ../../Make.inc

LIBDIR=../../lib
INCDIR=../../include
MODDIR=../../modules
HERE=..
OBJS=psb_s_prec_type_impl.o psb_d_prec_type_impl.o \
     psb_c_prec_type_impl.o psb_z_prec_type_impl.o \
     psb_d_diagprec_impl.o psb_d_bjacprec_impl.o psb_d_nullprec_impl.o \
     psb_dilu_fct.o psb_d_ilu0_fact.o psb_d_iluk_fact.o psb_d_ilut_fact.o \
     psb_dprecbld.o psb_dprecinit.o \
     psb_s_diagprec_impl.o psb_s_bjacprec_impl.o psb_s_nullprec_impl.o \
     psb_silu_fct.o psb_s_ilu0_fact.o psb_s_iluk_fact.o psb_s_ilut_fact.o \
     psb_sprecbld.o  psb_sprecinit.o \
     psb_c_diagprec_impl.o psb_c_bjacprec_impl.o psb_c_nullprec_impl.o \
     psb_cilu_fct.o psb_c_ilu0_fact.o psb_c_iluk_fact.o psb_c_ilut_fact.o \
     psb_cprecbld.o  psb_cprecinit.o \
     psb_z_diagprec_impl.o psb_z_bjacprec_impl.o psb_z_nullprec_impl.o \
     psb_zilu_fct.o psb_z_ilu0_fact.o psb_z_iluk_fact.o psb_z_ilut_fact.o \
     psb_zprecbld.o  psb_zprecinit.o \
     psb_c_sparsify.o psb_d_sparsify.o psb_s_sparsify.o psb_z_sparsify.o \
     psb_crwclip.o psb_drwclip.o psb_srwclip.o psb_zrwclip.o \
     psb_c_sp_drop.o psb_d_sp_drop.o psb_s_sp_drop.o psb_z_sp_drop.o \
     psb_dsparse_biconjg_llk_noth.o psb_dsparse_biconjg_llk.o \
     psb_dsparse_biconjg_mlk.o psb_dsparse_biconjg_s_ft_llk.o \
     psb_dsparse_biconjg_s_llk.o \
     psb_csparse_biconjg_llk_noth.o psb_csparse_biconjg_llk.o \
     psb_csparse_biconjg_mlk.o psb_csparse_biconjg_s_ft_llk.o \
     psb_csparse_biconjg_s_llk.o \
     psb_zsparse_biconjg_llk_noth.o psb_zsparse_biconjg_llk.o \
     psb_zsparse_biconjg_mlk.o psb_zsparse_biconjg_s_ft_llk.o \
     psb_zsparse_biconjg_s_llk.o \
     psb_ssparse_biconjg_llk_noth.o psb_ssparse_biconjg_llk.o \
     psb_ssparse_biconjg_mlk.o psb_ssparse_biconjg_s_ft_llk.o \
     psb_ssparse_biconjg_s_llk.o \
     psb_d_ainv_bld.o psb_c_ainv_bld.o psb_s_ainv_bld.o \
     psb_z_ainv_bld.o \
     psb_c_invt_fact.o psb_d_invt_fact.o psb_s_invt_fact.o psb_z_invt_fact.o\
     psb_c_invk_fact.o psb_d_invk_fact.o psb_s_invk_fact.o psb_z_invk_fact.o

LIBNAME=$(PRECLIBNAME)
COBJS=
FINCLUDES=$(FMFLAG).. $(FMFLAG)$(MODDIR)

lib:  $(OBJS)
	$(AR) $(HERE)/$(LIBNAME) $(OBJS)
	$(RANLIB) $(HERE)/$(LIBNAME)

veryclean: clean

clean:
	/bin/rm -f $(OBJS) $(LOCAL_MODS)
