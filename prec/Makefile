include ../Make.inc

LIBDIR=../lib
INCDIR=../include
MODDIR=../modules
HERE=.
MODOBJS=psb_prec_const_mod.o\
	psb_s_prec_type.o psb_d_prec_type.o psb_c_prec_type.o psb_z_prec_type.o \
	psb_s_base_prec_mod.o psb_d_base_prec_mod.o psb_c_base_prec_mod.o psb_z_base_prec_mod.o \
	psb_prec_type.o  \
	psb_prec_mod.o psb_s_prec_mod.o psb_d_prec_mod.o psb_c_prec_mod.o psb_z_prec_mod.o \
	psb_d_diagprec.o psb_d_nullprec.o psb_d_bjacprec.o psb_s_ilu_fact_mod.o \
	psb_s_diagprec.o psb_s_nullprec.o psb_s_bjacprec.o psb_d_ilu_fact_mod.o \
	psb_c_diagprec.o psb_c_nullprec.o psb_c_bjacprec.o psb_c_ilu_fact_mod.o \
	psb_z_diagprec.o psb_z_nullprec.o psb_z_bjacprec.o psb_z_ilu_fact_mod.o \
	psb_c_ainv_fact_mod.o psb_d_ainv_fact_mod.o psb_s_ainv_fact_mod.o psb_z_ainv_fact_mod.o \
	psb_c_ainv_tools_mod.o psb_d_ainv_tools_mod.o psb_s_ainv_tools_mod.o psb_z_ainv_tools_mod.o \
	psb_ainv_tools_mod.o \
	psb_biconjg_mod.o psb_c_biconjg_mod.o psb_d_biconjg_mod.o psb_s_biconjg_mod.o \
	psb_z_biconjg_mod.o \
	psb_c_invt_fact_mod.o psb_d_invt_fact_mod.o psb_s_invt_fact_mod.o \
	psb_z_invt_fact_mod.o\
	psb_c_invk_fact_mod.o psb_d_invk_fact_mod.o psb_s_invk_fact_mod.o \
	psb_z_invk_fact_mod.o

LIBNAME=$(PRECLIBNAME)
COBJS=
FINCLUDES=$(FMFLAG). $(FMFLAG)$(MODDIR)
OBJS=$(F90OBJS) $(COBJS) $(MPFOBJS) $(MODOBJS)

objs:  $(OBJS) impld
	/bin/cp -p $(CPUPDFLAG) *$(.mod) $(MODDIR)
lib:  objs
	$(MAKE) -C impl lib
	$(AR) $(HERE)/$(LIBNAME) $(OBJS)
	$(RANLIB) $(HERE)/$(LIBNAME)
	/bin/cp -p $(CPUPDFLAG) $(HERE)/$(LIBNAME) $(LIBDIR)

impld: $(OBJS)
	$(MAKE) -C impl objs

$(OBJS): $(MODDIR)/$(BASEMODNAME)$(.mod)


psb_s_base_prec_mod.o psb_d_base_prec_mod.o psb_c_base_prec_mod.o psb_z_base_prec_mod.o \
psb_s_ilu_fact_mod.o psb_d_ilu_fact_mod.o psb_c_ilu_fact_mod.o psb_z_ilu_fact_mod.o : psb_prec_const_mod.o
psb_s_prec_type.o: psb_s_base_prec_mod.o
psb_d_prec_type.o: psb_d_base_prec_mod.o
psb_c_prec_type.o: psb_c_base_prec_mod.o
psb_z_prec_type.o: psb_z_base_prec_mod.o
psb_s_prec_mod.o:  psb_prec_type.o
psb_d_prec_mod.o:  psb_prec_type.o
psb_c_prec_mod.o:  psb_prec_type.o
psb_z_prec_mod.o:  psb_prec_type.o
psb_prec_type.o: psb_s_prec_type.o psb_d_prec_type.o psb_c_prec_type.o psb_z_prec_type.o
psb_prec_mod.o: psb_s_prec_mod.o psb_d_prec_mod.o psb_c_prec_mod.o psb_z_prec_mod.o
psb_s_bjacprec.o psb_s_diagprec.o psb_s_nullprec.o: psb_prec_mod.o psb_s_base_prec_mod.o
psb_d_bjacprec.o psb_d_diagprec.o psb_d_nullprec.o: psb_prec_mod.o psb_d_base_prec_mod.o
psb_c_bjacprec.o psb_c_diagprec.o psb_c_nullprec.o: psb_prec_mod.o psb_c_base_prec_mod.o
psb_z_bjacprec.o psb_z_diagprec.o psb_z_nullprec.o: psb_prec_mod.o psb_z_base_prec_mod.o
psb_s_bjacprec.o: psb_s_ilu_fact_mod.o psb_s_ainv_fact_mod.o psb_s_invk_fact_mod.o psb_s_invt_fact_mod.o
psb_d_bjacprec.o: psb_d_ilu_fact_mod.o psb_d_ainv_fact_mod.o psb_d_invk_fact_mod.o psb_d_invt_fact_mod.o
psb_c_bjacprec.o: psb_c_ilu_fact_mod.o psb_c_ainv_fact_mod.o psb_c_invk_fact_mod.o psb_c_invt_fact_mod.o
psb_z_bjacprec.o: psb_z_ilu_fact_mod.o psb_z_ainv_fact_mod.o psb_z_invk_fact_mod.o psb_z_invt_fact_mod.o
psb_d_ainv_fact_mod.o: psb_prec_const_mod.o psb_ainv_tools_mod.o
psb_s_ainv_fact_mod.o: psb_prec_const_mod.o psb_ainv_tools_mod.o
psb_c_ainv_fact_mod.o: psb_prec_const_mod.o psb_ainv_tools_mod.o
psb_z_ainv_fact_mod.o: psb_prec_const_mod.o psb_ainv_tools_mod.o
psb_ainv_tools_mod.o: psb_c_ainv_tools_mod.o psb_d_ainv_tools_mod.o psb_s_ainv_tools_mod.o psb_z_ainv_tools_mod.o
psb_biconjg_mod.o: psb_c_biconjg_mod.o \
	psb_d_biconjg_mod.o psb_s_biconjg_mod.o psb_z_biconjg_mod.o
psb_c_biconjg_mod.o psb_d_biconjg_mod.o psb_s_biconjg_mod.o psb_z_biconjg_mod.o:  psb_prec_const_mod.o
psb_c_invt_fact_mod.o: psb_prec_const_mod.o psb_c_ilu_fact_mod.o
psb_d_invt_fact_mod.o: psb_prec_const_mod.o psb_d_ilu_fact_mod.o
psb_s_invt_fact_mod.o: psb_prec_const_mod.o psb_s_ilu_fact_mod.o
psb_z_invt_fact_mod.o: psb_prec_const_mod.o psb_z_ilu_fact_mod.o
psb_c_invk_fact_mod.o: psb_prec_const_mod.o psb_c_ilu_fact_mod.o
psb_d_invk_fact_mod.o: psb_prec_const_mod.o psb_d_ilu_fact_mod.o
psb_s_invk_fact_mod.o: psb_prec_const_mod.o psb_s_ilu_fact_mod.o
psb_z_invk_fact_mod.o: psb_prec_const_mod.o psb_z_ilu_fact_mod.o

veryclean: clean
	/bin/rm -f $(LIBNAME) *$(.mod)

iclean:
	cd impl && $(MAKE) clean
clean: iclean
	/bin/rm -f $(OBJS) $(LOCAL_MODS)
