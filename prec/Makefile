include ../Make.inc

LIBDIR=../lib
HERE=.
MODOBJS= psb_prec_type.o  psb_prec_mod.o psb_d_diagprec.o psb_d_nullprec.o \
	psb_d_bjacprec.o
F90OBJS=  psb_dilu_fct.o\
	psb_dprecbld.o psb_dprecset.o psb_dprecinit.o \
	psb_sbjac_bld.o psb_silu_fct.o\
	psb_sprecbld.o psb_sprecset.o psb_sprecinit.o \
	psb_sdiagsc_bld.o \
	psb_sprc_aply.o psb_sgprec_aply.o psb_sbjac_aply.o\
	psb_cbjac_bld.o psb_cilu_fct.o\
	psb_cprecbld.o psb_cprecset.o psb_cprecinit.o \
	psb_cdiagsc_bld.o  \
	psb_cprc_aply.o psb_cgprec_aply.o psb_cbjac_aply.o\
	psb_zbjac_bld.o psb_zilu_fct.o\
	psb_zprecbld.o psb_zprecset.o psb_zprecinit.o \
	psb_zdiagsc_bld.o  \
	psb_zprc_aply.o psb_zgprec_aply.o psb_zbjac_aply.o

LIBMOD=psb_prec_mod$(.mod)
LOCAL_MODS=$(MODOBJS:.o=$(.mod))
LIBNAME=$(PRECLIBNAME)
COBJS=
FINCLUDES=$(FMFLAG)$(LIBDIR) $(FMFLAG).
OBJS=$(F90OBJS) $(COBJS) $(MPFOBJS) $(MODOBJS)

lib:  $(OBJS) 
	$(AR) $(HERE)/$(LIBNAME) $(OBJS)
	$(RANLIB) $(HERE)/$(LIBNAME)
	/bin/cp -p $(HERE)/$(LIBNAME) $(LIBDIR)
	/bin/cp -p $(LIBMOD) $(LOCAL_MODS) $(LIBDIR)

$(OBJS): $(LIBDIR)/psb_base_mod$(.mod)

$(F90OBJS): $(MODOBJS)

psb_prec_mod.o:  psb_prec_type.o
psb_d_bjacprec.o psb_d_diagprec.o psb_d_nullprec.o:  psb_prec_type.o psb_prec_mod.o

veryclean: clean
	/bin/rm -f $(LIBNAME)

clean:
	/bin/rm -f $(OBJS) $(LOCAL_MODS)

