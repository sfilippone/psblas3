include ../../Make.inc


LIBDIR=../../lib/

MPFOBJS=psb_dilu_bld.o psb_dbldaggrmat.o psb_zilu_bld.o psb_zbldaggrmat.o 
F90OBJS=psb_dasmatbld.o psb_dslu_bld.o psb_dumf_bld.o psb_dilu_fct.o\
	psb_dmlprc_bld.o psb_dsp_renum.o\
	psb_dprecbld.o gps.o psb_dprecfree.o psb_dprecset.o \
	psb_dbaseprc_bld.o psb_ddiagsc_bld.o psb_dgenaggrmap.o \
	psb_dprc_aply.o psb_dmlprc_aply.o \
	psb_dbaseprc_aply.o psb_dbjac_aply.o\
	psb_zasmatbld.o psb_zslu_bld.o psb_zumf_bld.o psb_zilu_fct.o\
	psb_zmlprc_bld.o psb_zsp_renum.o\
	psb_zprecbld.o gps.o psb_zprecfree.o psb_zprecset.o \
	psb_zbaseprc_bld.o psb_zdiagsc_bld.o psb_zgenaggrmap.o \
	psb_zprc_aply.o psb_zmlprc_aply.o \
	psb_zbaseprc_aply.o psb_zbjac_aply.o\
	$(MPFOBJS) 

COBJS=psb_slu_impl.o psb_umf_impl.o psb_zslu_impl.o psb_zumf_impl.o
INCDIRS=-I. -I.. -I$(LIBDIR) 

OBJS=$(F90OBJS) $(COBJS) $(MPFOBJS)

lib:  mpobjs $(OBJS) 
	$(AR) $(LIBDIR)/$(LIBNAME) $(MPFOBJS) $(OBJS)
	$(RANLIB) $(LIBDIR)/$(LIBNAME)

mpobjs: 
	(make $(MPFOBJS) F90="$(MPF90)" F90COPT="$(F90COPT)")

veryclean: clean
	/bin/rm -f $(LIBNAME)

clean:
	/bin/rm -f $(OBJS) $(LOCAL_MODS)

veryclean: clean
