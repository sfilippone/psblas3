include ../../Make.inc


LIBDIR=../../lib/

MPFOBJS=psb_dilu_bld.o psb_dbldaggrmat.o 
F90OBJS=psb_dasmatbld.o psb_dslu_bld.o psb_dumf_bld.o psb_dilu_fct.o\
	psb_dmlprc_bld.o psb_dsp_renum.o\
	psb_dprec.o psb_dprecbld.o gps.o psb_dprecfree.o psb_dprecset.o \
	psb_dbaseprc_bld.o psb_dgenaggrmap.o $(MPFOBJS) 

COBJS=psb_slu_impl.o psb_umf_impl.o
INCDIRS=-I. -I.. -I$(LIBDIR) 

OBJS=$(F90OBJS) $(COBJS)

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
