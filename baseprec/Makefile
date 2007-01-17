include ../Make.inc

LIBDIR=../lib
HERE=.
MODOBJS= psb_prec_type.o  psb_prec_mod.o
F90OBJS= psb_dilu_bld.o psb_dilu_fct.o\
	psb_dsp_renum.o\
	psb_dprecbld.o psb_dprecfree.o psb_dprecset.o \
	psb_dbaseprc_bld.o psb_ddiagsc_bld.o \
	psb_dprc_aply.o \
	psb_dbaseprc_aply.o psb_dbjac_aply.o\
	psb_zilu_bld.o psb_zilu_fct.o\
	psb_zsp_renum.o\
	psb_zprecbld.o psb_zprecfree.o psb_zprecset.o \
	psb_zbaseprc_bld.o psb_zdiagsc_bld.o  \
	psb_zprc_aply.o psb_zbaseprc_aply.o psb_zbjac_aply.o

LIBMOD=psb_prec_mod$(.mod)
LOCAL_MODS=$(MODOBJS:.o=$(.mod))
LIBNAME=$(PRECLIBNAME)
COBJS=
INCDIRS=-I. -I$(LIBDIR) 
OBJS=$(F90OBJS) $(COBJS) $(MPFOBJS) $(MODOBJS)

lib:  $(OBJS) 
	$(AR) $(HERE)/$(LIBNAME) $(OBJS)
	$(RANLIB) $(HERE)/$(LIBNAME)
	/bin/cp -p $(HERE)/$(LIBNAME) $(LIBDIR)
	/bin/cp -p $(LIBMOD) $(LIBDIR)

$(F90OBJS): $(MODOBJS)
psb_prec_mod.o: psb_prec_type.o

veryclean: clean
	/bin/rm -f $(LIBNAME)

clean:
	/bin/rm -f $(OBJS) $(LOCAL_MODS)

