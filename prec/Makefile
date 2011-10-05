include ../Make.inc

LIBDIR=../lib
HERE=.
MODOBJS=psb_prec_const_mod.o\
	psb_s_prec_type.o psb_d_prec_type.o psb_c_prec_type.o psb_z_prec_type.o \
	psb_s_base_prec_mod.o psb_d_base_prec_mod.o psb_c_base_prec_mod.o psb_z_base_prec_mod.o \
	psb_prec_type.o  \
	psb_prec_mod.o psb_s_prec_mod.o psb_d_prec_mod.o psb_c_prec_mod.o psb_z_prec_mod.o \
	psb_d_diagprec.o psb_d_nullprec.o psb_d_bjacprec.o \
	psb_s_diagprec.o psb_s_nullprec.o psb_s_bjacprec.o \
	psb_c_diagprec.o psb_c_nullprec.o psb_c_bjacprec.o \
	psb_z_diagprec.o psb_z_nullprec.o psb_z_bjacprec.o 

F90OBJS=  psb_dilu_fct.o\
	psb_dprecbld.o psb_dprecset.o psb_dprecinit.o \
	psb_silu_fct.o\
	psb_sprecbld.o psb_sprecset.o psb_sprecinit.o \
	psb_cilu_fct.o\
	psb_cprecbld.o psb_cprecset.o psb_cprecinit.o \
	psb_zilu_fct.o\
	psb_zprecbld.o psb_zprecset.o psb_zprecinit.o 

LIBMOD=psb_prec_mod$(.mod)
LOCAL_MODS=$(MODOBJS:.o=$(.mod))
LIBNAME=$(PRECLIBNAME)
COBJS=
FINCLUDES=$(FMFLAG)$(LIBDIR) $(FMFLAG).
OBJS=$(F90OBJS) $(COBJS) $(MPFOBJS) $(MODOBJS)

lib:  $(LIBDIR)/$(LIBNAME) $(LIBDIR)/$(LIBMOD)

$(HERE)/$(LIBNAME): $(OBJS)
	$(AR) $(HERE)/$(LIBNAME) $(OBJS)
	$(RANLIB) $(HERE)/$(LIBNAME)

$(LIBDIR)/$(LIBNAME): $(HERE)/$(LIBNAME)
	/bin/cp -p $(CPUPDFLAG) $(HERE)/$(LIBNAME) $(LIBDIR)

$(LIBDIR)/$(LIBMOD):
	/bin/cp -p $(CPUPDFLAG) $(LIBMOD) $(LOCAL_MODS) $(LIBDIR)



$(OBJS): $(LIBDIR)/psb_base_mod$(.mod)

$(F90OBJS): $(MODOBJS)

psb_s_base_prec_mod.o psb_d_base_prec_mod.o psb_c_base_prec_mod.o psb_z_base_prec_mod.o: psb_prec_const_mod.o
psb_s_prec_type.o: psb_s_base_prec_mod.o
psb_d_prec_type.o: psb_d_base_prec_mod.o
psb_c_prec_type.o: psb_c_base_prec_mod.o
psb_z_prec_type.o: psb_z_base_prec_mod.o
psb_s_prec_mod.o:  psb_s_prec_type.o
psb_d_prec_mod.o:  psb_d_prec_type.o
psb_c_prec_mod.o:  psb_c_prec_type.o
psb_z_prec_mod.o:  psb_z_prec_type.o
psb_prec_type.o: psb_s_prec_type.o psb_d_prec_type.o psb_c_prec_type.o psb_z_prec_type.o
psb_prec_mod.o: psb_s_prec_mod.o psb_d_prec_mod.o psb_c_prec_mod.o psb_z_prec_mod.o 
psb_s_bjacprec.o psb_s_diagprec.o psb_s_nullprec.o: psb_prec_mod.o psb_s_base_prec_mod.o 
psb_d_bjacprec.o psb_d_diagprec.o psb_d_nullprec.o: psb_prec_mod.o psb_d_base_prec_mod.o 
psb_c_bjacprec.o psb_c_diagprec.o psb_c_nullprec.o: psb_prec_mod.o psb_c_base_prec_mod.o 
psb_z_bjacprec.o psb_z_diagprec.o psb_z_nullprec.o: psb_prec_mod.o psb_z_base_prec_mod.o 

veryclean: clean
	/bin/rm -f $(LIBNAME)

clean:
	/bin/rm -f $(OBJS) $(LOCAL_MODS)

