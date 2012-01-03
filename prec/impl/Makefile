include ../../Make.inc

LIBDIR=../../lib
HERE=..
OBJS=psb_s_prec_type_impl.o psb_d_prec_type_impl.o \
	psb_c_prec_type_impl.o psb_z_prec_type_impl.o 


F90OBJS=  psb_dilu_fct.o\
	psb_dprecbld.o psb_dprecset.o psb_dprecinit.o \
	psb_silu_fct.o\
	psb_sprecbld.o psb_sprecset.o psb_sprecinit.o \
	psb_cilu_fct.o\
	psb_cprecbld.o psb_cprecset.o psb_cprecinit.o \
	psb_zilu_fct.o\
	psb_zprecbld.o psb_zprecset.o psb_zprecinit.o 

LIBNAME=$(PRECLIBNAME)
COBJS=
FINCLUDES=$(FMFLAG).. $(FMFLAG)$(LIBDIR) 

lib:  $(OBJS) 
	$(AR) $(HERE)/$(LIBNAME) $(OBJS)
	$(RANLIB) $(HERE)/$(LIBNAME)

veryclean: clean

clean:
	/bin/rm -f $(OBJS) $(LOCAL_MODS)

