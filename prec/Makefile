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


LIBNAME=$(PRECLIBNAME)
COBJS=
FINCLUDES=$(FMFLAG). $(FMFLAG)$(LIBDIR) 
OBJS=$(F90OBJS) $(COBJS) $(MPFOBJS) $(MODOBJS)

lib:  $(OBJS) impld
	$(AR) $(HERE)/$(LIBNAME) $(OBJS)
	$(RANLIB) $(HERE)/$(LIBNAME)
	/bin/cp -p $(CPUPDFLAG) $(HERE)/$(LIBNAME) $(LIBDIR)
	/bin/cp -p $(CPUPDFLAG) *$(.mod) $(LIBDIR)

impld: $(OBJS)
	cd impl && $(MAKE)

$(OBJS): $(LIBDIR)/$(BASEMODNAME)$(.mod)


psb_s_base_prec_mod.o psb_d_base_prec_mod.o psb_c_base_prec_mod.o psb_z_base_prec_mod.o: psb_prec_const_mod.o
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

veryclean: clean
	/bin/rm -f $(LIBNAME)

iclean:
	cd impl && $(MAKE) clean
clean: iclean
	/bin/rm -f $(OBJS) $(LOCAL_MODS)

