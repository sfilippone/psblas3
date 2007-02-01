include ../../Make.inc

MODULES = psb_realloc_mod.o psb_string_mod.o psb_spmat_type.o \
	psb_desc_type.o psb_spsb_mod.o \
	psb_serial_mod.o psb_tools_mod.o \
	psb_error_mod.o \
	psb_const_mod.o \
	psb_comm_mod.o psb_psblas_mod.o psi_mod.o \
	psb_check_mod.o blacs_env.o psb_gps_mod.o

#	psb_methd_mod.o psb_prec_type.o  psb_prec_mod.o 

LIBMOD=psb_base_mod$(.mod)
MPFOBJS=psb_penv_mod.o 
OBJS = error.o psb_base_mod.o $(EXTRA_COBJS)
LIBDIR = ..
INCDIRS = -I . 


psb_realloc_mod.o : psb_error_mod.o 
psb_spmat_type.o : psb_realloc_mod.o psb_error_mod.o psb_const_mod.o psb_string_mod.o
psb_error_mod.o: psb_const_mod.o
psb_penv_mod.o: psb_const_mod.o psb_error_mod.o psb_realloc_mod.o
psi_mod.o: psb_penv_mod.o psb_error_mod.o psb_desc_type.o
psb_desc_type.o: psb_const_mod.o psb_error_mod.o psb_penv_mod.o
psb_check_mod.o: psb_desc_type.o
psb_methd_mod.o: psb_serial_mod.o psb_desc_type.o psb_prec_type.o
psb_tools_mod.o: psb_spmat_type.o psb_desc_type.o psi_mod.o psb_gps_mod.o
psb_gps_mod.o:	psb_realloc_mod.o
psb_base_mod.o: $(MODULES) $(MPFOBJS)


lib:  mpfobjs $(MODULES)  $(OBJS)
	$(AR) $(LIBDIR)/$(LIBNAME) $(MODULES) $(OBJS) $(MPFOBJS)
	$(RANLIB) $(LIBDIR)/$(LIBNAME)
	/bin/cp -p $(LIBMOD) ./parts.fh $(LIBDIR)
	/bin/cp -p *$(.mod) $(LIBDIR)


mpfobjs: 
	(make $(MPFOBJS) F90="$(MPF90)" FC="$(MPF90)" FCOPT="$(F90COPT)")

clean:
	/bin/rm -f $(MODULES) $(OBJS) $(MPFOBJS) *$(.mod)

veryclean: clean

