include ../../Make.inc

MODULES = psb_realloc_mod.o psb_string_mod.o psb_spmat_type.o \
	psb_desc_type.o psb_spsb_mod.o \
	psb_blacs_mod.o psb_serial_mod.o psb_tools_mod.o \
	psb_prec_type.o psb_error_mod.o psb_prec_mod.o \
	psb_methd_mod.o psb_const_mod.o \
	psb_comm_mod.o psb_psblas_mod.o psi_mod.o \
	psb_sparse_mod.o psb_check_mod.o psb_all_mod.o

OBJS = error.o

INCDIRS = -I ../../lib
LIBDIR = ../../lib

psb_realloc_mod.o : psb_error_mod.o 
psb_spmat_type.o : psb_realloc_mod.o psb_const_mod.o
psb_error_mod.o: psb_const_mod.o
psb_const_mod.f90: psb_const.fh
psb_blacs_mod.o	: psb_const_mod.o psb_error_mod.o

lib: $(MODULES) $(OBJS)
	$(AR) $(LIBDIR)/$(LIBNAME) $(MODULES) $(OBJS)
	$(RANLIB) $(LIBDIR)/$(LIBNAME)
	cp *$(.mod) ./psb_const.fh ./parts.fh ../../lib


clean:
	/bin/rm -f $(MODULES) $(OBJS) *$(.mod)

veryclean: clean

