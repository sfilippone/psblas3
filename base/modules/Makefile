include ../../Make.inc

BASIC_MODS= psb_const_mod.o psb_error_mod.o psb_realloc_mod.o
UTIL_MODS = psb_string_mod.o \
	psb_desc_type.o psb_sort_mod.o psb_serial_mod.o \
	psb_base_tools_mod.o psb_s_tools_mod.o  psb_d_tools_mod.o\
	psb_c_tools_mod.o  psb_z_tools_mod.o  psb_tools_mod.o \
	psb_penv_mod.o psi_comm_buffers_mod.o psi_penv_mod.o psi_bcast_mod.o \
	psi_reduce_mod.o psi_p2p_mod.o psb_error_impl.o \
	psb_linmap_type_mod.o psb_linmap_mod.o psb_comm_mod.o\
	psb_s_psblas_mod.o psb_c_psblas_mod.o \
	psb_d_psblas_mod.o psb_z_psblas_mod.o psb_psblas_mod.o \
	psi_serial_mod.o psi_mod.o  psb_ip_reord_mod.o\
	psb_check_mod.o psb_gps_mod.o psb_hash_mod.o\
	psb_base_mat_mod.o  psb_mat_mod.o\
	psb_s_base_mat_mod.o psb_s_csr_mat_mod.o psb_s_csc_mat_mod.o psb_s_mat_mod.o \
	psb_d_base_mat_mod.o psb_d_csr_mat_mod.o psb_d_csc_mat_mod.o psb_d_mat_mod.o \
	psb_c_base_mat_mod.o psb_c_csr_mat_mod.o psb_c_csc_mat_mod.o psb_c_mat_mod.o \
	psb_z_base_mat_mod.o psb_z_csr_mat_mod.o psb_z_csc_mat_mod.o psb_z_mat_mod.o 


MODULES=$(BASIC_MODS) $(UTIL_MODS)

LIBMOD=psb_sparse_mod$(.mod)
OBJS = error.o psb_sparse_mod.o $(EXTRA_COBJS) cutil.o
LIBDIR=..
CINCLUDES=-I. 
FINCLUDES=$(FMFLAG)$(LIBDIR) $(FMFLAG). $(FIFLAG).


lib:  $(BASIC_MODS) penvmod  $(UTIL_MODS)  $(OBJS)  $(LIBMOD)
	$(AR) $(LIBDIR)/$(LIBNAME) $(MODULES) $(OBJS) $(MPFOBJS)
	$(RANLIB) $(LIBDIR)/$(LIBNAME)
	/bin/cp -p $(LIBMOD) $(LIBDIR)
	/bin/cp -p *$(.mod) $(LIBDIR)


psi_penv_mod.o: psi_comm_buffers_mod.o psb_const_mod.o psb_realloc_mod.o 
psi_bcast_mod.o psi_reduce_mod.o psi_p2p_mod.o penvmod.o: psi_penv_mod.o	
psb_penv_mod.o: psi_bcast_mod.o psi_reduce_mod.o psi_p2p_mod.o 

psb_base_mat_mod.o: psb_string_mod.o psb_sort_mod.o psb_ip_reord_mod.o\
	psb_error_mod.o psi_serial_mod.o
psb_s_base_mat_mod.o psb_d_base_mat_mod.o psb_c_base_mat_mod.o psb_z_base_mat_mod.o: psb_base_mat_mod.o
psb_s_mat_mod.o: psb_s_base_mat_mod.o psb_s_csr_mat_mod.o psb_s_csc_mat_mod.o
psb_d_mat_mod.o: psb_d_base_mat_mod.o psb_d_csr_mat_mod.o psb_d_csc_mat_mod.o
psb_c_mat_mod.o: psb_c_base_mat_mod.o psb_c_csr_mat_mod.o psb_c_csc_mat_mod.o
psb_z_mat_mod.o: psb_z_base_mat_mod.o psb_z_csr_mat_mod.o psb_z_csc_mat_mod.o
psb_s_csc_mat_mod.o psb_s_csr_mat_mod.o: psb_s_base_mat_mod.o 
psb_d_csc_mat_mod.o psb_d_csr_mat_mod.o: psb_d_base_mat_mod.o 
psb_c_csc_mat_mod.o psb_c_csr_mat_mod.o: psb_c_base_mat_mod.o 
psb_z_csc_mat_mod.o psb_z_csr_mat_mod.o: psb_z_base_mat_mod.o 
psb_mat_mod.o: psb_s_mat_mod.o psb_d_mat_mod.o psb_c_mat_mod.o psb_z_mat_mod.o
psb_realloc_mod.o : psb_error_mod.o 
psb_spmat_type.o : psb_realloc_mod.o psb_error_mod.o psb_const_mod.o psb_string_mod.o psb_sort_mod.o
psb_error_mod.o: psb_const_mod.o
psb_ip_reord_mod.o: psb_const_mod.o
psb_blacs_mod.o: psb_const_mod.o
psi_serial_mod.o: psb_const_mod.o psb_realloc_mod.o
psi_mod.o: psb_penv_mod.o psb_error_mod.o psb_desc_type.o psb_const_mod.o psi_serial_mod.o psb_serial_mod.o
psb_desc_type.o: psb_const_mod.o psb_error_mod.o psb_penv_mod.o psb_realloc_mod.o  psb_hash_mod.o
psb_linmap_mod.o: psb_linmap_type_mod.o psb_mat_mod.o 
psb_linmap_type_mod.o: psb_desc_type.o psb_error_mod.o psb_serial_mod.o psb_comm_mod.o psb_mat_mod.o
psb_comm_mod.o: psb_desc_type.o psb_mat_mod.o
psb_check_mod.o: psb_desc_type.o
psb_serial_mod.o: psb_mat_mod.o psb_string_mod.o psb_sort_mod.o psi_serial_mod.o
psb_sort_mod.o: psb_error_mod.o psb_realloc_mod.o psb_const_mod.o
psb_tools_mod.o: 	psb_base_tools_mod.o psb_s_tools_mod.o  psb_d_tools_mod.o\
	psb_c_tools_mod.o  psb_z_tools_mod.o  
psb_base_tools_mod.o psb_s_tools_mod.o  psb_d_tools_mod.o psb_c_tools_mod.o  psb_z_tools_mod.o: psb_desc_type.o psi_mod.o psb_gps_mod.o psb_mat_mod.o
psb_psblas_mod.o: psb_mat_mod.o psb_desc_type.o psb_s_psblas_mod.o psb_c_psblas_mod.o \
 psb_d_psblas_mod.o psb_z_psblas_mod.o

psb_gps_mod.o:	psb_realloc_mod.o
psb_hash_mod.o: psb_const_mod.o psb_realloc_mod.o

psb_sparse_mod.o: $(MODULES)

newmods: $(BASIC_MODS)
	(cd ../newserial; make lib LIBNAME=$(LIBNAME))


penvmod:
	(make psb_penv_mod.o F90COPT="$(F90COPT) $(EXTRA_OPT)")



clean:
	/bin/rm -f $(MODULES) $(OBJS) $(MPFOBJS) *$(.mod)

veryclean: clean

