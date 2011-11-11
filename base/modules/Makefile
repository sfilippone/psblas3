include ../../Make.inc

BASIC_MODS= psb_const_mod.o psb_error_mod.o psb_realloc_mod.o
UTIL_MODS = psb_string_mod.o psb_desc_const_mod.o psb_indx_map_mod.o\
        psb_gen_block_map_mod.o psb_list_map_mod.o psb_repl_map_mod.o\
        psb_glist_map_mod.o psb_hash_map_mod.o \
	psb_desc_type.o psb_sort_mod.o psb_serial_mod.o \
	psb_base_tools_mod.o psb_s_tools_mod.o  psb_d_tools_mod.o\
	psb_c_tools_mod.o  psb_z_tools_mod.o  psb_tools_mod.o \
	psb_penv_mod.o psi_comm_buffers_mod.o psi_penv_mod.o psi_bcast_mod.o \
	psi_reduce_mod.o psi_p2p_mod.o psb_error_impl.o \
	psb_linmap_type_mod.o psb_linmap_mod.o \
	psb_s_linmap_mod.o psb_d_linmap_mod.o psb_c_linmap_mod.o psb_z_linmap_mod.o \
	psb_comm_mod.o psb_i_comm_mod.o psb_s_comm_mod.o psb_d_comm_mod.o\
	psb_c_comm_mod.o psb_z_comm_mod.o \
	psb_d_base_vect_mod.o psb_d_vect_mod.o\
	psb_s_base_vect_mod.o psb_s_vect_mod.o\
	psb_c_base_vect_mod.o psb_c_vect_mod.o\
	psb_z_base_vect_mod.o psb_z_vect_mod.o\
	psb_vect_mod.o\
	psb_s_psblas_mod.o psb_c_psblas_mod.o \
	psb_d_psblas_mod.o psb_z_psblas_mod.o psb_psblas_mod.o \
	psi_serial_mod.o \
	psi_mod.o psi_i_mod.o psi_s_mod.o psi_d_mod.o psi_c_mod.o psi_z_mod.o\
	psb_ip_reord_mod.o\
	psb_check_mod.o psb_hash_mod.o\
	psb_base_mat_mod.o  psb_mat_mod.o\
	psb_s_base_mat_mod.o psb_s_csr_mat_mod.o psb_s_csc_mat_mod.o psb_s_mat_mod.o \
	psb_d_base_mat_mod.o psb_d_csr_mat_mod.o psb_d_csc_mat_mod.o psb_d_mat_mod.o \
	psb_c_base_mat_mod.o psb_c_csr_mat_mod.o psb_c_csc_mat_mod.o psb_c_mat_mod.o \
	psb_z_base_mat_mod.o psb_z_csr_mat_mod.o psb_z_csc_mat_mod.o psb_z_mat_mod.o 


MODULES=$(BASIC_MODS) $(UTIL_MODS)

LIBMOD=psb_base_mod$(.mod)
OBJS = error.o psb_base_mod.o $(EXTRA_COBJS) cutil.o
LIBDIR=..
CINCLUDES=-I. 
FINCLUDES=$(FMFLAG)$(LIBDIR) $(FMFLAG). $(FIFLAG).


lib:  psb_penv_mod.o  $(UTIL_MODS)  $(OBJS)  $(LIBDIR)/$(LIBNAME) 
	/bin/cp -p $(CPUPDFLAG) $(LIBMOD) $(LIBDIR)
	/bin/cp -p $(CPUPDFLAG) *$(.mod) $(LIBDIR)

$(LIBDIR)/$(LIBNAME): $(MODULES) $(OBJS) $(MPFOBJS)
	$(AR) $(LIBDIR)/$(LIBNAME) $(MODULES) $(OBJS) $(MPFOBJS)
	$(RANLIB) $(LIBDIR)/$(LIBNAME)


psb_error_mod.o: psb_const_mod.o
psb_realloc_mod.o:  psb_error_mod.o
$(UTILS_MODS): $(BASIC_MODS)


psi_comm_buffers_mod.o: psb_const_mod.o
psi_penv_mod.o: psi_comm_buffers_mod.o 
psi_bcast_mod.o psi_reduce_mod.o psi_p2p_mod.o: psi_penv_mod.o	
psb_penv_mod.o: psi_bcast_mod.o psi_reduce_mod.o psi_p2p_mod.o 

psb_base_mat_mod.o: psi_serial_mod.o
psb_s_base_mat_mod.o psb_d_base_mat_mod.o psb_c_base_mat_mod.o psb_z_base_mat_mod.o: psb_base_mat_mod.o
psb_s_base_mat_mod.o: psb_s_base_vect_mod.o
psb_d_base_mat_mod.o: psb_d_base_vect_mod.o
psb_c_base_mat_mod.o: psb_c_base_vect_mod.o
psb_z_base_mat_mod.o: psb_z_base_vect_mod.o
psb_c_base_vect_mod.o psb_s_base_vect_mod.o psb_d_base_vect_mod.o psb_z_base_vect_mod.o: psi_serial_mod.o psb_realloc_mod.o
psb_s_mat_mod.o: psb_s_base_mat_mod.o psb_s_csr_mat_mod.o psb_s_csc_mat_mod.o psb_s_vect_mod.o
psb_d_mat_mod.o: psb_d_base_mat_mod.o psb_d_csr_mat_mod.o psb_d_csc_mat_mod.o psb_d_vect_mod.o
psb_c_mat_mod.o: psb_c_base_mat_mod.o psb_c_csr_mat_mod.o psb_c_csc_mat_mod.o psb_c_vect_mod.o
psb_z_mat_mod.o: psb_z_base_mat_mod.o psb_z_csr_mat_mod.o psb_z_csc_mat_mod.o psb_z_vect_mod.o
psb_s_csc_mat_mod.o psb_s_csr_mat_mod.o: psb_s_base_mat_mod.o 
psb_d_csc_mat_mod.o psb_d_csr_mat_mod.o: psb_d_base_mat_mod.o 
psb_c_csc_mat_mod.o psb_c_csr_mat_mod.o: psb_c_base_mat_mod.o 
psb_z_csc_mat_mod.o psb_z_csr_mat_mod.o: psb_z_base_mat_mod.o 
psb_mat_mod.o: psb_vect_mod.o psb_s_mat_mod.o psb_d_mat_mod.o psb_c_mat_mod.o psb_z_mat_mod.o
error.o psb_realloc_mod.o : psb_error_mod.o 
psb_error_impl.o: psb_penv_mod.o
psb_spmat_type.o: psb_string_mod.o psb_sort_mod.o
psi_i_mod.o: psb_desc_type.o 
psi_s_mod.o: psb_desc_type.o psb_s_base_vect_mod.o
psi_d_mod.o: psb_desc_type.o psb_d_base_vect_mod.o
psi_c_mod.o: psb_desc_type.o psb_c_base_vect_mod.o
psi_z_mod.o: psb_desc_type.o psb_z_base_vect_mod.o
psi_mod.o: psb_penv_mod.o psb_desc_type.o psi_serial_mod.o psb_serial_mod.o\
	psi_i_mod.o psi_s_mod.o psi_d_mod.o psi_c_mod.o psi_z_mod.o
psb_desc_type.o: psb_penv_mod.o psb_realloc_mod.o\
	  psb_hash_mod.o psb_hash_map_mod.o psb_list_map_mod.o \
	  psb_repl_map_mod.o psb_gen_block_map_mod.o psb_desc_const_mod.o\
	  psb_indx_map_mod.o 
psb_indx_map_mod.o: psb_desc_const_mod.o psb_error_mod.o psb_penv_mod.o
psb_hash_map_mod.o psb_list_map_mod.o psb_repl_map_mod.o psb_gen_block_map_mod.o:\
    psb_indx_map_mod.o psb_desc_const_mod.o  \
    psb_sort_mod.o psb_penv_mod.o 
psb_glist_map_mod.o: psb_list_map_mod.o
psb_hash_map_mod.o: psb_hash_mod.o psb_sort_mod.o
psb_hash_mod.o: psb_realloc_mod.o psb_const_mod.o
psb_linmap_mod.o: psb_s_linmap_mod.o psb_d_linmap_mod.o psb_c_linmap_mod.o psb_z_linmap_mod.o 
psb_s_linmap_mod.o: psb_linmap_type_mod.o psb_mat_mod.o psb_s_vect_mod.o
psb_d_linmap_mod.o: psb_linmap_type_mod.o psb_mat_mod.o psb_d_vect_mod.o
psb_c_linmap_mod.o: psb_linmap_type_mod.o psb_mat_mod.o psb_c_vect_mod.o
psb_z_linmap_mod.o: psb_linmap_type_mod.o psb_mat_mod.o psb_z_vect_mod.o
psb_linmap_type_mod.o: psb_desc_type.o psb_serial_mod.o psb_comm_mod.o psb_mat_mod.o
psb_comm_mod.o: psb_desc_type.o psb_mat_mod.o 
psb_check_mod.o: psb_desc_type.o
psb_serial_mod.o: psb_mat_mod.o psb_string_mod.o psb_sort_mod.o psi_serial_mod.o
psb_s_vect_mod.o: psb_s_base_vect_mod.o
psb_d_vect_mod.o: psb_d_base_vect_mod.o
psb_c_vect_mod.o: psb_c_base_vect_mod.o
psb_z_vect_mod.o: psb_z_base_vect_mod.o
psb_tools_mod.o: psb_base_tools_mod.o psb_s_tools_mod.o  psb_d_tools_mod.o\
	psb_c_tools_mod.o  psb_z_tools_mod.o  
psb_base_tools_mod.o psb_s_tools_mod.o psb_d_tools_mod.o psb_c_tools_mod.o psb_z_tools_mod.o: psb_desc_type.o psi_mod.o  psb_mat_mod.o
psb_s_tools_mod.o: psb_s_vect_mod.o 
psb_d_tools_mod.o: psb_d_vect_mod.o 
psb_c_tools_mod.o: psb_c_vect_mod.o 
psb_z_tools_mod.o: psb_z_vect_mod.o 

psb_s_psblas_mod.o: psb_s_vect_mod.o psb_s_mat_mod.o
psb_d_psblas_mod.o: psb_d_vect_mod.o psb_d_mat_mod.o
psb_c_psblas_mod.o: psb_c_vect_mod.o psb_c_mat_mod.o
psb_z_psblas_mod.o: psb_z_vect_mod.o psb_z_mat_mod.o
psb_psblas_mod.o: psb_s_psblas_mod.o psb_c_psblas_mod.o  psb_d_psblas_mod.o psb_z_psblas_mod.o
psb_s_psblas_mod.o psb_c_psblas_mod.o psb_d_psblas_mod.o psb_z_psblas_mod.o: psb_mat_mod.o psb_desc_type.o
psb_vect_mod.o: psb_d_vect_mod.o psb_s_vect_mod.o psb_c_vect_mod.o psb_z_vect_mod.o
psb_comm_mod.o: psb_i_comm_mod.o psb_s_comm_mod.o psb_d_comm_mod.o psb_c_comm_mod.o psb_z_comm_mod.o 
psb_i_comm_mod.o: psb_desc_type.o 
psb_s_comm_mod.o: psb_s_vect_mod.o psb_desc_type.o psb_mat_mod.o 
psb_d_comm_mod.o: psb_d_vect_mod.o psb_desc_type.o psb_mat_mod.o 
psb_c_comm_mod.o: psb_c_vect_mod.o psb_desc_type.o psb_mat_mod.o 
psb_z_comm_mod.o: psb_z_vect_mod.o psb_desc_type.o psb_mat_mod.o 

psb_base_mod.o: $(MODULES)

psb_penv_mod.o: $(BASIC_MODS) psb_penv_mod.F90
	$(F90) -c psb_penv_mod.F90  $(FINCLUDES) $(FDEFINES) $(F90COPT) $(EXTRA_OPT)



clean:
	/bin/rm -f $(MODULES) $(OBJS) $(MPFOBJS) *$(.mod)

veryclean: clean

