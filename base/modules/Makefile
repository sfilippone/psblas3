include ../../Make.inc

BASIC_MODS= psb_const_mod.o psb_error_mod.o psb_realloc_mod.o
COMMINT=psi_comm_buffers_mod.o psi_penv_mod.o psi_bcast_mod.o psi_reduce_mod.o psi_p2p_mod.o
UTIL_MODS = aux/psb_string_mod.o desc/psb_desc_const_mod.o desc/psb_indx_map_mod.o\
        desc/psb_gen_block_map_mod.o desc/psb_list_map_mod.o desc/psb_repl_map_mod.o\
        desc/psb_glist_map_mod.o desc/psb_hash_map_mod.o \
	desc/psb_desc_mod.o aux/psb_sort_mod.o \
	serial/psb_s_serial_mod.o serial/psb_d_serial_mod.o serial/psb_c_serial_mod.o serial/psb_z_serial_mod.o \
	serial/psb_serial_mod.o \
	tools/psb_cd_tools_mod.o tools/psb_i_tools_mod.o tools/psb_s_tools_mod.o  tools/psb_d_tools_mod.o\
	tools/psb_c_tools_mod.o  tools/psb_z_tools_mod.o  tools/psb_tools_mod.o \
	psb_penv_mod.o $(COMMINT) psb_error_impl.o \
	comm/psb_base_linmap_mod.o comm/psb_linmap_mod.o \
	comm/psb_s_linmap_mod.o comm/psb_d_linmap_mod.o comm/psb_c_linmap_mod.o comm/psb_z_linmap_mod.o \
	comm/psb_comm_mod.o comm/psb_i_comm_mod.o comm/psb_s_comm_mod.o comm/psb_d_comm_mod.o\
	comm/psb_c_comm_mod.o comm/psb_z_comm_mod.o \
	serial/psb_i_base_vect_mod.o serial/psb_i_vect_mod.o\
	serial/psb_d_base_vect_mod.o serial/psb_d_vect_mod.o\
	serial/psb_s_base_vect_mod.o serial/psb_s_vect_mod.o\
	serial/psb_c_base_vect_mod.o serial/psb_c_vect_mod.o\
	serial/psb_z_base_vect_mod.o serial/psb_z_vect_mod.o\
	serial/psb_vect_mod.o\
	psblas/psb_s_psblas_mod.o psblas/psb_c_psblas_mod.o \
	psblas/psb_d_psblas_mod.o psblas/psb_z_psblas_mod.o psblas/psb_psblas_mod.o \
	aux/psi_serial_mod.o aux/psi_i_serial_mod.o \
	aux/psi_s_serial_mod.o aux/psi_d_serial_mod.o aux/psi_c_serial_mod.o aux/psi_z_serial_mod.o \
	psi_mod.o psi_i_mod.o psi_s_mod.o psi_d_mod.o psi_c_mod.o psi_z_mod.o\
	aux/psb_ip_reord_mod.o\
	aux/psb_i_sort_mod.o aux/psb_s_sort_mod.o aux/psb_d_sort_mod.o \
	aux/psb_c_sort_mod.o aux/psb_z_sort_mod.o \
	psb_check_mod.o aux/psb_hash_mod.o\
	serial/psb_base_mat_mod.o  serial/psb_mat_mod.o\
	serial/psb_s_base_mat_mod.o serial/psb_s_csr_mat_mod.o serial/psb_s_csc_mat_mod.o serial/psb_s_mat_mod.o \
	serial/psb_d_base_mat_mod.o serial/psb_d_csr_mat_mod.o serial/psb_d_csc_mat_mod.o serial/psb_d_mat_mod.o \
	serial/psb_c_base_mat_mod.o serial/psb_c_csr_mat_mod.o serial/psb_c_csc_mat_mod.o serial/psb_c_mat_mod.o \
	serial/psb_z_base_mat_mod.o serial/psb_z_csr_mat_mod.o serial/psb_z_csc_mat_mod.o serial/psb_z_mat_mod.o 


MODULES=$(BASIC_MODS) $(UTIL_MODS)
OBJS = error.o psb_base_mod.o $(EXTRA_COBJS) cutil.o 
LIBDIR=..
CINCLUDES=-I. 
FINCLUDES=$(FMFLAG)$(LIBDIR) $(FMFLAG). $(FIFLAG).


lib:  $(LIBDIR)/$(LIBNAME) 
	/bin/cp -p $(CPUPDFLAG) *$(.mod) $(LIBDIR)

$(LIBDIR)/$(LIBNAME): $(MODULES) $(OBJS) $(MPFOBJS)
	$(AR) $(LIBDIR)/$(LIBNAME) $(MODULES) $(OBJS) $(MPFOBJS)
	$(RANLIB) $(LIBDIR)/$(LIBNAME)


psb_error_mod.o: psb_const_mod.o
psb_realloc_mod.o:  psb_error_mod.o
$(UTIL_MODS): $(BASIC_MODS)

psi_penv_mod.o: psi_comm_buffers_mod.o 
psi_bcast_mod.o psi_reduce_mod.o psi_p2p_mod.o: psi_penv_mod.o	


aux/psb_string_mod.o desc/psb_desc_const_mod.o psi_comm_buffers_mod.o: psb_const_mod.o
aux/psb_hash_mod.o: psb_realloc_mod.o psb_const_mod.o
aux/psb_i_sort_mod.o aux/psb_s_sort_mod.o aux/psb_d_sort_mod.o aux/psb_c_sort_mod.o aux/psb_z_sort_mod.o \
aux/psb_ip_reord_mod.o aux/psi_serial_mod.o aux/psb_sort_mod.o: $(BASIC_MODS)
aux/psb_sort_mod.o: aux/psb_i_sort_mod.o aux/psb_s_sort_mod.o aux/psb_d_sort_mod.o \
	aux/psb_c_sort_mod.o aux/psb_z_sort_mod.o aux/psb_ip_reord_mod.o aux/psi_serial_mod.o
aux/psi_serial_mod.o:  aux/psi_i_serial_mod.o \
	aux/psi_s_serial_mod.o aux/psi_d_serial_mod.o aux/psi_c_serial_mod.o aux/psi_z_serial_mod.o 
aux/psi_i_serial_mod.o aux/psi_s_serial_mod.o aux/psi_d_serial_mod.o aux/psi_c_serial_mod.o aux/psi_z_serial_mod.o: psb_const_mod.o

serial/psb_base_mat_mod.o: aux/psi_serial_mod.o
serial/psb_s_base_mat_mod.o serial/psb_d_base_mat_mod.o serial/psb_c_base_mat_mod.o serial/psb_z_base_mat_mod.o: serial/psb_base_mat_mod.o
serial/psb_s_base_mat_mod.o: serial/psb_s_base_vect_mod.o
serial/psb_d_base_mat_mod.o: serial/psb_d_base_vect_mod.o 
serial/psb_c_base_mat_mod.o: serial/psb_c_base_vect_mod.o
serial/psb_z_base_mat_mod.o: serial/psb_z_base_vect_mod.o
serial/psb_c_base_vect_mod.o serial/psb_s_base_vect_mod.o serial/psb_d_base_vect_mod.o serial/psb_z_base_vect_mod.o: serial/psb_i_base_vect_mod.o 
serial/psb_i_base_vect_mod.o serial/psb_c_base_vect_mod.o serial/psb_s_base_vect_mod.o serial/psb_d_base_vect_mod.o serial/psb_z_base_vect_mod.o: aux/psi_serial_mod.o psb_realloc_mod.o
serial/psb_s_mat_mod.o: serial/psb_s_base_mat_mod.o serial/psb_s_csr_mat_mod.o serial/psb_s_csc_mat_mod.o serial/psb_s_vect_mod.o
serial/psb_d_mat_mod.o: serial/psb_d_base_mat_mod.o serial/psb_d_csr_mat_mod.o serial/psb_d_csc_mat_mod.o serial/psb_d_vect_mod.o serial/psb_i_vect_mod.o
serial/psb_c_mat_mod.o: serial/psb_c_base_mat_mod.o serial/psb_c_csr_mat_mod.o serial/psb_c_csc_mat_mod.o serial/psb_c_vect_mod.o
serial/psb_z_mat_mod.o: serial/psb_z_base_mat_mod.o serial/psb_z_csr_mat_mod.o serial/psb_z_csc_mat_mod.o serial/psb_z_vect_mod.o
serial/psb_s_csc_mat_mod.o serial/psb_s_csr_mat_mod.o: serial/psb_s_base_mat_mod.o 
serial/psb_d_csc_mat_mod.o serial/psb_d_csr_mat_mod.o: serial/psb_d_base_mat_mod.o 
serial/psb_c_csc_mat_mod.o serial/psb_c_csr_mat_mod.o: serial/psb_c_base_mat_mod.o 
serial/psb_z_csc_mat_mod.o serial/psb_z_csr_mat_mod.o: serial/psb_z_base_mat_mod.o 
serial/psb_mat_mod.o: serial/psb_vect_mod.o serial/psb_s_mat_mod.o serial/psb_d_mat_mod.o serial/psb_c_mat_mod.o serial/psb_z_mat_mod.o
serial/psb_serial_mod.o: serial/psb_s_serial_mod.o serial/psb_d_serial_mod.o serial/psb_c_serial_mod.o serial/psb_z_serial_mod.o 
serial/psb_i_vect_mod.o: serial/psb_i_base_vect_mod.o
serial/psb_s_vect_mod.o: serial/psb_s_base_vect_mod.o serial/psb_i_vect_mod.o
serial/psb_d_vect_mod.o: serial/psb_d_base_vect_mod.o serial/psb_i_vect_mod.o
serial/psb_c_vect_mod.o: serial/psb_c_base_vect_mod.o serial/psb_i_vect_mod.o
serial/psb_z_vect_mod.o: serial/psb_z_base_vect_mod.o serial/psb_i_vect_mod.o
serial/psb_s_serial_mod.o serial/psb_d_serial_mod.o serial/psb_c_serial_mod.o serial/psb_z_serial_mod.o: serial/psb_mat_mod.o aux/psb_string_mod.o aux/psb_sort_mod.o aux/psi_serial_mod.o
serial/psb_vect_mod.o: serial/psb_i_vect_mod.o serial/psb_d_vect_mod.o serial/psb_s_vect_mod.o serial/psb_c_vect_mod.o serial/psb_z_vect_mod.o

error.o psb_realloc_mod.o: psb_error_mod.o 
psb_error_impl.o: psb_penv_mod.o
psb_spmat_type.o: aux/psb_string_mod.o aux/psb_sort_mod.o

desc/psb_desc_mod.o: psb_penv_mod.o psb_realloc_mod.o\
	  aux/psb_hash_mod.o desc/psb_hash_map_mod.o desc/psb_list_map_mod.o \
	  desc/psb_repl_map_mod.o desc/psb_gen_block_map_mod.o desc/psb_desc_const_mod.o\
	  desc/psb_indx_map_mod.o serial/psb_i_vect_mod.o
psi_i_mod.o: desc/psb_desc_mod.o serial/psb_i_vect_mod.o
psi_s_mod.o: desc/psb_desc_mod.o serial/psb_s_vect_mod.o
psi_d_mod.o: desc/psb_desc_mod.o serial/psb_d_vect_mod.o
psi_c_mod.o: desc/psb_desc_mod.o serial/psb_c_vect_mod.o
psi_z_mod.o: desc/psb_desc_mod.o serial/psb_z_vect_mod.o
psi_mod.o: psb_penv_mod.o desc/psb_desc_mod.o aux/psi_serial_mod.o serial/psb_serial_mod.o\
	psi_i_mod.o psi_s_mod.o psi_d_mod.o psi_c_mod.o psi_z_mod.o

desc/psb_indx_map_mod.o: desc/psb_desc_const_mod.o psb_error_mod.o psb_penv_mod.o
desc/psb_hash_map_mod.o desc/psb_list_map_mod.o desc/psb_repl_map_mod.o desc/psb_gen_block_map_mod.o:\
    desc/psb_indx_map_mod.o desc/psb_desc_const_mod.o  \
    aux/psb_sort_mod.o psb_penv_mod.o 
desc/psb_glist_map_mod.o: desc/psb_list_map_mod.o
desc/psb_hash_map_mod.o: aux/psb_hash_mod.o aux/psb_sort_mod.o
desc/psb_gen_block_map_mod.o: aux/psb_hash_mod.o

psb_check_mod.o: desc/psb_desc_mod.o

comm/psb_linmap_mod.o: comm/psb_s_linmap_mod.o comm/psb_d_linmap_mod.o comm/psb_c_linmap_mod.o comm/psb_z_linmap_mod.o 
comm/psb_s_linmap_mod.o: comm/psb_base_linmap_mod.o serial/psb_s_mat_mod.o serial/psb_s_vect_mod.o
comm/psb_d_linmap_mod.o: comm/psb_base_linmap_mod.o serial/psb_d_mat_mod.o serial/psb_d_vect_mod.o
comm/psb_c_linmap_mod.o: comm/psb_base_linmap_mod.o serial/psb_c_mat_mod.o serial/psb_c_vect_mod.o
comm/psb_z_linmap_mod.o: comm/psb_base_linmap_mod.o serial/psb_z_mat_mod.o serial/psb_z_vect_mod.o
comm/psb_base_linmap_mod.o: desc/psb_desc_mod.o serial/psb_serial_mod.o comm/psb_comm_mod.o 
comm/psb_comm_mod.o: desc/psb_desc_mod.o serial/psb_mat_mod.o 
comm/psb_comm_mod.o: comm/psb_i_comm_mod.o comm/psb_s_comm_mod.o comm/psb_d_comm_mod.o comm/psb_c_comm_mod.o comm/psb_z_comm_mod.o 
comm/psb_i_comm_mod.o: serial/psb_i_vect_mod.o desc/psb_desc_mod.o 
comm/psb_s_comm_mod.o: serial/psb_s_vect_mod.o desc/psb_desc_mod.o serial/psb_mat_mod.o 
comm/psb_d_comm_mod.o: serial/psb_d_vect_mod.o desc/psb_desc_mod.o serial/psb_mat_mod.o 
comm/psb_c_comm_mod.o: serial/psb_c_vect_mod.o desc/psb_desc_mod.o serial/psb_mat_mod.o 
comm/psb_z_comm_mod.o: serial/psb_z_vect_mod.o desc/psb_desc_mod.o serial/psb_mat_mod.o 

tools/psb_tools_mod.o: tools/psb_cd_tools_mod.o tools/psb_s_tools_mod.o  tools/psb_d_tools_mod.o\
	tools/psb_i_tools_mod.o tools/psb_c_tools_mod.o  tools/psb_z_tools_mod.o  
tools/psb_cd_tools_mod.o tools/psb_i_tools_mod.o tools/psb_s_tools_mod.o tools/psb_d_tools_mod.o tools/psb_c_tools_mod.o tools/psb_z_tools_mod.o: desc/psb_desc_mod.o psi_mod.o  serial/psb_mat_mod.o
tools/psb_i_tools_mod.o: serial/psb_i_vect_mod.o 
tools/psb_s_tools_mod.o: serial/psb_s_vect_mod.o 
tools/psb_d_tools_mod.o: serial/psb_d_vect_mod.o 
tools/psb_c_tools_mod.o: serial/psb_c_vect_mod.o 
tools/psb_z_tools_mod.o: serial/psb_z_vect_mod.o 

psblas/psb_s_psblas_mod.o: serial/psb_s_vect_mod.o serial/psb_s_mat_mod.o
psblas/psb_d_psblas_mod.o: serial/psb_d_vect_mod.o serial/psb_d_mat_mod.o
psblas/psb_c_psblas_mod.o: serial/psb_c_vect_mod.o serial/psb_c_mat_mod.o
psblas/psb_z_psblas_mod.o: serial/psb_z_vect_mod.o serial/psb_z_mat_mod.o
psblas/psb_psblas_mod.o: psblas/psb_s_psblas_mod.o psblas/psb_c_psblas_mod.o  psblas/psb_d_psblas_mod.o psblas/psb_z_psblas_mod.o
psblas/psb_s_psblas_mod.o psblas/psb_c_psblas_mod.o psblas/psb_d_psblas_mod.o psblas/psb_z_psblas_mod.o: serial/psb_mat_mod.o desc/psb_desc_mod.o


psb_base_mod.o: $(MODULES)

psi_penv_mod.o: psi_penv_mod.F90 $(BASIC_MODS) 
	$(FC)   $(FINCLUDES) $(FDEFINES) $(FCOPT) $(EXTRA_OPT) -c $< -o $@

psb_penv_mod.o: psb_penv_mod.F90 $(COMMINT) $(BASIC_MODS)
	$(FC)   $(FINCLUDES) $(FDEFINES) $(FCOPT) $(EXTRA_OPT) -c $< -o $@

psi_comm_buffers_mod.o: psi_comm_buffers_mod.F90 $(BASIC_MODS)
	$(FC)   $(FINCLUDES) $(FDEFINES) $(FCOPT) $(EXTRA_OPT) -c $< -o $@

psi_p2p_mod.o: psi_p2p_mod.F90 $(BASIC_MODS) 
	$(FC)   $(FINCLUDES) $(FDEFINES) $(FCOPT) $(EXTRA_OPT) -c $< -o $@

psi_bcast_mod.o: psi_bcast_mod.F90 $(BASIC_MODS) 
	$(FC)   $(FINCLUDES) $(FDEFINES) $(FCOPT) $(EXTRA_OPT) -c $< -o $@

psi_reduce_mod.o: psi_reduce_mod.F90 $(BASIC_MODS)
	$(FC)   $(FINCLUDES) $(FDEFINES) $(FCOPT) $(EXTRA_OPT) -c $< -o $@

clean:
	/bin/rm -f $(MODULES) $(OBJS) $(MPFOBJS) *$(.mod)

veryclean: clean

