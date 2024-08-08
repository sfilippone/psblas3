include ../../Make.inc
LIBDIR=../../lib
INCDIR=../../include
MODDIR=../../modules
UP=..
#
# Compilers and such
#
#CCOPT= -g
FINCLUDES=$(FMFLAG).. $(FMFLAG)$(MODDIR) $(FMFLAG)$(INCDIR) $(FIFLAG).. 
LIBNAME=libpsb_openacc.a

OBJS=   \
psb_s_oacc_csr_vect_mv.o \
psb_s_oacc_csr_inner_vect_sv.o \
psb_s_oacc_csr_csmm.o \
psb_s_oacc_csr_csmv.o \
psb_s_oacc_csr_scals.o \
psb_s_oacc_csr_scal.o \
psb_s_oacc_csr_allocate_mnnz.o \
psb_s_oacc_csr_reallocate_nz.o \
psb_s_oacc_csr_cp_from_coo.o \
psb_s_oacc_csr_cp_from_fmt.o \
psb_s_oacc_csr_mv_from_coo.o \
psb_s_oacc_csr_mv_from_fmt.o \
psb_s_oacc_csr_mold.o \
psb_s_oacc_mlt_v_2.o \
psb_s_oacc_mlt_v.o \
psb_d_oacc_csr_vect_mv.o \
psb_d_oacc_csr_inner_vect_sv.o \
psb_d_oacc_csr_csmm.o \
psb_d_oacc_csr_csmv.o \
psb_d_oacc_csr_scals.o \
psb_d_oacc_csr_scal.o \
psb_d_oacc_csr_allocate_mnnz.o \
psb_d_oacc_csr_reallocate_nz.o \
psb_d_oacc_csr_cp_from_coo.o \
psb_d_oacc_csr_cp_from_fmt.o \
psb_d_oacc_csr_mv_from_coo.o \
psb_d_oacc_csr_mv_from_fmt.o \
psb_d_oacc_csr_mold.o \
psb_d_oacc_mlt_v_2.o \
psb_d_oacc_mlt_v.o \
psb_c_oacc_csr_vect_mv.o \
psb_c_oacc_csr_inner_vect_sv.o \
psb_c_oacc_csr_csmm.o \
psb_c_oacc_csr_csmv.o \
psb_c_oacc_csr_scals.o \
psb_c_oacc_csr_scal.o \
psb_c_oacc_csr_allocate_mnnz.o \
psb_c_oacc_csr_reallocate_nz.o \
psb_c_oacc_csr_cp_from_coo.o \
psb_c_oacc_csr_cp_from_fmt.o \
psb_c_oacc_csr_mv_from_coo.o \
psb_c_oacc_csr_mv_from_fmt.o \
psb_c_oacc_csr_mold.o \
psb_c_oacc_mlt_v_2.o \
psb_c_oacc_mlt_v.o \
psb_z_oacc_csr_vect_mv.o \
psb_z_oacc_csr_inner_vect_sv.o \
psb_z_oacc_csr_csmm.o \
psb_z_oacc_csr_csmv.o \
psb_z_oacc_csr_scals.o \
psb_z_oacc_csr_scal.o \
psb_z_oacc_csr_allocate_mnnz.o \
psb_z_oacc_csr_reallocate_nz.o \
psb_z_oacc_csr_cp_from_coo.o \
psb_z_oacc_csr_cp_from_fmt.o \
psb_z_oacc_csr_mv_from_coo.o \
psb_z_oacc_csr_mv_from_fmt.o \
psb_z_oacc_csr_mold.o \
psb_z_oacc_mlt_v_2.o \
psb_z_oacc_mlt_v.o \
psb_s_oacc_ell_vect_mv.o \
psb_s_oacc_ell_inner_vect_sv.o \
psb_s_oacc_ell_csmm.o \
psb_s_oacc_ell_csmv.o \
psb_s_oacc_ell_scals.o \
psb_s_oacc_ell_scal.o \
psb_s_oacc_ell_reallocate_nz.o \
psb_s_oacc_ell_allocate_mnnz.o \
psb_s_oacc_ell_cp_from_coo.o \
psb_s_oacc_ell_cp_from_fmt.o \
psb_s_oacc_ell_mv_from_coo.o \
psb_s_oacc_ell_mv_from_fmt.o \
psb_s_oacc_ell_mold.o \
psb_s_oacc_hll_mold.o \
psb_s_oacc_hll_mv_from_fmt.o \
psb_s_oacc_hll_mv_from_coo.o \
psb_s_oacc_hll_cp_from_fmt.o \
psb_s_oacc_hll_cp_from_coo.o \
psb_s_oacc_hll_allocate_mnnz.o \
psb_s_oacc_hll_reallocate_nz.o \
psb_s_oacc_hll_scal.o \
psb_s_oacc_hll_scals.o \
psb_s_oacc_hll_csmv.o \
psb_s_oacc_hll_csmm.o \
psb_s_oacc_hll_inner_vect_sv.o \
psb_s_oacc_hll_vect_mv.o \
psb_d_oacc_ell_vect_mv.o \
psb_d_oacc_ell_inner_vect_sv.o \
psb_d_oacc_ell_csmm.o \
psb_d_oacc_ell_csmv.o \
psb_d_oacc_ell_scals.o \
psb_d_oacc_ell_scal.o \
psb_d_oacc_ell_reallocate_nz.o \
psb_d_oacc_ell_allocate_mnnz.o \
psb_d_oacc_ell_cp_from_coo.o \
psb_d_oacc_ell_cp_from_fmt.o \
psb_d_oacc_ell_mv_from_coo.o \
psb_d_oacc_ell_mv_from_fmt.o \
psb_d_oacc_ell_mold.o \
psb_d_oacc_hll_mold.o \
psb_d_oacc_hll_mv_from_fmt.o \
psb_d_oacc_hll_mv_from_coo.o \
psb_d_oacc_hll_cp_from_fmt.o \
psb_d_oacc_hll_cp_from_coo.o \
psb_d_oacc_hll_allocate_mnnz.o \
psb_d_oacc_hll_reallocate_nz.o \
psb_d_oacc_hll_scal.o \
psb_d_oacc_hll_scals.o \
psb_d_oacc_hll_csmv.o \
psb_d_oacc_hll_csmm.o \
psb_d_oacc_hll_inner_vect_sv.o \
psb_d_oacc_hll_vect_mv.o \
psb_c_oacc_ell_vect_mv.o \
psb_c_oacc_ell_inner_vect_sv.o \
psb_c_oacc_ell_csmm.o \
psb_c_oacc_ell_csmv.o \
psb_c_oacc_ell_scals.o \
psb_c_oacc_ell_scal.o \
psb_c_oacc_ell_reallocate_nz.o \
psb_c_oacc_ell_allocate_mnnz.o \
psb_c_oacc_ell_cp_from_coo.o \
psb_c_oacc_ell_cp_from_fmt.o \
psb_c_oacc_ell_mv_from_coo.o \
psb_c_oacc_ell_mv_from_fmt.o \
psb_c_oacc_ell_mold.o \
psb_c_oacc_hll_mold.o \
psb_c_oacc_hll_mv_from_fmt.o \
psb_c_oacc_hll_mv_from_coo.o \
psb_c_oacc_hll_cp_from_fmt.o \
psb_c_oacc_hll_cp_from_coo.o \
psb_c_oacc_hll_allocate_mnnz.o \
psb_c_oacc_hll_reallocate_nz.o \
psb_c_oacc_hll_scal.o \
psb_c_oacc_hll_scals.o \
psb_c_oacc_hll_csmv.o \
psb_c_oacc_hll_csmm.o \
psb_c_oacc_hll_inner_vect_sv.o \
psb_c_oacc_hll_vect_mv.o \
psb_z_oacc_ell_vect_mv.o \
psb_z_oacc_ell_inner_vect_sv.o \
psb_z_oacc_ell_csmm.o \
psb_z_oacc_ell_csmv.o \
psb_z_oacc_ell_scals.o \
psb_z_oacc_ell_scal.o \
psb_z_oacc_ell_reallocate_nz.o \
psb_z_oacc_ell_allocate_mnnz.o \
psb_z_oacc_ell_cp_from_coo.o \
psb_z_oacc_ell_cp_from_fmt.o \
psb_z_oacc_ell_mv_from_coo.o \
psb_z_oacc_ell_mv_from_fmt.o \
psb_z_oacc_ell_mold.o \
psb_z_oacc_hll_mold.o \
psb_z_oacc_hll_mv_from_fmt.o \
psb_z_oacc_hll_mv_from_coo.o \
psb_z_oacc_hll_cp_from_fmt.o \
psb_z_oacc_hll_cp_from_coo.o \
psb_z_oacc_hll_allocate_mnnz.o \
psb_z_oacc_hll_reallocate_nz.o \
psb_z_oacc_hll_scal.o \
psb_z_oacc_hll_scals.o \
psb_z_oacc_hll_csmv.o \
psb_z_oacc_hll_csmm.o \
psb_z_oacc_hll_inner_vect_sv.o \
psb_z_oacc_hll_vect_mv.o 

objs: $(OBJS)

lib: objs
	ar cur ../$(LIBNAME) $(OBJS)


#psb_d_oacc_csr_vect_mv.o psb_d_oacc_csr_inner_vect_sv.o \
	psb_d_oacc_csr_csmm.o psb_d_oacc_csr_csmv.o psb_d_oacc_csr_scals.o \
	psb_d_oacc_csr_scal.o psb_d_oacc_csr_allocate_mnnz.o \
	psb_d_oacc_csr_reallocate_nz.o psb_d_oacc_csr_cp_from_coo.o \
	psb_d_oacc_csr_cp_from_fmt.o psb_d_oacc_csr_mv_from_coo.o \
	psb_d_oacc_csr_mv_from_fmt.o psb_d_oacc_csr_mold.o: $(UP)/psb_d_oacc_csr_mat_mod.o

#psb_d_oacc_mlt_v_2.o psb_d_oacc_mlt_v.o: $(UP)/psb_d_oacc_vect_mod.o


clean:
	/bin/rm -f $(OBJS) *.smod

.c.o:
	$(CC) $(CCOPT) $(CCOPENACC) $(CINCLUDES) $(CDEFINES) -c $< -o $@
.f90.o:
	$(FC) $(FCOPT) $(FCOPENACC) $(FINCLUDES) -c $< -o $@
.F90.o:
	$(FC) $(FCOPT) $(FCOPENACC) $(FINCLUDES) $(FDEFINES) -c $< -o $@
.cpp.o:
	$(CXX) $(CXXOPT) $(CXXOPENACC) $(CXXINCLUDES) $(CXXDEFINES) -c $< -o $@
