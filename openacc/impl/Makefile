include ../../Make.inc
LIBDIR=../../lib
INCDIR=../../include
MODDIR=../../modules
#
# Compilers and such
#
#CCOPT= -g
FINCLUDES=$(FMFLAG).. $(FMFLAG)$(MODDIR) $(FMFLAG)$(INCDIR) $(FIFLAG).. 
LIBNAME=libpsb_openacc.a

OBJS= psb_d_oacc_csr_vect_mv.o psb_d_oacc_csr_inner_vect_sv.o \
	psb_d_oacc_csr_csmm.o psb_d_oacc_csr_csmv.o psb_d_oacc_csr_scals.o \
	psb_d_oacc_csr_scal.o psb_d_oacc_csr_allocate_mnnz.o \
	psb_d_oacc_csr_reallocate_nz.o psb_d_oacc_csr_cp_from_coo.o \
	psb_d_oacc_csr_cp_from_fmt.o psb_d_oacc_csr_mv_from_coo \
	psb_d_oacc_csr_mv_from_fmt psb_d_oacc_csr_mold.o \
	psb_d_oacc_mlt_v_2.o psb_d_oacc_mlt_v.o

objs: $(OBJS)

lib: objs
	ar cur ../$(LIBNAME) $(OBJS)

clean:
	/bin/rm -f $(OBJS)

.c.o:
	$(CC) $(CCOPT) $(CCOPENACC) $(CINCLUDES) $(CDEFINES) -c $< -o $@
.f90.o:
	$(FC) $(FCOPT) $(FCOPENACC) $(FINCLUDES) -c $< -o $@
.F90.o:
	$(FC) $(FCOPT) $(FCOPENACC) $(FINCLUDES) $(FDEFINES) -c $< -o $@
.cpp.o:
	$(CXX) $(CXXOPT) $(CXXOPENACC) $(CXXINCLUDES) $(CXXDEFINES) -c $< -o $@
