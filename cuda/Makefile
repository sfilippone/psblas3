include ../Make.inc
#
# Libraries used
#
LIBDIR=../lib
INCDIR=../include
MODDIR=../modules
PSBLAS_LIB= -lpsb_util  -lpsb_base
#-lpsb_util -lpsb_krylov -lpsb_prec -lpsb_base 
LDLIBS=$(PSBLDLIBS)
#
# Compilers and such
#
#CCOPT= -g
FINCLUDES=$(FMFLAG).   $(FMFLAG)$(INCDIR)  $(FMFLAG)$(MODDIR) $(FIFLAG).
CINCLUDES=$(SPGPU_INCLUDES) $(CUDA_INCLUDES) -I$(INCDIR)
LIBNAME=libpsb_cuda.a


FOBJS=cusparse_mod.o base_cusparse_mod.o \
 s_cusparse_mod.o  d_cusparse_mod.o  c_cusparse_mod.o  z_cusparse_mod.o \
	psb_vectordev_mod.o core_mod.o \
	psb_s_vectordev_mod.o psb_d_vectordev_mod.o psb_i_vectordev_mod.o\
	psb_c_vectordev_mod.o psb_z_vectordev_mod.o psb_base_vectordev_mod.o \
        elldev_mod.o hlldev_mod.o diagdev_mod.o hdiagdev_mod.o \
	psb_i_cuda_vect_mod.o \
	psb_d_cuda_vect_mod.o  psb_s_cuda_vect_mod.o\
	psb_z_cuda_vect_mod.o  psb_c_cuda_vect_mod.o\
        psb_d_cuda_elg_mat_mod.o   psb_d_cuda_hlg_mat_mod.o \
	psb_d_cuda_hybg_mat_mod.o  psb_d_cuda_csrg_mat_mod.o\
        psb_s_cuda_elg_mat_mod.o   psb_s_cuda_hlg_mat_mod.o \
	psb_s_cuda_hybg_mat_mod.o  psb_s_cuda_csrg_mat_mod.o\
        psb_c_cuda_elg_mat_mod.o   psb_c_cuda_hlg_mat_mod.o \
	psb_c_cuda_hybg_mat_mod.o  psb_c_cuda_csrg_mat_mod.o\
        psb_z_cuda_elg_mat_mod.o   psb_z_cuda_hlg_mat_mod.o \
	psb_z_cuda_hybg_mat_mod.o  psb_z_cuda_csrg_mat_mod.o\
	psb_cuda_env_mod.o psb_cuda_mod.o \
	psb_d_cuda_diag_mat_mod.o\
	psb_d_cuda_hdiag_mat_mod.o psb_s_cuda_hdiag_mat_mod.o\
	psb_s_cuda_dnsg_mat_mod.o psb_d_cuda_dnsg_mat_mod.o \
	psb_c_cuda_dnsg_mat_mod.o psb_z_cuda_dnsg_mat_mod.o \
	dnsdev_mod.o

COBJS= elldev.o hlldev.o diagdev.o hdiagdev.o vectordev.o  ivectordev.o dnsdev.o\
	svectordev.o dvectordev.o cvectordev.o zvectordev.o cuda_util.o \
	fcusparse.o scusparse.o dcusparse.o ccusparse.o zcusparse.o 

OBJS=$(COBJS) $(FOBJS)

lib: objs ilib cudalib spgpulib
	ar cur $(LIBNAME) $(OBJS)
	/bin/cp -p $(LIBNAME) $(LIBDIR)


objs: spgpuinc $(OBJS) iobjs cudaobjs spgpuobjs
	/bin/cp -p *$(.mod) $(MODDIR)
	/bin/cp -p *.h $(INCDIR)

spgpuinc:
	$(MAKE) -C spgpu includes	
spgpuobjs:
	$(MAKE) -C spgpu objs
spgpulib:
	$(MAKE) -C spgpu lib


dnsdev_mod.o hlldev_mod.o elldev_mod.o psb_base_vectordev_mod.o: core_mod.o
psb_d_cuda_vect_mod.o psb_s_cuda_vect_mod.o psb_z_cuda_vect_mod.o psb_c_cuda_vect_mod.o: psb_i_cuda_vect_mod.o
psb_i_cuda_vect_mod.o : psb_vectordev_mod.o psb_cuda_env_mod.o
cusparse_mod.o:  s_cusparse_mod.o  d_cusparse_mod.o  c_cusparse_mod.o  z_cusparse_mod.o 
s_cusparse_mod.o  d_cusparse_mod.o  c_cusparse_mod.o  z_cusparse_mod.o : base_cusparse_mod.o  
psb_d_cuda_hlg_mat_mod.o: hlldev_mod.o psb_d_cuda_vect_mod.o  psb_cuda_env_mod.o
psb_d_cuda_elg_mat_mod.o: elldev_mod.o psb_d_cuda_vect_mod.o
psb_d_cuda_diag_mat_mod.o: diagdev_mod.o psb_d_cuda_vect_mod.o
psb_d_cuda_hdiag_mat_mod.o: hdiagdev_mod.o psb_d_cuda_vect_mod.o
psb_s_cuda_dnsg_mat_mod.o: dnsdev_mod.o psb_s_cuda_vect_mod.o
psb_d_cuda_dnsg_mat_mod.o: dnsdev_mod.o psb_d_cuda_vect_mod.o
psb_c_cuda_dnsg_mat_mod.o: dnsdev_mod.o psb_c_cuda_vect_mod.o
psb_z_cuda_dnsg_mat_mod.o: dnsdev_mod.o psb_z_cuda_vect_mod.o
psb_s_cuda_hlg_mat_mod.o: hlldev_mod.o psb_s_cuda_vect_mod.o  psb_cuda_env_mod.o
psb_s_cuda_elg_mat_mod.o: elldev_mod.o psb_s_cuda_vect_mod.o 
psb_s_cuda_diag_mat_mod.o: diagdev_mod.o psb_s_cuda_vect_mod.o
psb_s_cuda_hdiag_mat_mod.o: hdiagdev_mod.o psb_s_cuda_vect_mod.o 
psb_s_cuda_csrg_mat_mod.o psb_s_cuda_hybg_mat_mod.o: cusparse_mod.o psb_vectordev_mod.o
psb_d_cuda_csrg_mat_mod.o psb_d_cuda_hybg_mat_mod.o: cusparse_mod.o psb_vectordev_mod.o
psb_z_cuda_hlg_mat_mod.o: hlldev_mod.o psb_z_cuda_vect_mod.o  psb_cuda_env_mod.o
psb_z_cuda_elg_mat_mod.o: elldev_mod.o psb_z_cuda_vect_mod.o 
psb_c_cuda_hlg_mat_mod.o: hlldev_mod.o psb_c_cuda_vect_mod.o  psb_cuda_env_mod.o
psb_c_cuda_elg_mat_mod.o: elldev_mod.o psb_c_cuda_vect_mod.o 
psb_c_cuda_csrg_mat_mod.o psb_c_cuda_hybg_mat_mod.o: cusparse_mod.o psb_vectordev_mod.o
psb_z_cuda_csrg_mat_mod.o psb_z_cuda_hybg_mat_mod.o: cusparse_mod.o psb_vectordev_mod.o
psb_vectordev_mod.o: psb_s_vectordev_mod.o psb_d_vectordev_mod.o psb_c_vectordev_mod.o psb_z_vectordev_mod.o psb_i_vectordev_mod.o
psb_i_vectordev_mod.o psb_s_vectordev_mod.o psb_d_vectordev_mod.o psb_c_vectordev_mod.o psb_z_vectordev_mod.o: psb_base_vectordev_mod.o
vectordev.o: cuda_util.o   vectordev.h
elldev.o: elldev.c
dnsdev.o: dnsdev.c
fcusparse.h elldev.c: elldev.h vectordev.h
fcusparse.o scusparse.o dcusparse.o ccusparse.o zcusparse.o : fcusparse.h
fcusparse.o scusparse.o dcusparse.o ccusparse.o zcusparse.o : fcusparse_fct.h
svectordev.o:  svectordev.h vectordev.h
dvectordev.o:  dvectordev.h vectordev.h
cvectordev.o:  cvectordev.h vectordev.h
zvectordev.o:  zvectordev.h vectordev.h
psb_cuda_env_mod.o: base_cusparse_mod.o
psb_cuda_mod.o:  psb_cuda_env_mod.o psb_i_cuda_vect_mod.o\
	psb_d_cuda_vect_mod.o  psb_s_cuda_vect_mod.o\
	psb_z_cuda_vect_mod.o  psb_c_cuda_vect_mod.o\
        psb_d_cuda_elg_mat_mod.o   psb_d_cuda_hlg_mat_mod.o \
	psb_d_cuda_hybg_mat_mod.o  psb_d_cuda_csrg_mat_mod.o\
        psb_s_cuda_elg_mat_mod.o   psb_s_cuda_hlg_mat_mod.o \
	psb_s_cuda_hybg_mat_mod.o  psb_s_cuda_csrg_mat_mod.o\
        psb_c_cuda_elg_mat_mod.o   psb_c_cuda_hlg_mat_mod.o \
	psb_c_cuda_hybg_mat_mod.o  psb_c_cuda_csrg_mat_mod.o\
        psb_z_cuda_elg_mat_mod.o   psb_z_cuda_hlg_mat_mod.o \
	psb_z_cuda_hybg_mat_mod.o  psb_z_cuda_csrg_mat_mod.o\
	psb_d_cuda_diag_mat_mod.o \
	psb_d_cuda_hdiag_mat_mod.o psb_s_cuda_hdiag_mat_mod.o\
	psb_s_cuda_dnsg_mat_mod.o psb_d_cuda_dnsg_mat_mod.o \
	psb_c_cuda_dnsg_mat_mod.o psb_z_cuda_dnsg_mat_mod.o 

iobjs: $(FOBJS)
	$(MAKE) -C impl objs
cudaobjs: $(FOBJS)
	$(MAKE) -C CUDA objs

ilib: objs
	$(MAKE) -C impl lib LIBNAME=$(LIBNAME)
cudalib: objs ilib
	$(MAKE) -C CUDA lib LIBNAME=$(LIBNAME)

clean: cclean iclean cudaclean spgpuclean
	/bin/rm -f  $(FOBJS) *$(.mod) *.a

cclean: 
	/bin/rm -f  $(COBJS) 
iclean:
	$(MAKE) -C impl clean
cudaclean:
	$(MAKE) -C CUDA clean
spgpuclean:
	$(MAKE) -C spgpu clean

veryclean:  clean
