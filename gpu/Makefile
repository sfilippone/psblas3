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
CINCLUDES=$(SPGPU_INCLUDES) $(CUDA_INCLUDES)
LIBNAME=libpsb_gpu.a


FOBJS=cusparse_mod.o base_cusparse_mod.o \
 s_cusparse_mod.o  d_cusparse_mod.o  c_cusparse_mod.o  z_cusparse_mod.o \
	psb_vectordev_mod.o core_mod.o \
	psb_s_vectordev_mod.o psb_d_vectordev_mod.o psb_i_vectordev_mod.o\
	psb_c_vectordev_mod.o psb_z_vectordev_mod.o psb_base_vectordev_mod.o \
        elldev_mod.o hlldev_mod.o diagdev_mod.o hdiagdev_mod.o \
	psb_i_gpu_vect_mod.o \
	psb_d_gpu_vect_mod.o  psb_s_gpu_vect_mod.o\
	psb_z_gpu_vect_mod.o  psb_c_gpu_vect_mod.o\
        psb_d_elg_mat_mod.o   psb_d_hlg_mat_mod.o \
	psb_d_hybg_mat_mod.o  psb_d_csrg_mat_mod.o\
        psb_s_elg_mat_mod.o   psb_s_hlg_mat_mod.o \
	psb_s_hybg_mat_mod.o  psb_s_csrg_mat_mod.o\
        psb_c_elg_mat_mod.o   psb_c_hlg_mat_mod.o \
	psb_c_hybg_mat_mod.o  psb_c_csrg_mat_mod.o\
        psb_z_elg_mat_mod.o   psb_z_hlg_mat_mod.o \
	psb_z_hybg_mat_mod.o  psb_z_csrg_mat_mod.o\
	psb_gpu_env_mod.o psb_gpu_mod.o \
	psb_d_diag_mat_mod.o\
	psb_d_hdiag_mat_mod.o psb_s_hdiag_mat_mod.o\
	psb_s_dnsg_mat_mod.o psb_d_dnsg_mat_mod.o \
	psb_c_dnsg_mat_mod.o psb_z_dnsg_mat_mod.o \
	dnsdev_mod.o

COBJS= elldev.o hlldev.o diagdev.o hdiagdev.o vectordev.o  ivectordev.o dnsdev.o\
	svectordev.o dvectordev.o cvectordev.o zvectordev.o cuda_util.o \
	fcusparse.o scusparse.o dcusparse.o ccusparse.o zcusparse.o 

OBJS=$(COBJS) $(FOBJS)

lib: objs

objs: $(OBJS) iobjs cudaobjs
	/bin/cp -p *$(.mod) $(MODDIR)
	/bin/cp -p *.h $(INCDIR)

lib: ilib cudalib
	ar cur $(LIBNAME) $(OBJS)
	/bin/cp -p $(LIBNAME) $(LIBDIR)

dnsdev_mod.o hlldev_mod.o elldev_mod.o psb_base_vectordev_mod.o: core_mod.o
psb_d_gpu_vect_mod.o psb_s_gpu_vect_mod.o psb_z_gpu_vect_mod.o psb_c_gpu_vect_mod.o: psb_i_gpu_vect_mod.o
psb_i_gpu_vect_mod.o : psb_vectordev_mod.o psb_gpu_env_mod.o
cusparse_mod.o:  s_cusparse_mod.o  d_cusparse_mod.o  c_cusparse_mod.o  z_cusparse_mod.o 
s_cusparse_mod.o  d_cusparse_mod.o  c_cusparse_mod.o  z_cusparse_mod.o : base_cusparse_mod.o  
psb_d_hlg_mat_mod.o: hlldev_mod.o psb_d_gpu_vect_mod.o  psb_gpu_env_mod.o
psb_d_elg_mat_mod.o: elldev_mod.o psb_d_gpu_vect_mod.o
psb_d_diag_mat_mod.o: diagdev_mod.o psb_d_gpu_vect_mod.o
psb_d_hdiag_mat_mod.o: hdiagdev_mod.o psb_d_gpu_vect_mod.o
psb_s_dnsg_mat_mod.o: dnsdev_mod.o psb_s_gpu_vect_mod.o
psb_d_dnsg_mat_mod.o: dnsdev_mod.o psb_d_gpu_vect_mod.o
psb_c_dnsg_mat_mod.o: dnsdev_mod.o psb_c_gpu_vect_mod.o
psb_z_dnsg_mat_mod.o: dnsdev_mod.o psb_z_gpu_vect_mod.o
psb_s_hlg_mat_mod.o: hlldev_mod.o psb_s_gpu_vect_mod.o  psb_gpu_env_mod.o
psb_s_elg_mat_mod.o: elldev_mod.o psb_s_gpu_vect_mod.o 
psb_s_diag_mat_mod.o: diagdev_mod.o psb_s_gpu_vect_mod.o
psb_s_hdiag_mat_mod.o: hdiagdev_mod.o psb_s_gpu_vect_mod.o 
psb_s_csrg_mat_mod.o psb_s_hybg_mat_mod.o: cusparse_mod.o psb_vectordev_mod.o
psb_d_csrg_mat_mod.o psb_d_hybg_mat_mod.o: cusparse_mod.o psb_vectordev_mod.o
psb_z_hlg_mat_mod.o: hlldev_mod.o psb_z_gpu_vect_mod.o  psb_gpu_env_mod.o
psb_z_elg_mat_mod.o: elldev_mod.o psb_z_gpu_vect_mod.o 
psb_c_hlg_mat_mod.o: hlldev_mod.o psb_c_gpu_vect_mod.o  psb_gpu_env_mod.o
psb_c_elg_mat_mod.o: elldev_mod.o psb_c_gpu_vect_mod.o 
psb_c_csrg_mat_mod.o psb_c_hybg_mat_mod.o: cusparse_mod.o psb_vectordev_mod.o
psb_z_csrg_mat_mod.o psb_z_hybg_mat_mod.o: cusparse_mod.o psb_vectordev_mod.o
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
psb_gpu_env_mod.o: base_cusparse_mod.o
psb_gpu_mod.o:  psb_gpu_env_mod.o psb_i_gpu_vect_mod.o\
	psb_d_gpu_vect_mod.o  psb_s_gpu_vect_mod.o\
	psb_z_gpu_vect_mod.o  psb_c_gpu_vect_mod.o\
        psb_d_elg_mat_mod.o   psb_d_hlg_mat_mod.o \
	psb_d_hybg_mat_mod.o  psb_d_csrg_mat_mod.o\
        psb_s_elg_mat_mod.o   psb_s_hlg_mat_mod.o \
	psb_s_hybg_mat_mod.o  psb_s_csrg_mat_mod.o\
        psb_c_elg_mat_mod.o   psb_c_hlg_mat_mod.o \
	psb_c_hybg_mat_mod.o  psb_c_csrg_mat_mod.o\
        psb_z_elg_mat_mod.o   psb_z_hlg_mat_mod.o \
	psb_z_hybg_mat_mod.o  psb_z_csrg_mat_mod.o\
	psb_d_diag_mat_mod.o \
	psb_d_hdiag_mat_mod.o psb_s_hdiag_mat_mod.o\
	psb_s_dnsg_mat_mod.o psb_d_dnsg_mat_mod.o \
	psb_c_dnsg_mat_mod.o psb_z_dnsg_mat_mod.o 

iobjs: $(FOBJS)
	$(MAKE) -C impl objs
cudaobjs: $(FOBJS)
	$(MAKE) -C CUDA objs

ilib: objs
	$(MAKE) -C impl lib LIBNAME=$(LIBNAME)
cudalib: objs
	$(MAKE) -C CUDA lib LIBNAME=$(LIBNAME)

clean: cclean iclean cudaclean
	/bin/rm -f  $(FOBJS) *$(.mod) *.a

cclean: 
	/bin/rm -f  $(COBJS) 
iclean:
	$(MAKE) -C impl clean
cudaclean:
	$(MAKE) -C CUDA clean

veryclean:  clean
