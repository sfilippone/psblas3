TOP=../../..
include $(TOP)/Make.inc
#
# Libraries used
#
LIBDIR=$(TOP)/lib
INCDIR=$(TOP)/include
MODDIR=$(TOP)/modules

LIBNAME=libpsb_cuda.a

CINCLUDES=-I$(INCDIR)
CU_INCLUDES=-I$(INCDIR) $(CUDA_INCLUDES)

OBJS= cabs.o camax.o casum.o caxpby.o ccolspan.o cmlt.o caxy.o cdot.o \
        cgath.o cnrm2.o cscal.o cscat.o csetscal.o cupd_xyz.o cxyzw.o \
      dabs.o damax.o dasum.o daxpby.o dcolspan.o dmlt.o daxy.o ddot.o \
        dgath.o dnrm2.o dscal.o dscat.o dsetscal.o dupd_xyz.o dxyzw.o \
      sabs.o samax.o sasum.o saxpby.o scolspan.o smlt.o saxy.o sdot.o \
        sgath.o snrm2.o sscal.o sscat.o ssetscal.o supd_xyz.o sxyzw.o \
      zabs.o zamax.o zasum.o zaxpby.o zcolspan.o zmlt.o zaxy.o zdot.o \
        zgath.o znrm2.o zscal.o zscat.o zsetscal.o zupd_xyz.o zxyzw.o \
      igath.o iscat.o isetscal.o \
      ell_ccsput.o ell_dcsput.o ell_scsput.o ell_zcsput.o \
      ell_cspmv.o ell_dspmv.o ell_sspmv.o ell_zspmv.o \
      dia_cspmv.o dia_dspmv.o dia_sspmv.o dia_zspmv.o \
      hdia_cspmv.o hdia_dspmv.o hdia_sspmv.o hdia_zspmv.o \
      hell_cspmv.o hell_dspmv.o hell_sspmv.o hell_zspmv.o \

objs: $(OBJS)

lib: objs
	$(AR) ../../$(LIBNAME) $(OBJS)

clean:
	/bin/rm -fr $(OBJS)

.cu.o:
	$(NVCC) $(CINCLUDES) $(CDEFINES) $(CUDEFINES) $(CU_INCLUDES) -c $<

