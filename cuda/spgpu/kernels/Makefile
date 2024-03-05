TOP=../../..
include $(TOP)/Make.inc
#
# Libraries used
#
LIBDIR=$(TOP)/lib
INCDIR=$(TOP)/include
MODDIR=$(TOP)/modules
UP=..
LIBNAME=$(UP)/libspgpu.a
CINCLUDES=-I$(INCDIR)

OBJS=cabs.o  camax.o  casum.o  caxpby.o  caxy.o  cdot.o  cgath.o \
  cnrm2.o cscal.o cscat.o csetscal.o cabgdxyz.o\
  dabs.o damax.o dasum.o daxpby.o daxy.o ddot.o dgath.o dabgdxyz.o\
  dia_cspmv.o dia_dspmv.o dia_sspmv.o dia_zspmv.o dnrm2.o \
  dscal.o dscat.o dsetscal.o ell_ccsput.o ell_cspmv.o \
  ell_dcsput.o ell_dspmv.o ell_scsput.o ell_sspmv.o ell_zcsput.o ell_zspmv.o \
  hdia_cspmv.o hdia_dspmv.o hdia_sspmv.o hdia_zspmv.o hell_cspmv.o hell_dspmv.o \
  hell_sspmv.o hell_zspmv.o igath.o iscat.o isetscal.o sabs.o samax.o sasum.o \
  saxpby.o saxy.o sdot.o sgath.o snrm2.o sscal.o sscat.o ssetscal.o zabs.o zamax.o sabgdxyz.o\
  zasum.o zaxpby.o zaxy.o zdot.o zgath.o znrm2.o zscal.o zscat.o zsetscal.o zabgdxyz.o \
  sxyzw.o cxyzw.o dxyzw.o zxyzw.o

objs: $(OBJS)
lib: objs
	ar cur $(UP)/$(LIBNAME) $(OBJS)


clean:
	/bin/rm -fr $(OBJS)
.cu.o:
	$(NVCC)   $(CINCLUDES) $(CDEFINES) $(CUDEFINES) $(CUDA_INCLUDES) -c $<

