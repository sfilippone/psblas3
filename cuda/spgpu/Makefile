TOP=../..
include $(TOP)/Make.inc
#
# Libraries used
#
UP=..
LIBDIR=$(TOP)/lib
INCDIR=$(TOP)/include
MODDIR=$(TOP)/modules
LIBNAME=$(UP)/libpsb_cuda.a

OBJS=coo.o core.o dia.o ell.o hdia.o hell.o 
CU_INCLUDES=-I$(INCDIR) $(CUDA_INCLUDES)

all: incl objs

objs: incl $(OBJS) iobjs

lib: objs iobjs ilib
	$(AR) $(LIBNAME) $(OBJS)
#	/bin/cp -p $(LIBNAME) $(LIBDIR)

iobjs:
	$(MAKE) -C kernels objs
ilib:
	$(MAKE) -C kernels lib LIBNAME=$(LIBNAME)

incl:
	/bin/cp -p *.h $(INCDIR)

clean: iclean
	/bin/rm -fr $(OBJS) $(LIBNAME)
iclean:
	$(MAKE) -C kernels clean 
.c.o:
	$(CC) $(CCOPT) $(CINCLUDES) $(CDEFINES) $(CU_INCLUDES) -c $<
.cpp.o:
	$(CXX) $(CXXOPT) $(CXXINCLUDES) $(CXXDEFINES) $(CU_INCLUDES) -c $< -o $@

.cu.o:
	$(NVCC)   $(CINCLUDES) $(CDEFINES) $(CUDEFINES) $(CU_INCLUDES) -c $<
