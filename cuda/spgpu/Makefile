TOP=../..
include $(TOP)/Make.inc
#
# Libraries used
#
LIBDIR=$(TOP)/lib
INCDIR=$(TOP)/include
MODDIR=$(TOP)/modules
LIBNAME=libspgpu.a

OBJS=coo.o core.o dia.o ell.o hdia.o hell.o 

all: includes objs

objs: $(OBJS) iobjs

lib: objs iobjs ilib
	ar cur $(LIBNAME) $(OBJS)
	/bin/cp -p $(LIBNAME) $(LIBDIR)

iobjs:
	$(MAKE) -C kernels objs
ilib:
	$(MAKE) -C kernels lib LIBNAME=$(LIBNAME)

includes:
	/bin/cp -p *.h $(INCDIR)

clean: iclean
	/bin/rm -fr $(OBJS) $(LIBNAME)
iclean:
	$(MAKE) -C kernels clean 
.c.o:
	$(CC) $(CCOPT) $(CINCLUDES) $(CDEFINES) $(CUDA_INCLUDES) -c $<
.cpp.o:
	$(CXX) $(CXXOPT) $(CXXINCLUDES) $(CXXDEFINES) $(CUDA_INCLUDES) -c $< -o $@

.cu.o:
	$(NVCC)   $(CINCLUDES) $(CDEFINES) $(CUDEFINES) $(CUDA_INCLUDES) -c $<
