INSTALLDIR=../..
INCDIR=$(INSTALLDIR)/include
MODDIR=$(INSTALLDIR)/modules/
include $(INCDIR)/Make.inc.psblas
#
# Libraries used
LIBDIR=$(INSTALLDIR)/lib
# BLASLIBDIR=/opt/share/sdk/intel/nvidia_hpc_sdk/Linux_x86_64/24.3/compilers/lib
# GKLIBDIR=/home/jalmerol/GKlib/build/Linux-x86_64
PSBLAS_LIB= -L$(LIBDIR) -lpsb_util -lpsb_linsolve -lpsb_prec -lpsb_base
LDLIBS= $(PSBLDLIBS)
#
# Compilers and such
#
CCOPT= -g
FINCLUDES=$(FMFLAG)$(MODDIR) $(FMFLAG).


EXEDIR=./runs

all: runsd psb_d_pde3d psb_s_pde3d psb_d_pde2d psb_s_pde2d psb_d_pde_nest_psblas psb_d_pde_nest_krylov

runsd:
	(if test ! -d runs ; then mkdir runs; fi)

psb_d_pde3d: psb_d_pde3d.o 
	$(FLINK)  psb_d_pde3d.o -o psb_d_pde3d $(PSBLAS_LIB) $(LDLIBS)
	/bin/mv psb_d_pde3d $(EXEDIR)

psb_s_pde3d: psb_s_pde3d.o 
	$(FLINK) psb_s_pde3d.o  -o psb_s_pde3d $(PSBLAS_LIB) $(LDLIBS)
	/bin/mv psb_s_pde3d $(EXEDIR)

psb_d_pde2d: psb_d_pde2d.o 
	$(FLINK)  psb_d_pde2d.o -o psb_d_pde2d $(PSBLAS_LIB) $(LDLIBS)
	/bin/mv psb_d_pde2d $(EXEDIR)


psb_s_pde2d: psb_s_pde2d.o 
	$(FLINK) psb_s_pde2d.o  -o psb_s_pde2d $(PSBLAS_LIB) $(LDLIBS)
	/bin/mv psb_s_pde2d $(EXEDIR)

psb_d_pde_nest_psblas: psb_d_pde_nest_psblas.o
	$(FLINK) psb_d_pde_nest_psblas.o -o psb_d_pde_nest_psblas $(PSBLAS_LIB) $(LDLIBS)
	/bin/mv psb_d_pde_nest_psblas $(EXEDIR)

psb_d_pde_nest_krylov: psb_d_pde_nest_krylov.o
	$(FLINK) psb_d_pde_nest_krylov.o -o psb_d_pde_nest_krylov $(PSBLAS_LIB) $(LDLIBS)
	/bin/mv psb_d_pde_nest_krylov $(EXEDIR)

clean: 
	/bin/rm -f psb_d_pde3d.o psb_d_oacc_pde3d.o psb_s_pde3d.o  psb_d_pde2d.o psb_s_pde2d.o psb_d_pde_nest_psblas.o psb_d_pde_nest_krylov.o *$(.mod) \
	$(EXEDIR)/psb_d_pde3d $(EXEDIR)/psb_s_pde3d $(EXEDIR)/psb_d_pde2d $(EXEDIR)/psb_s_pde2d $(EXEDIR)/psb_d_pde_nest_psblas $(EXEDIR)/psb_d_pde_nest_krylov
verycleanlib: 
	(cd ../..; make veryclean)
lib:
	(cd ../../; make library)



