INSTALLDIR=../..
INCDIR=$(INSTALLDIR)/include
MODDIR=$(INSTALLDIR)/modules/
include $(INCDIR)/Make.inc.psblas
#
# Libraries used
LIBDIR=$(INSTALLDIR)/lib
PSBLAS_LIB= -L$(LIBDIR) -lpsb_util -lpsb_krylov -lpsb_prec -lpsb_base 
LDLIBS=$(PSBLDLIBS)
#
# Compilers and such
#
CCOPT= -g
FINCLUDES=$(FMFLAG)$(MODDIR) $(FMFLAG).

LINKOPT=
EXEDIR=./runs

all: psb_d_pde3d psb_s_pde3d psb_d_pde2d psb_s_pde2d

psb_d_pde3d.o: psb_d_scsr_mat_mod.o 
psb_d_pde3d: psb_d_pde3d.o psb_d_scsr_mat_mod.o psb_d_scsr_impl.o 
	$(FLINK) $(LINKOPT) psb_d_pde3d.o psb_d_scsr_mat_mod.o psb_d_scsr_impl.o -o psb_d_pde3d $(PSBLAS_LIB) $(LDLIBS)
	/bin/mv psb_d_pde3d $(EXEDIR)


psb_s_pde3d: psb_s_pde3d.o 
	$(FLINK) $(LINKOPT) psb_s_pde3d.o  -o psb_s_pde3d $(PSBLAS_LIB) $(LDLIBS)
	/bin/mv psb_s_pde3d $(EXEDIR)

psb_d_pde2d: psb_d_pde2d.o 
	$(FLINK) $(LINKOPT)  psb_d_pde2d.o -o psb_d_pde2d $(PSBLAS_LIB) $(LDLIBS)
	/bin/mv psb_d_pde2d $(EXEDIR)


psb_s_pde2d: psb_s_pde2d.o 
	$(FLINK) $(LINKOPT) psb_s_pde2d.o  -o psb_s_pde2d $(PSBLAS_LIB) $(LDLIBS)
	/bin/mv psb_s_pde2d $(EXEDIR)


clean: 
	/bin/rm -f psb_d_pde3d.o psb_s_pde3d.o  psb_d_pde2d.o psb_s_pde2d.o  *$(.mod) \
	$(EXEDIR)/psb_d_pde3d $(EXEDIR)/psb_s_pde3d $(EXEDIR)/psb_d_pde2d $(EXEDIR)/psb_s_pde2d 
verycleanlib: 
	(cd ../..; make veryclean)
lib:
	(cd ../../; make library)



