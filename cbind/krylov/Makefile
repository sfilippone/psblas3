TOP=../..
include $(TOP)/Make.inc
LIBDIR=$(TOP)/lib
INCLUDEDIR=$(TOP)/include
MODDIR=$(TOP)/modules
HERE=..

FINCLUDES=$(FMFLAG). $(FMFLAG)$(HERE) $(FMFLAG)$(MODDIR) 
CINCLUDES=-I. -I$(HERE) -I$(INCLUDEDIR)

OBJS=psb_base_krylov_cbind_mod.o \
psb_skrylov_cbind_mod.o psb_dkrylov_cbind_mod.o psb_ckrylov_cbind_mod.o psb_zkrylov_cbind_mod.o 
CMOD=psb_krylov_cbind.h


LIBNAME=$(CKRYLOVLIBNAME)	


lib:  $(OBJS) $(CMOD)
	$(AR) $(HERE)/$(LIBNAME) $(OBJS)
	$(RANLIB) $(HERE)/$(LIBNAME)
	/bin/cp -p $(HERE)/$(LIBNAME) $(LIBDIR)
	/bin/cp -p *$(.mod) $(CMOD) $(HERE)

psb_skrylov_cbind_mod.o psb_dkrylov_cbind_mod.o psb_ckrylov_cbind_mod.o psb_zkrylov_cbind_mod.o: psb_base_krylov_cbind_mod.o
veryclean: clean
	/bin/rm -f $(HERE)/$(LIBNAME)

clean:
	/bin/rm -f $(OBJS) *$(.mod)

veryclean: clean
