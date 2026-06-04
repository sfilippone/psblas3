TOP=../..
include $(TOP)/Make.inc
LIBDIR=$(TOP)/lib
INCLUDEDIR=$(TOP)/include
MODDIR=$(TOP)/modules
HERE=..

FINCLUDES=$(FMFLAG). $(FMFLAG)$(HERE) $(FMFLAG)$(MODDIR) 
CINCLUDES=-I. -I$(HERE) -I$(INCLUDEDIR)

OBJS=psb_base_linsolve_cbind_mod.o \
psb_slinsolve_cbind_mod.o psb_dlinsolve_cbind_mod.o psb_clinsolve_cbind_mod.o psb_zlinsolve_cbind_mod.o 
CMOD=psb_linsolve_cbind.h


LIBNAME=$(CLINSOLVELIBNAME)	


objs:  $(OBJS) $(CMOD)
	/bin/cp -p *$(.mod) $(CMOD) $(HERE)
lib:  objs
	$(AR) $(HERE)/$(LIBNAME) $(OBJS)
	$(RANLIB) $(HERE)/$(LIBNAME)
	/bin/cp -p $(HERE)/$(LIBNAME) $(LIBDIR)

psb_slinsolve_cbind_mod.o psb_dlinsolve_cbind_mod.o psb_clinsolve_cbind_mod.o psb_zlinsolve_cbind_mod.o: psb_base_linsolve_cbind_mod.o
veryclean: clean
	/bin/rm -f $(HERE)/$(LIBNAME)

clean:
	/bin/rm -f $(OBJS) *$(.mod)

veryclean: clean
