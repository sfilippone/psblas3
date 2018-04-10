TOP=../..
include $(TOP)/Make.inc
LIBDIR=$(TOP)/lib
INCLUDEDIR=$(TOP)/include
MODDIR=$(TOP)/modules
HERE=..

FINCLUDES=$(FMFLAG). $(FMFLAG)$(HERE) $(FMFLAG)$(MODDIR) 
CINCLUDES=-I. -I$(HERE) -I$(INCLUDEDIR)

OBJS=psb_prec_cbind_mod.o \
psb_sprec_cbind_mod.o psb_dprec_cbind_mod.o psb_cprec_cbind_mod.o psb_zprec_cbind_mod.o \
psb_c_sprec.o psb_c_dprec.o psb_c_cprec.o psb_c_zprec.o
CMOD=psb_prec_cbind.h psb_c_sprec.h psb_c_dprec.h psb_c_cprec.h psb_c_zprec.h



LIBNAME=$(CPRECLIBNAME)	


lib:  $(OBJS) $(CMOD)
	$(AR) $(HERE)/$(LIBNAME) $(OBJS)
	$(RANLIB) $(HERE)/$(LIBNAME)
	/bin/cp -p $(HERE)/$(LIBNAME) $(LIBDIR)
	/bin/cp -p  *$(.mod) $(CMOD) $(HERE)

psb_prec_cbind_mod.o: psb_sprec_cbind_mod.o psb_dprec_cbind_mod.o psb_cprec_cbind_mod.o psb_zprec_cbind_mod.o
veryclean: clean
	/bin/rm -f $(HERE)/$(LIBNAME)

clean:
	/bin/rm -f $(OBJS)  *$(.mod)

veryclean: clean
