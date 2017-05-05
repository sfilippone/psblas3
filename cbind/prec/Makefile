TOP=../..
include $(TOP)/Make.inc
LIBDIR=$(TOP)/lib
INCLUDEDIR=$(TOP)/include
HERE=..

FINCLUDES=$(FMFLAG). $(FMFLAG)$(HERE) $(FMFLAG)$(INCLUDEDIR) 
CINCLUDES=-I. -I$(HERE) -I$(INCLUDEDIR)

OBJS=psb_prec_cbind_mod.o \
psb_sprec_cbind_mod.o psb_dprec_cbind_mod.o psb_cprec_cbind_mod.o psb_zprec_cbind_mod.o \
psb_c_dprec.o
CMOD=psb_prec_cbind.h psb_c_dprec.h


LIBMOD=psb_prec_cbind_mod$(.mod)  \
psb_sprec_cbind_mod$(.mod) psb_dprec_cbind_mod$(.mod) \
psb_cprec_cbind_mod$(.mod) psb_zprec_cbind_mod$(.mod) 

LOCAL_MODS=$(LIBMOD)
LIBNAME=$(CPRECLIBNAME)	


lib:  $(OBJS) $(CMOD)
	$(AR) $(HERE)/$(LIBNAME) $(OBJS)
	$(RANLIB) $(HERE)/$(LIBNAME)
	/bin/cp -p $(HERE)/$(LIBNAME) $(LIBDIR)
	/bin/cp -p $(LIBMOD) $(CMOD) $(HERE)
	/bin/cp -p $(LIBMOD) $(CMOD) $(INCLUDEDIR)

psb_prec_cbind_mod.o: psb_sprec_cbind_mod.o psb_dprec_cbind_mod.o psb_cprec_cbind_mod.o psb_zprec_cbind_mod.o
veryclean: clean
	/bin/rm -f $(HERE)/$(LIBNAME)

clean:
	/bin/rm -f $(OBJS) $(LOCAL_MODS)

veryclean: clean
