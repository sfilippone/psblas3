TOP=../..
include $(TOP)/Make.inc
LIBDIR=$(TOP)/lib
INCLUDEDIR=$(TOP)/include
MODDIR=$(TOP)/modules
HERE=..

FINCLUDES=$(FMFLAG). $(FMFLAG)$(HERE) $(FMFLAG)$(MODDIR)
CINCLUDES=-I. -I$(HERE) -I$(INCLUDEDIR)

OBJS=psb_util_cbind_mod.o \
psb_c_util_cbind_mod.o \
psb_d_util_cbind_mod.o \
psb_s_util_cbind_mod.o \
psb_z_util_cbind_mod.o
CMOD=psb_util_cbind.h psb_c_cutil.h psb_c_zutil.h psb_c_dutil.h psb_c_sutil.h

LIBNAME=$(CUTILLIBNAME)


lib:  $(OBJS) $(CMOD)
	$(AR) $(HERE)/$(LIBNAME) $(OBJS)
	$(RANLIB) $(HERE)/$(LIBNAME)
	/bin/cp -p $(HERE)/$(LIBNAME) $(LIBDIR)
	/bin/cp -p *$(.mod) $(CMOD) $(HERE)

psb_util_cbind_mod.o: psb_c_util_cbind_mod.o psb_d_util_cbind_mod.o psb_s_util_cbind_mod.o psb_z_util_cbind_mod.o
veryclean: clean
	/bin/rm -f $(HERE)/$(LIBNAME)

clean:
	/bin/rm -f $(OBJS) *$(.mod)

veryclean: clean
