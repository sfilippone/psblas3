TOP=../..
include $(TOP)/Make.inc
LIBDIR=$(TOP)lib
INCLUDEDIR=$(TOP)/include
MODDIR=$(TOP)/modules
HERE=..

FINCLUDES=$(FMFLAG). $(FMFLAG)$(HERE) $(FMFLAG)$(MODDIR)
CINCLUDES=-I. -I$(HERE) -I$(INCLUDEDIR)

FOBJS= psb_c_util_cbind_mod.o \
		  psb_d_util_cbind_mod.o \
			psb_s_util_cbind_mod.o \
			psb_z_util_cbind_mod.o

# COBJS= psb_c_sutil.o psb_c_dutil.o psb_c_cutil.o psb_c_zutil.o
CMOD= psb_util_cbind.h psb_c_cutil.h psb_c_zutil.h psb_c_dutil.h psb_c_sutil.h
OBJS=$(FOBJS) $(COBJS)

LIBNAME=$(CBINDLIBNAME)


lib:  $(OBJS) $(CMOD)
	$(AR) $(HERE)/$(LIBNAME) $(OBJS)
	$(RANLIB) $(HERE)/$(LIBNAME)
	/bin/cp -p  *$(.mod) $(CMOD) $(HERE)

# $(COBJS): $(CMOD)

psb_util_cbind_mod.o: psb_c_util_cbind_mod.o \
		  psb_d_util_cbind_mod.o \
			psb_s_util_cbind_mod.o \
			psb_z_util_cbind_mod.o

veryclean: clean
	/bin/rm -f $(HERE)/$(LIBNAME)

clean:
	/bin/rm -f $(OBJS) *$(.mod)

veryclean: clean
