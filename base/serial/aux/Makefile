include ../../../Make.inc
#
# The object files
#

FOBJS  = psi_acx_mod.o psi_alcx_mod.o psi_lcx_mod.o psi_tscx.o\
	isr.o isrx.o iasr.o iasrx.o msort_up.o msort_dw.o\
	imsr.o imsrx.o imsru.o iamsort_up.o iamsort_dw.o idot.o inrm2.o\
	dsr.o dsrx.o dasr.o dasrx.o dmsr.o dmsrx.o \
	dmsort_up.o dmsort_dw.o   damsort_up.o damsort_dw.o \
	ssr.o ssrx.o sasr.o sasrx.o smsr.o smsrx.o \
	smsort_up.o smsort_dw.o   samsort_up.o samsort_dw.o \
	clcmp_mod.o clsr.o clsrx.o  \
	calcmp_mod.o calsr.o calsrx.o \
	cacmp_mod.o casr.o casrx.o camsr.o camsrx.o camsort_up.o camsort_dw.o\
	zlcmp_mod.o zlsr.o zlsrx.o  \
	zalcmp_mod.o zalsr.o zalsrx.o \
	zacmp_mod.o zasr.o zasrx.o zamsr.o zamsrx.o zamsort_up.o zamsort_dw.o


OBJS=$(FOBJS)

#
# Where the library should go, and how it is called. 
# 
LIBDIR=../..
INCDIR=../..
MODDIR=../../modules
FINCLUDES=$(FMFLAG). $(FMFLAG)$(MODDIR) $(FMFLAG)$(INCDIR) 

#
# No change should be needed below 
#

default: lib

lib: $(OBJS)
	$(AR) $(LIBDIR)/$(LIBNAME) $(OBJS)
	$(RANLIB) $(LIBDIR)/$(LIBNAME)

clsr.o  clsrx.o:  clcmp_mod.o
calsr.o calsrx.o: calcmp_mod.o
camsort_up.o camsort_dw.o casr.o  casrx.o:  cacmp_mod.o
zlsr.o  zlsrx.o:  zlcmp_mod.o
zalsr.o zalsrx.o: zalcmp_mod.o
zamsort_up.o zamsort_dw.o zasr.o  zasrx.o:  zacmp_mod.o

clean: 
	/bin/rm -f $(OBJS) *$(.mod)

veryclean: clean

