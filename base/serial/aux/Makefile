include ../../../Make.inc
#
# The object files
#

FOBJS  = isr.o isrx.o iasr.o iasrx.o msort_up.o msort_dw.o\
	isaperm.o ibsrch.o issrch.o  imsr.o imsrx.o imsru.o\
	dsr.o dsrx.o dasr.o dasrx.o dmsr.o dmsrx.o  dmsort_up.o dmsort_dw.o \
	zlcmp_mod.o zlsr.o zlsrx.o  \
	zalcmp_mod.o zalsr.o zalsrx.o \
	zacmp_mod.o zasr.o zasrx.o zamsr.o zamsrx.o zamsort_up.o zamsort_dw.o


OBJS=$(FOBJS)

#
# Where the library should go, and how it is called. 
# Note that we are regenerating most of libsparker.a on the fly. 
SPARKERDIR=..
LIBDIR=../..
MODDIR=../../modules
FINCLUDES=$(FMFLAG)$(LIBDIR) $(FMFLAG)$(MODDIR) $(FMFLAG)$(SPARKERDIR)  $(FMFLAG). 

#
# No change should be needed below 
#

default: lib

lib: $(OBJS)
	$(AR) $(LIBDIR)/$(LIBNAME) $(OBJS)
	$(RANLIB) $(LIBDIR)/$(LIBNAME)
zlsr.o  zlsrx.o:  zlcmp_mod.o
zalsr.o zalsrx.o: zalcmp_mod.o
zasr.o  zasrx.o:  zacmp_mod.o

clean: 
	/bin/rm -f $(OBJS) *$(.mod)

veryclean: clean

