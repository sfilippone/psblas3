include ../../../Make.inc

#
# The object files
#
FOBJS = daxpby.o dcsmm.o dcsnmi.o dcsrp.o dcssm.o \
	dgelp.o dlpupd.o dswmm.o dswprt.o \
	dswsm.o smmp.o dcsrws.o \
	zcsnmi.o  zaxpby.o zcsmm.o zcssm.o zswmm.o zswsm.o\
	zcsrws.o zgelp.o zlpupd.o

OBJS=$(FOBJS)

#
# Where the library should go, and how it is called. 
# Note that we are regenerating most of libsparker.a on the fly. 
LIBDIR=../../../lib
#LIBNAME=libsparker.a	
LIBFILE=$(LIBDIR)/$(LIBNAME)
INCDIRS=-I. -I$(LIBDIR)

#
# No change should be needed below 
#


default: lib

lib: $(OBJS)
	$(AR) $(LIBDIR)/$(LIBNAME) $(OBJS)
	$(RANLIB) $(LIBDIR)/$(LIBNAME)

clean: cleanobjs

veryclean: cleanobjs

cleanobjs:
	/bin/rm -f $(OBJS) 

