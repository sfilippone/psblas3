include ../../../Make.inc

#
# The object files
#
FOBJS = dcooprt.o  dcoonrmi.o dcoomm.o dcoomv.o dcoosm.o dcoosv.o dcoorws.o\
	zcoomm.o zcoomv.o zcoonrmi.o zcoorws.o zcoosm.o zcoosv.o


OBJS=$(FOBJS)

#
# Where the library should go, and how it is called. 
# Note that we are regenerating most of libsparker.a on the fly. 
#LIBNAME=libsparker.a	
SPARKERDIR=..
LIBDIR = ../..
MODDIR = ../../modules
INCDIRS = -I $(LIBDIR) -I $(MODDIR) -I . -I$(SPARKERDIR)
LIBFILE=$(LIBDIR)/$(LIBNAME)

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

