include ../../../Make.inc

#
# The object files
#
FOBJS = ccoonrmi.o ccoomm.o ccoomv.o ccoosm.o ccoosv.o ccoorws.o\
	zcoomm.o zcoomv.o zcoonrmi.o zcoorws.o zcoosm.o zcoosv.o

#	dcoonrmi.o dcoomm.o dcoomv.o dcoosm.o dcoosv.o dcoorws.o\
#scoonrmi.o scoomm.o scoomv.o scoosm.o scoosv.o scoorws.o\


OBJS=$(FOBJS)

#
# Where the library should go, and how it is called. 
# Note that we are regenerating most of libsparker.a on the fly. 
#LIBNAME=libsparker.a	
SPARKERDIR=..
LIBDIR=../..
MODDIR=../../modules
FINCLUDES=$(FMFLAG)$(LIBDIR) $(FMFLAG)$(MODDIR) $(FMFLAG)$(SPARKERDIR)  $(FMFLAG). 
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

