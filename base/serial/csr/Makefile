include ../../../Make.inc

#
# The object files
#

FOBJS = dcsrck.o dcsrmm.o dcsrsm.o dcsrmv.o dcsrsv.o dcrnrmi.o  \
        dcsrprt.o dcsrmv4.o dcsrmv2.o dcsrmv3.o dcsrrws.o\
	zcrnrmi.o zcsrmm.o zcsrrws.o zcsrsm.o zsrmv.o zsrsv.o zcsrck.o

OBJS=$(FOBJS)

#
# Where the library should go, and how it is called. 
# Note that we are regenerating most of libsparker.a on the fly. 
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

