include ../../../Make.inc
#
# The object files
#

FOBJS = djadmm.o djadmv.o djadsm.o djadsv.o djdnrmi.o djadnr.o djadprt.o\
	djadmv2.o djadmv3.o djadmv4.o djadrws.o djdrws.o

OBJS=$(FOBJS)

#
# Where the library should go, and how it is called. 
# Note that we are regenerating most of libsparker.a on the fly. 
SPARKERDIR=..
LIBDIR = ../..
MODDIR = ../../modules
INCDIRS = -I $(LIBDIR) -I $(MODDIR) -I . -I$(SPARKERDIR)
#LIBNAME=libsparker.a	
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

