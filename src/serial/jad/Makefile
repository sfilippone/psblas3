include ../../../Make.inc
#
# The object files
#

FOBJS = dcojdupd.o djadmm.o djadmv.o djadsm.o djadsv.o djdnrmi.o djadnr.o djadprt.o\
	djadmv2.o djadmv3.o djadmv4.o

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

clean: cleanobjs

veryclean: cleanobjs

cleanobjs:
	/bin/rm -f $(OBJS) 

