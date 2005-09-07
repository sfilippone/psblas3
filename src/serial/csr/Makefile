include ../../../../Make.inc

#
# The object files
#

FOBJS = dcsrck.o dcsrmm.o dcsrsm.o dcsrmv.o dcsrsv.o dcrnrmi.o  \
        dcrcrupd.o dcocrupd.o dcsrprt.o dcsrmv4.o dcsrmv2.o dcsrmv3.o\
	zcsrck.o zcrnrmi.o zcsrmm.o zsrmv.o zcsrsm.o zsrsv.o \
	zcrcrupd.o zcocrupd.o zcsrprt.o

OBJS=$(FOBJS)

#
# Where the library should go, and how it is called. 
# Note that we are regenerating most of libsparker.a on the fly. 
#LIBDIR=../../LIB
#LIBNAME=libsparker.a	
LIBFILE=$(LIBDIR)/$(LIBNAME)
SPARKERDIR=..
INCDIRS=-I. -I$(SPARKERDIR) -I$(LIBDIR)

#
# No change should be needed below 
#


default: lib

lib: $(OBJS)
	$(AR) $(LIBFILE) $(OBJS) 
	$(RANLIB) $(LIBFILE)

$(OBJS): $(SPARKERDIR)/sparker.fh


clean: cleanobjs

veryclean: cleanobjs

cleanobjs:
	/bin/rm -f $(OBJS) 

