include ../../../../Make.inc

#
# The object files
#

#FOBJS = dcsrck.o dcsrmm.o dcsrsm.o dsrmv.o dsrsv.o dcrnrmi.o dcrrs.o \
#        dcrupdate.o dcooprt.o
FOBJS =  dcooprt.o  dcoonrmi.o dcoomm.o dcoomv.o dcoosm.o dcoosv.o\
	 zcoonrmi.o zcoomm.o zcoomv.o zcoosm.o zcoosv.o zcooprt.o

#zcsrck.o zcrnrmi.o zcsrmm.o zsrmv.o zcsrsm.o zsrsv.o  

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

