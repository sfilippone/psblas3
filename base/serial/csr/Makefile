include ../../../Make.inc

#
# The object files
#

FOBJS = dcsrck.o dcsrmm.o dcsrsm.o dcsrmv.o dcsrsv.o dcrnrmi.o  \
        dcsrmv4.o dcsrmv2.o dcsrmv3.o dcsrrws.o\
	scsrmm.o scsrmv.o scsrmv2.o scsrmv3.o scsrmv4.o scsrsm.o scsrsv.o\
	scrnrmi.o \
	ccrnrmi.o ccsrmm.o ccsrrws.o ccsrsm.o csrmv.o csrsv.o ccsrck.o\
	zcrnrmi.o zcsrmm.o zcsrrws.o zcsrsm.o zsrmv.o zsrsv.o zcsrck.o

OBJS=$(FOBJS)

#
# Where the library should go, and how it is called. 
# Note that we are regenerating most of libsparker.a on the fly. 
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

