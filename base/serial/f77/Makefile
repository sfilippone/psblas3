include ../../../Make.inc

#
# The object files
#
FOBJS =  daxpby.o 	saxpby.o \
	caxpby.o  zaxpby.o smmp.o


OBJS=$(FOBJS)

#
# Where the library should go, and how it is called. 
# Note that we are regenerating most of libsparker.a on the fly. 
SPARKERDIR=..
LIBDIR=../..
MODDIR=../../modules
FINCLUDES=$(FMFLAG)$(LIBDIR) $(FMFLAG)$(MODDIR) $(FMFLAG)$(SPARKERDIR)  $(FMFLAG). 
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

