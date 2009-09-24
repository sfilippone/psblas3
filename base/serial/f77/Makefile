include ../../../Make.inc

#
# The object files
#
FOBJS =  daxpby.o 	saxpby.o \
	sgelp.o\
	caxpby.o clpupd.o ccsmm.o cswmm.o ccsnmi.o ccsrws.o\
	cswsm.o ccssm.o cgelp.o\
	zcsnmi.o  zaxpby.o zcsmm.o zcssm.o zswmm.o zswsm.o\
	zcsrws.o zgelp.o zlpupd.o

#dcsmm.o dcsnmi.o dcsrp.o dcssm.o \
#	dgelp.o dlpupd.o dswmm.o \
#	dswsm.o smmp.o dcsrws.o \
#slpupd.o scsmm.o sswmm.o scsnmi.o scsrws.o\
#sswsm.o scssm.o 


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

