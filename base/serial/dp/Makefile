include ../../../Make.inc
#
# The object files
#

XOBJS = scrjd.o   ccrjd.o zcrjd.o
FOBJS = partition.o dgblock.o  dvtfg.o \
	check_dim.o \
	Max_nnzero.o \
	gen_block.o\
	scrco.o scrcr.o scocr.o scoco.o sgindex.o sgind_tri.o\
	ccoco.o ccocr.o ccrco.o ccrcr.o cgindex.o cgind_tri.o\
	zcoco.o zcocr.o zcrco.o zcrcr.o zgindex.o zgind_tri.o\
	$(XOBJS)

#dcoco.o dcocr.o dcrco.o djdcox.o djdco.o dgind_tri.o \
#dcrcr.o 
#dgindex.o djadrp.o djadrp1.o dcsrrp.o dcsrp1.o 
#dcrjd.o

# 
# dgind_tri.o
# 


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


lib: xobjs $(FOBJS)
	$(AR) $(LIBDIR)/$(LIBNAME) $(OBJS)
	$(RANLIB) $(LIBDIR)/$(LIBNAME)

xobjs:
	(make $(XOBJS)  FCOPT="$(FCOPT) $(EXTRA_OPT)")

clean: cleanobjs

veryclean: cleanobjs

cleanobjs:
	/bin/rm -f $(OBJS) 


