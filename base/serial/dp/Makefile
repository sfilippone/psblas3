include ../../../Make.inc
#
# The object files
#

FOBJS = dcrcr.o dcrjd.o dgblock.o partition.o  \
	dgindex.o djadrp.o djadrp1.o dcsrrp.o dcsrp1.o check_dim.o \
	Max_nnzero.o dcoco.o dcocr.o dcrco.o djdcox.o djdco.o dvtfg.o dgind_tri.o \
	gen_block.o reordvn.o\
	zcoco.o zcocr.o zcrco.o zcrcr.o zcrjd.o zgindex.o zgind_tri.o

# 
# dgind_tri.o
# 


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


lib: $(FOBJS)
	$(AR) $(LIBDIR)/$(LIBNAME) $(OBJS)
	$(RANLIB) $(LIBDIR)/$(LIBNAME)

clean: cleanobjs

veryclean: cleanobjs

cleanobjs:
	/bin/rm -f $(OBJS) 


